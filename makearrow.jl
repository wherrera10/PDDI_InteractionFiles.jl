#!/usr/bin/env julia

# build_pddi_arrow.jl
#
# Builds a single Arrow table combining three PDDI (potential drug-drug
# interaction) sources into a common schema:
#
#     Name_1      :: String   drug name
#     Name_2      :: String   drug name
#     is_severe   :: Bool     true = serious/high-priority, false = minor/precaution
#     description :: String   free-text description of the interaction
#     source      :: String   which dataset this row came from (handy for
#                              debugging / trust-weighting; drop it later
#                              with `select!(df, Not(:source))` if you don't
#                              want it in the final table)
#
# Verified against the actual files in dbmi-pitt/public-PDDI-analysis:
#
#   1) FrenchDDI => PDDI-Datasets/FrenchNatnlFormulary/frenchDDI.csv
#        tab-delimited, 9 columns, HAS a header, HAS a real is_severecode
#        in `niveau`: CI (contraindicated), AD (avoid if possible),
#        PE (precaution), PC (take into account). CI/AD -> serious (true),
#        PE/PC -> minor (false). Translated to English as frenchDDI_en.csv
#
#   2) ONC High-Priority => PDDI-Datasets/ONC-High-Priority/ONC_High_Priority_Mapped.csv
#        '$'-delimited, NO header, NO is_severecolumn, NO free-text
#        description. 5 fields per line because each record ends with a
#        trailing '$' (giving one empty trailing field). Because this file
#        *is* the FDA/ONC high-priority interaction list, every row is
#        treated as is_severe= true by definition.
#
#   3) db_drug_interactions.csv => the common Kaggle/DrugBank export.
#        This source has no explicit is_severefield, so is_severeis set to
#        missing if not derivable from a keyword in the description text.
#
# Usage: Run in the arrow directory, which contains the input CSV files, as:
#   julia build_pddi_arrow.jl
#
# Packages needed for installation: Arrow, CSV, DataFrames, RxNav, StatsBase
# Once processed, RxNav needs only Arrow.jl and DataFrames.jl for access.

using Arrow, CSV, DataFrames, RxNav, Serialization, StatsBase

# File pathnames for input CSVs and output Arrow file
const FRENCHDDI_EN_CSV = "frenchDDI_en.csv"
const ONC_HIGH_CSV = "ONC_High_Priority_Mapped.csv"
const DRUGBANK_CSV = "db_drug_interactions.csv"
const OUTPUT_ARROW = "../pddi_interactions.arrow"

# Translation cache
const FRENCH_TO_ENGLISH_DRUG = Ref(Dict(
    "paracétamol" => "acetaminophen",
    "warfarine" => "warfarin",
    "ibuprofène" => "ibuprofen",
    "amoxicilline" => "amoxicillin",
    "penicilline" => "penicillin",
))

const serialized_cache_name = "french_to_english_drug_cache.jls"

"""
    french2englishdrugname(frenchtext::String)::String

Check a French drug name and change to English via cache lookup or translation.
"""
function french2englishdrugname(frenchtext)
    if haskey(FRENCH_TO_ENGLISH_DRUG[], frenchtext)
        return FRENCH_TO_ENGLISH_DRUG[][frenchtext]
    else
        eng = RxNav.getSpellingSuggestions(String(frenchtext))
        if isnothing(eng) || isempty(eng)
            FRENCH_TO_ENGLISH_DRUG[][frenchtext] = frenchtext
            return frenchtext
        end
        FRENCH_TO_ENGLISH_DRUG[][frenchtext] = eng[1]
        return eng[1]
    end
end

"""
    load_frenchdb(path)::DataFrame

Loads the FrenchDB PDDI file. Confirmed header (tab-delimited, quoted with "):
mol2, mol1, prota2, prota1, description_interaction, mecanisme, niveau, DB1, DB2

IMPORTANT PREPROCESSING NOTE: 

The FrenchDB CSV was originally in French and was preprocessed to use English translations
of the French text before use here. This means that the FrenchDDI_en.csv file should already 
contain English translations of the French text before the data is loaded with this function. 

The translation was done via Google Drive / Sheets by uploading the FrenchDDI.csv file
and adding new columns called "mol1_en", "mol2_en", "description_interaction_en" and 
"mecanisme_en". These columns were created to store the English translations. 
Google Translate was used to create these translations of the mol1, mol2, 
description_interaction, and mecanisme columns by placing 
=GOOGLETRANSLATE(FRENCH_COLUMN_LETTER_NUMBER, "fr", "en")
into the corresponding new columns, where FRENCH_COLUMN_LETTER_NUMBER was per-cell letter 
and row number of the original French text cell (e.g., A2 for the mol2 column). 
We then downloaded the translated CSV file as "frenchDDI_en.csv" for use with this function.

There were many French medication terms, such as 'warfarine', which failed to be 
translated correctly by Google Translate. These are often obviously translatable to English
to someone with a good pharmaceutical vocabulary even when they were not translated by 
Google's vocabulary. Because we need consistent spelling for automated searches, French
spellings have been when possible changed to English ones using the function 
`french2englishdrugname`, which handles such cases by checking its cache and falling back
to RxNav's speller to obtain a translation if RxNav's speller has such a suggestion.

SO, AN IMPORTANT WARNING: translation errors often occur and are possible here as well.
"""
function load_frenchdb(path::AbstractString)
    # deserialize the cache if it exists
    if isfile(serialized_cache_name)
        FRENCH_TO_ENGLISH_DRUG[] = deserialize(serialized_cache_name)
    end

    raw = CSV.read(path, DataFrame; quotechar = '"', missingstring = "")

    # Confirmed is_severecodes in the :niveau column: CI, AD, PE, PC
    serious_codes = Set(["CI", "AD"]) # contraindicated / avoid if possible
    minor_codes = Set(["PE", "PC"])   # precaution / take into account

    is_severe = map(raw.niveau) do code
        c = strip(uppercase(coalesce(code, "")))
        if c in serious_codes
            true
        elseif c in minor_codes
            false
        else
            @warn "FrenchDB: unrecognized niveau code, defaulting to false" code
            false
        end
    end

    desc = map(raw.description_interaction_en, raw.mecanisme_en) do d, m
        d = strip(coalesce(d, ""))
        m = strip(coalesce(m, ""))
        isempty(d) ? m : d
    end

    df = DataFrame(
        Name_1 = titlecase.(french2englishdrugname.(strip.(coalesce.(raw.mol1_en, "")))),
        Name_2 = titlecase.(french2englishdrugname.(strip.(coalesce.(raw.mol2_en, "")))),
        is_severe = is_severe,
        description = String.(desc),
        source = fill("FrenchDDI, translated, originally via Github", nrow(raw)),
    )
    # ensure the cache is saved
    serialize(serialized_cache_name, FRENCH_TO_ENGLISH_DRUG[])
    return df
end

"""
    load_onc_highdb(path)::DataFrame

Loads the ONC High-Priority mapped file. Confirmed format: '\$'-delimited,
NO header, 5 fields per row (Name1, DBID1, Name2, DBID2, "") because each
line ends with a trailing '\$'. No description or is_severecolumns exist in
the source -- every row is a high-priority (serious) interaction by
definition, and we synthesize a short description.
"""
function load_onc_highdb(path::AbstractString)
    raw = CSV.read(
        path,
        DataFrame;
        delim = '$',
        header = [:Name1, :DBID1, :Name2, :DBID2, :Blank],
        missingstring = "",
    )

    n1 = titlecase.(strip.(coalesce.(raw.Name1, "")))
    n2 = titlecase.(strip.(coalesce.(raw.Name2, "")))

    DataFrame(
        Name_1 = n1,
        Name_2 = n2,
        is_severe = fill(true, nrow(raw)),
        description = [
            "ONC High-Priority potential drug-drug interaction between $(a) and $(b)."
            for (a, b) in zip(n1, n2)
        ],
        source = fill("ONC-High-Priority via Github", nrow(raw)),
    )
end

# Used to try to identify serious interactions in DrugBank csv from the description column
const SERIOUS_KEYWORDS = [
    "contraindicat",
    "serious",
    "severe",
    "major",
    "life-threatening",
    "life threatening",
    "fatal",
    "toxicity",
    "increased risk of death",
]

"""
    load_drugbank_sentences(path)::DataFrame

Loads the DrugBank-derived sentence-description file (commonly distributed
as `db_drug_interactions.csv`). This file was obtained from 
https://www.kaggle.com/datasets/mghobashy/drug-drug-interactions/data
on 25 September 2026. It may need to be updated if the source file changes.
"""
function load_drugbank_sentences(path::AbstractString)
    raw = CSV.read(path, DataFrame; missingstring = "")
    c1, c2, cd = Symbol("Drug 1"), Symbol("Drug 2"), Symbol("Interaction Description")

    desc = String.(strip.(coalesce.(raw[!, cd], "")))
    is_severe = map(desc) do d
        dl = lowercase(d)
        any(kw -> occursin(kw, dl), SERIOUS_KEYWORDS)
    end

    DataFrame(
        Name_1 = titlecase.(strip.(coalesce.(raw[!, c1], ""))),
        Name_2 = titlecase.(strip.(coalesce.(raw[!, c2], ""))),
        is_severe = is_severe,
        description = desc,
        source = fill("DrugBank via Kaggle", nrow(raw)),
    )
end

""" Build combined table from all csv, then output combined data to Arrow file """
function build_combined_arrow_table()::DataFrame
    tables = DataFrame[]

    if isfile(FRENCHDDI_EN_CSV)
        push!(tables, load_frenchdb(FRENCHDDI_EN_CSV))
    else
        @warn "Skipping FrenchDB, file not found" FRENCHDDI_EN_CSV
    end

    if isfile(ONC_HIGH_CSV)
        push!(tables, load_onc_highdb(ONC_HIGH_CSV))
    else
        @warn "Skipping ONC High-Priority, file not found" ONC_HIGH_CSV
    end

    if isfile(DRUGBANK_CSV)
        push!(tables, load_drugbank_sentences(DRUGBANK_CSV))
    else
        @warn "Skipping db_drug_interactions.csv, file not found" DRUGBANK_CSV
    end

    isempty(tables) &&
        error("No input files found -- check the paths at the top of the script.")

    combined = reduce(vcat, tables)

    # Normalized (titlecaseuppercased, trimmed) helper columns used only for fast,
    # case-insensitive lookups. Keep the originals for display.
    combined.Name_1 = titlecase.(strip.(combined.Name_1))
    combined.Name_2 = titlecase.(strip.(combined.Name_2))

    Arrow.write(OUTPUT_ARROW, combined)
    println("Wrote $(nrow(combined)) rows to $OUTPUT_ARROW")
    println(combined.source |> countmap)

    return combined
end

"""
    interactions_for(df, drug; severeonly=false)::DataFrame

Returns all interactions involving `drug` (case-insensitive, matches either
Name_1 or Name_2). Set `severeonly=true` to get only the serious /
high (is_severe == true) rows.
"""
function interactions_for(df::DataFrame, drug::AbstractString; severeonly::Bool = false)
    target = titlecase(strip(drug))
    mask = (df.Name_1 .== target) .| (df.Name_2 .== target)
    if severeonly
        mask = mask .& df.is_severe
    end
    return df[mask, [:Name_1, :Name_2, :is_severe, :description, :source]]
end

"""
    high_severity(df)::DataFrame

Returns every row flagged as serious/high-is_severe regardless of drug.
"""
high_severity(df::DataFrame) = df[
    df.is_severe, [:Name_1, :Name_2, :is_severe, :description, :source]
]



# Build the combined Arrow table and demonstrate example queries
df = build_combined_arrow_table()

# example usage -- replace "Warfarin" with a drug you know is in your data
example_drug = "Warfarin"
println("\nAll interactions for $example_drug:")
println(interactions_for(df, example_drug))

println("\nHigh (is_severe) only interactions for $example_drug:")
println(interactions_for(df, example_drug; severeonly = true))

println("\nTotal high (is_severe) rows in whole table: ", nrow(high_severity(df)))
