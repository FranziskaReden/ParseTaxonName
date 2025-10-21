# NCBI Taxon ID Finder

This repository provides a set of Python scripts to **resolve organism names to NCBI Taxonomy IDs** and retrieve their **taxonomic lineages**.  
It was developed to support downstream analyses in **phylogenomics**, **biodiversity informatics**, and **taxonomic data integration**, where sample or genome labels often need to be linked to standardized NCBI identifiers.

The tool is designed to work with **incomplete, inconsistent, or ambiguous names**, using a biologically informed name-matching strategy that minimizes false positives in the dense NCBI taxonomy database.

---

## Overview

The workflow accepts a list of taxon names and searches for corresponding entries in the NCBI taxonomy.  
It performs the search in multiple stages:

1. **Strict search** – finds exact matches using a pre-indexed dictionary of NCBI taxon names.  
2. **Relaxed search** – identifies near matches using substring filtering and fuzzy similarity scoring.  
3. **Lenient search** – iteratively reduces complex names (e.g., those including strain identifiers or qualifiers) to identify possible matches at higher ranks.  

The search strategy ensures that at least the **genus component** of a name matches exactly to NCBI taxonomy entries, preventing false positives caused by overlapping or homonymous names.

---

## Repository Structure

| File | Description |
|------|-------------|
| **find_taxon_ids.py** | Main driver script. Handles input, runs multi-stage name matching, manages multiprocessing, and writes results and checkpoints. |
| **get_lineage.py** | Retrieves full NCBI lineages for matched Taxonomy IDs and appends them to the results file. |
| **ncbi_tax.py** | Loads and preprocesses NCBI taxonomy data (`names.dmp`, `nodes.dmp`), builds internal indices by starting letter, and flags duplicate taxon names. |
| **utils.py** | Contains helper functions for name cleanup, I/O handling, checkpointing, and text normalization. |

---

## Matching Strategy

Initial tests using regular expressions were found to produce a large number of false positives due to overlapping and similar taxon names in the NCBI taxonomy.  
For this reason, the final implementation **avoids regex** and instead uses a controlled combination of:

- **Substring filtering** via `pandas.Series.str.contains()` (non-regex mode, case-insensitive).  
- **Fuzzy similarity scoring** using `rapidfuzz.ratio()` to quantify the match between candidate names and query variants.  
- **Biologically informed filtering** based on genus-level exact matching.  

This approach significantly reduces ambiguity while maintaining good runtime performance, aided by pre-indexing taxa by their initial letter and applying multiprocessing for large datasets.

---

## Dependencies

- Python ≥ 3.8  
- `pandas`  
- `rapidfuzz`  
- `tqdm`  

Install the dependencies using conda and the yaml file provided in this repository:

```bash
conda env create -f environment.yaml
```

## Usage

### Command-line example

```
python parse_taxon_name.py \
  -n 'Homo sp' 'Mus musculuss' 'Acanthotrema' 'Escherichia imaginarus' 'Unicorn' \
  --mode lenient \
  -l minimal \
  --db /PATH/TO/DB/
```

### Arguments

| Argument            | Description                                                                                                                                                                                                                                                                               |
| ------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| -n / --name         | One or more taxon names to resolve. Provide multiple names separated by spaces and enclose names containing spaces in quotes. Example: -n 'Homo sp' 'Mus musculuss'                                                                                                                       |
| --mode              | Matching strictness. Options:<br>• strict — exact text matches only (fastest, lowest recall)<br>• relaxed — includes substring and partial matches<br>• lenient — uses fuzzy matching and name reduction (genus/core name) to increase recall; recommended for noisy or incomplete inputs |
| -l / --lineage_mode | Lineage output format. Options:<br>• minimal — outputs only the main taxonomic ranks (compact format)<br>•reduced — outputs all ranks that are unique (omits clades with rank 'clade') <br>• full — includes all intermediate taxonomic levels                                                                                                                              |
| --redo              | Forces reprocessing of all names, overwriting existing checkpoints and cached results                                                                                                                                                                                                     |
| --cores             | Number of CPU cores to use for multiprocessing. Accepts an integer. Default is 1 core.                                                                                                                                                                         |
| --prefix            | Prefix for output files. All results will be written using this prefix (see Output Files section)                                                                                                                                                                                         |
| --score             | Minimum fuzzy similarity threshold (numeric). Candidates with a score below this value are ignored                                                                                                                                                                                        |
| --quiet             | Suppress or reduce console progress output                                                                                                                                                                                                                                                |
| --update-db | Updates the local NCBI taxonomy database to the latest version (requires internet access) |

**Note:** The first time the script is run, it will automatically download the NCBI taxonomy database if it is not present. Internet access is required for the initial download or for updates using the --update-db flag. The downloaded data is then processed (filtering of homonyms and indexing), which will take some time. However, this step is only required once or every time you wish to update the database.



### Example output

The output will look something like this: 
```
Reading in /PATH/TO/DB/taxa_names_sorted.tsv file...
Loaded 5 names. Of those 5 are unique.

Starting exact match search for 5 names...
name    tax_id  name_txt        name_class      strict_score    relaxed_score   reduced_name    no_number_name  min_name        time(s) comment
Homo sp None    None    None    0       0       None    None    None    0.0     None
Mus musculuss   None    None    None    0       0       None    None    None    0.0     None
Acanthotrema    None    None    None    0       0       None    None    None    0.01331 HOMONYM - multiple entries found: 378736, 1415158
Unicorn 1498384 Unicorn nan     100.0   0       None    None    None    0.00021 None
Escherichia imaginarus  None    None    None    0       0       None    None    None    0.0     None
Checkpoint saved: 1 names processed.

Starting lenient search for 4 names using 1 cores...
name    tax_id  name_txt        name_class      strict_score    relaxed_score   reduced_name    no_number_name  min_name        time(s) comment
Homo sp 2813599 Homo sp.        nan     93.333  100.0   Homo sp.        None    None    0.05315 None
Mus musculuss   10090   Mus musculus    nan     96.0    96.0    Mus musculuss   None    Mus     0.08457 None
Acanthotrema    None    None    None    0       0       Acanthotrema    None    None    0.05449 HOMONYM - multiple entries found: 378736, 1415158
Escherichia imaginarus  561     Escherichia     nan     66.667  100.0   Escherichia imaginarus  None    Escherichia     0.06989 None
Checkpoint saved: 4 names processed.

Matched names written to tax_ids.tsv. Failed names to tax_ids_failed.txt.
Reading in /home/frareden/.ncbi_tax/nodes.tsv file...

20/10/2025 16:41:49: Retrieving lineages (minimal)....
1498384 1       root:1;cellular organisms:131567;Eukaryota:2759;Metazoa:33208;Arthropoda:6656;Arachnida:6854;Araneae:6893;Oonopidae:1498386;Unicorn:1498384;
2813599 1       root:1;cellular organisms:131567;Eukaryota:2759;Metazoa:33208;Chordata:7711;Mammalia:40674;Primates:9443;Hominidae:9604;Homo:9605;Homo sp.:2813599;
10090   1       root:1;cellular organisms:131567;Eukaryota:2759;Metazoa:33208;Chordata:7711;Mammalia:40674;Rodentia:9989;Muridae:10066;Mus:10088;Mus musculus:10090;
561     1       root:1;cellular organisms:131567;Bacteria:2;Pseudomonadati:3379134;Pseudomonadota:1224;Gammaproteobacteria:1236;Enterobacterales:91347;Enterobacteriaceae:543;Escherichia:561;

20/10/2025 16:41:49: Results were written into lineage.tsv file.
```

### Interpretation of Results

- Exact (strict) search is attempted first; names with no exact match are passed to relaxed and lenient stages.
- Lenient matching applies a sequence of name reductions (tidying, removing strain numbers, shaving words) and requires the genus/core part to match while using fuzzy scoring to correct spelling and handle abbreviations.
- Homonym cases (multiple entries for the same name) are flagged with a HOMONYM comment and the candidate TaxIDs are reported for manual inspection.

**Example outcomes from the provided input:**

- Mus musculuss → matched to Mus musculus (TaxID 10090) by lenient search (fuzzy correction).
- Homo sp → matched to Homo sp. (TaxID 2813599); lenient mode  replaced sp to sp. (as found in the NCBI taxonomy database) to find an exact match.
- Unicorn → exact match found in the database (TaxID 1498384).
- Acanthotrema → flagged as HOMONYM because multiple NCBI entries match the reduced form; the script lists candidate TaxIDs (requires manual curation).
- Escherichia imaginarus → was reduced to the genus core Escherichia as no species with this name was found (TaxID 561).

**Lineage retrieval:**

After matching TaxIDs, the script reads the local nodes data and writes lineage entries in the chosen format (minimal). The lineage shows the hierarchical path from root to the resolved taxon.

---

### Output Files

**`<prefix>tax_ids.tsv`**

- Tab-delimited table of matched names and statistics. Columns include: original_name, tax_id, name_txt, name_class, strict_score, relaxed_score, reduced_name, no_number_name, min_name, time(s), comment.

**`<prefix>tax_ids_failed.txt`**

- Plain list of input names for which no reliable match was found. Useful for manual curation or rerunning with adjusted parameters.

**`lineage.tsv`**

- Contains tax_id, rank, and the lineage string (semicolon-separated rank:name pairs). Format depends on lineage_mode (minimal vs full).




