# Business Logic API - README

## Overview
This script provides a set of functions to interact with a genomic database, enabling users to retrieve, search, and analyze genetic data by various parameters. The primary goal of the API is to interface with the database to support bioinformatics applications, offering access to gene, protein, codon, and restriction enzyme information. The script consists of two major categories of functions:

1. **Search Functions**: These functions bridge the front-end with the database layer and return results directly from the database without modification.
2. **Calculation Functions**: These functions perform computations over the genetic data, either across the entire database or for specific genes, and return the results after processing.

The API interacts with an external `dbapi` module and optionally uses configuration settings from a `config` module.

## Requirements
- Python 3.x
- Modules required: `dbapi`, `config`
  
Ensure that the paths to `dbapi` and `config` are correctly set in the script, especially when running from different directories.

## Functions

### 1. **Search Functions** (Bridge between Front-End and Database)
These functions simply retrieve data from the database layer and return it to the front end.

#### `getAllEntries()`
Fetches all available entries in the database.
- **Returns:** List of all gene entries from the database.

#### `searchGeneID(geneID)`
Searches the database using a specific gene ID.
- **Input:** `geneID` (str) – the gene ID to search for.
- **Returns:** List of matching gene details, or an empty list if no matches.

#### `searchProteinProduct(proteinProductString)`
Searches the database by protein product name.
- **Input:** `proteinProductString` (str) – the protein product name.
- **Returns:** List of matching entries or an empty list if no matches.

#### `searchAccession(accessionCode)`
Searches the database using a gene accession code.
- **Input:** `accessionCode` (str) – the accession code.
- **Returns:** List of matching entries or an empty list if no matches.

#### `searchChromLocation(chromosomeLocation)`
Queries the database by chromosome location.
- **Input:** `chromosomeLocation` (str) – chromosome location.
- **Returns:** List of matching gene entries or an empty list if no matches.

### 2. **Calculation Functions** (Perform Calculations on the Data)
These functions process genetic data either across the whole database or for a specific gene and return the results after calculation.

#### a) **Functions that Perform Calculations on the Entire Database**

##### `getGeneralCodonFreq()`
Calculates codon frequency data for various amino acids across chromosome 9 and returns the results.
- **Returns:** A dictionary where the keys are amino acid codes (e.g., `Ala(A)`), and the values are lists of:
  - Codons used by the amino acid.
  - Corresponding percentage usage for each codon.
  - Corresponding frequencies for each codon.
  
  Example return:
  ```python
  {'Ala(A)': [('gct', 'gcc', 'gca', 'gcg'), [0.5, 1.35, 2.82, 1.75], [0.08, 0.21, 0.44, 0.27]]}
  ```

#### b) **Functions that Perform Calculations on a Single Gene**

These functions take a gene's accession code and compute results based on the gene-specific data retrieved from the database.

##### `getDNAseqAndCodingRegions(accession)`
Generates the entire DNA sequence of a gene and identifies its exons (coding regions).
- **Input:** `accession` (str) – the gene accession code.
- **Returns:** A list containing the gene sequence, with coding exons as lists of strings and non-coding sequences as plain strings.
  
  Example return:
  ```python
  ['cggttaagc', ['atgacggggctggc'], 'ccggttacgta']
  ```

##### `alignNucleotide_to_aminoacid(accession)`
Aligns each amino acid in the protein sequence with its corresponding codons in the gene's coding sequence.
- **Input:** `accession` (str) – the gene accession code.
- **Returns:** A list of strings, where each string contains the amino acid and its corresponding codon, separated by a colon.
  
  Example return:
  ```python
  ['M:atg', 'G:ggt', 'F:ttt', 'L:tta']
  ```

##### `RE_sites_list_cat(accession)`
Identifies restriction enzyme sites for a gene and categorizes their availability for cloning with five predefined enzymes: EcoRI, BamHI, BsuMI, HindIII, and StuI.
- **Input:** `accession` (str) – the gene accession code.
- **Returns:** Two dictionaries:
  1. **Enzyme Availability Dictionary**: Categorizes enzymes based on their applicability for gene cloning (e.g., enzymes cutting at the 5' end, 3' end, or not applicable).
  2. **Site Indices Dictionary**: Lists the positions of each restriction enzyme's cut sites in the gene sequence.
  
  Example return:
  ```python
  {'available_at_5p': [], 'available_at_3p': [], 'not_applicable': ['EcoRI', 'BamHI', 'BsuMI', 'HindIII', 'StuI']}
  {'EcoRI': [543, 1484, 1432, 6191, 9597, 9870], 'BamHI': [8128, 1073, 3004], 'BsuMI': [6043], 'HindIII': [533, 734, 8953, 11430, 12179], 'StuI': []}
  ```

##### `input_seq_RE_sites(accession, your_choice)`
Works similarly to `RE_sites_list_cat()` but allows the user to provide a custom restriction enzyme sequence.
- **Input:** 
  - `accession` (str) – the gene accession code.
  - `your_choice` (str) – the custom restriction enzyme sequence.
- **Returns:** Results similar to `RE_sites_list_cat()`, but for the custom restriction enzyme.

##### `geneCodonUsage(accession)`
Calculates codon usage statistics, including the percentage usage, frequency, and ratio of codon usage for the gene compared to chromosome 9.
- **Input:** `accession` (str) – the gene accession code.
- **Returns:** A dictionary with amino acid codes as keys (e.g., `Ala(A)`), and each value consists of:
  1. List of codons used by the amino acid.
  2. List of corresponding percentage usage for each codon.
  3. List of corresponding frequencies for each codon.
  - Additionally, a list showing whether the codon usage differs significantly (more than a twofold difference).
  
  Example return:
  ```python
  {'Ala(A)': [('gct', 'gcc', 'gca', 'gcg'), [0.97, 0.24, 0.97, 0.0], [0.44, 0.11, 0.44, 0.0]]}
  ```

## Usage

To use any function, ensure that the `dbapi` and `config` modules are available and correctly configured. These functions can be integrated into larger bioinformatics applications or used as backend services for web tools.

Example usage:
```python
from business_logic_api import *

# Fetch all database entries
entries = getAllEntries()

# Search by gene ID
gene_results = searchGeneID("AB065492")

# Retrieve codon frequency data across chromosome 9
codon_freq = getGeneralCodonFreq()

# Identify restriction enzyme sites for a gene
re_sites = RE_sites_list_cat("AB065492")
```

## Notes
- This API is designed to support CGI and front-end applications and can be adapted as needed.
- Ensure valid access to database functions and handle exceptions where required.
