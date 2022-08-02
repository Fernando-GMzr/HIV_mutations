

# HIV_mutations_scape


# HIV mutation classifier

Part of this script was used to obtain the results of HIV protein mutation results from treated and untreated donors published in:
*https://www.jci.org/articles/view/154422*

The script search motifs related to epitopes described for hiv proteins (accordly the information of " DB")
in a aminoacidic multifasta of HIV sequences. HIV-1 sequences are classified according to their genetic patterns within each T cell subsets. For this, in each sequence is identified the code of donnor, and used to filter the dataset of Donnor HLA information, consequently, each HLA pattern is searched in the sequence.\Once time that motif is found, a function evaluate patterns related to gene expression (relation of motif position with start and stop codons positions). Thus, each secuence if classified into six possibilites. The categories are: 1) when the pattern is detected in sequences with correct orf ("correct_orf"), 2) the pattern is detected after the reading frame (after the stop codon, "after_orf"), 3) the pattern is found in a sequence without a stop codon ("without_stop"), 4) sequences without methionine (start codon), with the pattern before the stop codon ("without_met_Pattern_before_stp"), 5) sequences without methionine with the pattern after the stop codon ("without_met_Pattern_after_stp") and 6) patterns detected in sequences without start and stop codons.\Afterly, this will write in a output dataframe whose information show the type of motif epitope identified, the HLA, header of secuence, position detected and the category of expression of each secuences. The patterns motifs related to HLA donnor not found, are counting and writted in the dataframe for future statistic calculations.


## Installation and getting started

To run the script it is necessary to have the modules installed:

**Biopython modules:**

- SeqIO
- Seq
- SeqUtils

**CSV modules:**

- csv
- pandas

**others:**

- itertools
- re

## Structure of data input

*The input format maintains the following structure:*

| Patient_ID | allele     | HLA  | Variant_Epitope | Protein | Epitope_WT |
|------------|------------|------|-----------------|---------|------------|
| 1756       | A.allele.1 | A*03 | KIRLRPGGQ       | Gag     | KIRLRPGGK  |
| 1756       | A.allele.1 | A*03 | KIRLRPGGQ       | Gag     | RLRPGGKKK  |
| 1756       | A.allele.1 | A*03 | KIRLRPGGQ       | Gag     | RLRPGGKKKY |

### Columns description


**Columns:** data related to Donnor, HLA and Aminoacidic motivs related to wild and variant epitopes of protein.

**Patient_ID:** Patient code used to filter the sequences to be identified.

**HLA:** HLA  of donnor

**Variant_epitope:** motifs that belong to mutated epitopes

**Epitope_WT:** motifs that belong to "Wildtype" epitopes

**Protein:** Proteins related to motivs described.

*The metadata information of proteins, motifs epitope and hla was obtain from "HIV database", where the original database was transformed to filter and match with the inmunogenic charachteristic of each Donnor*.



### Running

In a directory containing fasta files and database in csv, run:

```python
HIV_search_patterns.py "input_HLA_donnor.csv" "output_mutation_found.csv"
```

## Structure of data output

| Epitope_type | pattern_recover | hla  | position on the sequence | Genes | pattern    | sequences    | position | origin | Donor id |
|--------------|-----------------|------|--------------------------|-------|------------|--------------|----------|--------|----------|
| wild_variant | yes             | B*35 | silent                   | pol   | TVLDVGDAY  | 1292-17-N-23 | 262.0    | N      | 1292     |
| wild_variant | yes             | B*35 | silent                   | pol   | VPLDEDFRKY | 1292-17-N-23 | 273.0    | N      | 1292     |


## Features


* Add other possibilities of input formats.
* Configure HLA-related conditions to process mutation detection.
* Adding specific classifiers according to mutations and proteins
* Add graphical output of calculated statistics
