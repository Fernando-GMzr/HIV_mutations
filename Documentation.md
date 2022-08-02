
## Documentation

## script [flowchart](https://github.com/Fernando-GMzr/HIV_mutations/blob/master/fluxogram.png)

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

*This dataset comes from merging the information of the HLA donors with the dataset of epitope patterns of the "los alamos" database.*

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
* Add options of detection in mutations related do HLA.
* Adding more classifiers
* Add graphical output of calculated statistics
