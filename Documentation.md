
## Documentation

## Script [flowchart](https://github.com/Fernando-GMzr/HIV_mutations/blob/master/fluxogram.png)

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

## Preparation data input

* Dataset retrieve and process from "Los alamos"

 The datasets corresponding to:

* Protein epitope variant and related HLA

* Protein epitope Wild type and related HLA

* Dataset HLA information related to donnor


**Were combined  as shown in the figure:**


![Figure 1](https://github.com/Fernando-GMzr/HIV_mutations/blob/master/Figure1.png)

* The dataset of HLA donnor was matched  with information of epitopes variant and wild , generating two dataset.

**Donnor-HLA** match to **HLA-variant epitope** and  **Donnor-HLA** match to **HLA-Wild type epitope:**

**The process was carried out as shown in figure 2:**

![Figure 2](https://github.com/Fernando-GMzr/HIV_mutations/blob/master/Figure2.png)

## Structure of data input

* The **Donnor_HLA_WT_epitope** and **Donnor_HLA_Variant_epitope** was join in one dataset.

**The input format maintains the following structure:**

| Patient_ID | allele     | HLA  | Variant_Epitope | Protein | Epitope_WT |
|------------|------------|------|-----------------|---------|------------|
| 1756       | A.allele.1 | A*03 | KIRLRPGGQ       | Gag     | KIRLRPGGK  |
| 1756       | A.allele.1 | A*03 | KIRLRPGGQ       | Gag     | RLRPGGKKK  |
| 1756       | A.allele.1 | A*03 | KIRLRPGGQ       | Gag     | RLRPGGKKKY |

*This dataset comes from merging the information of the HLA donors with the dataset of epitope patterns of the "los alamos" database.*

*  *Each dataset with information on the donor, hla and wild and variant epitope pattern were combined in a single dataset (see available script join_variant_wild_DS.R") in order to gather the information in a single input. Although the combination generated epitope patterns duplicates, the python scripts filters patterns in order to analyze each pattern corresponding to HLA and donor only once.*

### Columns description


**Columns:** Data of HLA and Aminoacidic motivs related to wild and variant epitopes of protein.


* **Patient_ID:** Patient code used to filter the sequences to be identified.

* **HLA:** HLA  of donnor

* **Variant_epitope:** motifs that belong to mutated epitopes

* **Epitope_WT:** motifs that belong to "Wildtype" epitopes

* **Protein:** Proteins related to motivs described.

*The metadata information of proteins, motifs epitope and hla was obtain from "HIV database" (https://www.hiv.lanl.gov/content/immunology/index.html), where the original database was transformed to filter and match inmunogenic characteristic of each Donnor*.



### Running

In a directory containing fasta files and csv, run in terminal:

``` python 
python HIV_count_patterns.py "input_donnor_Fusion_epitope_info.csv" "output_mutation_found.csv" 
```
*If it is required to evaluate the expressed pol sequences (according to the conditions of the analyzed gag sequences of the same molecule), the arguments must be added in the command line:*

``` python 
python HIV_count_patterns.py "input_donnor_Fusion_epitope_info.csv" "output_patten_found.csv" "pol_express"  "output_gag_pol_expression.csv"
```
* **"pol_express":** adding this argument calls a new function to evaluate the expression of pol depending on the conditions of gag.


* **"output_gag_pol_expression.csv":** this "*.csv" file is a modification of the first output where information about the pol expression is added

## Structure of data output

| Epitope_type | pattern_recover | hla  | position on the sequence | Genes | pattern    | sequences    | position | origin | Donor id |
|--------------|-----------------|------|--------------------------|-------|------------|--------------|----------|--------|----------|
| wild_variant | yes             | B*35 | silent                   | pol   | TVLDVGDAY  | 1292-17-N-23 | 262.0    | N      | 1292     |
| wild_variant | yes             | B*35 | silent                   | pol   | VPLDEDFRKY | 1292-17-N-23 | 273.0    | N      | 1292     |


## Features


* Add other possibilities of input formats.
* Add options of detection to mutations related do HLA.
* Add graphical output of calculated statistics
