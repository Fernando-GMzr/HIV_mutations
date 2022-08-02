

# HIV_mutations_scape


# HIV mutation classifier

Part of this script was used to obtain the results of HIV protein mutation results from treated and untreated donors published in:
*https://www.jci.org/articles/view/154422*

The script search motifs related to epitopes described for hiv proteins (accordly the information of "Los alamos Data Base")
in a aminoacidic multifasta of HIV sequences. HIV-1 sequences are classified according to their genetic patterns within each T cell subsets. For this, in each sequence is identified the code of donnor, and used to filter the dataset of Donnor HLA information, consequently, each HLA pattern is searched in the sequence.<br /> Once time that motif is found, a function evaluate patterns related to gene expression (relation of motif position with start and stop codons positions). Thus, each secuence if classified into six possibilites. The categories are: **1)** when the pattern is detected in sequences with correct orf ("correct_orf"), **2)** the pattern is detected after the reading frame (after the stop codon, *"after_orf"*), **3)** the pattern is found in a sequence without a stop codon ("without_stop"), **4)** sequences without methionine (start codon), with the pattern before the stop codon (*"without_met_Pattern_before_stp"*), **5)** sequences without methionine with the pattern after the stop codon (*"without_met_Pattern_after_stp"*) and **6)** patterns detected in sequences without start and stop codons (* "without_start_stop"*).<br /> Afterly, this will write in a output dataframe whose information show the type of motif epitope identified, the HLA, header of secuence, position detected and the category of expression of each secuences. The patterns motifs related to HLA donnor not found, are counting and writted in the dataframe for future statistic calculations.

## script [flowchart](https://github.com/Fernando-GMzr/HIV_mutations/blob/master/fluxogram.png)

![flowchart](https://github.com/Fernando-GMzr/HIV_mutations/blob/master/fluxogram.png)


