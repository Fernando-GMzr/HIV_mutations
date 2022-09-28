

# HIV mutation classifier

Part of this script was used to obtain the results of HIV protein mutation results from treated and untreated donors published in:
*https://www.jci.org/articles/view/154422*

This scripts was developed to identify and quantify CTL/CD8+ T-cell epitopes and  escape variants in multifasta sequences of HIV-1 proteins (Gag, Pol and Nef). In addition, potentially silent and non-translated CTL epitopes can be also quantified using this script. Databases from Los Alamos were used to identify experimentally validated HIV-1 CTL epitopes and escape variants. In addition, the information contained in the sequence ID (e.g., donor ID, time point, cell type, HIV protein) can be also identified and classified.<br /> Since epitopes and escape variants are restricted to specific HLA haplotypes, the HLA type of each donor is considered. Therefore, the donor HLA type is matched with the corresponding epitopes (or escape variants) restricted to the donor HLA type and filtered. Next, the filtered epitopes or escape variants (from now patterns) are identified in all the sequences from the corresponding HIV protein and donor. 
After the patterns are identified, they are classified based on their potential to be translated. For this analysis the presence or absence and the position of start and stop codons within the open reading frame (ORF) of the protein sequence analysed is considered. Thus, each identified pattern can be classified as follows:<br /> 

    **1)** Translated (translated): the identified pattern is contained in an ORF sequence that may be translated - i.e. start codon at the beginning of the ORF and a stop codon after the pattern -.
    **2)** Translated but without stop codon (transl_no_stop_codon): The ORF contains the corresponding start codon but does not contain any stop codon. 
    **3)** Silenced by premature stop codon (silent_prem_stop_codon): The pattern may not be translated due to the presence of a premature stop codon. 
    **4)** Silenced by the absence of a start codon (silent_no_start_codon): The pattern may not be translated due to the absence of a start codon in the ORF.  
    **5)** Silenced by both the absence of a start codon and the presence of a premature stop codon (silent_no_start_codon_prem_stop_codon): The pattern may not be translated due to the absence of a start codon and the presence of a premature stop codon.  
    **6)** Silenced by the absence of a start codon (silent_no_start_codon_no_stop_codon): The pattern may not be translated due to the absence of a start codon and the ORF does not contain any stop codon.  

<br /> Then, all the information obtained is classified in a table or dataframe showing the pattern sequence, the category of each identified pattern as well as the sequence header,  HLA type, the sequence position of the pattern. The script also quantified as “no patterns” if the corresponding donor HLA type does not match with any HLA restricted epitope or no pattern is identified in the analysed sequence,

<br />  As the Pol ORF does not contain start codon, since its translation occurs through a frameshifting during the translation of the polyprotein protein Gag-Pol, patterns identified in Pol sequences were classified as translated if the corresponding Gag sequence contains a start codon and no premature stop codon were detected before the start of the Pol ORF. 

<br /> This pipeline was used to analyse sequences from naïve and memory CD4+ T-cells from individual living with HIV under antiretroviral therapy. The results obtained were published in Duette et al, 2022, Journal of Clinical Investigation  (https://www.jci.org/articles/view/154422)




## script [flowchart](https://github.com/Fernando-GMzr/HIV_mutations/blob/master/fluxogram.png)

![flowchart](https://github.com/Fernando-GMzr/HIV_mutations/blob/master/fluxogram.png)


