#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 17:43:18 2022

@author: fernando
"""
import pandas as pd
import sys
import os   

current= os.getcwd() 
os.chdir(current)
#os.chdir('/home/fernando/Documentos/Epitopes_pipi/redesign_inmuno_hiv')
#os.listdir()


#dataset = pd.read_csv('3_octubre_detect.csv') 
#output= pd.read_csv(csv_polexpression) 
#dataset = pd.read_csv('patterns_3sptember.csv')
#csv_polexpresion=patterns_match_3oct.csv'
#dataset.columns
#len(dataset)

''' The following function filters the sequences of the gag and pol genes, according to the conditions that allow the expression of pol. '''
def pol_express(dataset, output):
    g=["gag"]
    #Genes corresponding to gag are filtered from the dataset
    gag = dataset[dataset.Genes.isin(g)]
    gag['Genes']
    p=['pol']
    #Genes corresponding to pol are filtered from the dataset
    pol= dataset[dataset.Genes.isin(p)]
    len(pol)
    pol.columns
    
    # the sequences that meet the conditions "translated",'transl_no_stop_codon' are filtered out.
    
    filg= ["translated",'transl_no_stop_codon']
    
    gagF =gag[gag.pattern_pos_seq.isin(filg)]
    len(gagF)
    pol.columns
    # the sequences that meet the conditions 'silent_no_start_codon', "silent_no_start_codon_no_stop_codon" are filtered out, 
    
    filP = ['silent_no_start_codon', "silent_no_start_codon_no_stop_codon"]
    polF =pol[pol.pattern_pos_seq.isin(filP)]
    
    len(gagF)
    len(polF)
    
    #Fusion of previously filtered GAG protein and POL protein
    #A dataset is created containing only the names of sequences that mapped to gag and pol filtered:
    df3=pd.merge(gagF,polF, on='sequences', how='left', indicator=True)
    len(df3)
    len(df3.columns)
    #df3 = df3.iloc[:,[1]].drop_duplicates()
    len(df3)
    pd.unique(df3['sequences'])
    
    
    #Duplicate rows are deleted:
    
    
    df4=df3.drop_duplicates(keep='first')
    len(df4)
    both=['both']
    df4= df4[df4._merge.isin(both)]
    len(pd.unique(df4['sequences']))
    
    
    
    #df4 = df4[["sequences","Genes_y"]]
    df4=df4[['Genes_y','sequences']].drop_duplicates()
    len(df4)
    df4=df4.rename(columns={"Genes_y":"Genes"})
    df5=df4.assign(expression_pol="yes") ####
    len(df5)
    
    df_pol_Expresion= pd.merge(dataset,df5, on=['sequences',"Genes"], how='left', indicator= False)
    len(df_pol_Expresion)
    #df_pol_Expresion.iloc[:,[2,]]
    df_pol_Expresion.to_csv(output, index=False)
    return print(f' The process has ended , the file {output} has been saved')

#pol_express(dataset, "patterns_match_3oct.csv")
#len(df_pol_Expresion)
#len(df5[df5.expression_GagNef.isin(["yes"])])
#ex = ['yes']
#len(df_pol_Expresion[df_pol_Expresion.expression_pol.isin(ex)])

