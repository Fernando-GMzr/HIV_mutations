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
#s.listdir()


#ataset = pd.read_csv(csv_output) 
#output= pd.read_csv(csv_polexpression) 
#dataset = pd.read_csv('patterns_3sptember.csv')
#csv_polexpresion=patterns_match_3sptember.csv'
#dataset.columns
#len(dataset)
def pol_express(dataset, output):
    g=["gag"]
    #Se filtran del dataset los genes correspondientes a pol
    gag = dataset[dataset.Genes.isin(g)]
    gag['Genes']
    p=['pol']
    #Se filtran del dataset los genes correspondientes a gag' 
    pol= dataset[dataset.Genes.isin(p)]
    len(pol)
    pol.columns
    
    # the sequences that meet the conditions "Pattern_in_orf",'without_stop' are filtered out.
    
    filg= ["translated",'transl_no_stop_codon']
    
    gagF =gag[gag.pattern_pos_seq.isin(filg)]
    len(gagF)
    pol.columns
    # the sequences that meet the conditions 'without_start_Pattern_before_stop', "without_start-stop" are filtered out, 
    
    filP = ['translated','silent_prem_stop_codon']
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
#len(df_pol_Expresion)
#len(df5[df5.expression_GagNef.isin(["yes"])])
#ex = ['yes']
#len(df_pol_Expresion[df_pol_Expresion.expression_pol.isin(ex)])

