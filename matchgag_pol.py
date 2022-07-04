#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 17:43:18 2022

@author: fernando
"""
import pandas as pd
import sys
import os   

os.chdir('/home/fernando/Documentos/Epitopes_pipi/redesign_inmuno_hiv')
dataset = pd.read_csv('patterns_1jul.csv')

g=["gag"]
gag = dataset[dataset.Genes.isin(g)]
gag['Genes']
p=['pol']

pol= dataset[dataset.Genes.isin(p)]
pol.columns
filg= ["correct orf",'without stop']
gagF =gag[gag.pattern_pos_seq.isin(filg)]
pol.columns
filP = ['without_met_Pattern_before_stp',"without start-stop"]
polF =pol[pol.pattern_pos_seq.isin(filP)]

len(gagF)
df3=pd.merge(gagF,polF, on='sequences', how='left', indicator= True)
df3.columns
pd.unique(df3['sequences',''])
df4=df3.drop_duplicates(keep='first')
both=['both']
df4= df4[df4._merge.isin(both)]
df4 = df4[["sequences","Genes_y"]]
df4=df4.rename(columns={"Genes_y":"Genes"})
df5=df4.assign(expression="yes")
dataset.columns
df5.columns
len(df5)
df_pol_Expresion= pd.merge(dataset,df5, on=['sequences',"Genes"], how='left', indicator= False)
df_pol_Expresion= df_pol_Expresion.drop_duplicates()
df_pol_Expresion.to_csv("patterns_matched3jul.csv", index=False)
len(df_pol_Expresion)
