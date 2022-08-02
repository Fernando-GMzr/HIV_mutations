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
dataset = pd.read_csv('patterns_1agoust.csv')
dataset.columnsm
g=["gag"]
gag = dataset[dataset.Genes.isin(g)]
gag['Genes']
p=['pol']

pol= dataset[dataset.Genes.isin(p)]
len(pol)
pol.columns
filg= ["correct orf",'without stop']
gagF =gag[gag.pattern_pos_seq.isin(filg)]
len(gagF)
pol.columns
filP = ['without_met_Pattern_before_stp',"without start-stop"]
polF =pol[pol.pattern_pos_seq.isin(filP)]

len(gagF)
len(polF)

'''Fusion de la proteina GAG y la proteina POL previamente filtradas'''
df3=pd.merge(gagF,polF, on='sequences', how='left', indicator=True)
len(df3)
len(df3.columns)
#df3 = df3.iloc[:,[1]].drop_duplicates()
len(df3)
pd.unique(df3['sequences'])
len(pd.unique(df4['sequences']))

df4=df3.drop_duplicates(keep='first')
len(df4)
both=['both']
df4= df4[df4._merge.isin(both)]


'''Creo un dataset conteniendo unicamente los nombres de secuencias que mapearon para gag y pol'''
#df4 = df4[["sequences","Genes_y"]]
df4=df4[['Genes_y','sequences']].drop_duplicates()
df4=df4.rename(columns={"Genes_y":"Genes"})
df5=df4.assign(expression_pol="yes") ####
#dataset["pattern_pos_seq"].unique
#df5[['sequences', 'Genes', 'expression_GagNef']].nunique()
#df5 = df5[['sequences', 'Genes', 'expression_GagNef']].nunique()
#len(df5.drop_duplicates())
#len(dataset)
df5 = df5.drop_duplicates()

#df_pol_Expresion.drop_duplicates()
df_pol_Expresion= pd.merge(dataset,df5, on=['sequences',"Genes"], how='left', indicator= False)
len(df_pol_Expresion)
#df_pol_Expresion.iloc[:,[2,]]
df_pol_Expresion.to_csv("patterns_matched2ago.csv", index=False)
#len(df_pol_Expresion)
#len(df5[df5.expression_GagNef.isin(["yes"])])
