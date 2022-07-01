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
pol.columns
pol= dataset[dataset.Genes.isin(p)]
filg= ["correct orf",'without stop']
gagF =gag[gag.pattern_pos_seq.isin(filg)]
pol.columns
filP = ['without_met_Pattern_before_stp',"without start-stop"]
polF =pol[pol.pattern_pos_seq.isin(filP)]

len(gagF)
df3=pd.merge(gagF,polF, on='sequences', how='left', indicator= True)
df3.to_csv("patterns_matched1jul.csv")