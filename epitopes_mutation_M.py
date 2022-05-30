#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 04:23:01 2021

@author: bioinfo
"""
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 12:34:00 2021

@author: eares
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 18:43:16 2021

@author: eares
"""

#import matplotlib.pyplot as plt
#from matplotlib import pyplot as plt
#import seaborn as sns
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SeqUtils
import csv
import sys
import pandas as pd

import itertools
import re

os.getcwd()
os.listdir()
os.chdir('/home/fernando/Documentos/Epitopes_pipi')
os.getcwd()
os.listdir()
os.chdir('/home/bioinformatica/Documentos/Documentos/nuevos')
os.chdir('/home/fernando/Documentos/Epitopes_pipi/redesign_inmuno_hiv')
#Functions

def get_fastas(directorio):
    lista = []
    for root, dirs, files in os.walk(directorio): #camino por el directorio
        for file in files:
            if 'fas' in file: 
                lista.append(file)
                
    return lista

"""Detects aminoacidic patterns related to wild-type or escape mutations in amino acid sequences.
It takes with input a multifasta and a database containing the immunological information of each patient
and the sub sequences related to the wild-type or escape mutations."""

def hunt_mutations(db_mut, number_index,sequence):
    epitope_wild = db_mut["Epitope_WT"].iloc[number_index]
    epitope_scape = db_mut["Variant_Epitope"].iloc[number_index]
    hla_scape = db_mut["HLA"].iloc[number_index]
    position_wild = sequence.seq.find(str( epitope_wild))
    position_scape = sequence.seq.find(str( epitope_scape))
    position_asteris = sequence.seq.find(str('*'))
    position_metionina = sequence.seq.find(str('M'))
    
        
    return epitope_wild,epitope_scape,hla_scape, position_wild,position_scape,position_asteris,position_metionina

mutations.head()    

#hla = db_mut["HLA"].unique()
    #hla = hla.tolist()
    #HLA_df= db_mut[db_mut.HLA.isin(hla_scape)]
"""This function filters the detected amino acid immunogenic patterns that correspond to the HLA of the analyzed patient. """
def filters_mut(db_mut,Patient):#hla_scape,HLA_patient,position_scape,position_wild):
    HLA_patient_df = db_mut[db_mut.Patient_ID.astype(str)== str(Patient)]
    
    
    return HLA_patient_df

        

"""This function returns a dictionary containing a key that classifies the sequence according 
to the position of the amino acid pattern (type variant scape) with respect to the reading frame, and its value corresponds to a boolean.
"""
def write_scape(position_scape, hla_scape, position_asteris, position_metionina, HLA_patient_df, epitope_scape, count_pattern_sc):
    if hla_scape in HLA_patient_df["HLA"].tolist():
        if (position_scape != -1) and (len(HLA_patient_df) != 0):    
            if ( epitope_scape not in count_pattern_sc.keys()) or ( epitope_scape in count_pattern_sc.keys() and  hla_scape not in count_pattern_sc[epitope_scape]):
            
                #print(f'entro en scape {(position_scape != -1) and ( epitope_scape not in count_pattern_sc.keys()) or (epitope_scape in count_pattern_sc.keys() and (hla_scape in count_pattern_sc[epitope_scape] == False)) and (len(HLA_patient_df) != 0)}')
                if position_metionina== 0:
                    if position_asteris!= -1:
                        if position_scape < position_asteris:
                            
                            return dict(correct_orf = True)
                        else:
                            
                            return dict(after_orf= True)
                    else:
                            
                            return dict(without_stop=True)
                else: 
                    if position_asteris != -1:
                        if position_scape < position_asteris:
                            
                            return dict(silent = True)
                        else:
                            
                            return dict(without_met_Pattern_after_stp= True)
                    else:
                        
                        return dict(without_star_stop = True)
            else:
                return {}        
        else:
            return {} 
    else:
        return {}

"""This function returns a dictionary containing a key that classifies the sequence according 
to the position of the amino acid pattern (type variant wild) with respect to the reading frame, and its value corresponds to a boolean.
"""        
def write_wild(position_wild,hla_scape, position_asteris, position_metionina,   HLA_patient_df, epitope_wild, count_pattern_wd ):
    if hla_scape in HLA_patient_df["HLA"].tolist():
        if (position_wild != -1) and (len(HLA_patient_df) != 0):
            if ( epitope_wild not in count_pattern_wd.keys()) or ( epitope_wild in count_pattern_wd.keys() and  (hla_scape not in count_pattern_wd[epitope_wild])):
        
                
               #print(f'entro en wild { (position_wild != -1) and ( epitope_wild not in count_pattern_wd.keys()) or (epitope_wild in count_pattern_sc.keys() and (hla_scape in count_pattern_sc[epitope_wild] == False)) and (len(HLA_patient_df) != 0)}')                           
               if position_metionina== 0:
                   if position_asteris != -1:
                       if position_wild < position_asteris:
                           return dict(correct_orf= True)
                       else:
                           return dict(after_orf= True)
                   else:
                       return dict(without_stop=True)
               else:
                   if position_asteris!= -1:
                       if position_wild < position_asteris:
                           return dict(silent = True)
                       else:
                           return dict(without_met_Pattern_after_stp = True)
                   else:
                       return dict(without_star_stop = True)
            else:
                return {}    
        else:
            return {}
    else:
        return {}
                                              
""" the function write_rows takes the output (dictionary) of the functions write_scape and write_wild, 
evaluates the type of sequence and writes rows with the information sought.
"""
def write_rows(write_scape,write_wild,epitope_wild,epitope_scape,hla_scape, position_wild,position_scape,position_asteris,position_metionina,c):
            if write_scape == True:            
                if correct_orf in write_scape:
                    count_pattern_sc.append(epitope_scape)
                    df.loc[c, 'Epitope_type'] = 'scape_variant'
                    pos_raw_hla_df= HLA_patient_df['HLA'].values.tolist().index(hla_scape)
                    lista_alelos = HLA_patient_df['allele'].values.tolist()
                    df.loc[c, 'allele'] = lista_alelos[ pos_raw_hla_df]
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'correct orf'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_scape
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    return df.loc[c,] 
                if after_orf in write_scape:
                    count_pattern_sc.append(epitope_scape)
                    df.loc[c, 'Epitope_type'] = 'scape_variant'
                    pos_raw_hla_df= HLA_patient_df['HLA'].values.tolist().index(hla_scape)
                    lista_alelos = HLA_patient_df['allele'].values.tolist()
                    df.loc[c, 'allele'] = lista_alelos[ pos_raw_hla_df]
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'after orf'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_scape                
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)                
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    return df.loc[c,]         
                if without_stop in write_scape:
                    count_pattern_sc.append(epitope_scape)
                    df.loc[c, 'Epitope_type'] = 'scape_variant'
                    pos_raw_hla_df= HLA_patient_df['HLA'].values.tolist().index(hla_scape)
                    lista_alelos = HLA_patient_df['allele'].values.tolist()
                    df.loc[c, 'allele'] = lista_alelos[ pos_raw_hla_df]
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'without stop'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_scape                
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)                
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient                                
                    return df.loc[c,]
                if silent:
                    count_pattern_sc.append(epitope_scape)
                    df.loc[c, 'Epitope_type'] = 'scape_variant'
                    pos_raw_hla_df= HLA_patient_df['HLA'].values.tolist().index(hla_scape)
                    lista_alelos = HLA_patient_df['allele'].values.tolist()
                    df.loc[c, 'allele'] = lista_alelos[ pos_raw_hla_df]
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'silent'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_scape                    
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    return df.loc[c,]
                if without_met_Pattern_after_stp:
                    count_pattern_sc.append(epitope_scape)                                                
                    df.loc[c, 'Epitope_type'] = 'scape_variant'
                    pos_raw_hla_df= HLA_patient_df['HLA'].values.tolist().index(hla_scape)
                    lista_alelos = HLA_patient_df['allele'].values.tolist()
                    df.loc[c, 'allele'] = lista_alelos[ pos_raw_hla_df]
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'silent'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_scape                    
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient            
                    return df.loc[c,]
                if whitout_metionina_stop in write_scape:
                    count_pattern_sc.append(epitope_scape)
                    df.loc[c, 'Epitope_type'] = 'scape_variant'
                    pos_raw_hla_df= HLA_patient_df['HLA'].values.tolist().index(hla_scape)
                    lista_alelos = HLA_patient_df['allele'].values.tolist()
                    df.loc[c, 'allele'] = lista_alelos[ pos_raw_hla_df]
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'without Metionine-stop'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_scape
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    return df.loc[c,]
            else:
                #count_pattern_sc.append(epitope_scape)
                df.loc[c, 'Epitope_type'] = 'Not found scape_variant'
                
                df.loc[c, 'pattern_recover'] = 'Not found'
                df.loc[c,'hla'] = 'no'
                df.loc[c, 'position on the sequence'] = 'Not found'
                df.loc[c, 'Genes'] = gen
                #df.loc[c,'count'] = count1
                df.loc[c,'pattern'] = epitope_wild
                df.loc[c,'sequences'] = sequence.id 
                df.loc[c,'position'] = 'Not found'
                #df.loc[c,'positionb'] = positionb
                #df.loc[c,'overlap'] = count
                df.loc[c, "origin"] =   origin
                df.loc[c, 'Donor id'] =  Patient
                return df.loc[c,]
                
            if write_wild == True:             
                if correct_orf in write_wild:
                    count_pattern_wd.append(epitope_wild)
                    df.loc[c, 'Epitope_type'] = 'wild_variant'
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c,'hla'] = 'no'
                    df.loc[c, 'position on the sequence'] = 'correct orf'
                    df.loc[c, 'Genes'] = gen
                    #df.loc[c,'count'] = count1
                    df.loc[c,'pattern'] = epitope_wild
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_wild)
                    #df.loc[c,'positionb'] = positionb
                    #df.loc[c,'overlap'] = count
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    return df.loc[c,]
                if after_orf in write_wild:
                    count_pattern_wd.append(epitope_wild)
                    df.loc[c, 'Epitope_type'] = 'wild_variant'
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'after orf'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'pattern'] = epitope_wild
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_wild)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    return df.loc[c,]
                if without_stop in write_wild:
                    count_pattern_wd.append(epitope_wild)
                    df.loc[c, 'Epitope_type'] = 'wild_variant'
                    #df.loc[c, 'allele'] = lista_alelos[ pos_raw_hla_df]
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'without stop'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = "no" #hla_scape
                    df.loc[c,'pattern'] = epitope_wild                                        
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    return df.loc[c,]
                if silent in write_wild:    
                    count_pattern_wd.append(epitope_wild)
                    df.loc[c, 'Epitope_type'] = 'wild_variant'
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c,'hla'] = 'no'
                    df.loc[c, 'position on the sequence'] = 'silent'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'pattern'] = epitope_wild
                    df.loc[c,'hla'] = "no" #hla_scape
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_wild)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    return df.loc[c,]
                if without_met_Pattern_after_stp in write_wild:
                    count_pattern_wd.append(epitope_wild)
                    df.loc[c, 'Epitope_type'] = 'wild_variant'
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'without metionina after asterisc'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = "no" 
                    df.loc[c,'pattern'] = epitope_wild
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_wild)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    return df.loc[c,]           
                if whitout_metionina_stop in write_wild:
                   count_pattern_wd.append(epitope_wild)
                   df.loc[c, 'Epitope_type'] = 'wild_variant'
                   df.loc[c, 'pattern_recover'] = "yes"
                   df.loc[c, 'position on the sequence'] = 'without metionine-stop'
                   df.loc[c, 'Genes'] = gen
                   df.loc[c,'hla'] = "no" #hla_scape
                   df.loc[c,'pattern'] = epitope_wild                                        
                   df.loc[c,'sequences'] = sequence.id 
                   df.loc[c,'position'] = int(position_wild)
                   df.loc[c, "origin"] =   origin
                   df.loc[c, 'Donor id'] =  Patient            
                   return df.loc[c,]
            else:
                df.loc[c, 'Epitope_type'] = 'Not found wild_variant'
                df.loc[c, 'pattern_recover'] = 'Not found'
                df.loc[c,'hla'] = 'no'
                df.loc[c, 'position on the sequence'] = 'Not found'
                df.loc[c, 'Genes'] = gen
                #df.loc[c,'count'] = count1
                df.loc[c,'pattern'] = epitope_wild
                df.loc[c,'sequences'] = sequence.id 
                df.loc[c,'position'] = 'Not found'
                #df.loc[c,'positionb'] = positionb
                #df.loc[c,'overlap'] = count
                df.loc[c, "origin"] =   origin
                df.loc[c, 'Donor id'] =  Patient
                return df.loc[c,]
            
def dict_epitopes(dictionary, hla, epitope):   
    if epitope not in dictionary.keys():
        dictionary[epitope] = hla
    if epitope  in dictionary.keys() and hla not in dictionary[epitope]: 
        dictionary[epitope] = hla
    return dictionary                                      





#mut = open('ctl_variant_modified.csv')
#mutations = pd.read_csv(mut)    
#open de HLA information
mutation = open('epitopes_corregido_2.csv')
mutations = pd.read_csv(mutation)
mutations.head()
directorio = os.getcwd()
lista_fasta = get_fastas(directorio)
print(lista_fasta)
lista_fasta = lista_fasta[::-1]
lista_fasta = lista_fasta[0:3]
Gen = lista_fasta
df = pd.DataFrame()  
dicgenes = {}
c = 0
#count_pattern = []
for l, Gen in enumerate(lista_fasta):
    #print(str(gen))
    fasta = SeqIO.parse(Gen,'fasta') #parse

    print(f'entro con el gen {Gen}')
    for n, sequence in enumerate(fasta):
        #print(f'comenzando la sequencia {sequence} numero {n}')
        #dic = {}
        gen = Gen.split('_')[0] 
        #gen  = "gag"
        mutations['Protein'] = mutations['Protein'].str.lower()
        mutations_gen = mutations[mutations.Protein.isin([gen])]
        ids =  sequence.id
        split = ids.split('-')
        if len(split) == 3:
            Patient = ids.split('-')[0]
            origin = ids.split('-')[1]
            year_sample = "-"
        if len(split) == 4:
            Patient = ids.split('-')[0]
            Patient =str(Patient)
            year_sample = ids.split('-')[0]
            origin = ids.split('-')[2]
        #mutations_gen = mutations_gen.copy()
        #len(mutations_gen)
        count_pattern_sc = {}
        count_pattern_wd = {}
        HLA_patient_df = filters_mut(mutations_gen,Patient )
        #c = mutations_gen[mutations_gen.Patient_ID == 1292]#.unique()
        #c["HLA"].unique()
        for (index, row)  in enumerate(HLA_patient_df.iterrows()):
            ###function captar posiciones
                  
            epitope_wild,epitope_scape,hla_scape, position_wild,position_scape,position_asteris,position_metionina = hunt_mutations(HLA_patient_df,index,sequence)
            #hla_scape,Patient,position_scape,position_wild)
            #print(f'el paciente {Patient} tiene un df de {len(HLA_patient_df)}')
            dict_scape = write_scape(position_scape, hla_scape, position_asteris, position_metionina, HLA_patient_df,epitope_scape, count_pattern_sc)
            #print(f'diccionario scap {dict_scape} ')
            
                
            dict_wild = write_wild(position_wild,hla_scape, position_asteris, position_metionina,   HLA_patient_df, epitope_wild, count_pattern_wd)
            
            #write_rows(write_scape,write_wild,epitope_wild,epitope_scape,hla_scape, position_wild,position_scape,position_asteris,position_metionina,c)
            
            #print(f' diccionario wild {dict_wild}')
            if len(dict_scape) != 0:            
                if "correct_orf" in dict_scape.keys():
                    count_pattern_sc = dict_epitopes(count_pattern_sc, hla_scape, epitope_scape) #.append(hla_scape)
                    #count_pattern_sc.append(epitope_scape) 
                    df.loc[c, 'Epitope_type'] = 'scape_variant'
                    #pos_raw_hla_df= HLA_patient_df['HLA'].values.tolist().index(hla_scape)
                    #lista_alelos = HLA_patient_df['allele'].values.tolist()
                    #df.loc[c, 'allele'] = lista_alelos[ pos_raw_hla_df]
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'correct orf'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_scape
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    c += 1
                if "after_orf" in dict_scape.keys():
                    count_pattern_sc = dict_epitopes(count_pattern_sc, hla_scape, epitope_scape)
                    #count_pattern_sc.append(epitope_scape) 
                    df.loc[c, 'Epitope_type'] = 'scape_variant'
                    #os_raw_hla_df= HLA_patient_df['HLA'].values.tolist().index(hla_scape)
                    #lista_alelos = HLA_patient_df['allele'].values.tolist()
                    #df.loc[c, 'allele'] = lista_alelos[ pos_raw_hla_df]
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'after orf'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_scape                
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)                
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    c += 1
                if "without_stop" in dict_scape.keys():
                    count_pattern_sc = dict_epitopes(count_pattern_sc, hla_scape, epitope_scape)
                    #count_pattern_sc.append(epitope_scape) 
                    df.loc[c, 'Epitope_type'] = 'scape_variant'
                    #pos_raw_hla_df= HLA_patient_df['HLA'].values.tolist().index(hla_scape)
                    #lista_alelos = HLA_patient_df['allele'].values.tolist()
                    #df.loc[c, 'allele'] = lista_alelos[ pos_raw_hla_df]
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'without stop'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_scape                
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)                
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient                                
                    c += 1
                if "silent" in dict_scape.keys():
                    count_pattern_sc = dict_epitopes(count_pattern_sc, hla_scape, epitope_scape)
                    #count_pattern_sc.append(epitope_scape) 
                    df.loc[c, 'Epitope_type'] = 'scape_variant'
                    #pos_raw_hla_df= HLA_patient_df['HLA'].values.tolist().index(hla_scape)
                    #lista_alelos = HLA_patient_df['allele'].values.tolist()
                    #df.loc[c, 'allele'] = lista_alelos[ pos_raw_hla_df]
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'silent'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_scape                    
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    c += 1
                if "without_met_Pattern_after_stp" in dict_scape.keys():
                    count_pattern_sc = dict_epitopes(count_pattern_sc, hla_scape, epitope_scape)
                    #count_pattern_sc.append(epitope_scape)                                                 
                    df.loc[c, 'Epitope_type'] = 'scape_variant'
                    #pos_raw_hla_df= HLA_patient_df['HLA'].values.tolist().index(hla_scape)
                    #lista_alelos = HLA_patient_df['allele'].values.tolist()
                    #df.loc[c, 'allele'] = lista_alelos[ pos_raw_hla_df]
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'silent'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_scape                    
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient            
                    c += 1
                if "without_star_stop" in dict_scape.keys():
                    count_pattern_sc = dict_epitopes(count_pattern_sc, hla_scape, epitope_scape)
                    #count_pattern_sc.append(epitope_scape) 
                    df.loc[c, 'Epitope_type'] = 'scape_variant'
                    #pos_raw_hla_df= HLA_patient_df['HLA'].values.tolist().index(hla_scape)
                    #lista_alelos = HLA_patient_df['allele'].values.tolist()
                    #df.loc[c, 'allele'] = lista_alelos[ pos_raw_hla_df]
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'without Metionine-stop'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_scape
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    c += 1
                     
            if len(dict_scape) == 0:
                count_pattern_sc = dict_epitopes(count_pattern_sc, hla_scape, epitope_scape)
                #count_pattern_sc.append(epitope_scape)
                df.loc[c, 'Epitope_type'] = 'scape_variant'
                df.loc[c, 'pattern_recover'] = 'Not found scape_variant'
                df.loc[c,'hla'] = 'no'
                df.loc[c, 'position on the sequence'] = 'Not found'
                df.loc[c, 'Genes'] = gen
                #df.loc[c,'count'] = count1
                df.loc[c,'pattern'] = epitope_wild
                df.loc[c,'sequences'] = sequence.id 
                df.loc[c,'position'] = 'Not found'
                #df.loc[c,'positionb'] = positionb
                #df.loc[c,'overlap'] = count
                df.loc[c, "origin"] =   origin
                df.loc[c, 'Donor id'] =  Patient
                    
                    
                
            if len(dict_wild) != 0:             
                if "correct_orf" in dict_wild.keys():
                    count_pattern_wd = dict_epitopes(count_pattern_wd, hla_scape, epitope_wild)
                    #count_pattern_wd.append(epitope_wild)
                    df.loc[c, 'Epitope_type'] = 'wild_variant'
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c, 'position on the sequence'] = 'correct orf'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'pattern'] = epitope_wild
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_wild)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    c += 1
                if "after_orf" in dict_wild.keys():
                    count_pattern_wd = dict_epitopes(count_pattern_wd, hla_scape, epitope_wild)
                    #count_pattern_wd.append(epitope_wild)
                    df.loc[c, 'Epitope_type'] = 'wild_variant'
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'after orf'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_wild
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_wild)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    c += 1
                if "without_stop" in dict_wild.keys():
                    count_pattern_wd = dict_epitopes(count_pattern_wd, hla_scape, epitope_wild)
                    #count_pattern_wd.append(epitope_wild)
                    df.loc[c, 'Epitope_type'] = 'wild_variant'
                    #df.loc[c, 'allele'] = lista_alelos[ pos_raw_hla_df]
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'without stop'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_wild                                        
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                if "silent" in dict_wild.keys():    
                    count_pattern_wd = dict_epitopes(count_pattern_wd, hla_scape, epitope_wild)
                    #count_pattern_wd.append(epitope_wild)
                    df.loc[c, 'Epitope_type'] = 'wild_variant'
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c,'hla'] = 'no'
                    df.loc[c, 'position on the sequence'] = 'silent'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'pattern'] = epitope_wild
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_wild)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    c += 1
                if "without_met_Pattern_after_stp" in dict_wild.keys():
                    count_pattern_wd = dict_epitopes(count_pattern_wd, hla_scape, epitope_wild)
                    #count_pattern_wd.append(epitope_wild)
                    df.loc[c, 'Epitope_type'] = 'wild_variant'
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'without metionina Epitafter asterisc'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_wild
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_wild)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    c += 1                
                if "without_star_stop" in dict_wild.keys():
                   count_pattern_wd = dict_epitopes(count_pattern_wd, hla_scape, epitope_wild)
                   #count_pattern_wd.append(epitope_wild)
                   df.loc[c, 'Epitope_type'] = 'wild_variant'
                   df.loc[c, 'pattern_recover'] = "yes"
                   df.loc[c, 'position on the sequence'] = 'without metionine-stop'
                   df.loc[c, 'Genes'] = gen
                   df.loc[c,'hla'] = hla_scape
                   df.loc[c,'pattern'] = epitope_wild                                        
                   df.loc[c,'sequences'] = sequence.id 
                   df.loc[c,'position'] = int(position_wild)
                   df.loc[c, "origin"] =   origin
                   df.loc[c, 'Donor id'] =  Patient            
                   c += 1
            
            if len(dict_wild) == 0: 
                count_pattern_wd = dict_epitopes(count_pattern_wd, hla_scape, epitope_wild)
                #count_pattern_wd.append(epitope_wild)
                df.loc[c, 'Epitope_type'] = 'wild_variant'
                df.loc[c, 'pattern_recover'] = 'Not found wild_variant'
                df.loc[c,'hla'] = 'no'
                df.loc[c, 'position on the sequence'] = 'Not found'
                df.loc[c, 'Genes'] = gen
                #df.loc[c,'count'] = count1
                df.loc[c,'pattern'] = epitope_wild
                df.loc[c,'sequences'] = sequence.id 
                df.loc[c,'position'] = 'Not found'
                #df.loc[c,'positionb'] = positionb
                #df.loc[c,'overlap'] = count
                df.loc[c, "origin"] =   origin
                df.loc[c, 'Donor id'] =  Patient

                                       
            """if int(index + 1) == int(len(mutations_gen['Epitope_WT'])):          
                if len(count_pattern_sc)== 0 and len(count_pattern_wd)== 0:           
                    
                    #print(f'{c} nada en gen {gen} y paciente {Patient}')
                    df.loc[c, 'Epitope_type'] = 'no'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    df.loc[c, 'pattern_recover'] = 'no'
                    df.loc[c,'pattern'] = 'no'
                    #df.loc[c,'count'] = 'no'
                    #df.loc[c,'position'] = 'no'
                    #df.loc[c,'overlap'] = 'no'
                    df.loc[c,'hla'] = 'no'
                    c += 1"""
                    
print('grabando...')
df.to_csv('patterns_30may.csv')
                                        
                                    
                                    
        
#HLA_df= mutations_gen[mutations_gen['HLA']== hla_scape]
#HLA_patient_df = HLA_df[HLA_df.Patient_ID.astype(str)== str(Patient)]
print('grabando...')
df.to_csv('res_mutation_test_30May.csv')




#mejora de funciones


            
print('grabando...')
            
