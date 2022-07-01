#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
from  functions_hiv import *

os.chdir('/home/fernando/Documentos/Epitopes_pipi/redesign_inmuno_hiv')

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
                    df.loc[c, 'pattern_pos_seq'] = 'correct orf'
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
                    df.loc[c, 'pattern_pos_seq'] = 'after orf'
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
                    df.loc[c, 'pattern_pos_seq'] = 'without stop'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_scape                
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)                
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient                                
                    c += 1
                if "without_met_Pattern_before_stp" in dict_scape.keys():
                    count_pattern_sc = dict_epitopes(count_pattern_sc, hla_scape, epitope_scape)
                    #count_pattern_sc.append(epitope_scape) 
                    df.loc[c, 'Epitope_type'] = 'scape_variant'
                    #pos_raw_hla_df= HLA_patient_df['HLA'].values.tolist().index(hla_scape)
                    #lista_alelos = HLA_patient_df['allele'].values.tolist()
                    #df.loc[c, 'allele'] = lista_alelos[ pos_raw_hla_df]
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'pattern_pos_seq'] = 'without_met_Pattern_before_stp'
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
                    df.loc[c, 'pattern_pos_seq'] = "without_met_Pattern_after_stp"
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_scape                    
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient            
                    c += 1
                if "without_start_stop" in dict_scape.keys():
                    count_pattern_sc = dict_epitopes(count_pattern_sc, hla_scape, epitope_scape)
                    #count_pattern_sc.append(epitope_scape) 
                    df.loc[c, 'Epitope_type'] = 'scape_variant'
                    #pos_raw_hla_df= HLA_patient_df['HLA'].values.tolist().index(hla_scape)
                    #lista_alelos = HLA_patient_df['allele'].values.tolist()
                    #df.loc[c, 'allele'] = lista_alelos[ pos_raw_hla_df]
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'pattern_pos_seq'] = 'without start-stop'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_scape
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    c += 1
            
            '''if len(dict_scape) == 0:
                #count_pattern_sc = dict_epitopes(count_pattern_sc, hla_scape, epitope_scape)
                #count_pattern_sc.append(epitope_scape)
                df.loc[c, 'Epitope_type'] = 'scape_variant'
                df.loc[c, 'pattern_recover'] = 'Not found scape_variant'
                df.loc[c,'hla'] = 'no'
                df.loc[c, 'pattern_pos_seq'] = 'Not found'
                df.loc[c, 'Genes'] = gen
                #df.loc[c,'count'] = count1
                df.loc[c,'pattern'] = epitope_wild
                df.loc[c,'sequences'] = sequence.id 
                df.loc[c,'position'] = 'Not found'
                #df.loc[c,'positionb'] = positionb
                #df.loc[c,'overlap'] = count
int(len(mutations_gen['Epitope_WT']))                df.loc[c, "origin"] =   origin
                df.loc[c, 'Donor id'] =  Patient'''
            
            
            
            if len(dict_wild) != 0:             
                if "correct_orf" in dict_wild.keys():
                    count_pattern_wd = dict_epitopes(count_pattern_wd, hla_scape, epitope_wild)
                    #count_pattern_wd.append(epitope_wild)
                    df.loc[c, 'Epitope_type'] = 'wild_variant'
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c, 'pattern_pos_seq'] = 'correct orf'
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
                    df.loc[c, 'pattern_pos_seq'] = 'after orf'
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
                    df.loc[c, 'pattern_pos_seq'] = 'without stop'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_wild                                        
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                if "without_met_Pattern_before_stp" in dict_wild.keys():    
                    count_pattern_wd = dict_epitopes(count_pattern_wd, hla_scape, epitope_wild)
                    #count_pattern_wd.append(epitope_wild)
                    df.loc[c, 'Epitope_type'] = 'wild_variant'
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c,'hla'] = 'no'
                    df.loc[c, 'pattern_pos_seq'] = 'without_met_Pattern_before_stp'
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
                    df.loc[c, 'pattern_pos_seq'] = 'without_met_Pattern_after_stp'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_wild
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_wild)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    c += 1                
                if "without_start_stop" in dict_wild.keys():
                   count_pattern_wd = dict_epitopes(count_pattern_wd, hla_scape, epitope_wild)
                   #count_pattern_wd.append(epitope_wild)
                   df.loc[c, 'Epitope_type'] = 'wild_variant'
                   df.loc[c, 'pattern_recover'] = "yes"
                   df.loc[c, 'pattern_pos_seq'] = 'without start-stop'
                   df.loc[c, 'Genes'] = gen
                   df.loc[c,'hla'] = hla_scape
                   df.loc[c,'pattern'] = epitope_wild                                        
                   df.loc[c,'sequences'] = sequence.id 
                   df.loc[c,'position'] = int(position_wild)
                   df.loc[c, "origin"] =   origin
                   df.loc[c, 'Donor id'] =  Patient            
                   c += 1
            
            """if len(dict_wild) == 0: 
                #count_pattern_wd = dict_epitopes(count_pattern_wd, hla_scape, epitope_wild)
                #count_pattern_wd.append(epitope_wild)
                df.loc[c, 'Epitope_type'] = 'wild_variant'
                df.loc[c, 'pattern_recover'] = 'Not found wild_variant'
                df.loc[c,'hla'] = 'no'
                df.loc[c, 'pattern_pos_seq'] = 'Not found'
                df.loc[c, 'Genes'] = gen
                #df.loc[c,'count'] = count1
                df.loc[c,'pattern'] = epitope_wild
                df.loc[c,'sequences'] = sequence.id 
                df.loc[c,'position'] = 'Not found'
                #df.loc[c,'positionb'] = positionb
                #df.loc[c,'overlap'] = count
                df.loc[c, "origin"] =   origin
                df.loc[c, 'Donor id'] =  Patient"""
            
            
            if int(index + 1) == int(len(HLA_patient_df['Epitope_WT'])):          
                if len(count_pattern_sc)== 0:             
                    
                    #print(f'{c} nada en gen {gen} y paciente {Patient}')
                    df.loc[c, 'Epitope_type'] = 'not found vr scape'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    df.loc[c, 'pattern_recover'] = 'no'
                    df.loc[c,'pattern'] = 'no'
                    df.loc[c,'hla'] = 'no'
                    c += 1
                if len(count_pattern_wd)== 0:
                    df.loc[c, 'Epitope_type'] = 'not found vr wild'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    df.loc[c, 'pattern_recover'] = 'no'
                    df.loc[c,'pattern'] = 'no'
                    df.loc[c,'hla'] = 'no'
                    c += 1
df.to_csv('patterns_1jul.csv')
