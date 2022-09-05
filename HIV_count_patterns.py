#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
current= os.getcwd() 
os.chdir(current)
#os.chdir('/home/fernando/Documentos/Epitopes_pipi/redesign_inmuno_hiv')

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



if len(sys.argv) != 3:
    raise SystemExit(f'correct usage: {sys.argv[0]} ' 'input csv, output csv', f'arguments now: {sys.argv}')            
            
csv_input = sys.argv[1]
csv_output = sys.argv[2]
            
#mutation = open('epitopes_corregido_2.csv')
mutation = open(csv_input)
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

for l, Gen in enumerate(lista_fasta):
    #print(str(gen))
    fasta = SeqIO.parse(Gen,'fasta') #parse
    
    print(f' Analyzing the gene {Gen} ...')
    for n, sequence in enumerate(fasta):
        gen = Gen.split('_')[0] 
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
        
        count_pattern_sc = {}
        count_pattern_wd = {}
        HLA_patient_df = filters_mut(mutations_gen,Patient ) #filter HLA pattern related to donnor information
     
        for (index, row)  in enumerate(HLA_patient_df.iterrows()):
            ### Obtains the position of the patterns
            
            epitope_wild,epitope_scape,hla_scape, position_wild,position_scape,position_asteris,position_metionina = hunt_mutations(HLA_patient_df,index,sequence)
            
            dict_scape = write_scape(position_scape, hla_scape, position_asteris, position_metionina, HLA_patient_df,epitope_scape, count_pattern_sc)
           
            
            
            dict_wild = write_wild(position_wild,hla_scape, position_asteris, position_metionina,   HLA_patient_df, epitope_wild, count_pattern_wd)
            
            
            if len(dict_scape) != 0:            
                if "correct_orf" in dict_scape.keys():
                    count_pattern_sc = dict_epitopes(count_pattern_sc, hla_scape, epitope_scape) #.append(hla_scape)
                    df.loc[c, 'Epitope_type'] = 'scape_variant'
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
                    df.loc[c, 'Epitope_type'] = 'scape_variant'
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
                     
                    df.loc[c, 'Epitope_type'] = 'scape_variant'
                    
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
                    df.loc[c, 'Epitope_type'] = 'scape_variant'                    
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
                                                                   
                    df.loc[c, 'Epitope_type'] = 'scape_variant'                    
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
                    df.loc[c, 'Epitope_type'] = 'scape_variant'
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
                    df.loc[c, 'Epitope_type'] = 'wild_variant'
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c, 'pattern_pos_seq'] = 'Pattern_in_orf'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'pattern'] = epitope_wild
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_wild)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    c += 1
                if "after_orf" in dict_wild.keys():
                    count_pattern_wd = dict_epitopes(count_pattern_wd, hla_scape, epitope_wild)
                    df.loc[c, 'Epitope_type'] = 'wild_variant'
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'pattern_pos_seq'] = 'Pattern_after_orf'
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
                    df.loc[c, 'Epitope_type'] = 'wild_variant'
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'pattern_pos_seq'] = 'without_stop'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_wild                                        
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                if "without_met_Pattern_before_stp" in dict_wild.keys():    
                    count_pattern_wd = dict_epitopes(count_pattern_wd, hla_scape, epitope_wild)
                    df.loc[c, 'Epitope_type'] = 'wild_variant'
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c,'hla'] = 'no'
                    df.loc[c, 'pattern_pos_seq'] = 'without_start_Pattern_before_stop'
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
                    df.loc[c, 'Epitope_type'] = 'wild_variant'
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'pattern_pos_seq'] = 'without_start_Pattern_after_stop'
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
                   df.loc[c, 'Epitope_type'] = 'wild_variant'
                   df.loc[c, 'pattern_recover'] = "yes"
                   df.loc[c, 'pattern_pos_seq'] = 'without_start-stop'
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
                    
                    df.loc[c, 'Epitope_type'] = 'not_found_vr_scape'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    df.loc[c, 'pattern_recover'] = 'no'
                    df.loc[c,'pattern'] = 'no'
                    df.loc[c,'hla'] = 'no'
                    c += 1
                if len(count_pattern_wd)== 0:
                    df.loc[c, 'Epitope_type'] = 'not_found_vr_wild'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    df.loc[c, 'pattern_recover'] = 'no'
                    df.loc[c,'pattern'] = 'no'
                    df.loc[c,'hla'] = 'no'
                    c += 1
#df.to_csv('patterns_3sptember.csv', index=False)
df.to_csv(csv_outout , index=False)
sys.exit('The process has ended')