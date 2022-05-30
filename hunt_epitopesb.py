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
import os
import itertools
import re

os.getcwd()
os.listdir()
os.chdir('C:\\Users\\eares\\OneDrive\\Documentos\\HIV_scape_mut\\datos\\nuevos\\epitopes_pipi')
os.getcwd()
os.listdir()
os.chdir('/home/bioinformatica/Documentos/Documentos/nuevos')
os.chdir('/home/bioinfo/Documents/Epitopes_pipi')

os.chdir('/home/fernando/Documentos/Epitopes_pipi/redesign_inmuno_hiv')
###Functions
##obtain fasta sequences from directory specified
def get_fastas(directorio):
    lista = []
    for root, dirs, files in os.walk(directorio): #camino por el directorio
        for file in files:
            if 'fas' in file: 
                lista.append(file)
                
    return lista
## Detect a mutation, their type, their position, position of start and stop of ORF.
##Need a mutation database and index sequence to iterate in each cell to the WT and variant epitope columns.
##return variables with information if was detected signals of mutations wild and variants with their corresponding information of position anda origin.
def hunt_mutations(db_mut, number_index,sequence):
    epitope_wild = db_mut["Epitope_WT"].iloc[number_index]
    epitope_scape = db_mut["Variant_Epitope"].iloc[number_index]
    hla_scape = db_mut["HLA"].iloc[number_index]
    position_wild = sequence.seq.find(str( epitope_wild))
    position_scape = sequence.seq.find(str( epitope_scape))
    position_asteris = sequence.seq.find(str('*'))
    position_metionina = sequence.seq.find(str('M'))
    ids =  sequence.id
    split = ids.split('-')
    ##process of header to obtain information(code-cell origin and year) of each pacient:
    if len(split) == 3:
        Patient = ids.split('-')[0]
        origin = ids.split('-')[1]
        year_sample = "-"
    if len(split) == 4:
        Patient = ids.split('-')[0]
        Patient =str(Patient)
        year_sample = ids.split('-')[0]
        origin = ids.split('-')[2]
        
    return epitope_wild,epitope_scape,hla_scape, position_wild,position_scape,position_asteris,position_metionina,Patient, origin,year_sample

mutations.head()    
##This function filters the mutations that match with HLA related to mutation scape
def filters_mut_scape(db_mut,hla_scape,Patient):#,position_scape,position_wild):
    HLA_df= db_mut[db_mut['HLA']== hla_scape] #only stayed mutations equals to secuences with mutation scape
    HLA_patient_df = HLA_df[HLA_df.Patient_ID.astype(str)== str(Patient)] #the mutations pass to a new filter that only retain sequence related to pacient
    return HLA_patient_df

####This function filters the mutations that match with HLA related to mutation wild
def filters_mut(db_mut,hla_scape,Patient):#,position_scape,position_wild):
    HLA_df= db_mut[db_mut['HLA']== hla_scape]
    HLA_patient_df = HLA_df[HLA_df.Patient_ID.astype(str)== str(Patient)]
    return HLA_patient_df



 """ The followed functions classify the secuence from data information obtain previusly.
 """
def write_scape(epitope_scape,position_scape,count_pattern_sc,hla_scape, position_asteris, position_metionina,   HLA_patient_df):
    if (position_scape != -1) and (len(HLA_patient_df) != 0):
        #the folowed sentence prevent that a mutation 
        if ( epitope_scape not in count_pattern_sc.keys()) or ( epitope_scape in count_pattern_sc.keys() and  hla_scape not in count_pattern_sc[epitope_scape]):
            if position_metionina== 0:
                if position_asteris != -1:
                    if position_scape < position_asteris:                        
                        return dict(correct_orf= True)
                    else:                        
                        return dict(after_orf= True)
                else:                        
                        return dict(whitout_stop=True)
            else: 
                if position_asteris!= -1:
                    if position_scape < position_asteris:                        
                        return dict(without_metionina_B_stop = True)
                    else:                        
                        return dict(without_metionina_stop = True)
                else:                    
                    return dict(without_metionina_stop = True)     
        else:
            return dict()        
    else:
        return dict()
        
def write_wild(epitope_wild,position_wild,count_pattern_wd,hla_scape, position_asteris, position_metionina,   HLA_patient_df):
    if (position_wild != -1) and (len(HLA_patient_df) != 0):
        if ( epitope_wild not in count_pattern_wd.keys()) or ( epitope_wild in count_pattern_wd.keys() and  hla_scape not in count_pattern_wd[epitope_wild]):
            if position_metionina== 0:
                if position_asteris != -1:
                    if position_wild < position_asteris:                        
                        return dict(correct_orf= True)
                    else:                        
                        return dict(after_orf= True)
                else:                        
                        return dict(whitout_stop=True)
            else: 
                if position_asteris!= -1:
                    if position_wild < position_asteris:                        
                        return dict(without_metionina_B_stop = True)
                    else:                        
                        return dict(without_metionina_stop = True)
                else:                    
                    return dict(without_metionina_stop = True)     
        else:
            return dict()        
    else:
        return dict()
        
                                              
        
def write_rows(write_scape,write_wild):
            if write_scape == True:            
                if correct_orf in write_scape:
                    count_pattern_sc.append(epitope_scape)
                    df.loc[c, 'Epitope_recover'] = 'scape_variant'
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
                    df.loc[c, 'Epitope_recover'] = 'scape_variant'
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
                if whitout_stop in write_scape:
                    count_pattern_sc.append(epitope_scape)
                    df.loc[c, 'Epitope_recover'] = 'scape_variant'
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
                if whitout_metionina_B_asterisc:
                    count_pattern_sc.append(epitope_scape)
                    df.loc[c, 'Epitope_recover'] = 'scape_variant'
                    pos_raw_hla_df= HLA_patient_df['HLA'].values.tolist().index(hla_scape)
                    lista_alelos = HLA_patient_df['allele'].values.tolist()
                    df.loc[c, 'allele'] = lista_alelos[ pos_raw_hla_df]
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'without metionina_before asterisc'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_scape                    
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    return df.loc[c,]
                if whitout_metionina_A_asterisc:
                    count_pattern_sc.append(epitope_scape)                                                
                    df.loc[c, 'Epitope_recover'] = 'scape_variant'
                    pos_raw_hla_df= HLA_patient_df['HLA'].values.tolist().index(hla_scape)
                    lista_alelos = HLA_patient_df['allele'].values.tolist()
                    df.loc[c, 'allele'] = lista_alelos[ pos_raw_hla_df]
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'without metionina_before asterisc'
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
                    df.loc[c, 'Epitope_recover'] = 'scape_variant'
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
            if write_wild == True:             
                if correct_orf in write_wild:
                    count_pattern_wd.append(epitope_wild)
                    df.loc[c, 'Epitope_recover'] = 'wild_variant'
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
                    df.loc[c, 'Epitope_recover'] = 'wild_variant'
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'after orf'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'pattern'] = epitope_wild
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_wild)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    return df.loc[c,]
                if whitout_stop in write_wild:
                    count_pattern_wd.append(epitope_wild)
                    df.loc[c, 'Epitope_recover'] = 'wild_variant'
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
                if whitout_metionina_B_asterisc in write_wild:    
                    count_pattern_wd.append(epitope_wild)
                    df.loc[c, 'Epitope_recover'] = 'wild_variant'
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c,'hla'] = 'no'
                    df.loc[c, 'position on the sequence'] = 'without metionina_before asterisc'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'pattern'] = epitope_wild
                    df.loc[c,'hla'] = "no" #hla_scape
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_wild)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    return df.loc[c,]
                if whitout_metionina_A_asterisc in write_wild:
                    count_pattern_wd.append(epitope_wild)
                    df.loc[c, 'Epitope_recover'] = 'wild_variant'
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'without metionina Epitafter asterisc'
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
                   df.loc[c, 'Epitope_recover'] = 'wild_variant'
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

def dict_epitopes(dictionary, hla, epitope):   
    if epitope not in dictionary.keys():
        dictionary[epitope] = []
    else: 
        dictionary[epitope].append(hla)
#mutations = pd.read_csv(mut)    
    return dictionary                                      
    
#mut = open('ctl_variant_modified.csv')
#open de HLA information
#Open de epitopes info table.
mutation = open('epitopes_corregido_2.csv')
mutations = pd.read_csv(mutation)
mutations.head()
#Get the directory with fasta files 
directorio = os.getcwd()
lista_fasta = get_fastas(directorio)
print(lista_fasta)
lista_fasta = lista_fasta[::-1]
#Gen = lista_fasta
#Create dataset and dictionary to save data from protein analisis
df = pd.DataFrame()  
dicgenes = {}
c = 0
#count_pattern = []
for l, Gen in enumerate(lista_fasta, start= 1):
    #print(str(gen))
    fasta = SeqIO.parse(Gen,'fasta') #parse
    gen = Gen.split('_')[0] 
    mutations['Protein'] = mutations['Protein'].str.lower()
    mutations_gen = mutations[mutations.Protein.isin([gen])]
    mutations_gen = mutations_gen.copy()
    print(f'entro con el gen {Gen}')
    for n, sequence in enumerate(fasta):# , start= 2975):
        #print(f'comenzando la sequencia {sequence} numero {n}')
        #dic = {}
        
        len(mutations_gen)
        count_pattern_sc = {}
        count_pattern_wd = {}
        for (index, row)  in enumerate(mutations_gen.iterrows()):
            ###function captar posiciones
                  
            epitope_wild,epitope_scape,hla_scape,position_wild,position_scape,position_asteris,position_metionina,Patient, origin,year_sample = hunt_mutations(mutations_gen,index,sequence)
            HLA_patient_df = filters_mut(mutations_gen,hla_scape,Patient)#,position_scape,position_wild)
            #print(f'el paciente {Patient} tiene un df de {len(HLA_patient_df)}')
            dict_scape = write_scape(epitope_scape, position_scape, count_pattern_sc,hla_scape, position_asteris, position_metionina, HLA_patient_df)
            #print(f'diccionario scap {dict_scape} ')
            
                
            #dict_wild = write_wild(epitope_wild, position_wild, count_pattern_wd, hla_scape, position_asteris, position_metionina,   HLA_patient_df)
        
            #print(f' diccionario wild {dict_wild}')
            if len(dict_scape) != 0:            
                if "correct_orf" in dict_scape.keys():
                    dict_epitopes(count_pattern_sc, hla_scape, epitope_scape)
                        
                                  
                     #.append(hla_scape)
                    
                    df.loc[c, 'Epitope_recover'] = 'scape_variant'
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
                    c += 1
                if "after_orf" in dict_scape.keys():
                    dict_epitopes(count_pattern_sc, hla_scape, epitope_scape)
                    df.loc[c, 'Epitope_recover'] = 'scape_variant'
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
                    c += 1
                if "whitout_stop" in dict_scape.keys():
                    dict_epitopes(count_pattern_sc, hla_scape, epitope_scape)
                    df.loc[c, 'Epitope_recover'] = 'scape_variant'
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
                    c += 1
                if "whitout_metionina_B_asterisc" in dict_scape.keys():
                    dict_epitopes(count_pattern_sc, hla_scape, epitope_scape)
                    df.loc[c, 'Epitope_recover'] = 'scape_variant'
                    pos_raw_hla_df= HLA_patient_df['HLA'].values.tolist().index(hla_scape)
                    lista_alelos = HLA_patient_df['allele'].values.tolist()
                    df.loc[c, 'allele'] = lista_alelos[ pos_raw_hla_df]
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'without metionina_before asterisc'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_scape                    
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)
                    #df.loc[c,'positionb'] = positionb
                    #df.loc[c,'overlap'] = count
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    c += 1
                if "whitout_metionina_A_asterisc" in dict_scape.keys():
                    dict_epitopes(count_pattern_sc, hla_scape, epitope_scape)                                                
                    df.loc[c, 'Epitope_recover'] = 'scape_variant'
                    pos_raw_hla_df= HLA_patient_df['HLA'].values.tolist().index(hla_scape)
                    lista_alelos = HLA_patient_df['allele'].values.tolist()
                    df.loc[c, 'allele'] = lista_alelos[ pos_raw_hla_df]
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c, 'position on the sequence'] = 'without metionina_before asterisc'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'pattern'] = epitope_scape                    
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_scape)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient            
                    c += 1
                if "whitout_metionina_stop" in dict_scape.keys():
                    dict_epitopes(count_pattern_sc, hla_scape, epitope_scape)
                    df.loc[c, 'Epitope_recover'] = 'scape_variant'
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
                    #df.loc[c,'positionb'] = positionb
                    #df.loc[c,'overlap'] = count
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    c += 1
                     
                                        
                    
                    
           
            if int(index + 1) == int(len(mutations_gen['Epitope_WT'])):          
                if (len(count_pattern_sc) == 0) and (len(count_pattern_wd) == 0):           
                    
                    #print(f'{c} nada en gen {gen} y paciente {Patient}')
                    df.loc[c, 'Epitope_recover'] = 'no'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'sequences'] = sequence.id 
                    #df.loc[c,'positionb'] = positionb
                    #df.loc[c,'overlap'] = count
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    df.loc[c, 'pattern_recover'] = 'no'
                    df.loc[c,'pattern'] = 'no'
                    df.loc[c,'count'] = 'no'
                    df.loc[c,'position'] = 'no'
                    #df.loc[c,'positionb'] = positionb
                    df.loc[c,'overlap'] = 'no'
                    df.loc[c,'hla'] = 'no'
                    c += 1
                                        
                                    
                                    
        
#HLA_df= mutations_gen[mutations_gen['HLA']== hla_scape]
#HLA_patient_df = HLA_df[HLA_df.Patient_ID.astype(str)== str(Patient)]
print('grabando...')
df.to_csv('mutation_29may.csv')

######suspension temporal

     
            if len(dict_wild) != 0:             
                if "correct_orf" in dict_wild.keys():
                    dict_epitopes(count_pattern_wd, hla_scape, epitope_wild)
                    df.loc[c, 'Epitope_recover'] = 'wild_variant'
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c,'hla'] = hla_scape
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
                    c += 1
                if "after_orf" in dict_wild.keys():
                    dict_epitopes(count_pattern_wd, hla_scape, epitope_wild)
                    df.loc[c, 'Epitope_recover'] = 'wild_variant'
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
                if "whitout_stop" in dict_wild.keys():
                    dict_epitopes(count_pattern_wd, hla_scape, epitope_wild)
                    df.loc[c, 'Epitope_recover'] = 'wild_variant'
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
                if "whitout_metionina_B_asterisc" in dict_wild.keys():    
                    dict_epitopes(count_pattern_wd, hla_scape, epitope_wild)
                    df.loc[c, 'Epitope_recover'] = 'wild_variant'
                    df.loc[c, 'pattern_recover'] = "yes"
                    df.loc[c,'hla'] = 'no'
                    df.loc[c, 'position on the sequence'] = 'without metionina_before asterisc'
                    df.loc[c, 'Genes'] = gen
                    df.loc[c,'pattern'] = epitope_wild
                    df.loc[c,'hla'] = hla_scape
                    df.loc[c,'sequences'] = sequence.id 
                    df.loc[c,'position'] = int(position_wild)
                    df.loc[c, "origin"] =   origin
                    df.loc[c, 'Donor id'] =  Patient
                    c += 1
                if "whitout_metionina_A_asterisc" in dict_wild.keys():
                    dict_epitopes(count_pattern_wd, hla_scape, epitope_wild)
                    df.loc[c, 'Epitope_recover'] = 'wild_variant'
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
                if "without_metionina_stop" in dict_wild.keys():
                   dict_epitopes(count_pattern_wd, hla_scape, epitope_wild)
                   df.loc[c, 'Epitope_recover'] = 'wild_variant'
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
                             
                                       

#mejora de funciones


            
print('grabando...')
            
df.to_csv('mutation_crrctn_21M2.csv')
