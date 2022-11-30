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

#os.getcwd()
#os.listdir()
#os.chdir('/home/fernando/Documentos/Epitopes_pipi')
#os.getcwd()
#os.listdir()
#os.chdir('/home/bioinformatica/Documentos/Documentos/nuevos')
#os.chdir('/home/fernando/Documentos/Epitopes_pipi/redesign_inmuno_hiv')
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
    position_metionina = sequence.seq.find(str('M'),start=0, end=0)
    
        
    return epitope_wild,epitope_scape,hla_scape, position_wild,position_scape,position_asteris,position_metionina



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
                if position_metionina == 0:
                    if position_asteris!= -1:
                        if position_scape < position_asteris:
                            
                            return dict(translated = True)
                        else:
                            
                            return dict(after_orf= True)
                if position_metionina != -1 and position_asteris == -1:
                         
                        return dict(without_stop=True)  
                            
                else: 
                    if position_asteris != -1:
                        if position_scape < position_asteris:
                            
                            return dict(without_met_Pattern_before_stp = True)
                        else:
                            
                            return dict(without_met_Pattern_after_stp= True)
                    else:
                        
                        return dict(without_start_stop = True)
            else:
                return {}        
        else:
            return {} 
    else:
        return {}

"""This function returns a dictionary containing a key that classifies the sequence according 
to the position of the amino acid pattern (type variant wild) with respect to the reading frame, and its value corresponds to a boolean.
"""        
#


def write_wild(position_wild,hla_scape, position_asteris, position_metionina,   HLA_patient_df, epitope_wild, count_pattern_wd ):
    if hla_scape in HLA_patient_df["HLA"].tolist():
        if (position_wild != -1) and (len(HLA_patient_df) != 0):
            if ( epitope_wild not in count_pattern_wd.keys()) or ( epitope_wild in count_pattern_wd.keys() and  (hla_scape not in count_pattern_wd[epitope_wild])):
               if position_metionina == 0 and position_asteris == -1:
                   return dict(without_stop=True)
                
               #print(f'entro en wild { (position_wild != -1) and ( epitope_wild not in count_pattern_wd.keys()) or (epitope_wild in count_pattern_sc.keys() and (hla_scape in count_pattern_sc[epitope_wild] == False)) and (len(HLA_patient_df) != 0)}')                           
               if position_metionina == 0:
                   if position_asteris != -1:
                       if position_wild < position_asteris:
                           return dict(translated= True)
                       else:
                           return dict(after_orf= True)
               
               else:
                   if position_asteris!= -1:
                       if position_wild < position_asteris:
                           return dict(without_met_Pattern_before_stp = True)
                       else:
                           return dict(without_met_Pattern_after_stp = True)
                   else:
                       return dict(without_start_stop = True)
            else:
                return {}    
        else:
            return {}
    else:
        return {}
                      

def dict_epitopes(dictionary, hla, epitope):   
    if epitope not in dictionary.keys():
        dictionary[epitope] = hla
    if epitope  in dictionary.keys() and (hla not in dictionary[epitope]): 
        dictionary[epitope] = hla
    return dictionary                                      

