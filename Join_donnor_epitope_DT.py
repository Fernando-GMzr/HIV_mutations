#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 09:45:27 2022

@author: fernando
This script combine dataset information of donnor-HLA with dataset information
"""

import pandas as pd
import csv
import os



directorio = os.getcwd()
os.chdir('/home/fernando/Documentos/Epitopes_pipi/redesign_inmuno_hiv')
os.listdir('/home/fernando/Documentos/Epitopes_pipi/redesign_inmuno_hiv')
Donnor_dataset = open('input_HLA_donnor.csv')
epitope_WT = open('epitopeWT.csv')
epitope_VR = open('epitopeVARIANT.csv')
Donnor = pd.read_csv(Donnor_dataset)
EP_VR= pd.read_csv(epitope_VR)
EP_WT= pd.read_csv(epitope_WT)
Donnor.head()

Donnor_ep_VR=pd.merge(Donnor,EP_VR, on='HLA', how='left')
Donnor_ep_WT=pd.merge(Donnor,EP_WT, on='HLA', how='left')

Donnor_VR_WT = pd.merge(Donnor_ep_VR, Donnor_ep_WT, on=('Patient_ID','HLA', 'Protein'), how='left')
Donnor_VR_WT.to_csv('fusion_donnor_epitop2.csv', index=False)
