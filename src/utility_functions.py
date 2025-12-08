#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
UTILITY_FUNCTIONS.PY - Helper Functions for Copy Number Analysis
===============================================================================
Purpose:
    Provides core utility functions for copy number discordance calculations,
    sample name processing, cancer type mapping, and visualization helpers.
    
Key Functions:
    - get_disc_list(): Calculate bin-level discordance between two samples
    - perc_genome_disc(): Calculate percentage of genome discordant
    - arms_by_pipe(): Calculate number of arms discordant from pipeline output
    - get_arms_disc(): Calculate arm discordance from bins table
    - discordance_comparison(): Compare two samples (bins + arms)
    - discordance_comparison_bins(): Compare two samples (bins-based only)
    - map_cancer_type(): Standardize cancer type names across cohorts
    - get_passage(): Extract passage number from sample name
    - Other helper functions

Dependencies:
    - constants.py (ARM_COORDINATES, PERCENT_CUTOFF)
    - pandas, numpy, matplotlib, seaborn, re

Inputs:
    - None (imported by other scripts)

Outputs:
    - None (functions called by other scripts)

Usage:
    from utility_functions import perc_genome_disc, map_cancer_type

Author: Linoy
Created: Tue Mar 19 10:21:57 2024
===============================================================================
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re

from constants import *

def get_disc_list(sample_bins1, sample_bins2, discordance_threshold1, discordance_threshold2):
    '''
    Given two bin lists of different samples and two thresholds (top and bottom), returns a list with
    a  discordance value (0,1,NA) for each bin in the two input dictionaries. 
    Values of "0" mean the two samples were not discordant at that location, 
    values of "1" means there was a CN gain from pair1 to pair2, and values of "-1" means
    there was a CN loss from pair1 to pair2. 
    "NA" means there was no data for that location in the two dictionaries.
    '''
    is_discordant_list = [] 
    for i in range(len(sample_bins1)):
        logr_pair1 = sample_bins1[i] 
        logr_pair2 = sample_bins2[i]
        if pd.isna(logr_pair1) or pd.isna(logr_pair2):
            is_discordant_list.append("NA")
        elif (float(logr_pair1) >= 0.3) and ((float(logr_pair2) < discordance_threshold1) or ((float(logr_pair2) < 0.3) and (float(logr_pair1) - float(logr_pair2) >= discordance_threshold2))):
            is_discordant_list.append("1")
        elif (float(logr_pair2) >= 0.3) and ((float(logr_pair1) < discordance_threshold1) or ((float(logr_pair1) < 0.3) and (float(logr_pair2) - float(logr_pair1) >= discordance_threshold2))):
            is_discordant_list.append("1")
        elif (float(logr_pair1) <= -0.3) and ((float(logr_pair2) > -discordance_threshold1) or ((float(logr_pair2) > -0.3) and (float(logr_pair1) - float(logr_pair2) <= -discordance_threshold2))):
            is_discordant_list.append("-1")
        elif (float(logr_pair2) <= -0.3) and ((float(logr_pair1) > -discordance_threshold1) or ((float(logr_pair1) > -0.3) and (float(logr_pair2) - float(logr_pair1) <= -discordance_threshold2))):
            is_discordant_list.append("-1")
        else:
            is_discordant_list.append("0")
    return is_discordant_list

def perc_genome_disc (bins1,bins2): 
    disc_bins = get_disc_list(bins1, bins2, 0.1, 0.3)
    num_bins_disc = disc_bins.count("1") + disc_bins.count("-1")
    perc_genome_disc = num_bins_disc /(num_bins_disc + disc_bins.count("0"))*100 
    return (perc_genome_disc)


def discordance_comparison (sample1, sample2, bins_table, arms_table):
    if sample1 == "NA" or sample2 == "NA":
        return (np.nan,np.nan)
    else:
        sample1_bins = bins_table.loc[sample1]
        sample2_bins = bins_table.loc[sample2]
        sample1_arms = arms_table.loc[sample1]
        sample2_arms = arms_table.loc[sample2]
        perc_genome_disco = perc_genome_disc(sample1_bins, sample2_bins)
        num_arms_disc = arms_by_pipe(sample1_arms, sample2_arms)
        return (perc_genome_disco, num_arms_disc)

def map_cancer_type(df):
    cancer_type_aliases = {
        'Adrenocortical': ['acc'],
        'Bladder': ['bladder cancer', 'urothelial carcinoma', 'bladder', 'blca'],
        'Breast': ['breast cancer', 'brca', 'breast'],
        'Cholangiocarcinoma': ['chol'],
        'Colorectal': ['colorectal cancer', 'colorectal', 'crc', 'coad'],
        'Esophageal': ['esophageal adenocarcinoma', 'esca', 'esophageal'],
        'Glioblastoma': ['gbm','glioblastoma'],
        'Glioma': ['lgg'],
        'Gastric': ['gastric cancer', 'gastric'],
        'Head and Neck': ['headneck', 'head and neck', 'hnsc'],
        'HCC': ['hcc', 'lihc'],
        'Leukemia': ['laml'],
        'Leiomyosarcoma': ['leiomyosarcoma'],
        'Liver': ['liver cancer', 'liver'],
        'Lung': ['luad', 'lung', 'lusc', 'otherlung'],
        'Lymphoma': ['dlbc'],
        'Melanoma': ['melanoma', 'skcm', 'uvm'],
        'Mesothelioma': ['meso'],
        'Multiple Cancers': ['multiple_cancers'],
        'Neuroendocrine': ['neuroendocrine tumor'],
        'Ovarian': ['serous carcinoma of the ovary', 'ov'],
        'Other Tumor': ['othertumor'],
        'Pancreatic': ['pancreas', 'pancreatic', 'pancreatic adenocarcinoma', 'paad'],
        'Prostate': ['prostate cancer', 'prad', 'prostate'],
        'Rectal': ['read'],
        'Renal': ['renal', 'renal cancer', 'kirc', 'kirp', 'kich'],
        'Salivary': ['salivary gland cancer', 'salivary gland','salivary'],
        'Sarcoma': ['sarcoma', 'sarc'],
        'Skin': ['skin', 'skin cancer', 'Epithelial cancer', 'ec'],
        'Stomach': ['stad'],
        'Testicular': ['tgct'],
        'Thymoma': ['thym'],
        'Thyroid': ['thca'],
        'Uterine': ['endometrial adenocarcinoma', 'endometrial cancer', 'uterine', 'ucs'],}

    reverse_map = {}
    for canonical_name, aliases in cancer_type_aliases.items():
        for alias in aliases:
            reverse_map[alias.lower()] = canonical_name

    df['Cancer_type'] = df['Cancer_type'].str.lower()
    df['Cancer_type'] = df['Cancer_type'].apply(lambda x: reverse_map.get(x, x))

    return df


def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct * total / 100.0))
        return '{p:.2f}% ({v:d})'.format(p=pct, v=val) if pct >= 3 else ''
    return my_autopct

def clean_cohort_name(file_name):
    return os.path.splitext(file_name)[0].replace("_logr", "")

def arms_by_pipe (sample1_arms, sample2_arms):
    ''' This function calculates the number of arms discordant between two samples using the arm input
    derived directly from the pipe
    Decide which threshold to use for discordance because currently the output is 0'''
    arms_disc_list = []
    num_arms_no_data = 0 
    
    for i in range (46):
        sample1_arm = sample1_arms.iloc[i]
        sample2_arm = sample2_arms.iloc[i]
        
        if (float(sample1_arm) >= 0.3) and ((float(sample2_arm) < 0.1) or ((float(sample2_arm) < 0.3) and (float(sample1_arm) - float(sample2_arm) >= 0.3))):
            arms_disc_list.append("-1")
        elif (float(sample2_arm) >= 0.3) and ((float(sample1_arm) < 0.1) or ((float(sample1_arm) < 0.3) and (float(sample2_arm) - float(sample1_arm) >= 0.3))):
            arms_disc_list.append("1")
        elif (float(sample1_arm) <= -0.3) and ((float(sample2_arm) > -0.1) or ((float(sample2_arm) > -0.3) and (float(sample1_arm) - float(sample2_arm) <= -0.3))):
            arms_disc_list.append("1")
        elif (float(sample2_arm) <= -0.3) and ((float(sample1_arm) > -0.1) or ((float(sample1_arm) > -0.3) and (float(sample2_arm) - float(sample1_arm) <= -0.3))):
            arms_disc_list.append("-1")

        else:
            arms_disc_list.append("0")
         
    num_gains = arms_disc_list.count("1")
    num_losses = arms_disc_list.count("-1")
    num_neutral = arms_disc_list.count("0")
    num_arms_disc = num_gains + num_losses
    return num_arms_disc

def get_passage(sample_name):
    """
    Given a string sample_name, returns an integer with the passage information
    contained in the sample name (PT, P0, P1, etc...).
    """
    regex_results = re.match(r'([A-Z0-9a-z_-]+).(P[T0-9]+)', sample_name)
    if regex_results:
        passage = regex_results.groups()[1]  
        return int(passage[1:])  
    return np.nan  

def get_full_sample_name_from_bins(sample_name, bins_table):
    sample_name = sample_name.strip()
    if '_T1' in sample_name or '_O1' in sample_name:
        return sample_name
    
    matching_samples = [index for index in bins_table.index if sample_name in index]
    if matching_samples:
        return matching_samples[0]
    
    
def needs_correction(bins_table, arms_table, matches_table):
    '''Check if sample name correction in the arms table is needed based on 
    samples that exist in bins_table but not in arms_table'''
    bins_sample_names = bins_table.index
    tumor_samples = matches_table['Tumor'].unique()
    organoid_samples = [row for col in matches_table.columns if col.startswith('Organoid') for row in matches_table[col].unique()]
    all_samples_to_check = set(tumor_samples) | set(organoid_samples)
    
    missing_in_arms = [sample for sample in all_samples_to_check if sample in bins_sample_names and sample not in arms_table.index]
    if missing_in_arms:
        return True
    return False

def get_arms_disc(sample1_bins, sample2_bins, disc_threshold1, disc_threshold2):
    '''
    Calculates the arm discordance based on bins table and not arms table.
    Calculates the arm discordance per 2 samples and not per cohort.
    Given two bin lists of different samples and two thresholds (top and bottom), returns the number of arms discordant 
    and the percent genome discordant between the two samples
    Values of "0" mean the two samples were not discordant at that location, 
    "1" means there was a CN gain from pair1 to pair2, and values of "-1" means there was a CN loss from pair1 to pair2. 
    "NA" means there was no data for that location in the two dictionaries.
    '''
    num_arms_discordant = 0
    num_arms_no_data = 0 

    for i in range(len(ARM_COORDINATES) - 1):
        start = ARM_COORDINATES[i]
        stop = ARM_COORDINATES[i + 1]

        chr_arm = get_disc_list(sample1_bins, sample2_bins, 0.1, 0.3)[start - 1:stop - 1]
        num_bins_chr_arm = len(chr_arm)
        num_gains = chr_arm.count("1")
        num_losses = chr_arm.count("-1")
        num_neutral = chr_arm.count("0")
        
        if num_gains + num_losses + num_neutral == 0:
            num_arms_no_data += 1
        
        if (num_gains / num_bins_chr_arm) >= (PERCENT_CUTOFF / 100) or (num_losses / num_bins_chr_arm) >= (PERCENT_CUTOFF / 100):
            num_arms_discordant += 1

    return num_arms_discordant 


def discordance_comparison_bins(sample1_bins, sample2_bins):
    perc_genome_disco = perc_genome_disc(sample1_bins, sample2_bins)
    num_arms_disc = get_arms_disc(sample1_bins, sample2_bins, 0.1, 0.3)
    return perc_genome_disco, num_arms_disc