#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
03_PDO_PDX_ANALYSIS.PY - PDX Sample Analysis and Integration with PDO Data
===============================================================================
Execution Order: 3

Purpose:
    Analyzes PDX (Patient-Derived Xenograft) samples from Hoge et al. 2018 data
    and integrates them with PDO analysis results for comparative studies.
    
    Part 1: Creates tumor-PDX discordance table from Hoge data
    Part 2: Reshapes to sample-level (one row per PDX)
    Part 3: Merges PDX and PDO sample tables for unified analysis

Dependencies:
    - utility_functions.py
    - constants.py
    - all_sample_matches_reshaped.csv (from 02_PDO_sample_analysis.py)

Inputs:
    - PDX_data/1mb_bins_logr_tables/*.csv (Hoge et al. PDX bins data)
    - PDX_data/pairings_files/*.csv (PDX-tumor pairings)
    - all_sample_matches_reshaped.csv

Outputs:
    - Hoge_data_PDX_disc_table.csv 
      (tumor-level: multiple PDXs per tumor, includes PDX1-PDX2 comparisons)
    - Hoge_data_PDX_disc_table_single_PDXs.csv 
      (sample-level: one row per PDX sample)
    - merged_pdx_and_organoid_disc_table.csv 
      (combined PDO+PDX sample-level table for comparative analysis)

Key Metrics:
    - PT_PDX1_gen_disc, PT_PDX2_gen_disc: Tumor-PDX genome discordance
    - PT_PDX1_arms_disc, PT_PDX2_arms_disc: Tumor-PDX arm discordance
    - PDX1_PDX2_gen_disc: Between-passage PDX discordance
    - Passage_difference: Passage gap between PDX1 and PDX2

Usage:
    python 03_PDO_PDX_analysis.py

Author: Linoy
===============================================================================
"""

import os
import sys
import pandas as pd
import numpy as np
import re


#root = r"/workspaces/Organoid_Stability"
#sys.path.append(os.path.abspath(root))

os.chdir(r"C:\Users\linoy\OneDrive\Desktop\New folder\Organoid_Stability")
root = r"C:\Users\linoy\OneDrive\Desktop\New folder\Organoid_Stability"

from utility_functions import *
from constants import *

bins_tables_path = os.path.join(root, "PDX_data", "1mb_bins_logr_tables")
pairings_files_path = os.path.join(root, "PDX_data", "pairings_files")
bins_tables_dict = {}  
pairings_dict = {}  

#####################################
## Part 1 ##
'''Creating a tumor-PDX discordance table using the Hoge et al data from 2018\
and the functions used for this work'''

pdx_disc_table = pd.DataFrame()

for bins_file_name in os.listdir(bins_tables_path):
    cohort_name = clean_cohort_name(bins_file_name)
    pairings_file_name = next((f for f in os.listdir(pairings_files_path) if cohort_name in f), None)

    if pairings_file_name:
        bins_table = pd.read_csv(os.path.join(bins_tables_path, f"{cohort_name}_logr.csv"), header=None, index_col=0)
        
        pairings_file = pd.read_csv(os.path.join(pairings_files_path, f"{cohort_name}_with_pt.csv"))
        pairings_file = pairings_file.dropna(subset=["pt"])
        if "pt" not in pairings_file.columns or ("pair1" not in pairings_file.columns and "pair2" not in pairings_file.columns): continue
    
        for i, row in pairings_file.iterrows():
            tumor = row["pt"] 
            passage_dif = np.nan
            pdx2, pdx2_p = np.nan, np.nan
            
            if row["pair1"] == tumor or "pair1" not in pairings_file.columns:
                if pairings_file["pt"].value_counts()[tumor] > 1: continue
                pdx1 = row["pair2"]
                pdx1_p = get_passage(pdx1)            
            else:
                pdx1 = row["pair1"]
                pdx1_p = get_passage(pdx1)
                pdx2 = row["pair2"]
                pdx2_p = get_passage(pdx2)            
            
            if pdx1_p is not np.nan and pdx2_p is not np.nan:
                passage_dif = pdx2_p - pdx1_p  
                
            if tumor in bins_table.index:
                tumor_bins = bins_table.loc[tumor]
            else:
                tumor_bins = np.nan
            
            if pdx1 in bins_table.index:
                pdx1_bins = bins_table.loc[pdx1]
                tumor_bins_reset = tumor_bins.reset_index(drop=True)
                pdx1_bins_reset = pdx1_bins.reset_index(drop=True)
                pt_pdx1_genome_disc, pt_pdx1_arms_disc = discordance_comparison_bins(tumor_bins_reset, pdx1_bins_reset)
            else:
                pdx1_bins = pt_pdx1_genome_disc = pt_pdx1_arms_disc = np.nan
            
            if pdx2 != np.nan and pdx2 in bins_table.index:
                pdx2_bins = bins_table.loc[pdx2]
                pdx2_bins_reset = pdx2_bins.reset_index(drop=True)
                pt_pdx2_genome_disc, pt_pdx2_arms_disc = discordance_comparison_bins(tumor_bins_reset, pdx2_bins_reset)
                pdx1_pdx2_genome_disc, pdx1_pdx2_arms_disc = discordance_comparison_bins(pdx1_bins_reset, pdx2_bins_reset)
               
                avg_pt_pdx_gen_disc = (pt_pdx1_genome_disc + pt_pdx2_genome_disc) / 2
                avg_pt_pdx_arms_disc = (pt_pdx1_arms_disc + pt_pdx2_arms_disc) / 2
            else:
                pdx2_bins = pt_pdx2_genome_disc = pt_pdx2_arms_disc = np.nan
                pdx1_pdx2_genome_disc = pdx1_pdx2_arms_disc = np.nan                
                avg_pt_pdx_gen_disc = avg_pt_pdx_arms_disc = np.nan
    
            row_data = {
                "Cohort": cohort_name,
                "Tumor": tumor,
                "PDX1": pdx1, "PDX1_passage": pdx1_p,
                "PDX2": pdx2, "PDX2_passage": pdx2_p,
                "PT_PDX1_gen_disc": pt_pdx1_genome_disc,
                "PT_PDX2_gen_disc": pt_pdx2_genome_disc,
                "PDX1_PDX2_gen_disc": pdx1_pdx2_genome_disc,
                "Average_PT_PDX_genome_disc": avg_pt_pdx_gen_disc,
                "PT_PDX1_arms_disc": pt_pdx1_arms_disc,
                "PT_PDX2_arms_disc": pt_pdx2_arms_disc,
                "PDX1_PDX2_arms_disc": pdx1_pdx2_arms_disc,
                "Average_PT_PDX_arms_disc": avg_pt_pdx_arms_disc,
                "Passage_difference": passage_dif}
            pdx_disc_table = pd.concat([pdx_disc_table, pd.DataFrame(row_data, index=[0])])

pdx_disc_table.to_csv(os.path.join(root, "tables","Hoge_data_PDX_disc_table.csv"), index=False)

############################################
## Part 2 ##
'''Reshaping the same table so that each row is a tumor-pdx comparison (each tumor can have multiple rows).'''

reshaped_pdx_disc_table = pd.DataFrame()
for i, row in pdx_disc_table.iterrows():
    cohort_name = row['Cohort']
    tumor = row['Tumor']
    
    cancer_type = cohort_name.split('_')[-1]
    
    pdx1_data = row[['PDX1', 'PDX1_passage', 'PT_PDX1_gen_disc', 'PT_PDX1_arms_disc']]
    pdx1_data = pdx1_data.rename({
        'PDX1': 'PDX', 'PDX1_passage': 'Passage', 
        'PT_PDX1_gen_disc': 'Genome_Discordance', 
        'PT_PDX1_arms_disc': 'Arms_Discordance',})
    
    pdx1_data['Tumor'] = tumor
    pdx1_data['Cohort_Name'] = cohort_name
    pdx1_data['Cancer_type'] = cancer_type
    pdx1_data['Comparison'] = 'PDX1'  
    
    reshaped_pdx_disc_table = pd.concat([reshaped_pdx_disc_table, pdx1_data.to_frame().T], ignore_index=True)
    
    pdx2_bins = row.get('PDX2', np.nan) 
    
    if pd.isna(pdx2_bins):  
        continue
    else:
        pdx2_data = row[['PDX2', 'PDX2_passage', 'PT_PDX2_gen_disc', 'PT_PDX2_arms_disc']]
        pdx2_data = pdx2_data.rename({
            'PDX2': 'PDX', 'PDX2_passage': 'Passage', 
            'PT_PDX2_gen_disc': 'Genome_Discordance', 
            'PT_PDX2_arms_disc': 'Arms_Discordance'})

        pdx2_data['Tumor'] = tumor
        pdx2_data['Cohort_Name'] = cohort_name
        pdx2_data['Cancer_type'] = cancer_type
        pdx2_data['Comparison'] = 'PDX2'  
        
        reshaped_pdx_disc_table = pd.concat([reshaped_pdx_disc_table, pdx2_data.to_frame().T], ignore_index=True)

reshaped_pdx_disc_table['Cancer_type'] = map_cancer_type(reshaped_pdx_disc_table)['Cancer_type']

cols_order = ['Cohort_Name', 'Tumor', 'Cancer_type', 'Comparison', 'PDX', 'Passage', 'Genome_Discordance', 'Arms_Discordance']
reshaped_pdx_disc_table = reshaped_pdx_disc_table[cols_order]

reshaped_pdx_disc_table.to_csv(os.path.join(root, "tables", "Hoge_data_PDX_disc_table_single_PDXs.csv"), index=False)


############################################
## Part 3 ##
'''Merging the per-sample PDX table with the per-sample organoids table which
also includes matched PDXs from the same cohorts'''

organoid_sample_matches = pd.read_csv(os.path.join(root, "tables", "all_sample_matches_reshaped.csv"))
reshaped_pdx_disc_table = pd.read_csv(os.path.join(root, "tables", "Hoge_data_PDX_disc_table_single_PDXs.csv"))


merged_rows = []
for _, row in reshaped_pdx_disc_table.iterrows():
    merged_rows.append({
        'Cohort': row.get('Cohort_Name'),
        'Cancer_type': row.get('Cancer_type'),
        'Tumor': row.get('Tumor'),
        'Sample': row.get('PDX'),
        'Model': 'PDX',
        'Passage': row.get('Passage'),
        'Percent_Genome_Discordance': row.get('Genome_Discordance'),
        'Num_arms_disc': row.get('Arms_Discordance')})

for _, row in organoid_sample_matches.iterrows():
    if pd.notna(row.get('Tumor')):
        merged_rows.append({
            'Cohort': row.get('Cohort_Name'),
            'Cancer_type': row.get('Cancer_type'),
            'Tumor': row.get('Tumor'),
            'Sample': row.get('Sample'),
            'Model': row.get('Model'),
            'Passage': row.get('Passage'),
            'Percent_Genome_Discordance': row.get('Percent_genome_discordance'),
            'Num_arms_disc': row.get('Num_arms_disc')})

merged_pdo_pdx = pd.DataFrame(merged_rows)
merged_pdo_pdx = merged_pdo_pdx[
    ['Cohort', 'Cancer_type', 'Tumor', 'Sample', 'Model', 'Passage', 'Percent_Genome_Discordance', 'Num_arms_disc']]
output_path = os.path.join(root, "tables", "merged_pdx_and_organoid_disc_table.csv")
merged_pdo_pdx.to_csv(output_path, index=False)
