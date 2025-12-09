#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
01_UPDATE_METADATA.PY - Update Metadata Summary Statistics
===============================================================================
Execution Order: 1 

Purpose:
    Updates the metadata_organoids.csv file with current statistics from all
    cohorts including patient counts, tumor counts, organoid counts, PDX counts,
    and passage information. Generates a summary row with totals across all cohorts.

Dependencies:
    - utility_functions.py
    - constants.py

Inputs:
    - metadata_organoids.csv (existing metadata file)
    - {cohort_index}_{cohort_name}_matches_table.csv files from each cohort folder

Outputs:
    - metadata_organoids.csv (updated with current statistics)

Key Statistics Calculated:
    - Number of unique patients per cohort
    - Number of tumor samples
    - Number of organoid samples (including passages)
    - Number of PDX, 2D, PDXO, spheroid, and normal tissue samples
    - Counts of organoid passages and multiple tumor sections

Usage:
    python 01_update_metadata.py

Author: Linoy Raz
Contributors: Haia Khoury 
===============================================================================
"""

import os
import sys
import pandas as pd
import numpy as np
import re

# Auto-detect repository root (2 levels up from this script)
# Script location: repo_root/src/postprocessing/01_update_metadata.py
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(script_dir))

# Add postprocessing directory to path for imports
sys.path.insert(0, script_dir)

from utility_functions import *
from constants import *

# Set root to repository root for data access
root = repo_root

metadata_filepath = os.path.join(root, "data", "metadata_organoids.csv")
metadata = pd.read_csv(metadata_filepath, header=0)

summary_index = metadata.index[-1]
last_row_cohort_name = metadata.loc[summary_index, 'cohort_name']

cohort_metadata = metadata[metadata['cohort_name'] != last_row_cohort_name].copy()

models = ["Normal_tissue", "PDX", "2D", "PDXO", "spheroids"]
additionals = ["organoid_passages", "multiple_tumor_sections", "multiple_pdxs"]

PATIENT_ID_COL = 'patient_id'

cohort_metadata.rename(columns={"Tumor": "Patients"}, inplace=True)
cohort_metadata.insert(cohort_metadata.columns.get_loc("Patients") + 1, "Tumor", 0)

for idx, row in cohort_metadata.iterrows():
    cohort_index = row["cohort_index"]
    cohort_name = row["cohort_name"]

    cohort_folder = os.path.join(root, "results", "Cohort_data", f"{int(cohort_index)}_{cohort_name}")
    matches_table_path = os.path.join(cohort_folder, f"{int(cohort_index)}_{cohort_name}_matches_table.csv")

    matches_table = pd.read_csv(matches_table_path)

    if PATIENT_ID_COL in matches_table.columns:
        unique_count_column = PATIENT_ID_COL
    else:
        unique_count_column = "Tumor"

    id_series = matches_table[unique_count_column].dropna().astype(str)
    id_series = id_series.str.replace(r'[^\x20-\x7E]+', ' ', regex=True)
    id_series = id_series.str.strip().str.lower()
    id_series = id_series[id_series != '']
    cohort_metadata.loc[idx, "Patients"] = len(id_series.unique())

    cohort_metadata.loc[idx, "Tumor"] = len(matches_table)

    organoid_cols = [col for col in matches_table.columns if 'organoid' in col.lower()]
    if organoid_cols:
        total_organoid_count = matches_table[organoid_cols].notna().sum().sum()
    else:
        total_organoid_count = 0
    cohort_metadata.loc[idx, "Organoid"] = total_organoid_count

    for model in models:
        model_column = [col for col in matches_table.columns if col.lower() == model.lower()]
        cohort_metadata.loc[idx, model] = matches_table[model_column[0]].count() if model_column else 0

    for add_col in additionals:
        prefix = add_col.split("_")[0].lower()
        if add_col == "organoid_passages":
            passage_cols = [col for col in matches_table.columns if col.lower().startswith("organoid_pass")]
            cohort_metadata.loc[idx, add_col] = matches_table[passage_cols].count().sum() if passage_cols else 0
        else:
            model_columns = [col for col in matches_table.columns if col.lower().startswith(prefix) and col.lower() != prefix]
            cohort_metadata.loc[idx, add_col] = matches_table[model_columns].count().sum() if model_columns else 0

summary_row = cohort_metadata.select_dtypes(include=[np.number]).sum()
summary_row["cohort_name"] = "Total Summary"
summary_row["cohort_index"] = 999
summary_row["Cancer_type"] = "All cancers"

metadata.update(cohort_metadata)
metadata = pd.concat([cohort_metadata, pd.DataFrame([summary_row])], ignore_index=True)

metadata.to_csv(os.path.join(root, "data", "metadata_organoids.csv"), index=False)
