#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
04_DESCRIPTIVE_PLOTS.PY - Descriptive Statistics and Sample Distribution Plots
===============================================================================
Execution Order: 4 (First plotting script)

Purpose:
    Creates descriptive visualizations showing the distribution of samples,
    models, and cancer types in the dataset. Generates publication-quality
    figures for cohort characterization.
    
    Part 1: Pie charts of cancer type distribution (by sample and by model)
    Part 2: Histograms of PDO/PDX passage distributions
    Part 3: Bar charts of sample matching statistics (PT-PDO, PT-PDO-PDX, etc.)

Dependencies:
    - utility_functions.py
    - constants.py
    - all_sample_matches_reshaped.csv (from 02_PDO_sample_analysis.py)
    - all_sample_matches_table.csv (from 02_PDO_sample_analysis.py)
    - merged_pdx_and_organoid_disc_table.csv (from 03_PDO_PDX_analysis.py)
    - metadata_organoids.csv (from 01_update_metadata.py)

Inputs:
    - all_sample_matches_reshaped.csv
    - all_sample_matches_table.csv
    - merged_pdx_and_organoid_disc_table.csv
    - metadata_organoids.csv

Outputs (plots/ directory):
    - Sample_distribution_by_cancer_type.pdf
    - Model_distribution_by_cancer_type.pdf
    - Organoid_passage_distribution_seaborn_clean.pdf
    - PDX_passage_distribution_all_seaborn_clean.pdf
    - sample_matching_distribution.pdf

Usage:
    python 04_descriptive_plots.py

Author: Linoy
===============================================================================
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from collections import Counter
import matplotlib as mpl
import seaborn as sns

''' This code uses:
    * the all_sample_matches_reshaped table that shows all the 
matches by sample (row = sample = organoid/pdx) + discordance.
    * the all_sample_matches_table
    * the metadata table
    It only includes PDXs collected for this study.''' 


#root = r"/workspaces/Organoid_Stability"
#sys.path.append(os.path.abspath(root))

os.chdir(r"C:\Users\linoy\OneDrive\Desktop\New folder\Organoid_Stability")
root = r"C:\Users\linoy\OneDrive\Desktop\New folder\Organoid_Stability"

from utility_functions import *
from constants import *

all_sample_matches_reshaped = pd.read_csv(os.path.join(root, "tables", "all_sample_matches_reshaped.csv"))
all_samples = pd.read_csv(os.path.join(root, "tables", "all_sample_matches_table.csv"))

#######################################
## Part 1 ##
'''Pie plots showing cancer type distribution by sample (single organoid = sample)
   or by model (single patient = model)'''

sample_df = all_sample_matches_reshaped
organoid_sample_df = sample_df[sample_df['Model'] == 'Organoid']
sample_distribution = organoid_sample_df['Cancer_type'].value_counts()

sample_distribution = sample_distribution.sort_values(ascending=True)

sample_colors = [CUSTOM_COLORS.get(cancer_type, 'grey') for cancer_type in list(sample_distribution.index)]
plt.figure(figsize=(8, 8))
patches, texts, autotexts = plt.pie(sample_distribution, autopct=make_autopct(sample_distribution), startangle=90, labeldistance=1.1,
                                    wedgeprops={"alpha": 0.6}, colors=sample_colors)
plt.title('Distribution of PDO samples by cancer type', y=1.09, fontsize=16)

for i, (percentage, label) in enumerate(zip(sample_distribution, sample_distribution.index)):
    angle = patches[i].theta1 + 0.5 * (patches[i].theta2 - patches[i].theta1)
    radius = 1 * patches[i].r
    x = radius * np.cos(np.deg2rad(angle))
    y = radius * np.sin(np.deg2rad(angle))
    if percentage >= 2:
        rotation = angle + 180 if angle > 90 and angle < 270 else angle
        plt.text(x, y, f'{label}', horizontalalignment='center', verticalalignment='center', rotation=rotation, rotation_mode='anchor')
    else:
        line_length = 1.2 * patches[i].r
        x_line = line_length * np.cos(np.deg2rad(angle))
        y_line = line_length * np.sin(np.deg2rad(angle)) - (i / 15 - 0.17)
        plt.plot([0, x_line], [0, y_line - 0.2], color='black', linestyle='-', linewidth=0.5)
        rotation = 0 if angle < 180 else 180
        plt.text(x_line + 0.1, y_line, f'{label}', horizontalalignment='center', verticalalignment='center', rotation=rotation, rotation_mode='anchor')

total_samples = sample_distribution.sum()
plt.text(0, -1.3, f'Total samples: {total_samples}', horizontalalignment='center', verticalalignment='center', fontsize=12)

plt.gca().tick_params(axis='both', which='both', bottom=True, top=False, left=True, right=False, direction='out', length=6, width=0.5)

plt.savefig(os.path.join(root, "plots", "Sample_distribution_by_cancer_type.pdf"), dpi=300)
plt.show()

model_df = all_samples
model_distribution = model_df['Cancer_type'].value_counts()

model_distribution = model_distribution.sort_values(ascending=True)

model_colors = [CUSTOM_COLORS.get(cancer_type, 'grey') for cancer_type in list(model_distribution.index)]
plt.figure(figsize=(8, 8))
patches, texts, autotexts = plt.pie(model_distribution, autopct=make_autopct(model_distribution), startangle=90, labeldistance=1.1,
                                    wedgeprops={"alpha": 0.6}, colors=model_colors)
plt.title('Distribution of PDO models by cancer type', y=1.09, fontsize=16)

for i, (percentage, label) in enumerate(zip(model_distribution, model_distribution.index)):
    angle = patches[i].theta1 + 0.5 * (patches[i].theta2 - patches[i].theta1)
    radius = 1 * patches[i].r
    x = radius * np.cos(np.deg2rad(angle))
    y = radius * np.sin(np.deg2rad(angle))
    if percentage >= 2:
        rotation = angle + 180 if angle > 90 and angle < 270 else angle
        plt.text(x, y, f'{label}', horizontalalignment='center', verticalalignment='center', rotation=rotation, rotation_mode='anchor')
    else:
        line_length = 1.2 * patches[i].r
        x_line = line_length * np.cos(np.deg2rad(angle))
        y_line = line_length * np.sin(np.deg2rad(angle)) - (i / 15 - 0.17)
        plt.plot([0, x_line], [0, y_line - 0.2], color='black', linestyle='-', linewidth=0.5)
        rotation = 0 if angle < 180 else 180
        plt.text(x_line + 0.1, y_line, f'{label}', horizontalalignment='center', verticalalignment='center', rotation=rotation, rotation_mode='anchor')

total_models = model_distribution.sum()
plt.text(0, -1.2, f'Total models: {total_models}', horizontalalignment='center', verticalalignment='center', fontsize=12)

plt.gca().tick_params(axis='both', which='both', bottom=True, top=False, left=True, right=False, direction='out', length=6, width=0.5)

plt.savefig(os.path.join(root, "plots", "Model_distribution_by_cancer_type.pdf"), dpi=300)
plt.show()

########################
## Part 2 ##
'''Distribution histogram of orgnaoid/PDX sample passages
** currently this doesn't include PDXs collected by me, only Hoge PDXs'''

sns.set_theme(style="white")

merged_df = pd.read_csv(os.path.join(root, "tables", "merged_pdx_and_organoid_disc_table.csv"))
merged_df['Passage'] = pd.to_numeric(merged_df['Passage'], errors='coerce')
all_sample_matches_reshaped['Passage'] = pd.to_numeric(all_sample_matches_reshaped['Passage'], errors='coerce')

#PDO
organoid_data = all_sample_matches_reshaped[all_sample_matches_reshaped['Model'] == 'Organoid']
organoid_data_no_passage = organoid_data[organoid_data['Passage'].isna()]
organoid_data_with_passage = organoid_data.dropna(subset=['Passage'])

organoid_median_passage = organoid_data_with_passage['Passage'].median()
print(f"Median Passage for Organoids: {organoid_median_passage}")

plt.figure(figsize=(12, 8))
ax = sns.countplot(
    data=organoid_data_with_passage, x='Passage', color='#A3D9B1', edgecolor='black', linewidth=0.5)

ax.set_facecolor('white')
plt.gcf().set_facecolor('white')

for spine in ax.spines.values():
    spine.set_linewidth(0.5)

ax.tick_params(axis='x', which='both', bottom=True, top=False, direction='out', length=6, width=0.5)
ax.tick_params(axis='y', which='both', left=True, right=False, direction='out', length=6, width=0.5)

xticks = sorted(organoid_data_with_passage['Passage'].unique())
plt.xticks(ticks=range(len(xticks)), labels=[int(x) for x in xticks], fontsize=18)
plt.yticks(fontsize=18)

plt.title('Distribution of Organoid Passages', fontsize=26)
plt.xlabel('Passage', fontsize=24)
plt.ylabel('Frequency', fontsize=24)

total_nan = len(organoid_data_no_passage)
total_all = len(organoid_data_with_passage) + total_nan
plt.text(0.5, -0.17, f'Total organoid samples with passage: {len(organoid_data_with_passage)}',
         ha='center', va='center', transform=plt.gca().transAxes, fontsize=18)
plt.text(0.5, -0.22, f'Total organoid Samples: {total_all}',
         ha='center', va='center', transform=plt.gca().transAxes, fontsize=18)

plt.tight_layout()
plt.savefig(os.path.join(root, "plots", "Organoid_passage_distribution_seaborn_clean.pdf"), dpi=300)
plt.show()

#PDX
pdx_data_with_passage = merged_df[merged_df['Model'] == 'PDX'].dropna(subset=['Passage'])
pdx_data_no_passage = merged_df[(merged_df['Model'] == 'PDX') & (merged_df['Passage'].isna())]

plt.figure(figsize=(8, 8))
ax = sns.countplot(
    data=pdx_data_with_passage, x='Passage', color='#C9A7DA', edgecolor='black', linewidth=0.5)

ax.set_facecolor('white')
plt.gcf().set_facecolor('white')

for spine in ax.spines.values():
    spine.set_linewidth(0.5)

ax.tick_params(axis='x', which='both', bottom=True, top=False, direction='out', length=6, width=0.5)
ax.tick_params(axis='y', which='both', left=True, right=False, direction='out', length=6, width=0.5)

xticks = sorted(pdx_data_with_passage['Passage'].unique())
plt.xticks(ticks=range(len(xticks)), labels=[int(x) for x in xticks], fontsize=18)
plt.yticks(fontsize=18)

plt.title('Distribution of PDX Passages \n(Unmatched models)', fontsize=22)
plt.xlabel('Passage', fontsize=22)
plt.ylabel('Frequency', fontsize=22)

total_nan = len(pdx_data_no_passage)
total_all = len(pdx_data_with_passage) + total_nan
plt.text(0.5, -0.17, f'Total PDX samples with passage: {len(pdx_data_with_passage)}',
         ha='center', va='center', transform=plt.gca().transAxes, fontsize=18)
plt.text(0.5, -0.22, f'Total PDX Samples: {total_all}',
         ha='center', va='center', transform=plt.gca().transAxes, fontsize=18)

plt.tight_layout()
plt.savefig(os.path.join(root, "plots", "PDX_passage_distribution_all_seaborn_clean.pdf"), dpi=300)
plt.show()

###############################################################
## Part 3 ##
'''Showing all the existing matches in the dataset'''

all_organoid_cols = [col for col in all_samples.columns if col.lower().startswith('organoid')]

pt_pdo_count = all_samples[all_organoid_cols].count().sum()

pt_pdo_pdx_count = all_samples[(all_samples['Tumor'].notna()) & (all_samples['Organoid'].notna()) & (all_samples['PDX'].notna())].shape[0]

df_filtered = all_samples[all_samples['Tumor'].notna()]
organoid_passage_cols = [col for col in df_filtered.columns if col.lower().startswith('organoid_pass')]
organoid_pass_count = df_filtered[organoid_passage_cols].count().sum()

counts = {
    'PT-PDO': pt_pdo_count,
    'PT-PDO-PDX': pt_pdo_pdx_count,
    'Organoid Passages': organoid_pass_count}

fig, ax = plt.subplots(figsize=(6, 7))

max_count = max(counts.values())
min_count = min(counts.values())

cmap = cm.Greens
norm = mpl.colors.Normalize(vmin=min_count, vmax=max_count)

colors = [cmap(norm(c)) for c in counts.values()]

bars = ax.bar(counts.keys(), counts.values(), color=colors, edgecolor='black')

for i, (name, count) in enumerate(counts.items()):
    ax.text(i, count + max_count * 0.01, str(count), ha='center', va='bottom', fontsize=12)

sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

cbar = fig.colorbar(sm, ax=ax, pad=0.02)
cbar.set_label("Count", fontsize=12)

ax.set_ylabel("Samples", fontsize = 14)
ax.set_title("Sample Matching Distribution", fontsize = 16)

ax.tick_params(axis='x', which='both', bottom=True, top=False, direction='out', length=6, width=0.5)
ax.tick_params(axis='y', which='both', left=True, right=False, direction='out', length=6, width=0.5)

plt.tight_layout()
plt.show()

################
''' Same but also including patients and tumor models'''

metadata_filepath = os.path.join(root, "metadata_organoids.csv")
metadata_df = pd.read_csv(metadata_filepath)

summary_row = metadata_df[metadata_df['cohort_name'] == 'Total Summary'].iloc[0]

patient_count = summary_row['Patients']
model_count = summary_row['Tumor']

all_organoid_cols = [col for col in all_samples.columns if col.lower().startswith('organoid')]

pt_pdo_count = all_samples[all_organoid_cols].count().sum()

pt_pdo_pdx_count = all_samples[(all_samples['Tumor'].notna()) & (all_samples['Organoid'].notna()) & (all_samples['PDX'].notna())].shape[0]

df_filtered = all_samples[all_samples['Tumor'].notna()]
organoid_passage_cols = [col for col in df_filtered.columns if col.lower().startswith('organoid_pass')]
organoid_pass_count = df_filtered[organoid_passage_cols].count().sum()

counts = {
    'Patients': patient_count,
    'Models (Tumors)': model_count,
    'Matched PT-PDO\nsamples': pt_pdo_count,
    'Matched PT-PDO-PDX\nsamples': pt_pdo_pdx_count,
    'Organoid passages\nsamples': organoid_pass_count}

fig, ax = plt.subplots(figsize=(8, 7))

max_count = max(counts.values())
min_count = min(counts.values())

cmap = cm.Greens
norm = mpl.colors.Normalize(vmin=min_count, vmax=max_count)

colors = [cmap(norm(c)) for c in counts.values()]

bars = ax.bar(counts.keys(), counts.values(), color=colors, edgecolor='black')

for i, (name, count) in enumerate(counts.items()):
    ax.text(i, count + max_count * 0.01, str(count), ha='center', va='bottom', fontsize=12)

sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

cbar = fig.colorbar(sm, ax=ax, pad=0.02)

ax.set_ylabel("Count", fontsize = 17)

plt.setp(ax.get_xticklabels(), rotation=45, ha='center', fontsize=14)
ax.tick_params(axis='x', which='both', bottom=True, top=False, direction='out', length=6, width=0.5)
ax.tick_params(axis='y', which='both', left=True, right=False, direction='out', length=6, width=0.5, labelsize=14)

for spine in ax.spines.values():
    spine.set_linewidth(0.5)

plt.tight_layout()
plt.savefig(os.path.join(root, "plots","sample_matching_distribution.pdf"), dpi=300)
plt.show()