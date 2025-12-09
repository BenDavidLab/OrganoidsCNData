#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
05_PDO_PLOTS.PY - PDO-Specific Analysis and Visualization
===============================================================================
Execution Order: 5

Purpose:
    Creates comprehensive visualizations for PDO (Patient-Derived Organoid)
    discordance analysis including cancer type comparisons, passage effects,
    and temporal stability assessments.
    
    Part 1: Boxplots of genome/arm discordance by cancer type (sample & model analysis)
    Part 2: Scatter plots of discordance vs passage number
    Part 3: Barplots comparing genome vs arm discordance by cohort
    Part 4: Waterfall and line plots showing temporal changes (earliest vs latest)

Dependencies:
    - utility_functions.py
    - constants.py
    - all_sample_matches_reshaped.csv (from 02_PDO_sample_analysis.py)
    - all_sample_matches_reshaped_average.csv (from 02_PDO_sample_analysis.py)
    - all_sample_matches_table.csv (from 02_PDO_sample_analysis.py)

Inputs:
    - all_sample_matches_reshaped.csv (sample-level data)
    - all_sample_matches_reshaped_average.csv (model-level averages)
    - all_sample_matches_table.csv (patient-level data)

Outputs (plots/ directory):
    - PT-PDO_samples_genome_disc_by_cancer_type_sns.pdf
    - PT-PDO_samples_arms_disc_by_cancer_type.pdf
    - PT-PDO_models_genome_disc_by_cancer_type.pdf
    - PT-PDO_models_arms_disc_by_cancer_type.pdf
    - PT-PDO_genome_disc_vs_passage.pdf
    - PT-PDO_arms_disc_vs_passage.pdf
    - Genome_vs_arm_disc_by_cohort.pdf
    - waterfall_genome.pdf
    - waterfall_arms.pdf
    - PT_PDO_genome_discordance_early_vs_latest.pdf
    - PT_PDO_arm_discordance_early_vs_latest.pdf

Statistical Tests:
    - One-sample t-test for temporal changes (waterfall plots)
    - Comparison of earliest vs latest passages

Usage:
    python 05_PDO_plots.py

Author: Linoy Raz
Contributors: Haia Khoury 
===============================================================================
"""

import os
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.stats import mannwhitneyu, ttest_1samp
from statannotations.Annotator import Annotator
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as mtick
from Code.utility_functions import *
from Code.constants import *

''' This code uses the all_sample_matches_reshaped table that shows all the organoid
matches we have by organoid (organoid = row), adds genome discordance and creates
different boxplots '''

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) #Assumes scripts are in project root

all_sample_matches_reshaped = pd.read_csv(os.path.join(root, "tables", "all_sample_matches_reshaped.csv"))
all_sample_matches_reshaped_average = pd.read_csv(os.path.join(root, "tables", "all_sample_matches_reshaped_average.csv"))
all_sample_matches_table = pd.read_csv(os.path.join(root, "tables", "all_sample_matches_table.csv"))

all_sample_matches_reshaped = map_cancer_type(all_sample_matches_reshaped)


##################################################
## Part 1 ##
''' Tumor-organoid %genome and arm discordance by cancer type (boxplots).
    * Each tumor-organoid pair is a data point (= each organoid=sample is
    a separate data point.
    * The plot includes multiple organoid from the same tumor
    * Only cohorts with >3 organoids are included'''

## Sample analysis ##
sns.set_theme(style="whitegrid")

filtered_data = all_sample_matches_reshaped[
    (all_sample_matches_reshaped['Model'] == "Organoid") &
    (~pd.isna(all_sample_matches_reshaped['Tumor']))]

def plot_discordance_sns(data, metric, title, ylabel, colors, filename):
    grouped = data.groupby('Cancer_type')
    medians = grouped[metric].median().sort_values()
    counts = grouped[metric].count().loc[medians.index]

    mean_of_means = grouped[metric].mean().mean()
    median_of_medians = grouped[metric].median().median()

    filtered_colors = {k: v for k, v in colors.items() if k in medians.index}

    plt.figure(figsize=(12, 8))
    sns.set_style("white")
    ax = plt.gca() 

    sns.boxplot(
        x='Cancer_type', y=metric, data=data, order=medians.index,
        palette=filtered_colors, showfliers=False, linewidth=0.5, ax=ax)
    
    sns.stripplot(
        x='Cancer_type', y=metric, data=data, order=medians.index,
        color='grey', size=6, jitter=True, dodge=True, ax=ax)

    plt.xlabel('Cancer type', fontsize=18, labelpad=0)
    plt.ylabel(ylabel, fontsize=18, labelpad=10)
    plt.title(title, fontsize=24, pad=15)
    if metric == 'Percent_genome_discordance':
        ax.yaxis.set_major_formatter(mtick.PercentFormatter())

    xticks = plt.xticks(
        ticks=range(len(counts)),
        labels=[f"{cancer}\nN={n}" for cancer, n in zip(counts.index, counts)],
        #labels=[f"{cancer}" for cancer, n in zip(counts.index, counts)],
        rotation=45,
        ha='center',
        fontsize=18)

    plt.yticks(fontsize=18)
    
    plt.tick_params(axis='x', which='both', bottom=True, top=False, length=6, width=0.5, direction='out')
    plt.tick_params(axis='y', which='both', left=True, right=False, length=6, width=0.5, direction='out')
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)

    for i, median_val in enumerate(medians):
        plt.text(
            i,
            median_val + 0.35,
            f"{median_val:.2f}" + ('%' if metric == 'Percent_genome_discordance' else ''),
            horizontalalignment='center',
            color='black',
            fontsize=16)

    y_max = data[metric].max()
    plt.text(-0.35, y_max * 0.93, 
        f"Mean of means: {mean_of_means:.2f}" + ('%' if metric == 'Percent_genome_discordance' else ''),
        fontsize=18, color='black', ha='left')
    plt.text(-0.35, y_max * 0.99, 
        f"Median of medians: {median_of_medians:.2f}" + ('%' if metric == 'Percent_genome_discordance' else ''),
        fontsize=18, color='black', ha='left')
    
    handles = [plt.Line2D([0], [0], marker='o', color='w', label=ct, 
                         markerfacecolor=filtered_colors[ct], markersize=10)
                for ct in filtered_colors]
    plt.legend(handles=handles, title='Cancer type', bbox_to_anchor=(1, 0.5),
               loc='center left', frameon=False, fontsize=14, title_fontsize=14)


    plt.tight_layout()
    plt.savefig(os.path.join(root, "plots", f"{filename}"), dpi=300)
    plt.show()

filtered_data_genome = filtered_data.groupby('Cancer_type').filter(
    lambda x: x['Percent_genome_discordance'].notna().sum() >= 3)
filtered_data_arms = filtered_data.groupby('Cancer_type').filter(
    lambda x: x['Num_arms_disc'].notna().sum() >= 3)

plot_discordance_sns(
    data=filtered_data_genome,
    metric='Percent_genome_discordance',
    title='PT-PDO genome discordance by cancer type (N>=3)\n(Sample analysis)',
    ylabel='Genome discordance',
    colors={k: v for k, v in CUSTOM_COLORS.items() if k in filtered_data_genome['Cancer_type'].unique()},
    filename='PT-PDO_samples_genome_disc_by_cancer_type_sns.pdf')

plot_discordance_sns(
    data=filtered_data_arms,
    metric='Num_arms_disc',
    title='PT-PDO arm discordance by cancer type (N>=3)\n(Sample analysis)',
    ylabel='Number of arms discordant',
    colors={k: v for k, v in CUSTOM_COLORS.items() if k in filtered_data_arms['Cancer_type'].unique()},
    filename='PT-PDO_samples_arms_disc_by_cancer_type.pdf')

######################
## Model analysis ##

sns.set_theme(style="whitegrid")
filtered_data = all_sample_matches_reshaped_average[
    (all_sample_matches_reshaped_average['Model'] == "Organoid") &
    (~pd.isna(all_sample_matches_reshaped_average['Tumor']))]

filtered_data = filtered_data.groupby('Cancer_type').filter(lambda x: len(x) >= 3)
grouped_data = filtered_data.groupby('Cancer_type')

median_genome_discordance = grouped_data['Percent_genome_discordance'].median()
mean_of_means = grouped_data['Percent_genome_discordance'].mean().mean()
median_of_medians = median_genome_discordance.median()
cancer_order_genome = median_genome_discordance.sort_values().index

median_arm_level_discordance = grouped_data['Num_arms_disc'].median()
mean_of_means_arm = grouped_data['Num_arms_disc'].mean().mean()
median_of_medians_arm = median_arm_level_discordance.median()
cancer_order_arm = median_arm_level_discordance.sort_values().index

sns.set_theme(style="white")
plt.rcParams.update({
    'axes.edgecolor': 'black',
    'axes.linewidth': 0.5,
    'axes.titlesize': 22,
    'axes.labelsize': 20,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'figure.titlesize': 24})

# Genome discordance
g1 = sns.catplot(
    data=filtered_data,
    kind="box",
    x="Cancer_type",
    y="Percent_genome_discordance",
    order=cancer_order_genome,
    palette=CUSTOM_COLORS,
    height=8,
    aspect=1.5,
    linewidth=0.5, 
    showfliers=False)

sns.stripplot(
    data=filtered_data,
    x="Cancer_type",
    y="Percent_genome_discordance",
    order=cancer_order_genome,
    color='grey',
    size=5,
    jitter=True,
    ax=g1.ax)

g1.set_axis_labels("Cancer type", "Genome discordance")
g1.ax.set_title("PT-PDO genome discordance by cancer type (N≥3)\n(Model analysis)", fontsize=22)

sample_counts = grouped_data.size().loc[cancer_order_genome]
g1.set_xticklabels([f"{ct}\nN={n}" for ct, n in zip(sample_counts.index, sample_counts)], rotation=45, ha='center')
#g1.set_xticklabels([f"{ct}" for ct, n in zip(sample_counts.index, sample_counts)], rotation=45, ha='center')

for i, median in enumerate(median_genome_discordance.loc[cancer_order_genome]):
    g1.ax.text(i, median + 0.5, f"{median:.2f}%", ha='center', color='black', fontsize=16)

g1.ax.text(-0.2, filtered_data['Percent_genome_discordance'].max() - 4, f"Mean of means: {mean_of_means:.2f}%", fontsize=18)
g1.ax.text(-0.2, filtered_data['Percent_genome_discordance'].max(), f"Median of medians: {median_of_medians:.2f}%", fontsize=18)
g1.ax.yaxis.set_major_formatter(mtick.PercentFormatter())

for spine in g1.ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

g1.ax.tick_params(axis='both', which='both', direction='out', length=6, width=0.5,bottom=True,
    left=True, top=False, right=False)

handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=CUSTOM_COLORS[ct], markersize=10) for ct in cancer_order_genome]
g1.ax.legend(handles, cancer_order_genome, title="Cancer type", bbox_to_anchor=(1.01, 1), loc='upper left', frameon=False)

g1.fig.tight_layout()
g1.savefig(os.path.join(root, "plots", "PT-PDO_models_genome_disc_by_cancer_type.pdf"), dpi=300)
plt.show()

# Arms discordance
g2 = sns.catplot(data=filtered_data, kind="box", x="Cancer_type", y="Num_arms_disc",
    order=cancer_order_arm, palette=CUSTOM_COLORS, height=8, aspect=1.5, linewidth=0.5,
    showfliers=False)

sns.stripplot(data=filtered_data, x="Cancer_type", y="Num_arms_disc", order=cancer_order_arm,
    color='grey', size=5, jitter=True, ax=g2.ax)

g2.set_axis_labels("Cancer type", "Number of arms discordant")
g2.ax.set_title("PT-PDO number of arms discordant by cancer type (N≥3)\n(Model analysis)", fontsize=22)

sample_counts_arm = grouped_data.size().loc[cancer_order_arm]
g2.set_xticklabels([f"{ct}\nN={n}" for ct, n in zip(sample_counts_arm.index, sample_counts_arm)], rotation=45, ha='right')
#g2.set_xticklabels([f"{ct}" for ct, n in zip(sample_counts_arm.index, sample_counts_arm)], rotation=45, ha='right')

for i, median_arm in enumerate(median_arm_level_discordance.loc[cancer_order_arm]):
    g2.ax.text(i, median_arm + 0.2, f"{median_arm:.1f}", ha='center', color='black', fontsize=16)


g2.ax.text(-0.2, filtered_data['Num_arms_disc'].max() - 1.5, f"Mean of means: {mean_of_means_arm:.2f}", fontsize=18)
g2.ax.text(-0.2, filtered_data['Num_arms_disc'].max(), f"Median of medians: {median_of_medians_arm:.2f}", fontsize=18)

for spine in g2.ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

g2.ax.tick_params(axis='both', which='both', direction='out', length=6, width=0.5,bottom=True,
    left=True, top=False, right=False)

handles_arm = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=CUSTOM_COLORS[ct], markersize=10) for ct in cancer_order_arm]
g2.ax.legend(handles_arm, cancer_order_arm, title="Cancer type", bbox_to_anchor=(1.01, 1), loc='upper left', frameon=False)

g2.fig.tight_layout()
g2.savefig(os.path.join(root, "plots", "PT-PDO_models_arms_disc_by_cancer_type.pdf"), dpi=300)
plt.show()

###########

print("\nP-values for Organoid vs PDX comparisons by cancer type - sample analysis) (one-sided t-test, Organoid < PDX):")
print("="*70)

for cancer_type in sorted_cohorts:
    group_data = df_filtered[df_filtered['Cancer_type'] == cancer_type]

    org_gen = group_data[group_data['Model'] == 'Organoid']['Percent_genome_discordance']
    pdx_gen = group_data[group_data['Model'] == 'PDX']['Percent_genome_discordance']

    org_arm = group_data[group_data['Model'] == 'Organoid']['Num_arms_disc']
    pdx_arm = group_data[group_data['Model'] == 'PDX']['Num_arms_disc']

    _, p_gen = ttest_ind(org_gen, pdx_gen, alternative='less')
    _, p_arm = ttest_ind(org_arm, pdx_arm, alternative='less')

    print(f"{cancer_type}:")
    print(f"  Genome discordance p-value: {p_gen:.3f}")
    print(f"  Arm discordance p-value:    {p_arm:.3f}")
    print("-" * 70)


##############################################
## Part 2 ##
''' Tumor-organoid %genome and discordance by cancer type (boxplots).
    * Each tumor-organoid pair is a data point (= each organoid=sample is
    a separate data point.
    * The plot includes multiple organoid from the same tumor
    * Only cohorts with >3 organoids and >3 PDXs
    * One sided t-test'''
    
cohort_counts = all_sample_matches_reshaped.groupby(['Cancer_type', 'Model']).size().unstack(fill_value=0)
cohorts_of_interest = cohort_counts[(cohort_counts['PDX'] >= 3) & (cohort_counts['Organoid'] >= 3)].index.tolist()
filtered_df = all_sample_matches_reshaped[(all_sample_matches_reshaped['Cancer_type'].isin(cohorts_of_interest)) & (~all_sample_matches_reshaped['Tumor'].isna())]

median_discordance_organoids = filtered_df[filtered_df['Model'] == 'Organoid'].groupby('Cancer_type')['Percent_genome_discordance'].median()
median_discordance_pdx = filtered_df[filtered_df['Model'] == 'PDX'].groupby('Cancer_type')['Percent_genome_discordance'].median()
median_arm_discordance_organoids = filtered_df[filtered_df['Model'] == 'Organoid'].groupby('Cancer_type')['Num_arms_disc'].median()
median_arm_discordance_pdx = filtered_df[filtered_df['Model'] == 'PDX'].groupby('Cancer_type')['Num_arms_disc'].median()
sorted_cancer_types = median_discordance_organoids.sort_values().index

sns.set_style("white")

plt.figure(figsize=(8, 8))
sns.boxplot(data=filtered_df, x='Cancer_type', y='Percent_genome_discordance', showfliers=False, hue='Model', order=sorted_cancer_types, palette={"Organoid": "#A3D9B1", "PDX": "#C9A7DA"}, linewidth=0.5)
sns.stripplot(data=filtered_df, x='Cancer_type', y='Percent_genome_discordance', hue='Model', dodge=True, alpha=0.5, order=sorted_cancer_types, color="grey", jitter=True)
plt.title('Tumor-model percent genome discordance\n(Matched models)', fontsize=20)
plt.xlabel('Cancer type', fontsize=18, labelpad=35)
plt.ylabel('Genome discordance', fontsize=18)
plt.legend([], [], frameon=False)

organoid_mean_of_means = filtered_df[filtered_df['Model'] == 'Organoid']['Percent_genome_discordance'].mean()
pdx_mean_of_means = filtered_df[filtered_df['Model'] == 'PDX']['Percent_genome_discordance'].mean()
organoid_median_of_medians = filtered_df[filtered_df['Model'] == 'Organoid'].groupby('Cancer_type')['Percent_genome_discordance'].median().median()
pdx_median_of_medians = filtered_df[filtered_df['Model'] == 'PDX'].groupby('Cancer_type')['Percent_genome_discordance'].median().median()

y_max = filtered_df['Percent_genome_discordance'].max()
plt.text(-0.5, y_max * 0.98, f'Median of medians (PDOs): {organoid_median_of_medians:.2f}%\nMedian of Medians (PDXs): {pdx_median_of_medians:.2f}%', color='black', fontsize=14, ha='left', va='top')
plt.text(-0.5, y_max * 1.08, f'Mean of means (PDOs): {organoid_mean_of_means:.2f}%\nMean of Means (PDXs): {pdx_mean_of_means:.2f}%', color='black', fontsize=14, ha='left', va='top')

pvalue_dict = {}

for i, cancer_type in enumerate(sorted_cancer_types):
    org_data = filtered_df[(filtered_df['Cancer_type'] == cancer_type) & (filtered_df['Model'] == 'Organoid')]['Percent_genome_discordance']
    pdx_data = filtered_df[(filtered_df['Cancer_type'] == cancer_type) & (filtered_df['Model'] == 'PDX')]['Percent_genome_discordance']
    t_stat, p_value = stats.ttest_ind(org_data, pdx_data, nan_policy='omit', alternative='less')
    pvalue_dict[cancer_type] = p_value

    if p_value < 0.001:
        stars = '***'
    elif p_value < 0.01:
        stars = '**'
    elif p_value < 0.05:
        stars = '*'
    else:
        stars = 'ns'

    y_max_local = filtered_df[filtered_df['Cancer_type'] == cancer_type]['Percent_genome_discordance'].max()
    y_pos = y_max_local * 1.05
    plt.plot([i - 0.2, i + 0.2], [y_pos, y_pos], color='black', lw=0.5)
    plt.text(i, y_pos + 0.03, stars, ha='center', va='bottom', fontsize=12, color='black')

ax = plt.gca()
ax.xaxis.set_tick_params(which='both', bottom=True, top=False, direction='out', length=6, width=0.5)
ax.set_ylim(bottom=-1)
ax.yaxis.tick_left()
ax.set_yticklabels([f"{int(t)}%" for t in ax.get_yticks()])

ax.tick_params(axis='y', direction='out', length=6, width=0.5, labelsize=14)
ax.tick_params(axis='x', direction='out', length=6, width=0.5, labelsize=14) 

for spine in ax.spines.values():
    spine.set_linewidth(0.5)

new_labels = []
for cancer_type in sorted_cancer_types:
    n_orgs = filtered_df[(filtered_df['Cancer_type'] == cancer_type) & (filtered_df['Model'] == 'Organoid')].shape[0]
    n_pdxs = filtered_df[(filtered_df['Cancer_type'] == cancer_type) & (filtered_df['Model'] == 'PDX')].shape[0]
    label = f"{cancer_type}\nN PDOs={n_orgs}\nN PDXs={n_pdxs}"
    label = f"{cancer_type}"
    new_labels.append(label)

plt.xticks(ticks=range(len(sorted_cancer_types)), labels=new_labels, ha='center', fontsize=14)
plt.xlabel('Cancer type', fontsize=18, labelpad=35)
plt.ylabel('Genome discordance', fontsize=18)
plt.title('Tumor-model genome discordance\n(Matched models)', fontsize=20)
plt.tight_layout(rect=[0, 0, 1, 1])
plt.savefig(os.path.join(root, "plots", "Matched_TO_TP_perc_genome_disc_by_cancer_type.pdf"), dpi=300)
plt.show()

for cancer_type, pval in pvalue_dict.items():
    print(f"Percent genome discordance comparison for {cancer_type}: p = {pval:.5f}")

plt.figure(figsize=(8, 8))
sns.boxplot(data=filtered_df, x='Cancer_type', y='Num_arms_disc', showfliers=False, hue='Model', order=sorted_cancer_types, palette={"Organoid": "#A3D9B1", "PDX": "#C9A7DA"}, linewidth=0.5)
sns.stripplot(data=filtered_df, x='Cancer_type', y='Num_arms_disc', hue='Model', dodge=True, alpha=0.5, order=sorted_cancer_types, color="darkgrey", jitter=True)
plt.title('Tumor-model arm discordance\n(Matched models)', fontsize=20)
plt.xlabel('Cancer type', fontsize=20, labelpad=10)
plt.ylabel('Number of arms discordant', fontsize=20)
plt.legend([], [], frameon=False)

organoid_mean_arms = filtered_df[filtered_df['Model'] == 'Organoid']['Num_arms_disc'].mean()
pdx_mean_arms = filtered_df[filtered_df['Model'] == 'PDX']['Num_arms_disc'].mean()
organoid_median_arms = filtered_df[filtered_df['Model'] == 'Organoid'].groupby('Cancer_type')['Num_arms_disc'].median().median()
pdx_median_arms = filtered_df[filtered_df['Model'] == 'PDX'].groupby('Cancer_type')['Num_arms_disc'].median().median()

y_max_arms = filtered_df['Num_arms_disc'].max()
plt.text(-0.5, y_max_arms * 0.98, f'Median of medians (PDOs): {organoid_median_arms:.2f}\nMedian of medians (PDXs): {pdx_median_arms:.2f}', color='black', fontsize=14, ha='left', va='top')
plt.text(-0.5, y_max_arms * 1.08, f'Mean of means (PDOs): {organoid_mean_arms:.2f}\nMean of means (PDXs): {pdx_mean_arms:.2f}', color='black', fontsize=14, ha='left', va='top')

pvalue_dict_arms = {}

for i, cancer_type in enumerate(sorted_cancer_types):
    org_data = filtered_df[(filtered_df['Cancer_type'] == cancer_type) & (filtered_df['Model'] == 'Organoid')]['Num_arms_disc']
    pdx_data = filtered_df[(filtered_df['Cancer_type'] == cancer_type) & (filtered_df['Model'] == 'PDX')]['Num_arms_disc']
    t_stat, p_value = stats.ttest_ind(org_data, pdx_data, nan_policy='omit', alternative='less')
    pvalue_dict_arms[cancer_type] = p_value

    if p_value < 0.001:
        stars = '***'
    elif p_value < 0.01:
        stars = '**'
    elif p_value < 0.05:
        stars = '*'
    else:
        stars = 'ns'

    y_max_local = filtered_df[filtered_df['Cancer_type'] == cancer_type]['Num_arms_disc'].max()
    y_pos = y_max_local * 1.05
    plt.plot([i - 0.2, i + 0.2], [y_pos, y_pos], color='black', lw=0.5)
    plt.text(i, y_pos + 0.03, stars, ha='center', va='bottom', fontsize=12, color='black')

ax = plt.gca()
ax.set_ylim(bottom=-1)
ax.yaxis.tick_left()

ax.tick_params(axis='y', direction='out', length=6, width=0.5, labelsize=14)
ax.tick_params(axis='x', direction='out', length=6, width=0.5, labelsize=14)
ax.xaxis.set_tick_params(which='both', bottom=True, top=False, direction='out', length=6, width=0.5)
 
for spine in ax.spines.values():
    spine.set_linewidth(0.5)

plt.xticks(ticks=range(len(sorted_cancer_types)), labels=new_labels, ha='center', fontsize=14)
plt.xlabel('Cancer type', fontsize=18, labelpad=15)
plt.ylabel('Number of arms discordant', fontsize=18)
plt.title('Tumor-model arm discordance\n(Matched models)', fontsize=20)
plt.tight_layout(rect=[0, 0, 1, 1])
plt.savefig(os.path.join(root, "plots", "Matched_TO_TP_num_arms_disc_by_cancer_type.pdf"), dpi=300)
plt.show()

for cancer_type, pval in pvalue_dict_arms.items():
    print(f"Arm discordance comparison for {cancer_type}: p = {pval:.5f}")

##########################################
## Part 3 ##
''' Tumor-organoid %genome and arm discordance by cohort (cancer type = color).
    * Each tumor-organoid pair is a data point (= each organoid=sample is
    a separate data point.
    * The plot includes multiple organoid from the same tumor
    * Only cohorts with >=3 organoids'''

filtered_data_genome = all_sample_matches_reshaped[
    (all_sample_matches_reshaped['Model'] == "Organoid") & 
    (~pd.isna(all_sample_matches_reshaped['Tumor'])) & 
    (all_sample_matches_reshaped["Percent_genome_discordance"] > 0)
].groupby(['Cohort_Name', 'Cancer_type']).filter(lambda x: len(x) >= 3)

filtered_data_genome['Mini_Cohort'] = filtered_data_genome['Cohort_Name'] + " | " + filtered_data_genome['Cancer_type']
median_genome_discordance = filtered_data_genome.groupby('Mini_Cohort')['Percent_genome_discordance'].median()
sorted_mini_cohorts_genome = median_genome_discordance.sort_values().index

color_mapping_genome = {mini_cohort: CUSTOM_COLORS.get(cancer, 'grey') 
                             for mini_cohort, cancer in zip(filtered_data_genome['Mini_Cohort'], filtered_data_genome['Cancer_type'])}

sns.set_theme(style="white")
os.makedirs(os.path.join(root, "plots"), exist_ok=True)

def plot_discordance(data, value_col, ylabel, title, output_name, decimal_places=2):
    df = data.loc[(data["Model"] == "Organoid") & (~pd.isna(data["Tumor"])) & (~pd.isna(data[value_col]))].copy()
    if value_col == "Percent_genome_discordance":
        df = df.loc[df[value_col] > 0]
    df["Mini_Cohort"] = df["Cohort_Name"] + " | " + df["Cancer_type"]
    df = df.groupby("Mini_Cohort").filter(lambda x: len(x) >= 3).reset_index(drop=True)
    median_vals = df.groupby("Mini_Cohort")[value_col].median()
    sorted_mini_cohorts = median_vals.sort_values().index
    color_mapping = {mc: CUSTOM_COLORS.get(ct, "grey") for mc, ct in zip(df["Mini_Cohort"], df["Cancer_type"])}
    
    g = sns.catplot(data=df, x="Mini_Cohort", y=value_col, kind="box", order=sorted_mini_cohorts, palette=color_mapping, showfliers=False, height=6, aspect=2, linewidth=0.5)
    sns.stripplot(data=df, x="Mini_Cohort", y=value_col, order=sorted_mini_cohorts, color="grey", size=4, jitter=True, alpha=0.6, ax=g.ax)

    for i, mini_cohort in enumerate(sorted_mini_cohorts):
        median_val = median_vals[mini_cohort]
        g.ax.text(i, median_val + 0.3, f"{median_val:.{decimal_places}f}" + ("%" if "Percent" in value_col else ""), 
        ha="center", va="bottom", fontsize=9, color="black", fontweight="medium")

    sample_counts = df["Mini_Cohort"].value_counts()
    tick_labels = [f"{mc}\nN = {sample_counts[mc]}" for mc in sorted_mini_cohorts]
    #tick_labels = [f"{mc}" for mc in sorted_mini_cohorts]
    g.ax.set_xticklabels(tick_labels, rotation=45, ha="right")

    mean_of_means = df.groupby("Mini_Cohort")[value_col].mean().mean()
    median_of_medians = df.groupby("Mini_Cohort")[value_col].median().median()
    y_max = df[value_col].max()

    is_genome_disc = (value_col == "Percent_genome_discordance")

    plt.text(-0.2, y_max * 0.90,
        f"Mean of means: {mean_of_means:.2f}{'%' if is_genome_disc else ''}",
        fontsize=12, color='black', ha='left')
    
    plt.text(-0.2, y_max * 0.98,
        f"Median of medians: {median_of_medians:.2f}{'%' if is_genome_disc else ''}",
        fontsize=12, color='black', ha='left')

    g.set_axis_labels("Mini cohort (Cohort | Cancer type)", ylabel)
    g.ax.set_xlabel(g.ax.get_xlabel(), fontsize=15, labelpad=-3)
    g.ax.set_ylabel(g.ax.get_ylabel(), fontsize=15, labelpad=5)
    
    g.ax.set_title(title, fontsize=18, pad=7, loc="center")

    g.ax.tick_params(axis="x", which="both", bottom=True, top=False, length=6, width=0.5, direction="out")
    g.ax.tick_params(axis="y", which="both", left=True, right=False, length=6, width=0.5, direction="out")
    for spine in g.ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.5)
    g.ax.grid(False)

    unique_cancers = df["Cancer_type"].unique()
    legend_handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=CUSTOM_COLORS.get(cancer, "grey"), markersize=10, label=cancer) for cancer in unique_cancers]
    g.ax.legend(handles=legend_handles, title="Cancer type", loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    if value_col == 'Percent_genome_discordance':
        g.ax.yaxis.set_major_formatter(mtick.PercentFormatter())
    g.fig.tight_layout(rect=[0, 0, 0.9, 0.93])
    g.fig.subplots_adjust(top=0.9)
    output_path = os.path.join(root, "plots", output_name)
    g.fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.show()

plot_discordance(data=all_sample_matches_reshaped, value_col="Percent_genome_discordance", ylabel="Genome discordance", title="PT-PDO genome discordance by cohort", output_name="TO_samples_genome_disc_by_mini_cohort.pdf", decimal_places=2)
plot_discordance(data=all_sample_matches_reshaped, value_col="Num_arms_disc", ylabel="Number of arms discordant", title="PT-PDO arm discordance by cohort", output_name="TO_samples_arms_disc_by_mini_cohort.pdf", decimal_places=1)

################################################
## Part 4 ##

'''Creating waterfall plots showing the delta in the genome/arm discordance between the earliest and
latest organoid we have for matches samples with more than one organoid regardless of passage.
The statistics is one-sided one sample ttest assuming that the delta will be positive.
Statistics results are in printouts'''

all_sample_matches_reshaped = pd.read_csv(os.path.join(root, "tables", "all_sample_matches_reshaped.csv"))
all_sample_matches_table = pd.read_csv(os.path.join(root, "tables", "all_sample_matches_table.csv"))

organoid_pass_cols = [col for col in all_sample_matches_table.columns if col.startswith("Organoid")]

def get_first_and_last_organoids(row):
    organoids = [row[c] for c in organoid_pass_cols if pd.notna(row[c])]
    if len(organoids) >= 2:
        return pd.Series({'Earliest_Organoid': organoids[0], 'Latest_Organoid': organoids[-1]})
    else:
        return pd.Series({'Earliest_Organoid': np.nan, 'Latest_Organoid': np.nan})

paired_samples = all_sample_matches_table.copy()
paired_samples[['Earliest_Organoid','Latest_Organoid']] = paired_samples.apply(get_first_and_last_organoids, axis=1)
paired_samples_latest = paired_samples.dropna(subset=['Earliest_Organoid','Latest_Organoid']).copy()

paired_samples_latest = paired_samples_latest.merge(
    all_sample_matches_reshaped[['Sample','Percent_genome_discordance','Num_arms_disc']],
    left_on='Earliest_Organoid', right_on='Sample', how='left'
).rename(columns={'Percent_genome_discordance':'Early_genome_disc',
                  'Num_arms_disc':'Early_num_arms_disc'}).drop(columns='Sample')

paired_samples_latest = paired_samples_latest.merge(
    all_sample_matches_reshaped[['Sample','Percent_genome_discordance','Num_arms_disc']],
    left_on='Latest_Organoid', right_on='Sample', how='left'
).rename(columns={'Percent_genome_discordance':'Late_genome_disc',
                  'Num_arms_disc':'Late_num_arms_disc'}).drop(columns='Sample')

paired_samples_latest = paired_samples_latest.dropna(subset=['Early_genome_disc','Late_genome_disc',
                                           'Early_num_arms_disc','Late_num_arms_disc'])

paired_samples_latest['Delta_genome'] = paired_samples_latest['Late_genome_disc'] - paired_samples_latest['Early_genome_disc']
paired_samples_latest['Delta_arms'] = paired_samples_latest['Late_num_arms_disc'] - paired_samples_latest['Early_num_arms_disc']

# Waterfall plot - genome discordance
waterfall_genome = paired_samples_latest.sort_values('Delta_genome').reset_index(drop=True)
plt.figure(figsize=(9,6))
colors = ['#DB3547' if x>0 else '#5BBA56' if x<0 else 'grey' for x in waterfall_genome['Delta_genome']]

plt.bar(range(len(waterfall_genome)), waterfall_genome['Delta_genome'], color=colors, width=0.8, linewidth=0.5, edgecolor='black')

plt.axhline(0,color='k',linewidth=0.5)
plt.ylabel('Genome discordance difference', fontsize = 18)
plt.xlabel('Number of models', fontsize = 18, labelpad = 8)
plt.title('Genome discordance difference - Earliest vs latest passage', fontsize = 20)
plt.xlim(-0.5,len(waterfall_genome)-0.5)

ax = plt.gca()
ax.xaxis.set_tick_params(which='both', bottom=True, top=False, direction='out', length=6, width=0.5, labelsize=14)
ax.yaxis.set_tick_params(which='both', left=True, right=False, direction='out', length=6, width=0.5, labelsize=14)

for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.5)
    spine.set_color('black')

ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=100, decimals=0))
output_folder = os.path.join(root, "plots")
os.makedirs(output_folder, exist_ok=True)
plt.savefig(os.path.join(output_folder, "waterfall_genome.pdf"), dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()

# Waterfall plot - arm discordance
waterfall_arms = paired_samples_latest.sort_values('Delta_arms').reset_index(drop=True)
plt.figure(figsize=(9,6))  
colors = ['#DB3547' if x>0 else '#5BBA56' if x<0 else 'grey' for x in waterfall_arms['Delta_arms']]

plt.bar(range(len(waterfall_arms)), waterfall_arms['Delta_arms'], color=colors, width=0.8, linewidth=0.5, edgecolor='black')

plt.axhline(0,color='k',linewidth=0.5)
plt.ylabel('Arm discordance difference', fontsize=18)
plt.xlabel('Number of models', fontsize=18, labelpad=8)
plt.title('Arm discordance difference - earliest vs latest passage', fontsize=20)
plt.xlim(-0.5,len(waterfall_arms)-0.5)

ax = plt.gca()
ax.xaxis.set_tick_params(which='both', bottom=True, top=False, direction='out', length=6, width=0.5, labelsize=14)
ax.yaxis.set_tick_params(which='both', left=True, right=False, direction='out', length=6, width=0.5, labelsize=14)

for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.5)
    spine.set_color('black')

ax.yaxis.set_major_formatter(mtick.ScalarFormatter())
plt.savefig(os.path.join(output_folder, "waterfall_arms.pdf"), dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()

p_genome_wf = ttest_1samp(waterfall_genome['Delta_genome'], 0, alternative='greater')[1]
p_arms_wf = ttest_1samp(waterfall_arms['Delta_arms'], 0, alternative='greater')[1]

print(f"Genome delta waterfall p-value: {p_genome_wf:.4e}")
print(f"Arm delta waterfall p-value: {p_arms_wf:.4e}")

#########

''' Creating matched lineplots showing the PT-PDO % genome/arm discordance between the earliest and
latest matched organoids we have regardless of passage.
This analysis includes only PDOs for which there is more than 1 passage'''


def plot_matched_lines(paired_df, early_col, late_col, ylabel, title, filename):
    plt.figure(figsize=(5,6))
    ax = plt.gca()
    for _, row in paired_df.iterrows():
        if row[late_col] > row[early_col]:
            color = '#DB3547'
        elif row[late_col] < row[early_col]:
            color = '#5BBA56'
        else:
            color = 'grey'
        plt.plot([0,1], [row[early_col], row[late_col]], marker='o', color=color, alpha=0.7, linewidth=0.5)
    plt.xticks([0,1], ['Earliest','Latest'], fontsize=16)
    plt.yticks(fontsize=14)
    plt.ylabel(ylabel, fontsize=16)
    plt.title(title, fontsize=16)
    
    ax.xaxis.set_tick_params(which='both', bottom=True, top=False, direction='out', length=6, width=0.5, labelsize=14)
    ax.yaxis.set_tick_params(which='both', left=True, right=False, direction='out', length=6, width=0.5, labelsize=14)

    for spine in ax.spines.values(): 
        spine.set_visible(True)
        spine.set_linewidth(0.5)
    plt.tight_layout()
    os.makedirs(output_folder, exist_ok=True)
    plt.savefig(os.path.join(root, "plots", f"{filename}.pdf"), dpi=300, bbox_inches='tight')
    plt.show()

import matplotlib as mpl

def plot_matched_lines_genome(paired_df, early_col, late_col, ylabel, title, filename):
    plt.figure(figsize=(5,6))
    ax = plt.gca()
    for _, row in paired_df.iterrows():
        if row[late_col] > row[early_col]:
            color = '#DB3547'
        elif row[late_col] < row[early_col]:
            color = '#5BBA56'
        else:
            color = 'grey'
        plt.plot([0,1], [row[early_col], row[late_col]], marker='o', color=color, alpha=0.7, linewidth=0.5)
    plt.xticks([0,1], ['Earliest','Latest'], fontsize=16)
    plt.yticks(fontsize=14)
    plt.ylabel("Genome discordance", fontsize=16)
    plt.title(title, fontsize=16)

    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(decimals=0))
    
    ax.xaxis.set_tick_params(which='both', bottom=True, top=False, direction='out', length=6, width=0.5, labelsize=14)
    ax.yaxis.set_tick_params(which='both', left=True, right=False, direction='out', length=6, width=0.5, labelsize=14)

    for spine in ax.spines.values(): 
        spine.set_visible(True)
        spine.set_linewidth(0.5)
    plt.tight_layout()
    os.makedirs(output_folder, exist_ok=True)
    plt.savefig(os.path.join(root, "plots", f"{filename}.pdf"), dpi=300, bbox_inches='tight')
    plt.show()

plot_matched_lines_genome(paired_samples_latest, 'Early_genome_disc', 'Late_genome_disc',
                          'Genome discordance', 'PT-PDO genome discordance\nEarliest vs latest passage',
                          'PT_PDO_genome_discordance_early_vs_latest')

plot_matched_lines(paired_samples_latest, 'Early_num_arms_disc', 'Late_num_arms_disc',
                   'Number of arms discordant', 'PT-PDO arm discordance\nEarliest vs latest passage',
                   'PT_PDO_arm_discordance_early_vs_latest')
