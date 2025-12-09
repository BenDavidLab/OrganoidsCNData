#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
06_PDO_PDX_PLOTS.PY - Comparative PDO vs PDX Analysis and Visualization
===============================================================================
Execution Order: 6 (Final plotting script)

Purpose:
    Creates comparative visualizations between PDO and PDX models, analyzing
    differences in genomic stability, passage effects, and cumulative discordance
    distributions. Provides statistical comparisons between model types.
    
    Part 1: Side-by-side boxplots comparing PDO vs PDX by cancer type (with p-values)
    Part 2: Regression analysis of discordance vs passage (with interaction terms)
    Part 3: Cumulative distribution (ECDF) plots comparing PDO vs PDX stability

Dependencies:
    - utility_functions.py
    - constants.py
    - all_sample_matches_reshaped.csv (from 02_PDO_sample_analysis.py)
    - merged_pdx_and_organoid_disc_table.csv (from 03_PDO_PDX_analysis.py)
    - hoge_pdx_disc_table_Linoy.csv (PDX data)

Inputs:
    - all_sample_matches_reshaped.csv
    - merged_pdx_and_organoid_disc_table.csv
    - PDX_data/hoge_pdx_disc_table_Linoy.csv

Outputs (plots/ directory):
    - Tumor-model_percent_genome_discordanceUnmatched_analysis_pval.pdf
    - Tumor-model_number_of_arms_discordantUnmatched_analysis_pval.pdf
    - PT_model_genome_discordance_vs_passage_with_ci_log.pdf
    - PT_model_arms_discordance_vs_passage_with_ci.pdf
    - genome_discordance_ecdf_mw.pdf
    - arms_discordance_ecdf_mw.pdf

Statistical Tests:
    - Independent t-tests (PDO vs PDX by cancer type)
    - Linear regression with interaction terms (passage effects)
    - Mann-Whitney U test (overall model comparison)
    - Kolmogorov-Smirnov test (distribution comparison)

Key Comparisons:
    - PDO vs PDX discordance levels across cancer types
    - Passage-dependent changes in both model types
    - Cumulative percentage of models exceeding stability thresholds

Usage:
    python 06_PDO_PDX_plots.py

Author: Linoy Raz
Contributors: Haia Khoury 
===============================================================================
"""

import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, linregress, ks_2samp, mannwhitneyu, t
import statsmodels.formula.api as smf
from sklearn.linear_model import LinearRegression
from matplotlib.ticker import PercentFormatter, MultipleLocator
from Code.utility_functions import *
from Code.constants import *

############## MAIN #####################
## Part 1 ##
''' Tumor-organoid %genome and arm discordance by cancer type (boxplots).
    This plot contains match organoids and PDXs (from current study) and also
    unmatched PDXs from the same cancer types (from Hoge)
    * Each tumor-organoid/tumor-pdx pair is a data point (= each organoid/pdx=
                                                          = sample = point)
    * Multiple PDOs from the same PT are individual data points
    * Only cohorts with >3 organoids'''

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) #Assumes scripts are in project root

all_sample_matches_reshaped = pd.read_csv(os.path.join(root, "tables", "all_sample_matches_reshaped.csv"))
hoge_PDX_table_linoy = pd.read_csv(os.path.join(root, "PDX_data", "hoge_pdx_disc_table_Linoy.csv"))
hoge_PDX_table_linoy = hoge_PDX_table_linoy.rename(columns={"cohort_name": "Cohort_Name"})

hoge_PT_PDX_only = hoge_PDX_table_linoy[hoge_PDX_table_linoy["comparison_type"] == "PT-PDX"]
hoge_PT_PDX_only.loc[:, "Model"] = "PDX"
hoge_PT_PDX_only.rename(columns={"pair2": "Sample", "pair2_passage": "Passage"}, inplace=True)

merged_df = pd.concat([all_sample_matches_reshaped, 
                       hoge_PT_PDX_only[['Cohort_Name','Cancer_type', 'Tumor', 'Model', 'Sample', "Passage", 'Percent_genome_discordance', 'Num_arms_disc']]], 
                       axis=0, ignore_index=True)
merged_df = merged_df.dropna(subset=['Tumor'])

cancer_types_with_organoid_data = merged_df[merged_df['Model'] == 'Organoid']['Cancer_type'].unique()
cancer_types_with_PDX_data = merged_df[merged_df['Model'] == 'PDX']['Cancer_type'].unique()
cancer_types_with_both = set(cancer_types_with_organoid_data).intersection(set(cancer_types_with_PDX_data))
df_filtered = merged_df[merged_df['Cancer_type'].isin(cancer_types_with_both)]

df_filtered = df_filtered.groupby('Cancer_type').filter(lambda x: (x['Model'] == 'Organoid').sum() >= 3 and (x['Model'] == 'PDX').sum() >= 3)

sorted_cohorts = df_filtered.groupby('Cancer_type')['Percent_genome_discordance'].median().sort_values(ascending=True).index

def plot_with_pvalues_and_sample_sizes(data, metric, title, y_label, save_dir, save_suffix):
    plt.figure(figsize=(14, 8))

    sorted_data = data[data['Cancer_type'].isin(sorted_cohorts)]

    sns.set_style("white")
    sns.set_style("ticks")
    plt.tick_params(axis='x', which='both', bottom=True, top=False)

    sns.boxplot(data=sorted_data, x='Cancer_type', y=metric, hue='Model', dodge=True, showfliers=False,
                order=sorted_cohorts, palette={'Organoid': '#A3D9B1', 'PDX':'#C9A7DA'})
    sns.stripplot(data=sorted_data, x='Cancer_type', y=metric, hue='Model', dodge=True, jitter=True,
                  alpha=0.5, color='gray', order=sorted_cohorts)
    plt.legend([], [], frameon=False)  

    y_height_offset = 0.05
    line_offset = 0.02

    new_labels = []
    for i, cancer_type in enumerate(sorted_cohorts):
        group_data = df_filtered[df_filtered['Cancer_type'] == cancer_type]
        organoid_data = group_data[group_data['Model'] == 'Organoid'][metric]
        pdx_data = group_data[group_data['Model'] == 'PDX'][metric]
        t_statistic, p_value = ttest_ind(organoid_data, pdx_data, alternative='less')

        max_value = group_data[metric].max()
        p_value_y = max_value + (y_height_offset * max_value)

        plt.plot([i - 0.2, i + 0.2], [p_value_y - line_offset, p_value_y - line_offset], color='black', lw=1)

        if p_value < 0.001:
            star_label = '***'
        elif p_value < 0.01:
            star_label = '**'
        elif p_value < 0.05:
            star_label = '*'
        else:
            star_label = 'ns'

        plt.text(i, p_value_y, star_label, fontsize=16, ha='center', va='bottom')

        organoid_n = len(organoid_data)
        pdx_n = len(pdx_data)
        new_labels.append(f"{cancer_type}")

    plt.xticks(ticks=[i for i in range(len(sorted_cohorts))], labels=new_labels, ha='center', fontsize=12, fontweight='normal')

    ax = plt.gca()
    ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=6))
    ax.yaxis.set_ticks_position('left')
    ax.tick_params(axis='both', which='major', length=6, width=0.5, labelsize=16)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)

    metric_values = sorted_data[metric]
    plt.ylim(metric_values.min() * 0.9, metric_values.max() * 1.1)

    plt.title(title, fontsize=24)
    plt.xlabel('Cancer type', fontsize=20)
    plt.ylabel(y_label, fontsize=20)
    if metric == 'Percent_genome_discordance':
        ax.yaxis.set_major_formatter(PercentFormatter(100))


    PDO_vals = sorted_data[sorted_data["Model"] == "Organoid"][metric].dropna()
    PDX_vals = sorted_data[sorted_data["Model"] == "PDX"][metric].dropna()

    ax.text(
        0.01, 0.89,
        f"Median of medians (PDOs): {PDO_vals.median():.2f}{'%' if metric=='Percent_genome_discordance' else ''}\n"
        f"Median of medians (PDXs): {PDX_vals.median():.2f}{'%' if metric=='Percent_genome_discordance' else ''}",
        transform=ax.transAxes,
        ha='left',
        va='top',
        fontsize=14,
        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    ax.text(
        0.01, 0.98,
        f"Mean of means (PDOs): {PDO_vals.mean():.2f}{'%' if metric=='Percent_genome_discordance' else ''}\n"
        f"Mean of means (PDXs): {PDX_vals.mean():.2f}{'%' if metric=='Percent_genome_discordance' else ''}",
        transform=ax.transAxes,
        ha='left',
        va='top',
        fontsize=14,
        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    plt.tight_layout()
    plt.subplots_adjust(left=0.12, bottom=0.25)

    save_title = title.replace(' ', '_').replace('\n', '').replace('\\', '_')
    save_path = os.path.join(save_dir, f"{save_title}_{save_suffix}.pdf")
    os.makedirs(save_dir, exist_ok=True)
    plt.savefig(save_path, dpi=300)
    plt.show()


plot_with_pvalues_and_sample_sizes(df_filtered,'Percent_genome_discordance',
    'Tumor-model percent genome discordance\nUnmatched analysis','Genome discordance',
    'plots','pval')

plot_with_pvalues_and_sample_sizes(df_filtered,'Num_arms_disc',
    'Tumor-model number of arms discordant\nUnmatched analysis','Number of arms discordant',
    'plots','pval')


#######################################
## Part 2 ##
'''Plotting the change in percent genome discordance or number of arms discordant
   over passages in matched + unmatched organoids and PDXs
Statistical test - linear regression for each model, with an interaction term'''

merged_pdo_pdx = pd.read_csv(os.path.join(root, "tables", "merged_pdx_and_organoid_disc_table.csv"))
merged_pdo_pdx = merged_pdo_pdx[pd.to_numeric(merged_pdo_pdx['Passage'], errors='coerce').notnull()]
merged_pdo_pdx['Passage_Num'] = merged_pdo_pdx['Passage'].astype(int)
merged_pdo_pdx['Passage_Cat'] = merged_pdo_pdx['Passage_Num'].astype(str)
palette = {'Organoid': '#A3D9B1', 'PDX': '#C9A7DA'}
hue_order = ['Organoid', 'PDX']
legend_labels = {'Organoid': 'PT-PDO', 'PDX': 'PT-PDX'}
trendline_colors = {'Organoid': '#7BAF89', 'PDX': '#A582B9'}

def get_regression_ci(X_fit, y_fit, X_pred, reg_model, alpha=0.05):
    n = len(y_fit)
    if n < 3:
        return None, None
    df = n-2
    X_mean = np.mean(X_fit)
    ss_x = np.sum((X_fit - X_mean)**2)
    y_pred_fit = reg_model.predict(X_fit)
    rss = np.sum((y_fit - y_pred_fit)**2)
    mse = rss / df
    t_crit = t.ppf(1-alpha / 2,df)
    X_pred_flat = X_pred.flatten()
    if ss_x == 0:
        se_mean_pred = np.sqrt(mse * (1/n))
    else:
        se_mean_pred = np.sqrt(mse * (1/n + (X_pred_flat - X_mean)**2 / ss_x))
    margin_of_error = t_crit * se_mean_pred
    y_pred = reg_model.predict(X_pred)
    ci_lower = y_pred.flatten() - margin_of_error
    ci_upper = y_pred.flatten() + margin_of_error
    return ci_lower, ci_upper

def plot_with_linear_regression(df, value_col, label, show_trendline_ci=True):
    plt.figure(figsize=(13, 7))
    results = {'mann_whitney_u_p': None, 'trendlines': {}}
    numerical_passages = df['Passage_Num'].unique()
    numerical_passages.sort()
    plot_order = [str(p) for p in numerical_passages]
    ax = sns.barplot(data=df, x='Passage_Cat', y=value_col, hue='Model', palette=palette, errorbar=('ci', 95), err_kws={'linewidth': 1.0}, edgecolor='black', linewidth=0.8, hue_order=hue_order, order=plot_order)
    for patch in ax.patches:
        patch.set_linewidth(0.5)
    ax.set_xlabel('Passage', fontsize=18)
    ax.set_ylabel(label, fontsize=18)
    ax.set_title(label + ' by passage', fontsize=22)
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    ax.tick_params(axis='x', which='both', direction='out', length=6, width=0.5, bottom=True, labelsize=16)
    ax.tick_params(axis='y', which='both', direction='out', length=6, width=0.5, left=True, labelsize=16)
    org_all = df.loc[df['Model']=='Organoid', value_col].dropna()
    pdx_all = df.loc[df['Model']=='PDX', value_col].dropna()
    stat, pval = mannwhitneyu(org_all, pdx_all, alternative='less')
    results['mann_whitney_u_p'] = pval
    if ax.patches:
        tallest_bar = max(ax.patches, key=lambda b: b.get_height())
        ax.set_ylim(0, tallest_bar.get_height() * 1.05)
    passage_mapping = {p_str: i for i, p_str in enumerate(plot_order)}
    ax.set_xlim(-0.5, len(plot_order)-0.5)
    if value_col == 'Percent_Genome_Discordance':
        from matplotlib.ticker import FuncFormatter
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{y:.0f}%'))

    for model in hue_order:
        df_model = df[df['Model'] == model].groupby('Passage_Num')[value_col].mean().reset_index()
        X_fit = df_model['Passage_Num'].values.reshape(-1, 1)
        y_fit = df_model[value_col].values
        X_pred = None
        x_plot_range = None
        if model == 'Organoid':
            mask = X_fit.flatten() <= 31
            X_limited = X_fit[mask]
            y_limited = y_fit[mask]
            if X_limited.size > 0:
                reg = LinearRegression().fit(X_limited, y_limited)
                X_pred = np.linspace(X_limited.min(), 31.0, 100).reshape(-1, 1)
                start_passage_str = str(X_limited.min().item())
                start_index = passage_mapping.get(start_passage_str, 0)
                try:
                    end_index = passage_mapping['31']
                except KeyError:
                    end_index = len(plot_order) - 1
                x_plot_range = np.linspace(start_index, end_index, 100)
            else:
                continue
        else:
            if X_fit.size > 0:
                reg = LinearRegression().fit(X_fit, y_fit)
                X_pred = np.linspace(X_fit.min(), X_fit.max(), 100).reshape(-1, 1)
                start_passage_str = str(X_fit.min().item())
                end_passage_str = str(X_fit.max().item())
                start_index = passage_mapping.get(start_passage_str, 0)
                end_index = passage_mapping.get(end_passage_str, len(plot_order) - 1)
                x_plot_range = np.linspace(start_index, end_index, 100)
            else:
                continue
        results['trendlines'][model] = {'slope': reg.coef_[0], 'intercept': reg.intercept_}
        if show_trendline_ci:
            X_data = X_limited if model == 'Organoid' else X_fit
            y_data = y_limited if model == 'Organoid' else y_fit
            ci_lower, ci_upper = get_regression_ci(X_data, y_data, X_pred, reg)
            if ci_lower is not None:
                ax.fill_between(x_plot_range, ci_lower, ci_upper, color=trendline_colors[model], alpha=0.2, label=f'{legend_labels[model]} 95% CI')
        y_pred = reg.predict(X_pred)
        ax.plot(x_plot_range, y_pred, color=trendline_colors[model], linewidth=2.5)
    final_handles = [plt.Rectangle((0,0), 1, 1, fc=palette[m], ec='black', lw=0.5) for m in hue_order]
    final_labels = [legend_labels[m] for m in hue_order]
    ax.legend(final_handles, final_labels, title='Model', fontsize=10, loc='upper left', bbox_to_anchor=(1.01, 1))
    plt.tight_layout(rect=[0,0,0.95,1])
    return results

plt.close('all')
genome_results = plot_with_linear_regression(merged_pdo_pdx, 'Percent_Genome_Discordance', 'Average genome discordance', show_trendline_ci=True)
plt.savefig(os.path.join(root,"plots", "average_genome_disc_over_passages.pdf"), dpi=300)
plt.show()

print(f"Percent genome discordance\nOne sided Mann-Whitney U test P-value = {genome_results['mann_whitney_u_p']:.10f}")
print(f"Trendline: y = slope * Passage + intercept")
for model, data in genome_results['trendlines'].items():
    label = legend_labels[model]
    slope = data['slope']
    intercept = data['intercept']
    slope_str = f"{slope:.4f}"
    intercept_str = f"{intercept:.4f}"
    if slope >= 0:
        equation = f"y = {slope_str} * Passage + {intercept_str}"
    else:
        equation = f"y = {slope_str} * Passage {intercept_str}"
    print(f"* {label}: {equation}")
print("\n" + "="*50 + "\n")
plt.close('all')
arms_results = plot_with_linear_regression(merged_pdo_pdx, 'Num_arms_disc', 'Average number of arms discordant', show_trendline_ci=True)
plt.savefig(os.path.join(root,"plots", "average_arm_disc_over_passages.pdf"), dpi=300)
plt.show()

print(f"Number of arms discordant\nOne-sided Mann-Whitney U test P-value = {arms_results['mann_whitney_u_p']:.10f}")
print(f"Trendline: y = slope * Passage + intercept")
for model, data in arms_results['trendlines'].items():
    label = legend_labels[model]
    slope = data['slope']
    intercept = data['intercept']
    slope_str = f"{slope:.4f}"
    intercept_str = f"{intercept:.4f}"
    if slope >= 0:
        equation = f"y = {slope_str} * Passage + {intercept_str}"
    else:
        equation = f"y = {slope_str} * Passage {intercept_str}"
    print(f"* {label}: {equation}")
print("\n" + "="*50 + "\n")
plt.close('all')
plot_with_linear_regression(merged_pdo_pdx, 'Percent_Genome_Discordance', 'Average percent genome discordance', show_trendline_ci=False)
plt.show()
plt.close('all')
plot_with_linear_regression(merged_pdo_pdx, 'Num_arms_disc', 'Average number of arms discordant', show_trendline_ci=False)
plt.show()

#######
'''printing the genome and arm discordance of each model at passage 0 and across all passages'''
first_passage = merged_pdo_pdx.groupby('Model')['Passage_Num'].min()

for model in hue_order:
    first_passage_num = first_passage[model]
    subset = merged_pdo_pdx[(merged_pdo_pdx['Model'] == model) & 
                            (merged_pdo_pdx['Passage_Num'] == first_passage_num)]
    
    avg_genome_disc = subset['Percent_Genome_Discordance'].mean()
    avg_arm_disc = subset['Num_arms_disc'].mean()
    
    label = legend_labels[model]
    print(f"{label} at first passage (Passage {first_passage_num}):")
    print(f"  Average percent genome discordance: {avg_genome_disc:.2f}%")
    print(f"  Average number of arms discordant: {avg_arm_disc:.2f}\n")
    
####### 
'''Printing average and median genome and arm discordance across all passages for each model'''

hue_order = ['Organoid', 'PDX']
legend_labels = {'Organoid': 'PT-PDO', 'PDX': 'PT-PDX'}

first_passage = merged_pdo_pdx.groupby('Model')['Passage_Num'].min()

for model in hue_order:
    label = legend_labels[model]
    
    first_passage_num = first_passage[model]
    subset_first = merged_pdo_pdx[
        (merged_pdo_pdx['Model'] == model) &
        (merged_pdo_pdx['Passage_Num'] == first_passage_num)]
    avg_genome_first = subset_first['Percent_Genome_Discordance'].mean()
    avg_arm_first = subset_first['Num_arms_disc'].mean()
    
    print(f"{label} at first passage (Passage {first_passage_num}):")
    print(f"  Average percent genome discordance: {avg_genome_first:.2f}%")
    print(f"  Average number of arms discordant: {avg_arm_first:.2f}\n")
    
    subset_all = merged_pdo_pdx[merged_pdo_pdx['Model'] == model]
    avg_genome_all = subset_all['Percent_Genome_Discordance'].mean()
    median_genome_all = subset_all['Percent_Genome_Discordance'].median()
    avg_arm_all = subset_all['Num_arms_disc'].mean()
    median_arm_all = subset_all['Num_arms_disc'].median()
    
    print(f"{label} across all passages:")
    print(f"  Average percent genome discordance: {avg_genome_all:.2f}%")
    print(f"  Median percent genome discordance: {median_genome_all:.2f}%")
    print(f"  Average number of arms discordant: {avg_arm_all:.2f}")
    print(f"  Median number of arms discordant: {median_arm_all:.2f}\n")
    print("="*50 + "\n")

###############################
## Part 3 ##
'''Plotting the cumulative percentage of samples (unmatched organoids and PDXs)
that have at least a specific change in their genome (% of genome / # of arms).
The comparison is with a KS statistical test*** '''

merged_df['Num_arms_disc'] = pd.to_numeric(merged_df['Num_arms_disc'], errors='coerce')
merged_df['Percent_genome_discordance'] = pd.to_numeric(merged_df['Percent_genome_discordance'], errors='coerce')

filtered_df = merged_df[merged_df['Model'].isin(['Organoid', 'PDX'])]

mw_stat_genome, mw_pvalue_genome = mannwhitneyu(
    filtered_df[filtered_df['Model'] == 'Organoid']['Percent_genome_discordance'].dropna(),
    filtered_df[filtered_df['Model'] == 'PDX']['Percent_genome_discordance'].dropna(),
    alternative='less')
print(f"Genome discordance one-sided MWU: U={mw_stat_genome:.3f}, p={mw_pvalue_genome:.10f}")

sns.set_theme(style="white", context="talk")

mw_stat_arms, mw_pvalue_arms = mannwhitneyu(
    filtered_df[filtered_df['Model'] == 'Organoid']['Num_arms_disc'].dropna(),
    filtered_df[filtered_df['Model'] == 'PDX']['Num_arms_disc'].dropna(),
    alternative='less')
print(f"Arm discordance one-sided MWU: U={mw_stat_arms:.3f}, p={mw_pvalue_arms:.10f}")

x_line = 25

plot_df = filtered_df.dropna(subset=['Percent_genome_discordance'])
pdo_values = plot_df[plot_df['Model'] == 'Organoid']['Percent_genome_discordance'].values
pdx_values = plot_df[plot_df['Model'] == 'PDX']['Percent_genome_discordance'].values

pdo_fraction = np.sum(pdo_values >= x_line) / len(pdo_values)
pdx_fraction = np.sum(pdx_values >= x_line) / len(pdx_values)

print(f"Cumulative percent of models with ≥25% genome discordant:")
print(f"PDOs: {pdo_fraction*100:.2f}%")
print(f"PDXs: {pdx_fraction*100:.2f}%")

plt.figure(figsize=(8, 6))
palette = {'Organoid': '#A3D9B1', 'PDX': '#C9A7DA'}
sns.ecdfplot(data=plot_df, x='Percent_genome_discordance', hue='Model', complementary=True, linewidth=2, palette=palette)
plt.title('Genome discordance across samples\n(Unmatched models)', fontsize=16)
plt.xlabel('Genome discordance', fontsize=16)
plt.ylabel('Cumulative percent of models', fontsize=16)
x_min = plot_df['Percent_genome_discordance'].min()
x_max = plot_df['Percent_genome_discordance'].max()
padding = 0.02 * (x_max - x_min)
plt.xlim(x_min - padding, 100)
plt.ylim(-0.05, 1.05)
plt.gca().xaxis.set_major_locator(MultipleLocator(25))
plt.gca().xaxis.set_major_formatter(PercentFormatter(100))
plt.gca().yaxis.set_major_locator(MultipleLocator(0.25))
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.tick_params(axis='x', which='major', direction='out', length=6, width=0.5, labelsize=14, bottom=True, top=False)
plt.tick_params(axis='y', which='major', direction='out', length=6, width=0.5, labelsize=14, left=True, right=False)
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(0.5)
leg = ax.get_legend()
if leg:
    leg.set_title('Model', prop={'size': 12})
    for text in leg.get_texts():
        if text.get_text() == 'Organoid':
            text.set_text('PDO')
        text.set_fontsize(12)
    for line in leg.get_lines():
        line.set_linewidth(3)
    leg.set_bbox_to_anchor((1, 1))
    leg.set_loc('upper left')
    leg.set_frame_on(False)
x_start = ax.get_xlim()[0]
y_bottom = ax.get_ylim()[0]
for model in palette.keys():
    y_values = np.sort(plot_df[plot_df['Model'] == model]['Percent_genome_discordance'].values)
    n = len(y_values)
    y = np.sum(y_values >= x_line) / n
    plt.vlines(x_line, ymin=y_bottom, ymax=y, linestyles='dashed', colors='grey', linewidth=1.5)
    plt.hlines(y, x_start, x_line, linestyles='dashed', colors='grey', linewidth=1.5)
plt.tight_layout()
plt.savefig('plots/genome_discordance_ecdf_mw.pdf', dpi=300)
plt.show()

plt.figure(figsize=(8, 6))
plot_df = filtered_df.dropna(subset=['Num_arms_disc'])
sns.ecdfplot(data=plot_df, x='Num_arms_disc', hue='Model', complementary=True, linewidth=2, palette=palette)
plt.title('Arm discordance across samples\n(Unmatched models)', fontsize=16)
plt.xlabel('Number of arms discordant', fontsize=16)
plt.ylabel('Cumulative percent of models', fontsize=16)
plt.xlim(plot_df['Num_arms_disc'].min() - 0.02*(plot_df['Num_arms_disc'].max() - plot_df['Num_arms_disc'].min()), 30)
plt.ylim(-0.05, 1.05)
plt.gca().xaxis.set_major_locator(MultipleLocator(5))
plt.gca().yaxis.set_major_locator(MultipleLocator(0.25))
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.tick_params(axis='x', which='major', direction='out', length=6, width=0.5, labelsize=14, bottom=True, top=False)
plt.tick_params(axis='y', which='major', direction='out', length=6, width=0.5, labelsize=14, left=True, right=False)
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(0.5)
leg = ax.get_legend()
if leg:
    leg.set_title('Model', prop={'size': 12})
    for text in leg.get_texts():
        if text.get_text() == 'Organoid':
            text.set_text('PDO')
        text.set_fontsize(12)
    for line in leg.get_lines():
        line.set_linewidth(3)
    leg.set_bbox_to_anchor((1, 1))
    leg.set_loc('upper left')
    leg.set_frame_on(False)
plt.tight_layout()
plt.savefig('plots/arms_discordance_ecdf_mw.pdf', dpi=300)
plt.show()


##########################

'''Calculating the cumulative percent of models that have >25% genome discordance in PDOs vs PDXs.'''

organoid_data = filtered_df[filtered_df['Model'] == 'Organoid']['Percent_genome_discordance'].dropna()
pdx_data = filtered_df[filtered_df['Model'] == 'PDX']['Percent_genome_discordance'].dropna()
x_value = 25

organoid_count_above_x = (organoid_data > x_value).sum()
pdx_count_above_x = (pdx_data > x_value).sum()

organoid_total_count = len(organoid_data)
pdx_total_count = len(pdx_data)

if organoid_total_count > 0:
    organoid_cumulative_percent = (organoid_count_above_x / organoid_total_count) * 100
else:
    organoid_cumulative_percent = 0

if pdx_total_count > 0:
    pdx_cumulative_percent = (pdx_count_above_x / pdx_total_count) * 100
else:
    pdx_cumulative_percent = 0

print(f"--- Cumulative Percent of Models at x={x_value}% Genome Discordance ---")
print(f"For Organoid models (PDO):")
print(f"  Total models analyzed: {organoid_total_count}")
print(f"  Models with discordance > {x_value}%: {organoid_count_above_x}")
print(f"  Cumulative percent of models (y) at x={x_value}%: {organoid_cumulative_percent:.2f}%")
print(f"\nFor PDX models:")
print(f"  Total models analyzed: {pdx_total_count}")
print(f"  Models with discordance > {x_value}%: {pdx_count_above_x}")
print(f"  Cumulative percent of models (y) at x={x_value}%: {pdx_cumulative_percent:.2f}%")

##########################
