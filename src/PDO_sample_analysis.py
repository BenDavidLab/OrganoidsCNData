import os
import sys
import pandas as pd
import numpy as np
import re
from collections import defaultdict

#root = r"/workspaces/Organoid_Stability"
#sys.path.append(os.path.abspath(root))

os.chdir(r"C:\Users\linoy\OneDrive\Desktop\New folder\Organoid_Stability")
root = r"C:\Users\linoy\OneDrive\Desktop\New folder\Organoid_Stability"

from Code.utility_functions import *
from Code.constants import *

metadata = pd.read_csv(os.path.join(root, "metadata_organoids.csv"), header=0)
metadata = metadata.drop(metadata.index[-1]).copy()

## Main ##
##############################
## Part 1 ##
''' Creating a matches table for all the samples by patient.
    Each patient is a row that can include multiple organoids, PDXs, normal
    tissue samples, etc.'''

all_sample_matches_table = pd.DataFrame()

for _, row in metadata.iterrows():
    cohort_name = row["cohort_name"]
    matches_table = pd.read_csv(os.path.join(root, "Cohort_data", f"{row['cohort_index']}_{cohort_name}", f"{row['cohort_index']}_{cohort_name}_matches_table.csv"))
    cancer_type = row["Cancer_type"]

    for j, row2 in matches_table.iterrows(): 
        tumor_columns = [col for col in matches_table.columns if col.startswith('Tumor')]         
        organoid_columns = [col for col in matches_table.columns if col.startswith('Organoid')]
        pdx_columns = [col for col in matches_table.columns if col.startswith('PDX')]        
        
        new_row = {
            "Cohort_Name": cohort_name,
            "Cancer_type": cancer_type,
            "Tumor": row2.get("Tumor", np.nan),
            "Organoid": row2.get("Organoid", np.nan),
            "Normal_tissue": row2.get("Normal_tissue", np.nan)}
        
        for tumor_col in tumor_columns:
            new_row[tumor_col] = row2.get(tumor_col, np.nan)
            
        if cancer_type == "multiple_cancers":
            if pd.notna(new_row["Organoid"]):  
                matches_table_row = matches_table[matches_table["Organoid"] == new_row["Organoid"]]
                if not matches_table_row.empty:
                    new_row["Cancer_type"] = matches_table_row["Cancer_type"].iloc[0]
                else:
                    new_row["Cancer_type"] = np.nan
            else:
                new_row["Cancer_type"] = np.nan


        new_row = pd.DataFrame([new_row])
        

        for organoid_col in organoid_columns:
            new_row[organoid_col] = row2.get(organoid_col, np.nan)
        for pdx_col in pdx_columns:
            new_row[pdx_col] = row2.get(pdx_col, np.nan)
        
      
        if pd.notna(new_row["Organoid"].iloc[0]) and any(new_row[key].iloc[0] is not np.nan for key in ["Tumor", "Normal_tissue"] + organoid_columns + pdx_columns):
            all_sample_matches_table = pd.concat([all_sample_matches_table, new_row], ignore_index=True)

all_sample_matches_table = map_cancer_type(all_sample_matches_table)
all_sample_matches_table.to_csv(os.path.join(root, "tables", "all_sample_matches_table.csv"), index=False)

######################################################################
## Part 2 ##
''' Reshaping the previous table and keeping only the tumor-organoid matches
    so that each row is an sample (a single tumor will appear across multiple rows).
    This includes multiple organoid/PDX passages from the same tumor.'''

tumor_columns = [col for col in all_sample_matches_table.columns if col.startswith('Tumor')]
organoid_columns = [col for col in all_sample_matches_table.columns if col.startswith('Organoid')]
pdx_columns = [col for col in all_sample_matches_table.columns if col.startswith('PDX')]

sample_to_tumors = defaultdict(set)
for _, row in all_sample_matches_table.iterrows():
    tumors_in_row = [str(row[col]).strip() for col in tumor_columns if pd.notna(row[col]) and str(row[col]).strip() != "nan"]
    for col in organoid_columns + pdx_columns:
        val = row.get(col)
        if pd.notna(val):
            sample_to_tumors[val].update(tumors_in_row)

reshaped_data = []
for _, row in all_sample_matches_table.iterrows():
    for col in organoid_columns:
        val = row.get(col)
        if pd.notna(val):
            sample_name = str(val).strip()
            tumor_for_this = sorted(sample_to_tumors.get(sample_name, []))
            num_tumor_sections = len(tumor_for_this)
            passage_match = re.search(r'P(\d+)$', sample_name)
            passage = passage_match.group(1) if passage_match else np.nan
            reshaped_data.append({
                "Cohort_Name": row["Cohort_Name"],
                "Cancer_type": row["Cancer_type"],
                "Normal_tissue": row["Normal_tissue"],
                "Tumor": ";".join(tumor_for_this),
                "Sample": sample_name,
                "Model": "Organoid",
                "Passage": passage,
                "Num_tumor_sections": num_tumor_sections})
    
    for col in pdx_columns:
        val = row.get(col)
        if pd.notna(val):
            sample_name = str(val).strip()
            tumor_for_this = sorted(sample_to_tumors.get(sample_name, []))
            num_tumor_sections = len(tumor_for_this)
            passage_match = re.search(r'P(\d+)$', sample_name)
            passage = passage_match.group(1) if passage_match else np.nan
            reshaped_data.append({
                "Cohort_Name": row["Cohort_Name"],
                "Cancer_type": row["Cancer_type"],
                "Normal_tissue": row["Normal_tissue"],
                "Tumor": ";".join(tumor_for_this),
                "Sample": sample_name,
                "Model": "PDX",
                "Passage": passage,
                "Num_tumor_sections": num_tumor_sections})

all_sample_matches_reshaped = pd.DataFrame(reshaped_data)

percent_genome_discordance_list = []
num_arms_disc_pipe_list = []
for _, row in all_sample_matches_reshaped.iterrows():
    cohort_name = row["Cohort_Name"]
    cohort_index = metadata.loc[metadata['cohort_name'] == cohort_name, 'cohort_index'].values[0]
    bins_table = pd.read_csv(os.path.join(root, "Cohort_data", f"{cohort_index}_{cohort_name}", f"{cohort_index}_{cohort_name}_mb_bins_normalized.csv"), index_col=0)
    arms_table = pd.read_csv(os.path.join(root, "Cohort_data", f"{cohort_index}_{cohort_name}", f"{cohort_index}_{cohort_name}_arm_level_normalized.csv"), index_col="Sample_ID", header=0)
    tumor_values = [t.strip() for t in str(row["Tumor"]).split(";") if t.strip() and t.strip() != "nan"]
    sample = row["Sample"]
    genome_discordances = []
    arm_discordances = []
    for tumor_section in tumor_values:
        if tumor_section in bins_table.index and sample in bins_table.index:
            genome_discordances.append(perc_genome_disc(bins_table.loc[tumor_section], bins_table.loc[sample]))
        if tumor_section in arms_table.index and sample in arms_table.index:
            arm_discordances.append(arms_by_pipe(arms_table.loc[tumor_section], arms_table.loc[sample]))
    avg_genome_disc = np.nanmean(genome_discordances) if genome_discordances else np.nan
    avg_arm_disc = np.nanmean(arm_discordances) if arm_discordances else np.nan
    percent_genome_discordance_list.append(avg_genome_disc)
    num_arms_disc_pipe_list.append(avg_arm_disc)

all_sample_matches_reshaped["Percent_genome_discordance"] = all_sample_discordance = pd.Series(percent_genome_discordance_list).round(4)
all_sample_matches_reshaped["Num_arms_disc"] = pd.Series(num_arms_disc_pipe_list).round(4)
all_sample_matches_reshaped.to_csv(os.path.join(root, "tables", "all_sample_matches_reshaped.csv"), index=False)


##########################################
## Part 3 ##
'''Creating an average parallel to all_sample_matches_reshaped, showing the average discordance 
values between tumors and models. 
* For tumors that had multiple models derived (e.g. organoid passages), the discordance will be 
the average of all models. 
* For tumors with singular models derived the values will not change.

This df does not address cases in which one patient has multiple tumors - each tumor is regarded as a 
separate entity.
'''

all_sample_matches_reshaped_average = all_sample_matches_reshaped.copy()

all_sample_matches_reshaped_average["Tumor_list"] = all_sample_matches_reshaped_average["Tumor"].str.split(";")
all_sample_matches_reshaped_average["Organoid_list"] = all_sample_matches_reshaped_average["Sample"].str.split(";")

all_sample_matches_reshaped_average["Passage"] = pd.to_numeric(all_sample_matches_reshaped_average["Passage"], errors="coerce")

average_rows = []

for _, row in all_sample_matches_reshaped_average.iterrows():
    tumors = [t.strip() for t in row["Tumor_list"] if t.strip() != ""]
    organoids = [o.strip() for o in row["Organoid_list"] if o.strip() != ""]

    percent_discordances = []
    num_arms_discordances = []

    for org_idx, org in enumerate(organoids):
        for tum_idx, tum in enumerate(tumors):
            percent_discordances.append(row["Percent_genome_discordance"])
            num_arms_discordances.append(row["Num_arms_disc"])

    avg_percent = np.nanmean(percent_discordances) if percent_discordances else np.nan
    avg_arms = np.nanmean(num_arms_discordances) if num_arms_discordances else np.nan

    average_rows.append({
        "Cohort_Name": row["Cohort_Name"],
        "Cancer_type": row["Cancer_type"],
        "Normal_tissue": row["Normal_tissue"],
        "Tumor": ";".join(tumors),
        "Sample": ";".join(organoids),
        "Model": row["Model"],
        "Percent_genome_discordance": round(avg_percent, 4),
        "Num_arms_disc": round(avg_arms, 4),
        "Passage": round(row["Passage"], 2),
        "Num_tumor_sections_used": len(tumors),
        "Num_samples_used": len(organoids)})

all_sample_matches_reshaped_average = pd.DataFrame(average_rows)
all_sample_matches_reshaped_average.to_csv(os.path.join(root, "tables", "all_sample_matches_reshaped_average.csv"), index=False)
