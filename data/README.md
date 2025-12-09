# Data Directory
This directory contains input data files required for the analysis pipeline.

### Files
- metadata_organoids.csv
Cohort-level metadata file that tracks all patient cohorts in the study.
Required Columns:

cohort_index (int): Unique cohort identifier    
cohort_name (str): Cohort name    
Cancer_type (str): Cancer type designation    
Patients (int): Number of patients    
Tumor (int): Number of tumor samples    
Organoid (int): Number of organoid samples    

Note: This file is updated by src/postprocessing/01_update_metadata.py to add calculated statistics.

- reference_files/
Contains genome reference files used for alignment and CNV calling in the pre-processing pipeline.


| Cohort             | Genome build used (based on original study/kit) |
|--------------------|---------------------------------------|
| 1-Abud_Lab         | hg38                                  |
| 2-Blandino_Lab     | hg19                                  |
| 3-De_Julio_Lab     | hg38                                  |
| 4-Furukawa_Lab     | hg19                                  |
| 5-Huss_Lab         | hg38                                  |
| 6-Khurana_Lab      | hg19                                  |
| 7-Leung_Lab        | hg19                                  |
| 8-Rubin_Lab        | hg19                                  |
| 9-Shen_Lab         | hg19                                  |
| 10-Yeung_Lab       | hg38                                  |
| 11-Heim_Lab        | hg19                                  |
| 12-Clevers_Lab     | hg19                                  |
