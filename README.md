# CNV Analysis from WES Data

## Repository Structure

This repository contains Copy Number Variation (CNV) calling results for multiple cohorts analyzed using whole-exome sequencing (WES).

### Cohorts 1-10 (Legacy)
- Analyzed using original pipeline (see `legacy/`)
- Results are final and preserved for reproducibility

### Cohorts 11-13 (Current Pipeline)
- Analyzed using improved pipeline with enhanced logging, validation, and error handling
- Pipeline is modular, documented, and production-ready (see `pipeline/`)
- **Methodologically equivalent** to legacy pipeline but with improved code quality

### Results
All cohort results follow the same structure:
- `mb_bins_normalized.csv`
- `arm_level_normalized.csv`
- `segments_of_all_samples.seg`

## Running the Pipeline
- Aligns FASTQ files using BWA-MEM2, handles multiple lanes, marks duplicates.
- Runs CNVkit to call copy number variations.
- Combines all individual SEG files into one master file.
- Generates 1Mb bins and arm-level summaries.

## Citation
[When paper is published]


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
| 13-Blandino2_Lab   | hg19                                  |
