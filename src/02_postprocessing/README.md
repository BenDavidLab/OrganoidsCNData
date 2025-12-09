# Post-Processing Pipeline: CNV Discordance Analysis

This directory contains scripts for analyzing copy number variation (CNV) discordance between primary tumors and patient-derived models (organoids and xenografts) and generating publication-quality visualizations.

### Pipeline Overview
Pre-processed CNV data (from preprocessing/) → [01_update_metadata.py] → Updated cohort statistics  
[02_PDO_sample_analysis.py] → Sample-level discordance tables    
[03_PDO_PDX_analysis.py] → Integrated PDO-PDX analysis    
[04_descriptive_plots.py] → Descriptive statistics figures    
[05_PDO_plots.py] → PDO-specific analysis plots    
[06_PDO_PDX_plots.py] → Comparative PDO vs PDX plots        
results/tables/ + results/plots/    

#### Requirements    
Python Dependencies - All dependencies are listed in requirements.txt:    

bash - Install from repository root    
pip install -r requirements.txt    

#### Module Files    
* constants.py:    

Chromosome arm coordinates    
Arm names list    
Discordance thresholds (PERCENT_CUTOFF = 75)    
Custom color schemes for cancer types    

* utility_functions.py:    

Core discordance calculation functions    
Cancer type mapping    
Sample name parsing    
Helper visualization functions    



