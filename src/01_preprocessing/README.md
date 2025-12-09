# Pre-Processing Pipeline: FASTQ → CNV    
This directory contains the pre-processing pipeline that converts raw sequencing data (FASTQ) to normalized copy number variation (CNV) calls at 1MB bin and chromosome arm resolution.    

### Pipeline Overview    
FASTQ files → [00_fastq_to_bam.py] → BAM files    
BAM files → [01_bam_to_cnv.py] → CNV calls    
CNV calls → [02_cnv_to_segments.py] → Segmented CNV    
Segmented CNV → [03_segments_to_bins.py] → 1MB binned data    
Segmented CNV → [04_bins_to_arm_level.py] → Arm-level aggregation    
results/Cohort_data/        

#### Requirements
Python Dependencies - All dependencies are listed in requirements.txt:

bash - Install from repository root    
pip install -r requirements.txt

