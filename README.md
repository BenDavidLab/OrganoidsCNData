# Organoid Stability Analysis Pipeline
Complete computational pipeline for analyzing copy number variation (CNV) discordance between primary tumors and patient-derived models (organoids and xenografts).

### Overview
This repository contains the complete analysis pipeline used in our manuscript, from raw sequencing data (FASTQ) to publication-quality figures. The pipeline consists of two main components:

Pre-Processing: FASTQ → BAM → CNV → Segments → Bins → Arm-level normalization
Post-Processing: CNV discordance analysis and visualization

### Quick Start
#### Prerequisites
1- Python 3.8 or higher
2- Required data:
    Genome reference files in data/reference_files/
    Cohort metadata in data/metadata_organoids.csv
3- For complete analysis: Raw sequencing data (FASTQ files)

### Documentation
- Pre-processing README: FASTQ → CNV pipeline details
- Post-processing README: Analysis and visualization details
- Data README: Data structure and requirements

### Citation
If you use this pipeline in your research, please cite:
[Your manuscript citation will go here]

### License
This project is licensed under the MIT License - see the LICENSE file for details.

### Contact 
For questions about the analysis or code:
[fill emails]

### Acknowledgments 
All contributors and collaborators listed in the manuscript [add paper url]
