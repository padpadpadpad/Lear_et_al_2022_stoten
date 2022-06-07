# Analyses and data of:

_Lear et al. (2022) Bacterial colonisation dynamics of household plastics in a coastal environment. Science of the Total Environment_. 

DOI of paper: [doi.org/10.1016/j.scitotenv.2022.156199](https://www.sciencedirect.com/science/article/pii/S004896972203296X)

ENA Accession number for sequencing: PRJEB5334

### Outline

This repository contains the final datasets, analyses and figures of the above-mentioned paper. It can recreate all figures and tables in both the main text and the supplement apart from Tables 1 and 2.

The code for running the bioinformatics pipelines for the amplicon 16S sequencing and the isolate whole genome sequencing are not included in the repository. These are available upon request. The raw sequencing data (for the 16S sequencing) and the genome assemblies (isolate WGS) are available on ENA here. They can be downloaded using [enaBrowserTools](https://github.com/enasequence/enaBrowserTools) (other tools are likely available).

### Feedback

- Please report any problems or bugs in the code in the [Issues](https://github.com/padpadpadpad/Lear_et_al_2022_stoten/issues) tab of the GitHub repository. Alternatively, please email _d.padfield@exeter.ac.uk_.

### Licensing

This code is licensed under GPL-3.

### Running the scripts and analyses

- The project can be `cloned` or for those not familiar with GitHub, a zip file of this project can be downloaded using the "Clone or download" button at the top right of this page.
- All of the scripts for the analyses can be found in `scripts/`.
- Before running any of the R scripts, open the R project file in the downloaded folder. [R projects](https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects) automatically assigns the root directory to the directory in which the project resides. Consequently all of the analyses should be runnable without altering paths. These are very easy to open using RStudio. 
- `scripts/make_map.R` makes the map shown in Figure 1.
- `scripts/analyse_abundance.R` analyses the abundance data of the location and colonisation experiments and recreates Figure 2, Figure 3, and Table S1.
- `scripts/extra_functions.R` contains extra functions used in the analysis of the amplicon 16S data.
- `scripts/prevalence_filtering.R` runs through the filtering steps used on the original raw output of the amplicon 16S analysis. The raw object is `phyloseq_16s.rds`.
- `scripts/amplicon_16s_diversity.R` analyses estimates of alpha diversiy and pielou's evenness from the amplicon 16s data and recreates Figure 4, Table S2, and Table S3.
- `scripts/amplicon_16s_analysis.R` analyses differences in community composition from the amplicon 16s data and recreates Figure 5, Figure S2, Figure S3, and Figure 4.
- `scripts/analyse_virulence.R` analyses the survival data and recreates Figure 6 and Table S5. For this script, installing the survival branch of **rstanarm** can be difficult so please see this [GitHub Issue](https://github.com/stan-dev/rstanarm/issues/426) for how I installed it successfully most recently.
