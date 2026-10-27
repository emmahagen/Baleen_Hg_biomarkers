# Longitudinal reconstruction of mercury exposure and multiple biomarkers of physiology in baleen

Hagen EE, Fernández Ajó A, Reed J, Rice JM, Wiley VA, Buck CL

Code and data associated with the manuscript of the same title. 

Corresponding Author: Emma E. Hagen; Status: Under review at *Environmental Pollution*. Preprint available at SSRI: http://dx.doi.org/10.2139/ssrn.7674826.

## Overview

This repository contains the R code and data used to analyze longitudinal total mercury (THg), hormone, and stable isotope profiles along baleen plates from five subadult eastern North Pacific gray whales (*Eschrichtius robustus*) stranded along the Oregon coast, USA (2019–2021).
These analyses evaluate the potential of longitudinal baleen profiles as a framework for multi-biomarker integration and the reconstruction of antemortem analyte trajectories. 
The provided code supports three study aims:
1. Characterize patterns of THg exposure over time.
2. Examine how THg dynamics relate to previously characterized endocrine and stable isotope profiles.
3. Develop and assess an analytical approach to integrate multiple biomarkers derived from baleen.

## Files

Files and contents:  
`scripts/Supplementary_Code_Hagen_et_al.R`: Statistical analyses and diagnostic and trend plots  
`data/Supplementary_Data_Hagen_et_al.csv`: Measurements and specimen information for 110 baleen samples  
`.gitignore`: Excludes local R session files and generated output  

The script evaluates longitudinal THg patterns using linear mixed-effects models and individual linear models, compares 32 candidate biomarker models by AIC, and fits additive models for THg and cortisol with and without Er_5. Er_5 was associated with acute mortality from predation; causes of death for Er_1–Er_4 were undetermined.

## Running the analysis

Analyses were conducted in R version 4.5.0 and nlme 3.1-168. The commands below install available package versions rather than recreating our original software environment.

Install the required packages once in R:

```r
install.packages(c("nlme", "ggplot2", "dplyr", "broom", "mgcv",
                   "itsadug", "gratia", "performance"),
                 repos = "https://cloud.r-project.org")
```

Use a current dplyr version with `reframe()` and a gratia version supporting `derivatives(..., select = ...)` and the `.derivative`, `.se`, and `cm` output columns used by the script.

Start a fresh R session and set the working directory to the repository folder containing this README. 
Run:

```r
set.seed(1)  # Fixed seed for stochastic diagnostic routines; not an original study seed.
source("scripts/Supplementary_Code_Hagen_et_al.R",
       echo = TRUE, max.deparse.length = Inf, encoding = "UTF-8")
sessionInfo()
```

Model summaries and comparisons are printed in the console; plots appear on the active graphics device. Models, tables, and plot objects remain in the R session. The script does not export manuscript tables or assemble all figures. 
To save the text results and plots as files, run this command in the terminal from the repository folder:

```sh
Rscript --vanilla -e 'set.seed(1); source("scripts/Supplementary_Code_Hagen_et_al.R", echo = TRUE, max.deparse.length = Inf, encoding = "UTF-8"); sessionInfo()' > analysis.log 2>&1
```

Per R's default settings, plots are saved to `Rplots.pdf`.

## Data

Each row represents one sampling position. Total sample counts for Er_1–Er_5 are 12, 26, 20, 25, and 27, respectively. The provided CSV contains no missing values or duplicate per-individual positions.

CSV columns and descriptions:

`Whale_ID`: Official specimen identifier  
`ID`: Individual identifier used in the analysis (Er_1–Er_5)  
`cm`: Position along the baleen plate in cm; values represent distance from the gumline (i.e., most recent growth)  
`THg_negs_to_zero_extract_eff_corrected`: Extraction-efficiency-corrected THg, with negative values set to zero (µg/g dry weight)  
`Strand_Date`: Stranding date, supplied as month/day/two-digit year  
`Sex`: Female or Male  
`BGR_mm_per_week`: Estimated baleen growth rate (mm/week)  
`TL_cm`: Total body length from snout to fluke notch (cm)  
`log_Cortisol`: Log1p-transformed cortisol concentration; original units ng/g  
`log_T3`: Log1p-transformed triiodothyronine (T3) concentration; original units ng/g  
`dN15`: Nitrogen stable isotope ratio (δ15N) (‰)  
`dC13`: Carbon stable isotope ratio (δ13C) (‰)  

The data include extraction-efficiency corrected THg and log1p-transformed hormone values, not raw instrument measurements. The script assumes cortisol values were transformed using log1p, calculates log1p(THg), and reverses baleen position within each whale using cm_countdown = max(cm) - cm. Sample positions start at 1 cm; there are no 0-cm samples. Er_1 lacks the most recent proximal baleen, as described in the manuscript.

## Citation

If you use this code or data, please cite the manuscript:

Hagen, E.E., et al. Longitudinal reconstruction of mercury exposure and multiple biomarkers of physiology in baleen. [Journal, year, DOI to be added upon publication]

## Reference for original hormone and stable isotope analyses

Fernández Ajó, A., Teixeira, C., de Mello, D. M. D., Dillon, D., Rice, J. M., Buck, C. L., Hunt, K. E., Rogers, M. C., & Torres, L. G. (2024). A longitudinal study of endocrinology and foraging ecology of subadult gray whales prior to death based on baleen analysis. General and Comparative Endocrinology, 352, 114492. https://doi.org/10.1016/j.ygcen.2024.114492

## Contact

Emma E. Hagen - eeh248@nau.edu, Department of Biological Sciences, Northern Arizona University
