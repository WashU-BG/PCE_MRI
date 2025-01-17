# PCE_MRI
Contains code for the paper "Prenatal cannabis exposure, the brain, and psychopathology during early adolescence"  [https://doi.org/10.1101/2023.09.19.23295792](https://doi.org/10.1038/s44220-024-00281-7) 

## 1. Setup data for analysis
1. clean_merge_v5.R
2. prenatal_mri_setup.R

## 2. Run analyses on cluster
1. pce_allbrain_regressions_5b.R - Associations btw PCE and MRI
2. pce_cbcl_regressions_5.R - Associations btw MRI and psychopathology
3. pce_mediation_analysis.R - Mediation analyses - PCE --> MRI --> psychopathology

## 3. Parse/plot
1. parse_prenatal_output.R - Compute corrected p-values
2. prenatal_mri_figure.R - Figueres 2 & 3
3. pce_manhattan_plot2.R - Figure 1
