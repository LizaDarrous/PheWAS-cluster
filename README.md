# PheWAS-cluster

:warning: :grey\_exclamation: This repository contains **personalised scripts** that were used to run the full analysis. A **working example** is provided instead with instructions below to help you replicate and run your own analysis :grey\_exclamation:

## Overview

The scripts in 'PheWAS-cluster' represent our approach of Phe-WAS based clustering of Mendelian Randomisation instruments (PWC-MR). 
Our approach was used to investigate the large causal effect of body mass index (BMI) on educational attainment (EDU) -0.19 [-0.22, -0.16], where we hypothesise that potential horizontal pleiotropic effects (emerging due to heritable confounders, dynastic effects, genetic subtypes of obesity and other pleiotropic mechanisms, see **panel a** in the figure below) are biassing BMI's effect on educational attainment.

The main steps of the PWC-MR method are (illustrated in **panel b** of the figure below):
- Instrument selection and PheWAS
- IV clustering 
- Enrichment analysis and cluster specific MR
<p align="center">
<img src="misc/DAG_flowgram.jpg" align="center" height=420/>
</p>

## Working Example - BMI

Please download the `working-example` folder and set it as your working directory in R. There are two subfolders you would use, and one file to download from [here](https://drive.google.com/file/d/1KIwu8z2gBr616ZyNgOSQ8tK73AlbDTxL/view?usp=sharing) (150 MBs):
- **data**: this subfolder contains all the data files needed to run the example scripts (main analysis of PWC-MR), these are:
    - `unstdBeta_df.csv` / `unstdSE_df.csv` / `tstat_df.csv` / `pval_df.csv` these data frames contain the effect, standard error, t-statistic (beta/SE), and p-values respectively for the 348 genome wide significant BMI SNPs across 408 traits (filtered for sample size > 50'000). For other trait analysis, these data frames must be obtained either from UKBB or from PhenoScanner with similar filtering criteria (sample size, no duplicate traits).
    - `trait_info_nfil.csv` contains information on the 408 traits including: trait, description, effective sample size, variable_type...
    - `fpaths_fil_nfil.txt` contains the file paths of the traits used, this file is needed simple to remove duplicate traits if they have multiple versions (UKBB artefact).
    - `sig-clumped-IVs_21001.csv` contains the genome-wide significant and clumped BMI SNPs/IVs used for the TwoSampleMR causal effect estimate using all SNPs and the various clustered SNPs. The columns needed/included are: SNP, variant, beta, se, pval.exposure, effect_allele, other_allele, chr, tstat, N.
    - `clumped-IVs_845.csv` contains the effects of the same BMI SNPs/IVs but for the outcome of interest, in this case it is **EDU** (Age completed full time education). The columns needed/included are: SNP, variant, beta, se, pval, effect_allele, other_allele, chr, tstat, N.
    
- **scripts**: this subfolder has 2 main scripts numbered in order of use. In both of these scripts, the `# variable set-up` section in the begining should be updated to include the correct path for the `working-example` directory in the variable `res_dir`.  
    The scripts are: 
    - `1_QC_filtering.R` This script reads in `unstdBeta_df.csv`, `unstdSE_df.csv`, `tstat_df.csv`, `pval_df.csv`, `trait_info_nfil.csv`, and `fpaths_fil_nfil.txt` from the **data** subfolder. It then proceeds to filter out traits that have NA effects, duplicate traits (specifically exposure), traits with an exposure-genetic correlation > 0.75 (can be changed).
    
      The remaining SNP-trait effect matrix is then standardised, and SNPs are further removed if they are more strongly associated with traits other than the exposure.
      
      Lastly, all the variables are saved into an '.RData' file in the main directory called `QCdata_21001.Rdata` to be used in the second script.
    - `2_Clumping.R` This script reads in the previously created `QCdata_21001.Rdata` as well as the SNP effects for the exposure and outcome traits to be used in TwoSampleMR; `sig-clumped-IVs_21001.csv` and `clumped-IVs_845.csv`.
      
      The script then normalises the absolute value of the SNPxTrait effect matrix by row (SNP), and proceeds to run K-means clustering on the matrix after determining what the best number of clusters are (ranging from 2 to 50) using the AIC score. Then, the SNPs in each cluster are used to estimate a causal effect estimate on the outcome, as well as the SNPs altogether.
    
      Lastly, an enrichment ratio is calculated for each trait across all the clusters, and then the top 10 enriched traits for each cluster are written into an output file. 
      
      _There are several outputs from this script including plots for the AIC score of cluster numbers ranging from 2 to 50, SNP allocation into various clusters, MR estimates for the various clusters in .csv and .pdf format. Top 10 enriched traits for each cluster are also output in both .csv and .pdf format._
    
 - `Hail_AllxAll.csv` is a data frame containing the genetic and phenotypic correlation of multiple UKBB traits, downloaded from Neale's lab [here](https://ukbb-rg.hail.is/rg_browser/) in October 2021. 

------------ 
This analysis was run entirely on _R version 4.2.2 (2022-10-31_ as well as _R version 4.1.3 (2022-03-10)_.  
Analysis run time takes on average 15 minutes.

Secondary analysis is provided as non-custom scripts (systematic confounder search using MR followed by MVMR), that can be modified to include proper paths for analysis. 
    
