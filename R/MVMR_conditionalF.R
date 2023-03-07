#### 5b) Running MVMR on traits that have survived the stepwise MVMR, have a significant uni-MR effect on outcome, 
#### and provide the main exposure a conditional F-stat > 8

## Repurposed from bGWAS package: https://github.com/n-mounier/bGWAS/blob/master/R/makeMR_ZMatrix.R and 
## https://github.com/n-mounier/bGWAS/blob/master/R/identify_StudiesMR.R

library(data.table)
library(tidyverse)
library("readxl")
library("ggrepel")


### Useful small function (re-used by other functions)
update_log <- function(log_obj, text, verbose=F){
  log_obj = c(log_obj, text)
  if(verbose) cat(text, sep = "")
  return(log_obj)
}

get_names <- function(Files, Z_matrices = "~/ZMatrices/") {
  Studies = readr::read_tsv(file.path(Z_matrices, "AvailableStudies.tsv"), progress = FALSE, col_types = readr::cols())
  if(any(!Files %in% Studies$File)) stop("Some files are not part of Prior GWASs")
  Studies %>%
    slice(match(Files, .data$File)) %>%
    pull(.data$Name)  -> Names
  return(Names)
}

# Function to automatically get subsetted matrix of instruments
get_Instruments <- function(Zmat, Zlim){
  # get exposures
  Zmat %>%
    names() %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    slice(6:(ncol(Zmat)-1)) %>%
    pull() -> exposures
  
  # subset : any row, with abs(exposure)>Zlimit
  Zmat %>%
    filter_at(vars(exposures), any_vars(abs(.data$.)>Zlim)) -> Zmat
  return(Zmat)
}

# Function to automatically generate the formula for linear model
generate_Formula <- function(outcome, study_names, with_intercept=F ) {
  formula = ifelse(with_intercept,
                   paste0('`',outcome,'`', ' ~  1 + ',paste(collapse=' + ', paste(sep='','`',study_names,'`'))),
                   paste0('`',outcome,'`', ' ~ -1 + ',paste(collapse=' + ', paste(sep='','`',study_names,'`'))))
  #paste(outcome, ' ~  1 + ',paste(collapse=' + ', paste(sep='','`',study_names,'`'))),
  #paste(outcome, ' ~ -1 + ',paste(collapse=' + ', paste(sep='','`',study_names,'`'))))
  return(formula)
}

# variable set up
Z_matrices = "/Users/ldarrous/Desktop/UNIL/TraitCorr/PheWAS/21001/MVMR_21001-845/" #directory where ZMatrices and AvailableStudies are found
ZMatrix = fread("/Users/ldarrous/Desktop/UNIL/TraitCorr/PheWAS/21001/MVMR_21001-845/ZMatrix_confounder.csv")
trait_info = fread("/Users/ldarrous/Desktop/UNIL/TraitCorr/PheWAS/21001/trait_info_nfil.csv")
MR_threshold = 5e-8 #1e-5
MR_ninstruments = 3
#MR_pruning_dist = 500
#MR_pruning_LD = 0
MR_shrinkage = 0.05
stepwise_threshold = NULL
#prior_shrinkage = NULL
sign_method = "p"
sign_thresh = 5e-8
use_permutations= FALSE
res_pruning_dist = 500
res_pruning_LD = 0
save_files = FALSE
verbose = FALSE
Zlimit = stats::qnorm(MR_threshold/2, lower.tail = F)
source("/Users/ldarrous/Desktop/UNIL/TraitCorr/PheWAS/bash_trial/bGWAS/identify_studiesMR.R")

# try once without BMI and then add it back in
ZMatrix1 = ZMatrix[,-"21001_irnt"]
MVMR_noBMI = identify_studiesMR(ZMatrix1, MR_shrinkage, MR_threshold, stepwise_threshold, Z_matrices, save_files=FALSE, verbose=TRUE)
MVMR_noBMI$studies

# make a sub ZMatrix
ZMatrix %>%
  select(c("rs","chrm","pos","alt",'ref',MVMR_noBMI$studies,"21001_irnt","845")) %>%  ## add in BMI and outcome
  get_Instruments(Zlim=Zlimit) -> ZMatrix_subset
### UPDATE Z-MATRIX
myFormula = generate_Formula("845", c(MVMR_noBMI$studies,"21001_irnt")) #generate_Formula("845", c(MVMR_noBMI$studies,"21001_irnt")) #1687

## RUN MODEL - with studies already included and BMI
model = stats::lm(data=ZMatrix_subset, formula = myFormula)
stats::lm(data=ZMatrix_subset, formula = myFormula) %>%
  summary %>%
  stats::coef() %>%
  as.data.frame() %>% # 1st create a data.frame (to get rownames)
  tibble::rownames_to_column("nm") %>% # then convert to tibble
  as_tibble()  -> coefs
## ADDED LIZA
coefs$nm = gsub('`','',coefs$nm)
coefs %>%
  set_names(c("study", "estimate", "std_error", "Tstat", "P")) -> coefs
  
cat("Excluding BMI: \n")
inner_join(MVMR_noBMI$coeffs, trait_info[,c(1,3)], by=c("study"="phenotype"))
cat("Including BMI: \n")
inner_join(coefs, trait_info[,c(1,3)], by=c("study"="phenotype"))

##### testing for conditional F statistic 

## running univ MR on each of the exposures 
## get the uni_MR from bi_MR ran in step #4
EDU_biMR = fread("/Users/ldarrous/Desktop/UNIL/TraitCorr/PheWAS/bash_trial/LocalEpi/localEpi_EDU_407.csv")

grand_mr = as.data.frame(matrix(NA, nrow = length(MVMR_noBMI$studies), ncol = 5))
colnames(grand_mr) = c("trait","axy","se","pval","method")

for(i in 1:length(MVMR_noBMI$studies)){
  #print("enter")
  exp = MVMR_noBMI$studies[i]
  grand_mr$trait[i] = exp
  tab = EDU_biMR[which(EDU_biMR$exposure==exp),]
  if(nrow(tab)>1){
    #print(exp)
    tab = tab[-(which(tab$method=="MR Egger")),]
    sig_tab = tab[order(tab$pval)[2],]
    grand_mr$axy[i] = sig_tab$b
    grand_mr$se[i] = sig_tab$se
    grand_mr$pval[i] = sig_tab$pval
    grand_mr$method[i] = sig_tab$method
  }else if(nrow(tab)==1){
    #print(i)
    grand_mr$axy[i] = tab$b
    grand_mr$se[i] = tab$se
    grand_mr$pval[i] = tab$pval
    grand_mr$method[i] = tab$method
  }
}

grand_mr_sort = grand_mr[order(grand_mr$pval),]
print("The univariate MR estimates in order of signifcance: ")
#grand_mr_sort
inner_join(grand_mr_sort, trait_info[,c(1,3)], by=c("trait"="phenotype"))

stdB_matrix = fread("/Users/ldarrous/Desktop/UNIL/TraitCorr/PheWAS/21001/MVMR_21001-845/clumped_matrixBsted.csv")
stdSE_matrix = fread("/Users/ldarrous/Desktop/UNIL/TraitCorr/PheWAS/21001/MVMR_21001-845/clumped_matrixSEstd.csv")
Pval_matrix = fread("/Users/ldarrous/Desktop/UNIL/TraitCorr/PheWAS/21001/MVMR_21001-845/clumped_matrixP.csv")

res <- Map(combn, list(c(1:nrow(grand_mr_sort))), seq_along(c(1:nrow(grand_mr_sort))), simplify = FALSE)
med_comb = unlist(res, recursive = FALSE)
Fstat_arr = c()
med_arr = c()
axy_arr = c()
se_arr = c()
pval_arr = c()
MVMR_list = list()
counter = 1
for(med in med_comb){

  ZMatrix_subset %>%
    select(c("rs","chrm","pos","alt",'ref', grand_mr_sort$trait[c(med)], "21001_irnt","845")) %>%
    get_Instruments(Zlim=Zlimit) -> ZMatrix_subset_fcond
  
  stdB_matrix1 = stdB_matrix[which(stdB_matrix$rs %in% ZMatrix_subset_fcond$rs),]
  stdSE_matrix1 = stdSE_matrix[which(stdSE_matrix$rs %in% ZMatrix_subset_fcond$rs),]
  Pval_matrix1 = Pval_matrix[which(Pval_matrix$rs %in% ZMatrix_subset_fcond$rs),]
  all(stdB_matrix1$rs==stdSE_matrix1$rs)
  all(Pval_matrix1$rs==stdB_matrix1$rs)
  
  stdB_matrix1 %>%
    select(-c("rs","chrm","pos","alt",'ref',"845")) %>%
    select(c(grand_mr_sort$trait[c(med)], "21001_irnt"))-> X
  
  stdSE_matrix1 %>%
    select(-c("rs","chrm","pos","alt",'ref',"845")) %>%
    select(c(grand_mr_sort$trait[c(med)], "21001_irnt"))-> SE
  
  K = ncol(X)-1 # number of mediators
  L = nrow(X)
  
  E = as.vector(unlist(X[,"21001_irnt"]))
  M = as.matrix(X[,1:(ncol(X)-1)])
  
  delta.reg = lm(E ~ -1 + M)
  delta = delta.reg$coefficients
  res = delta.reg$residuals
  
  delta = as.vector(c(-1,delta))
  
  var_IV =  as.matrix(SE**2) %*% (delta**2)
  
  Q = sum(1/(var_IV) * res**2)
  
  Fstat = Q/(L-K)
  Fstat_arr=c(Fstat_arr,Fstat)
  med_arr = c(med_arr,stringr::str_c(grand_mr_sort$trait[med], collapse="-"))
  
  ## caluclate MVMR for that med comb
  
  myFormula = generate_Formula("845", c(grand_mr_sort$trait[c(med)],"21001_irnt")) #generate_Formula("845", c(MVMR_noBMI$studies,"21001_irnt")) #1687
  
  ## RUN MODEL - with studies already included and BMI
  #model = stats::lm(data=ZMatrix_subset_fcond, formula = myFormula)
  stats::lm(data=ZMatrix_subset_fcond, formula = myFormula) %>%
    summary %>%
    stats::coef() %>%
    as.data.frame() %>% # 1st create a data.frame (to get rownames)
    tibble::rownames_to_column("nm") %>% # then convert to tibble
    as_tibble()  -> coefs
  ## ADDED LIZA
  coefs$nm = gsub('`','',coefs$nm)
  coefs %>%
    set_names(c("study", "estimate", "std_error", "Tstat", "P")) -> coefs
  
  res = inner_join(coefs, trait_info[,c(1,3)], by=c("study"="phenotype"))
  MVMR_list[[counter]]=res
  
  axy_arr=c(axy_arr,unlist(res[nrow(res),2]))
  se_arr=c(se_arr,unlist(res[nrow(res),3]))
  pval_arr=c(pval_arr,unlist(res[nrow(res),5]))
  
  counter = counter + 1
}

print("BMI F-stat caluclated for different combinations of exposure traits")
MVMRcondF_df = cbind.data.frame(med_arr,Fstat_arr,axy_arr,se_arr,pval_arr)
## view different EXP conditional F-stat with different possibe combinations of candidate confounders
MVMRcondF_df

## select from the data frame above the row which shows the best combination of confounders that give a conditional Fstat > 8
MVMR_list[[5]]
