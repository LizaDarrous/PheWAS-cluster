### 1) Retrieve IVs for exposure of interest and run PheWAS loclaly to obtain SNPxTraits association matrix
### also retrieves the outcome SNP effect for later MR analysis
library(tidyverse)
library(data.table)
library(dplyr)
library(tictoc)

# variable set-up
EXP_dir = "/data/sgg3/data/neale_files/both_sexes/body/21001_irnt.gwas.imputed_v3.both_sexes.tsv.gz"
EXP_pheno = "21001"
OUT_dir = "/data/sgg3/data/neale_files/both_sexes/lifestyle/845.gwas.imputed_v3.both_sexes.tsv.gz"
OUT_pheno = "845"

VARinfo_dir = "/data/sgg2/liza/SEM_Real/test/data/variants.tsv.bgz"
UKBB_info_both_dir = "/data/sgg3/liza/TraitCorr/pheWAS/bash_trial/phenotypes.both_sexes.tsv.bgz"
UKBB_summ_both_dir = "/data/sgg3/data/neale_files/both_sexes/"
fpaths_fil_dir = "/data/sgg3/liza/TraitCorr/pheWAS/data/fpaths_fil.txt"  ## need to create this prior
res_dir = paste0("/data/sgg3/liza/TraitCorr/pheWAS/",EXP_pheno)
bash_pheWAS = "/data/sgg3/liza/TraitCorr/pheWAS/scripts/rslurm_bash_pheWAS.sh"
system(paste0("mkdir ",res_dir))
system(paste0("mkdir ",res_dir,"/data"))

MR_pval=2*pnorm(-abs(5.45)) # MR significance thresholds
n_thresh = 50000 # sample size threshold

# read in exposure trait
X = fread(EXP_dir)

# read in VARinfo to get RSID
VARinfo = fread(cmd = paste0(" zcat < ",VARinfo_dir), sep="\t", header=TRUE) ## Variant info file from UKBB GWAS Imputed v3 - File Manifest Release 20180731

# add in data from VARinfo including rsid, alt, ref, chr, pos, by=variant - UKBB specific
X1 = dplyr::inner_join(X, VARinfo[,c(1:6)])
X1$chr = as.numeric(X1$chr)
X1 = X1[!is.na(X1$chr),] # keep only autosomal chromosomes
print(paste0("Dimension of raw autosomal SNP matrix: ", dim(X1)[1],", ",dim(X1)[2]))

# functions to use
# get SNPs significant over a certain p-value / Z-threshold
prune_X = function(zX,p_limit=1e-5){
  zX=zX
  z_limit=abs(qnorm(0.5*p_limit))
  ind_keep=which(abs(zX)>z_limit)
  ind_keep=unique(ind_keep)
  ind_keep=list(ind_keep)
  return(ind_keep)
}

# Taken from Jonathan Sulc to create bins that fit a maximum of 50k SNPs (max threshold for clumping)
snp_bin  =  function( snp_ranks, chunk_size = SNP_CHUNK_SIZE ){
  if (nrow( snp_ranks ) == 0) {
    return()
  }
  
  max_chr  =  snp_ranks$chr %>%
    table %>%
    cumsum %>%
    (function(x) x < chunk_size) %>%
    (function(x) names(x)[ max(which(x)) ] ) %>%
    as.numeric
  if (is.na( max_chr )) {
    max_chr = min( snp_ranks$chr )
  }
  
  bin = snp_ranks %>%
    filter( chr <= max_chr ) %>%
    list
  return( c( bin,
             snp_bin( snp_ranks[ snp_ranks$chr > max_chr, ],
                      chunk_size ) ) )
}

# prune SNPs before formatting into TwoSampleMR data to reduce lag
mr_ind=unlist(prune_X(X1$tstat,MR_pval))
if(length(mr_ind)==0){
  print("No IVs available for this threshold, please choose another and re-run")
  stop()
}
X1= X1[mr_ind,]
print(paste0("Dimension of signficant SNP matrix: ", dim(X1)[1],", ",dim(X1)[2]))

# data formatting - sig IVs and clumped
mr_dataX = cbind.data.frame(SNP = X1$rsid, variant = X1$variant, beta = X1$beta, se = X1$se, pval.exposure = X1$pval, effect_allele = X1$alt, other_allele = X1$ref, chr=X1$chr, tstat=X1$tstat, N=X1$n_complete_samples)

clump_bin = snp_bin(mr_dataX,50000)
exp_data = c()
for (x in 1:length(clump_bin)) {
  temp = mr_dataX[mr_dataX$SNP %in% clump_bin[[x]]$SNP,]
  temp1 = TwoSampleMR::clump_data(temp)
  exp_data=rbind(exp_data,temp1)
}
print(paste0("Dimension of signficant and pruned SNP matrix: ", dim(exp_data)[1],", ",dim(exp_data)[2]))
# save sig SNPs for X
write.csv(exp_data, paste0(res_dir,"/sig-clumped-IVs_", EXP_pheno, ".csv"), row.names = FALSE)

## since UKBB Neale files have the same order, we can use the same indices to pull out PheWA
# get row indices from VARinfo
SNP_ind = which(VARinfo$rsid %in% exp_data$SNP) #in the order of VARinfo
SNP_ind1 = c(1, SNP_ind + 1) # have to do this to get the header. bash also reads header as ind 1 so +1 to all other indices
write.table(SNP_ind1, paste0(res_dir,"/",EXP_pheno,"-IV_ind.txt"), row.names = F, col.names = F)

# before getting the pheWAS for all traits in directory, filter out traits with low N (speeds up)
fpaths_fil = fread(fpaths_fil_dir, header = FALSE)
fpaths_fil1 = unlist(fpaths_fil) %>%
  strsplit( "/" ) %>%
  sapply( tail, 1 ) %>%
  gsub(".gz", "",.)

nSNP = length(SNP_ind)
ntrait = length(fpaths_fil1)
trait_info = matrix(NA, nrow = ntrait, ncol = 9)
UKBB_info_both = data.table::fread(cmd = paste0("zcat < ",UKBB_info_both_dir))
## loop over files and get info needed
for(i in 1:ntrait){
  trt = strsplit(fpaths_fil1[i],"\\.")[[1]][1] #get trait code
  trait_info[i,1] = trt
  trait_info[i,2] = NA  #same as n_non_missing
  if(length(unlist(UKBB_info_both[which(UKBB_info_both$phenotype==trt),2:8]))>0){
    trait_info[i,3:9] = unlist(UKBB_info_both[which(UKBB_info_both$phenotype==trt),2:8])}
}

trait_info = as.data.frame(trait_info, stringsAsFactors = F)
cols.num <- c("V2","V6","V7","V8","V9")
trait_info[cols.num] <- sapply(trait_info[cols.num],as.numeric)

colnames(trait_info) = c("phenotype", "n_complete_samples", "description","variable_type",
                         "source","n_non_missing","n_missing", "n_controls","n_cases")
write.csv(trait_info, paste0(res_dir,"/trait_info_all.csv"), row.names = F)

# removing traits with N_eff < 50k
trait_info$n_eff = trait_info$n_non_missing
ncont_ind = which(!is.na(trait_info$n_cases))
trait_info$n_eff[ncont_ind] = (4*trait_info$n_cases[ncont_ind]*trait_info$n_controls[ncont_ind]) / (trait_info$n_cases[ncont_ind]+trait_info$n_controls[ncont_ind]) # (4*(control*case))/totN
small_n = which(trait_info$n_eff < n_thresh)
trait_info = trait_info[-small_n,]
print(paste0("Number of traits remaining after filtering for small sample size (<",n_thresh,"): ", dim(trait_info)[1]))
write.csv(trait_info, paste0(res_dir,"/trait_info_nfil.csv"), row.names = F)

fpaths_fil_nfil = fpaths_fil[-small_n]
write.table(fpaths_fil_nfil, file = paste0(res_dir,"/data/fpaths_fil_nfil.txt"), sep = "\t", row.names = FALSE, col.names = FALSE)


## rslurm the local PheWAS process - bash script
library(rslurm)

get_job_status2 <- function (slr_job)
{
  if (!(class(slr_job) == "slurm_job"))
    stop("input must be a slurm_job")
  stat <- suppressWarnings(system(paste("squeue -n", slr_job$jobname),
                                  intern = TRUE))
  if (length(stat) > 1) {
    res = "Job running or in queue."
    completed = FALSE
  }
  else {
    res = "Job completed or stopped."
    completed = TRUE
  }
  return(completed)
}

setwd(res_dir)
system(paste0("mkdir ", res_dir, "/pheWAS"))

rslurm_bash_pheWAS <- function(file_path, SNPind_path, out_path){
  system(paste("bash", bash_pheWAS, file_path, SNPind_path, out_path))
}

pars <- data.frame(paste0(UKBB_summ_both_dir,unlist(fpaths_fil_nfil)),
                   paste0(res_dir,"/",EXP_pheno,"-IV_ind.txt"),
                   paste0(res_dir,"/pheWAS"))
colnames(pars) = c("file_path", "SNPind_path", "out_path")
tic()
sjob = slurm_apply(f = rslurm_bash_pheWAS, params = pars, jobname = paste0("pheWAS-",EXP_pheno), nodes = nrow(pars), cpus_per_node = 1, global_objects = c("bash_pheWAS"), slurm_options = list(partition = "sgg", account = "sgg", mem = "5000"), submit = TRUE)
wait_counter = 0
while (wait_counter < 1) {
  wait_counter = 0
  if(get_job_status2(sjob)){
    wait_counter = wait_counter + 1
  } else{
    wait_counter = wait_counter
  }
}
toc()

# get data needed from PheWAS
fpaths_fil_nfil1 = unlist(fpaths_fil_nfil) %>%
  strsplit( "/" ) %>%
  sapply( tail, 1 ) %>%
  gsub(".gz", "",.)

ntrait = length(fpaths_fil_nfil1)

# make all the matrixes I want with dim SNP vs traits
unstdBeta_df = matrix(NA, nrow = nSNP, ncol = ntrait)
unstdSE_df = matrix(NA, nrow = nSNP, ncol = ntrait)
tstat_df = matrix(NA, nrow = nSNP, ncol = ntrait)
pval_df = matrix(NA, nrow = nSNP, ncol = ntrait)

## loop over files and get info needed
for(i in 1:ntrait){
  tmp = fread(paste0(res_dir,"/pheWAS/",fpaths_fil_nfil1[i]))
  if(all(VARinfo$variant[SNP_ind]==tmp$variant)){
    unstdBeta_df[,i] = tmp$beta
    unstdSE_df[,i] = tmp$se
    tstat_df[,i] = tmp$tstat
    pval_df[,i] = tmp$pval
  }
}

colnames(unstdBeta_df) = colnames(unstdSE_df) = colnames(tstat_df) = colnames(pval_df) = trait_info$phenotype
rownames(unstdBeta_df) = rownames(unstdSE_df) = rownames(tstat_df) = rownames(pval_df) = VARinfo$rsid[SNP_ind]

write.csv(unstdBeta_df, paste0(res_dir,"/unstdBeta_df.csv"), row.names = T)
write.csv(unstdSE_df, paste0(res_dir,"/unstdSE_df.csv"), row.names = T)
write.csv(tstat_df, paste0(res_dir,"/tstat_df.csv"), row.names = T)
write.csv(pval_df, paste0(res_dir,"/pval_df.csv"), row.names = T)

### save exp and outcome data
Y = fread(OUT_dir)
Y = fread(cmd = paste0(" zcat < ",OUT_dir), sep="\t", header=TRUE)
#Y$variant = gsub("_",":",Y$MarkerName)
# exp is already pruned
Y1 = dplyr::inner_join(exp_data[,c(1,2)], Y)
Y1 = dplyr::inner_join(Y1, VARinfo[,c(1:6)])

# Data set up - sig IVs and clumped
mr_dataY = cbind.data.frame(SNP = Y1$SNP, variant = Y1$variant, beta = Y1$beta, se = Y1$se, pval.outcome = Y1$pval, effect_allele = Y1$alt, other_allele = Y1$ref, chr=Y1$chr, tstat=Y1$tstat, N=Y1$n_complete_samples)

write.csv(mr_dataY, paste0(res_dir,"/clumped-IVs_",OUT_pheno,".csv"), row.names = FALSE)
