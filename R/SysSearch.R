### Systematic search for confounders
## needs 35 GBs of space on the server
library(data.table)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggforce)
library(ggrepel)
library(writexl)
library(TwoSampleMR)

#variables
EXP_pheno = "21001" #for grepl purposes, do a check to make sure not a lot of traits in trait_info start with this
OUT_pheno = "845"
EXP_dir = "/data/sgg3/data/neale_files/both_sexes/body/21001_irnt.gwas.imputed_v3.both_sexes.tsv.gz"
OUT_dir = "/data/sgg3/data/neale_files/both_sexes/lifestyle/845.gwas.imputed_v3.both_sexes.tsv.gz"
VARinfo_dir = "/data/sgg2/liza/SEM_Real/test/data/variants.tsv.bgz"
UKBB_info_both_dir = "/data/sgg3/liza/TraitCorr/pheWAS/bash_trial/phenotypes.both_sexes.tsv.bgz"
UKBB_root_dir = "/data/sgg3/data/neale_files/both_sexes/"
res_dir = paste0("/data/sgg3/liza/TraitCorr/pheWAS/",EXP_pheno)
fpaths_fil_dir = paste0(res_dir,"/data/fpaths_fil_nfil.txt")
EXP_irnt = TRUE  #if true change all grepl to paste0(EXP_pheno,"_irnt")
OUT_irnt = FALSE  #if true change all grepl to paste0(EXP_pheno,"_irnt")
hail_gcorr_dir = "/data/sgg3/liza/TraitCorr/pheWAS/data/Hail_AllxAll.csv"
out_gcorr_thresh = 0.75
bi_MR_thresh = 0.05
stepMVMR_R = "/data/sgg3/liza/TraitCorr/pheWAS/scripts/identify_studiesMR.R"

if(EXP_irnt){EXP_pheno2=paste0(EXP_pheno,"_irnt")}else{EXP_pheno2=EXP_pheno}
if(OUT_irnt){OUT_pheno2=paste0(OUT_pheno,"_irnt")}else{OUT_pheno2=OUT_pheno}

setwd(res_dir)
load(paste0(res_dir,"/QCdata_",EXP_pheno,".Rdata"))

## Load in results of bidirectional MR on BMI/EDU
mrEXP = fread(paste0(res_dir,"/biMR_", EXP_pheno, ".csv"))
mrOUT = fread(paste0(res_dir,"/biMR_", OUT_pheno, ".csv"))

print(paste0("Bi-directional EXP-MR calculated for ", length(unique(mrEXP$outcome))-1," traits with IVs significant at at least 5e-6."))
print(paste0("Bi-directional OUT-MR calculated for ", length(unique(mrOUT$outcome))-1," traits with IVs significant at at least 5e-6."))

## removing highly correlated traits with OUT (>0.75)  ## everything named 1 after this has no highly correlated edu traits
df_gcorr_OUT = hail_df[which(hail_df$ID1 == OUT_pheno | hail_df$ID2 == OUT_pheno),]
high_gcorr = which(abs(df_gcorr_OUT$rg)>=out_gcorr_thresh)
df_gcorr_OUT_high = df_gcorr_OUT[high_gcorr,] 
# HAIL doens't have _irnt, so compare against description
# also get rid of outcome in list of traits
trait_info_fil = trait_info_noEXP[-which(trait_info_noEXP$description %in% unique(df_gcorr_OUT_high$ph1)),] 

## make a table for all traits with their effect on EXP, EXP's effect on them, same for OUT(x2)
grand_mr = as.data.frame(matrix(NA, nrow = nrow(trait_info_fil), ncol = 18))
colnames(grand_mr) = c("Phenotype","Description","onEXP","onEXPse","onEXPpval","onEXPmethod","EXPon","EXPonse","EXPonpval","EXPonmethod","onOUT","onOUTse","onOUTpval","onOUTmethod","OUTon","OUTonse","OUTonpval","OUTonmethod")
grand_mr$Phenotype = trait_info_fil$phenotype
grand_mr$Description = trait_info_fil$description

## signifcant by at least 2 methods after removing MR-Egger
for(i in 1:nrow(grand_mr)){
  tab1 = mrEXP[which(mrEXP$exposure==grand_mr$Phenotype[i])]
  if(nrow(tab1)>1){
    tab1 = tab1[-(which(tab1$method=="MR Egger")),]
    sig_tab1 = tab1[order(tab1$pval)[2],]
    grand_mr$onEXP[i] = sig_tab1$b
    grand_mr$onEXPse[i] = sig_tab1$se
    grand_mr$onEXPpval[i] = sig_tab1$pval
    grand_mr$onEXPmethod[i] = sig_tab1$method
  }else if(nrow(tab1)==1){
    grand_mr$onEXP[i] = tab1$b
    grand_mr$onEXPse[i] = tab1$se
    grand_mr$onEXPpval[i] = tab1$pval
    grand_mr$onEXPmethod[i] = tab1$method
  }
  
  tab2 = mrEXP[which(mrEXP$outcome==grand_mr$Phenotype[i])]
  if(nrow(tab2)>1){
    tab2 = tab2[-(which(tab2$method=="MR Egger")),]
    sig_tab2 = tab2[order(tab2$pval)[2],]
    grand_mr$EXPon[i] = sig_tab2$b
    grand_mr$EXPonse[i] = sig_tab2$se
    grand_mr$EXPonpval[i] = sig_tab2$pval
    grand_mr$EXPonmethod[i] = sig_tab2$method
  }else if(nrow(tab2)==1){
    grand_mr$EXPon[i] = tab2$b
    grand_mr$EXPonse[i] = tab2$se
    grand_mr$EXPonpval[i] = tab2$pval
    grand_mr$EXPonmethod[i] = tab2$method
  }
  
  tab3 = mrOUT[which(mrOUT$exposure==grand_mr$Phenotype[i])]
  if(nrow(tab3)>1){
    tab3 = tab3[-(which(tab3$method=="MR Egger")),]
    sig_tab3 = tab3[order(tab3$pval)[2],]
    grand_mr$onOUT[i] = sig_tab3$b
    grand_mr$onOUTse[i] = sig_tab3$se
    grand_mr$onOUTpval[i] = sig_tab3$pval
    grand_mr$onOUTmethod[i] = sig_tab3$method
  }else if(nrow(tab3)==1){
    grand_mr$onOUT[i] = tab3$b
    grand_mr$onOUTse[i] = tab3$se
    grand_mr$onOUTpval[i] = tab3$pval
    grand_mr$onOUTmethod[i] = tab3$method
  }
  
  tab4 = mrOUT[which(mrOUT$outcome==grand_mr$Phenotype[i])]
  if(nrow(tab4)>1){
    tab4 = tab4[-(which(tab4$method=="MR Egger")),]
    sig_tab4 = tab4[order(tab4$pval)[2],]
    grand_mr$OUTon[i] = sig_tab4$b
    grand_mr$OUTonse[i] = sig_tab4$se
    grand_mr$OUTonpval[i] = sig_tab4$pval
    grand_mr$OUTonmethod[i] = sig_tab4$method
  }else if(nrow(tab4)==1){
    grand_mr$OUTon[i] = tab4$b
    grand_mr$OUTonse[i] = tab4$se
    grand_mr$OUTonpval[i] = tab4$pval
    grand_mr$OUTonmethod[i] = tab4$method
  }
}

table(c(grand_mr$EXPonmethod, grand_mr$onEXPmethod, grand_mr$OUTonmethod, grand_mr$onOUTmethod))
sum(table(c(grand_mr$EXPonmethod, grand_mr$onEXPmethod, grand_mr$OUTonmethod, grand_mr$onOUTmethod))
)

# not relying on absolute significance from p-value to categorise trait, instead dealing with more strongly significant directions to divide directly into 4 groups
EXP_tdiff = (abs(grand_mr$EXPon)-abs(grand_mr$onEXP))/sqrt(grand_mr$EXPonse^2 + grand_mr$onEXPse^2)
EXP_tpval = pnorm(EXP_tdiff)
onEXP = which(EXP_tpval<0.05)  # which have stronger onEXP
EXPon = which(EXP_tpval>0.95)  # which have stronger EXPon
EXPneut = which(EXP_tpval>=0.05 & EXP_tpval<=0.95)

OUT_tdiff = (abs(grand_mr$OUTon)-abs(grand_mr$onOUT))/sqrt(grand_mr$OUTonse^2 + grand_mr$onOUTse^2)
OUT_tpval = pnorm(OUT_tdiff)
onOUT = which(OUT_tpval<0.05)  # which have stronger onOUT
OUTon = which(OUT_tpval>0.95)  # which have stronger OUTon
OUTneut = which(OUT_tpval>=0.05 & OUT_tpval<=0.95)

colliderTraits = grand_mr[unique(Reduce(intersect, list(EXPon,OUTon))),]
confounderTraits = grand_mr[unique(Reduce(intersect, list(onEXP, onOUT))),]
mediatorTraits = grand_mr[unique(Reduce(intersect, list(EXPon,onOUT))),]

colliderTraits = colliderTraits[which(colliderTraits$EXPonpval<bi_MR_thresh & colliderTraits$OUTonpval<bi_MR_thresh),]
confounderTraits = confounderTraits[which(confounderTraits$onEXPpval<bi_MR_thresh & confounderTraits$onOUTpval<bi_MR_thresh),]
mediatorTraits = mediatorTraits[which(mediatorTraits$EXPonpval<bi_MR_thresh & mediatorTraits$onOUTpval<bi_MR_thresh),]

write_xlsx(list(
  "Collider" = as.data.frame(colliderTraits),
  "Confounder" = as.data.frame(confounderTraits),
  "Mediator" = as.data.frame(mediatorTraits)
), path = paste0(res_dir,"/",EXP_pheno,"-",OUT_pheno,"_MR_TraitSearch_3groups.xlsx"))

## testing stepwise-MVMR using bGWAS. Need to prepare files (Zmatrix)
MVMR_dir = paste0("/MVMR_",EXP_pheno,"-",OUT_pheno)
traits_MVMR = unique(confounderTraits$Phenotype)

# create AvailableStudies.tsv
AS = trait_info[trait_info$phenotype %in% c(traits_MVMR,EXP_pheno2,OUT_pheno2) ,c(1,3,5)]
AS %>% arrange(match(phenotype, c(traits_MVMR,EXP_pheno2,OUT_pheno2))) -> AS

AS_df = cbind.data.frame("File"=AS$phenotype, "Name"=AS$description, "ID"=1:nrow(AS), "Trait"=AS$description, "Consortium	Reference"=AS$source, "Download"=NA, "Remarks"=NA, "N_SNPs"=NA, "N_Instruments"=NA)

system(paste0("mkdir ", res_dir, MVMR_dir))
write.table(AS_df, paste0(res_dir, MVMR_dir, "/AvailableStudies.tsv"), row.names = FALSE, sep = "\t")

# get confounder fpaths
mr_fpaths = fread(paste0(res_dir,"/data/mr_fpaths.csv"))
mr_fpaths = rbind(mr_fpaths, strsplit(EXP_dir,UKBB_root_dir)[[1]][2], use.names=FALSE) # OUT already present

par_exps = stringr::str_match(mr_fpaths$V1, "/(.*?)\\.")[,2]
pars <- data.frame(par_exp = par_exps,
                   par_fpaths = mr_fpaths$V1)

mvmr_fpaths = pars[which(pars$par_exp %in% AS_df$File),]
mvmr_fpaths %>% arrange(match(par_exp, AS_df$File)) %>% select(par_fpaths) -> mvmr_fpaths #rearrange so that they match with AS
write.csv(mvmr_fpaths, paste0(res_dir,MVMR_dir,"/MVMRfpaths_conf.csv"), row.names = FALSE)

# 
#mvmr_fpaths = fread(paste0(res_dir,MVMR_dir,"/MVMRfpaths_conf.csv"))
VARinfo = fread(cmd = " zcat < /data/sgg2/liza/SEM_Real/test/data/variants.tsv.bgz", sep="\t", header=TRUE) ## Variant info file from UKBB GWAS Imputed v3 - File Manifest Release 20180731

## create Z matrix by looping over traits, obtaining their sig SNPs (Z/tstat) and performing an outer join
setwd(UKBB_root_dir)
library(stringr)
Z_thresh_gw = qnorm(0.5*(5E-8), lower.tail=FALSE)
Z_snps = c()
for(ind in 1:(nrow(mvmr_fpaths)-1)){ #no outcome trait
  print(ind)
  trait_nm = gsub(".gwas.imputed_v3.both_sexes.tsv.gz","",tail(str_split(mvmr_fpaths[ind],"\\/")[[1]],1))
  f1 = fread(unlist(mvmr_fpaths[ind]))
  Z_snps = c(Z_snps, unlist(f1[abs(f1$tstat)>Z_thresh_gw, 'variant']))
}
length(Z_snps)
#head(Z_snps)
Z_matrix = as.data.frame(unique(Z_snps))
colnames(Z_matrix) = "variant"

for(ind in 1:nrow(mvmr_fpaths)){ #adding in outcome trait
  print(ind)
  trait_nm = gsub(".gwas.imputed_v3.both_sexes.tsv.gz","",tail(str_split(mvmr_fpaths[ind],"\\/")[[1]],1))
  f1 = fread(unlist(mvmr_fpaths[ind]))
  Z_matrix = dplyr::inner_join(Z_matrix,f1[,c('variant','tstat')], by=c("variant"))
  colnames(Z_matrix)[ncol(Z_matrix)] = trait_nm
}

write.csv(Z_matrix, paste0(res_dir, MVMR_dir, "/raw_Zmatrix.csv"), row.names = FALSE)

## need to do a rank-based clumping
# this ranks the absolute Z value (descending order - minus), keeping NAs, and averaging the ranking of duplicate values
Z_matrix_noOut = Z_matrix[,-ncol(Z_matrix)]
rank_matrix = sapply(-(abs(Z_matrix_noOut[,-1])), rank, na.last="keep", ties.method = "average")
# replace NAs with large value, get row rank min and then normalise that
top_rank = max(rank_matrix, na.rm=T)*ncol(rank_matrix)
rank_matrix1 = replace(rank_matrix, is.na(rank_matrix), top_rank)
rank_min_rows = apply(rank_matrix1, 1, min)
rank_min_pvals = rank_min_rows/max(rank_min_rows)
dummy_df = cbind.data.frame("variant"=Z_matrix_noOut[,1], "pval"=rank_min_pvals)
dummy_df = inner_join(dummy_df, VARinfo[,c(1:6)])
dummy_df$chr = as.numeric(dummy_df$chr)
dummy_df = dummy_df[!is.na(dummy_df$chr),]
dim(dummy_df)

# Taken from Jonathan Sulc to create bins that fit a maximum of 50k SNPs (max threshold for clumping)
snp_bin  =  function( snp_ranks,
                      chunk_size = SNP_CHUNK_SIZE ){
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
    dplyr::filter( chr <= max_chr ) %>%
    list
  return( c( bin,
             snp_bin( snp_ranks[ snp_ranks$chr > max_chr, ],
                      chunk_size ) ) )
}
clump_bin = snp_bin(dummy_df,50000)

clumped_sig_snps = c()
for (x in 1:length(clump_bin)) {
  temp = dummy_df[dummy_df$rsid %in% clump_bin[[x]]$rsid,]
  temp1 = ieugwasr::ld_clump(temp, clump_r2 = 0.01, clump_kb = 5000) ##updated these to be different than BGWAS's
  clumped_sig_snps=rbind(clumped_sig_snps,temp1) # 693
}

Bmatrix_clumped = matrix(NA, nrow = nrow(clumped_sig_snps), ncol = dim(mvmr_fpaths)[1])
SEmatrix_clumped = matrix(NA, nrow = nrow(clumped_sig_snps), ncol = dim(mvmr_fpaths)[1])

for(ind in 1:nrow(mvmr_fpaths)){
  print(ind)
  trait_nm = gsub(".gwas.imputed_v3.both_sexes.tsv.gz","",tail(str_split(mvmr_fpaths[ind],"\\/")[[1]],1))
  f1 = fread(unlist(mvmr_fpaths[ind]))
  Bmatrix_clumped[,ind] = f1$tstat[which(f1$variant %in% clumped_sig_snps$variant)] / sqrt(f1$n_complete_samples[which(f1$variant %in% clumped_sig_snps$variant)])  # Z_matrix1
  SEmatrix_clumped[,ind] = 1/sqrt(f1$n_complete_samples[which(f1$variant %in% clumped_sig_snps$variant)])
  colnames(Bmatrix_clumped)[ncol(Bmatrix_clumped)] = trait_nm
  colnames(SEmatrix_clumped)[ncol(SEmatrix_clumped)] = trait_nm
}
var_ind = f1[f1$variant %in% clumped_sig_snps$variant, 1]
rsid_temp = inner_join(var_ind, clumped_sig_snps[,c(1,7)])

rownames(Bmatrix_clumped) = rownames(SEmatrix_clumped) = rsid_temp$rsid
colnames(Bmatrix_clumped) = colnames(SEmatrix_clumped) = colnames(Z_matrix)[-1]

write.csv(Bmatrix_clumped, paste0(res_dir, MVMR_dir, "/Bmatrix_clumped.csv"), row.names = TRUE)
write.csv(SEmatrix_clumped, paste0(res_dir, MVMR_dir, "/SEmatrix_clumped.csv"), row.names = TRUE)

# clump Zmatrix
Z_matrix_clumped = inner_join(Z_matrix, clumped_sig_snps, by="variant")
write.csv(Z_matrix_clumped, paste0(res_dir, MVMR_dir, "/Zmatrix_clumped.csv"), row.names = FALSE)

## Stepwise-MVMR using clumped SNPs from the Zmatrix - re-purposing bGWAS scripts
Z_matrices = paste0(res_dir,MVMR_dir)
#prior_studies = NULL
MR_threshold = 5e-8 #1e-5
MR_ninstruments = 3
#MR_pruning_dist = 500
#MR_pruning_LD = 0
MR_shrinkage = 0.05  #usually 1
stepwise_threshold = NULL
#prior_shrinkage = NULL
#sign_method = "p"
#sign_thresh = 5e-8
#use_permutations = FALSE
#res_pruning_dist = 500
#res_pruning_LD = 0
save_files = FALSE
verbose = FALSE

# Useful small functions (re-used by other functions)
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
  return(formula)
}

## Zmatrix MVMR creation
setwd(Z_matrices)
matrixZ_clumped =  Z_matrix_clumped #fread("Zmatrix_clumped.csv")
matrixZ_clumped %>% 
  select(c("rsid","chr","pos","alt",'ref',colnames(matrixZ_clumped)[2:(length(traits_MVMR)+3)])) %>% #with outcome at first
  rename("rs" = "rsid", "chrm"="chr")-> matrixZ_clumped
## remove HLA region
dim(matrixZ_clumped)
HLA_n_ind = which(!(matrixZ_clumped$chrm==6 & matrixZ_clumped$pos>=28.5e6 & matrixZ_clumped$pos<=33.5e6))
length(HLA_n_ind)
if(length(HLA_n_ind) > 0) {
  matrixZ_clumped = matrixZ_clumped[HLA_n_ind,]
}

Log = c()
# filter the Z matrix based on the traits kept
Zlimit = stats::qnorm(MR_threshold/2, lower.tail = F)
print("# Thresholding... \n")
# DO NOT USE THE LAST COLUMN - OUTCOME
matrixZ_clumped %>%
  filter_at(vars(-c(1:5,as.numeric(ncol(matrixZ_clumped)))),
            any_vars(abs(.data$.) > Zlimit)) -> matrixZ_fil

tmp = (paste0("There are ", dim(matrixZ_fil)[1]," IVs significant for at least one trait at the threshold specified."))
print(tmp)
Log = update_log(Log, tmp, verbose)

ZMatrix = matrixZ_fil
ZMatrix %>%
  select(-c(1:5,as.numeric(ncol(ZMatrix)))) %>%  #no outcome
  apply(MARGIN = 2, FUN = function(col){
    sum(abs(col)>Zlimit)>=MR_ninstruments})  -> StudiesToKeep

# remove studies with instruments left less than MR threshold
if(!all(StudiesToKeep)){
  ZMatrix  %>%
    select(-c(1:5,as.numeric(ncol(ZMatrix)))) %>%
    select(names(StudiesToKeep[!StudiesToKeep])) %>%
    colnames -> StudiesToRemove
  tmp = paste0(paste0(get_names(StudiesToRemove, Z_matrices), collapse=" - "), " : removed (less than ", MR_ninstruments, " instrument after thresholding) \n")
  print(tmp)
  Log = update_log(Log, tmp, verbose)
  
  if(save_files){
    Files_Info$status[Files_Info$File %in% colnames(ZMatrix[,-c(1:5)])[!StudiesToKeep]] =
      paste0("Excluded for MR: less than ",  MR_ninstruments, " strong instrument left after thresholding/pruning")
  }
  
  ZMatrix %>%
    select(-StudiesToRemove) -> ZMatrix
}

# further checking of the SNPs to remove SNPs associated with studies removed because only one SNP
ZMatrix %>%
  select(-c(1:5,as.numeric(ncol(ZMatrix)))) %>%
  apply(MARGIN = 1, FUN = function(row){
    any(abs(row)>Zlimit)}) -> SNPsToKeep  

NAllStudies = nrow(ZMatrix)
if(sum(SNPsToKeep) != NAllStudies){
  ZMatrix %>%
    filter(SNPsToKeep) -> ZMatrix
  tmp = paste0(format(nrow(ZMatrix), big.mark = ",", scientific = F), " SNPs left after removing studies with only one strong instrument \n")
  print(tmp)
  Log = update_log(Log, tmp, verbose)
}

push_extreme_zs_back_a_little_towards_zero <- function(d) { # Some z-scores are just too far from zero
  maxAllowed_z = abs(stats::qnorm(1e-300 / 2)) # p=1e-300 is the max allowed now, truncate z-scores accordingly
  names(d)[!names(d) %in% c("rs","chrm","pos","alt","ref")] -> studies_here
  for(n in studies_here) {
    #n = paste0('`',n,'`')
    d %>%
      mutate(!!n := case_when(
        abs(eval(parse(text=paste0('`',n,'`')))) > maxAllowed_z  ~  maxAllowed_z,
        abs(eval(parse(text=paste0('`',n,'`')))) < -maxAllowed_z ~ -maxAllowed_z,
        TRUE ~ eval(parse(text=paste0('`',n,'`'))))) -> d
  }
  return(d)
}
# Truncate Z-scores -- generally no Zscore is larger than maxAllowed_z = 37.06579
ZMatrix %>%
  push_extreme_zs_back_a_little_towards_zero() -> ZMatrix_trunc

# Steiger filtering test
max_Z = apply(abs(ZMatrix_trunc[,-c(1:5)]),1, which.max)
if(length(which(max_Z==ncol(ZMatrix_trunc[,-c(1:5)]))) > 0){
  ZMatrix_trunc = ZMatrix_trunc[-which(max_Z==ncol(ZMatrix_trunc[,-c(1:5)])),]
  tmp = (paste0("There are ", dim(ZMatrix_trunc)[1]," remaining significant IVs that have a larger standerdised effect on the exposures compared to the outcome."))
  print(tmp)
  Log = update_log(Log, tmp, verbose)
}

# try also with shrinkage for ZMatrix
# Set the z-scores to 0 for the regression if shrinkage
if(MR_shrinkage < 1.0) {
  names(ZMatrix_trunc)[!names(ZMatrix_trunc) %in% c('rs','chrm','pos','alt','ref', OUT_pheno2 )] -> Prior_study_names
  threshold = abs(stats::qnorm(MR_shrinkage/2))
  for(column_of_zs in Prior_study_names) { 
    ZMatrix_trunc %>%
      mutate(!!column_of_zs := case_when(
        abs(eval(parse(text=paste0('`',column_of_zs,'`')))) < threshold ~ 0,
        TRUE ~ eval(parse(text=paste0('`',column_of_zs,'`'))))) -> ZMatrix_shrunk
  }
  tmp = paste0("Applying shrinkage (threshold = ", MR_shrinkage, ") before performing MR. \n")
  print(tmp)
  Log = update_log(Log, tmp, verbose)
}

# save both but mostly working with non-shrunk
write.csv(ZMatrix_trunc, paste0(Z_matrices,"/ZMatrix_filtered.csv"), row.names = F)
write.csv(ZMatrix_shrunk, paste0(Z_matrices,"/ZMatrix_filtered-shrunk.csv"), row.names = F)


## running stepwise MVMR
source(stepMVMR_R)

# without EXP
ZMatrix_subset = as.data.frame(ZMatrix_trunc)[, -(which(colnames(ZMatrix_trunc)==EXP_pheno2))]
pheno_nm = EXP_pheno2
MVMR_noEXP = identify_studiesMR(ZMatrix_subset, MR_shrinkage, MR_threshold, stepwise_threshold, Z_matrices, save_files=FALSE, verbose=TRUE)
Log = c(Log, MVMR_noEXP$log_info)
tmp = "The remaining studies from the stepwise MVMR are:"
print(tmp)
Log = update_log(Log, tmp, verbose)
tmp = c(paste0(paste0("- ", get_names(MVMR_noEXP$studies, Z_matrices)), collapse="\n"))
print(tmp)
Log = update_log(Log, tmp, verbose)

# Adding in the EXP to the remianing studies
# make a sub ZMatrix
ZMatrix_trunc %>%
  select(c("rs","chrm","pos","alt",'ref',MVMR_noEXP$studies,EXP_pheno2,OUT_pheno2)) %>%  ## add in EXP and outcome
  get_Instruments(Zlim=Zlimit) -> ZMatrix_subset1
myFormula = generate_Formula(OUT_pheno2, c(MVMR_noEXP$studies,EXP_pheno2))
# run model - with studies already included and BMI
stats::lm(data=ZMatrix_subset1, formula = myFormula) %>%
  summary %>%
  stats::coef() %>%
  as.data.frame() %>% # 1st create a data.frame (to get rownames)
  tibble::rownames_to_column("nm") %>% # then convert to tibble
  as_tibble()  -> MVMR_EXP
# added naming
MVMR_EXP$nm = gsub('`','',MVMR_EXP$nm)
MVMR_EXP %>%
  set_names(c("study", "estimate", "std_error", "Tstat", "P")) -> MVMR_EXP

log_info = apply(as.array(Log), 1,function(x) gsub("\n", "", x, fixed=T))

write(log_info, paste0(pheno_nm,"_MVMR.log"))

  write_xlsx(list(
    "confounder-noEXP" = as.data.frame(inner_join(MVMR_noEXP$coeffs, trait_info[,c(1,3)], by=c("study"="phenotype"))),
    "confounder-addedEXP" = as.data.frame(inner_join(MVMR_EXP, trait_info[,c(1,3)], by=c("study"="phenotype")))
  ),
  path = paste0(Z_matrices,"/StepwiseMVMR.xlsx"))


## Testing for conditional F statistic with EXP
  
# checking univ MR on each of the exposures  
# get the uni_MR from bi_MR ran in stepwise MVMR
uniMR = fread(paste0(Z_matrices,"/",pheno_nm,"_uniMR.csv"))
uniMR = uniMR[uniMR$study %in% MVMR_noEXP$studies,]
uniMR = uniMR[order(uniMR$P),]
print("The univariate MR estimates in order of signifcance: ")
uniMR
  
# calculating conditional F stat
stdB_matrix = fread(paste0(Z_matrices,"/Bmatrix_clumped.csv"))
stdSE_matrix = fread(paste0(Z_matrices,"/SEmatrix_clumped.csv"))
colnames(stdB_matrix)[1]=colnames(stdSE_matrix)[1] = "rs"

stdB_matrix = stdB_matrix[which(stdB_matrix$rs %in% ZMatrix_subset1$rs),]
stdSE_matrix = stdSE_matrix[which(stdSE_matrix$rs %in% ZMatrix_subset1$rs),]

res <- Map(combn, list(c(1:nrow(uniMR))), seq_along(c(1:nrow(uniMR))), simplify = FALSE)
med_comb = unlist(res, recursive = FALSE)

Fstat_arr = c()
med_arr = c()
axy_arr = c()
se_arr = c()
pval_arr = c()
MVMR_list = list()
counter = 1
for(med in med_comb){
  med = unlist(med)
  ZMatrix_subset1 %>%
    select(c("rs","chrm","pos","alt",'ref',uniMR$study[c(med)],EXP_pheno2,OUT_pheno2)) %>%
    get_Instruments(Zlim=Zlimit) -> ZMatrix_subset_fcond
  
  stdB_matrix1 = stdB_matrix[which(stdB_matrix$rs %in% ZMatrix_subset_fcond$rs),]
  stdSE_matrix1 = stdSE_matrix[which(stdSE_matrix$rs %in% ZMatrix_subset_fcond$rs),]
  stdB_matrix1 = stdB_matrix1[complete.cases(stdB_matrix1),]  
  stdSE_matrix1 = stdSE_matrix1[complete.cases(stdSE_matrix1),]
  all(stdB_matrix1$rs==stdSE_matrix1$rs)

  stdB_matrix1 %>%
    select(-c("rs",OUT_pheno2)) %>%
    select(c(uniMR$study[c(med)], EXP_pheno2))-> X
  
  stdSE_matrix1 %>%
    select(-c("rs",OUT_pheno2)) %>%
    select(c(uniMR$study[c(med)], EXP_pheno2))-> SE
  
  K = ncol(X)-1 # number of mediators
  L = nrow(X)
  
  E = as.vector(unlist(X[,EXP_pheno2]))
  M = as.matrix(X[,1:(ncol(X)-1)])
  
  delta.reg = lm(E ~ -1 + M)
  delta = delta.reg$coefficients
  res = delta.reg$residuals
  
  delta = as.vector(c(-1,delta))
  
  var_IV =  as.matrix(SE**2) %*% (delta**2)
  
  Q = sum(1/(var_IV) * res**2)
  
  Fstat = Q/(L-K)
  Fstat_arr=c(Fstat_arr,Fstat)
  med_arr = c(med_arr,stringr::str_c(uniMR$study[med], collapse="-"))
  
  ## caluclate MVMR for that confounder comb
  myFormula = generate_Formula(OUT_pheno2, c(uniMR$study[c(med)],EXP_pheno2)) 
  # run model - with studies already included and primary exposure
  stats::lm(data=ZMatrix_subset_fcond, formula = myFormula) %>%
    summary %>%
    stats::coef() %>%
    as.data.frame() %>% # 1st create a data.frame (to get rownames)
    tibble::rownames_to_column("nm") %>% # then convert to tibble
    as_tibble()  -> coefs
  # added naming
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

print("Exposure F-stat caluclated for different combinations of confounder traits")
MVMRcondF_df = cbind.data.frame(med_arr,Fstat_arr,axy_arr,se_arr,pval_arr)
MVMR_F8 = MVMR_list[which(MVMRcondF_df$Fstat_arr>8)]
write.csv(MVMRcondF_df, "MVMR_res_condFstat.csv", row.names = FALSE)
MVMR_F8.df <- do.call("rbind", lapply(MVMR_F8, as.data.frame)) 
write.csv(MVMR_F8.df, "MVMR_res_condFstatGT8.csv", row.names = FALSE)
