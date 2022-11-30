## Systematic search for potential confounders/colliders/mediators acting on EXP-OUT
library(data.table)
library(TwoSampleMR)
library(lhcMR)
library(dplyr)
library(rslurm)

# variables
EXP_dir = "/data/sgg3/data/neale_files/both_sexes/body/21001_irnt.gwas.imputed_v3.both_sexes.tsv.gz"
EXP_pheno = "21001"
OUT_dir = "/data/sgg3/data/neale_files/both_sexes/lifestyle/845.gwas.imputed_v3.both_sexes.tsv.gz"
OUT_pheno = "845"
VARinfo_dir = "/data/sgg2/liza/SEM_Real/test/data/variants.tsv.bgz"
res_dir = paste0("/data/sgg3/liza/TraitCorr/pheWAS/",EXP_pheno)
fpaths_fil_dir = paste0(res_dir,"/data/fpaths_fil_nfil.txt")
EXP_irnt = TRUE  #if true change all grepl to paste0(EXP_pheno,"_irnt")
OUT_irnt = FALSE  #if true change all grepl to paste0(EXP_pheno,"_irnt")
bi_MR_dir = "/data/sgg3/liza/TraitCorr/pheWAS/bash_trial/bi_MR.R"
setwd(res_dir)

#load QC-filtered data
load(paste0(res_dir,"/QCdata_",EXP_pheno,".Rdata"))

## get file paths to run MR on the remaining candidate traits 
fpaths_fil_nfil =fread(fpaths_fil_dir, header=FALSE)
fpaths_fil1 = fpaths_fil_nfil %>% tidyr::separate(V1, c("folder","file"), sep="/", remove=F)

fpath_ind = c()
for(i in 1:nrow(trait_info_noEXP)){
  fpath_ind = c(fpath_ind, which(fpaths_fil1$file==paste0(trait_info_noEXP$phenotype[i],".gwas.imputed_v3.both_sexes.tsv.gz")))
}
length(fpath_ind)
mr_fpaths = fpaths_fil1[fpath_ind,1]
write.csv(mr_fpaths, paste0(res_dir,"/data/mr_fpaths.csv"), row.names = F)

# read in EXP trait
X = fread(EXP_dir)
# read in OUT trait
Y = fread(OUT_dir)
# read in VARinfo to get RSID
VARinfo = fread(cmd = paste0(" zcat < ",VARinfo_dir), sep="\t", header=TRUE) ## Variant info file from UKBB GWAS Imputed v3 - File Manifest Release 20180731

X = dplyr::inner_join(X, VARinfo[,c(1:6)]) #Joining, by = "variant"
Y = dplyr::inner_join(Y, VARinfo[,c(1:6)]) #Joining, by = "variant"
rm(VARinfo)

#QC
colnames(X) = toupper(colnames(X))
colnames(Y) = toupper(colnames(Y))
#remove chromosome X
X = X[-which(X$CHR=="X"),]
Y = Y[-which(Y$CHR=="X"),]

# functions to use/repeat
source(bi_MR_dir)

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

#get confounder fpaths
mr_fpaths = fread(paste0(res_dir,"/data/mr_fpaths.csv"))
par_exps = stringr::str_match(mr_fpaths$V1, "/(.*?)\\.")[,2]

####### once exposure
pars1 <- data.frame(par_exp = par_exps,
                    par_fpaths = mr_fpaths$V1,
                    par_out = EXP_pheno)

OUTdf = X
sjob1 = slurm_apply(f = bi_MR, params = pars1, jobname = paste0("biMR-",EXP_pheno), nodes = nrow(pars1), cpus_per_node = 1,
                    global_objects = c("OUTdf"),
                    slurm_options = list(partition = "sgg", account = "sgg", mem = "20000"),
                    submit = TRUE)
# Keep a loop open till the job is completely done.
wait_counter = 0
while (wait_counter < 1) {
  wait_counter = 0
  if(lhcMR:::get_job_status2(sjob1)){
    wait_counter = wait_counter + 1
  } else{
    wait_counter = wait_counter
  }
}

res1 = get_slurm_out(sjob1, outtype = 'table')
write.csv(res1, paste0("biMR_",EXP_pheno,".csv"), row.names = F)

######### once outcome
pars2 <- data.frame(par_exp = par_exps,
                    par_fpaths = mr_fpaths$V1,
                    par_out = OUT_pheno)

OUTdf = Y
sjob2 = slurm_apply(f = bi_MR, params = pars2, jobname = paste0("biMR-",OUT_pheno), nodes = nrow(pars2), cpus_per_node = 1,
                    global_objects = c("OUTdf"),
                    slurm_options = list(partition = "sgg", account = "sgg", mem = "20000"),
                    submit = TRUE)
# Keep a loop open till the job is completely done.
wait_counter = 0
while (wait_counter < 1) {
  wait_counter = 0
  if(lhcMR:::get_job_status2(sjob2)){
    wait_counter = wait_counter + 1
  } else{
    wait_counter = wait_counter
  }
}

res2 = get_slurm_out(sjob2, outtype = 'table')
write.csv(res2, paste0("biMR_",OUT_pheno,".csv"), row.names = F)
