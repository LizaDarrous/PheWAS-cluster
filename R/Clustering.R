### 3) Clustering

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
#library(xlsx)
library(TwoSampleMR)

# variables
EXP_pheno = "21001"
OUT_pheno = "845"
pheno_irnt = TRUE  #maybe if true change all grepl to paste0(EXP_pheno,"_irnt")
res_dir = paste0("/data/sgg3/liza/TraitCorr/pheWAS/",EXP_pheno)
setwd(res_dir)

# cochran's Qtest >>>>
q.meta.test = function(b_x, se_x){  #https://cjvanlissa.github.io/Doing-Meta-Analysis-in-R/heterogeneity-statistics.html
  se_x[which(se_x==0)]=1  ### check!
  w0 = se_x^-2
  w = w0/sum(w0)
  
  meta.bet_x = sum(b_x*w)
  meta.se_x = (sum(w^2*se_x^2))^0.5
  Q = sum( (b_x-meta.bet_x)^2 * w0 )
  df = length(b_x)-1
  Q.pval = 1 - pchisq(Q,df)
  return(list("Q"=Q, "Q.pval"=Q.pval))
}

#load QC-filtered data
load(paste0(res_dir,"/QCdata_",EXP_pheno,".Rdata"))

# k-means clustering data set up
stdBeta_df_noEXP_nosus = stdBeta_df_noEXP[-sus_SNP_ind,]
eff_df = abs(stdBeta_df_noEXP_nosus)
eff_dfs = t(scale(t(eff_df)))

# Fitting K-Means clustering Model 
kmeansIC = function(fit){
  #https://stackoverflow.com/questions/15839774/how-to-calculate-bic-for-k-means-clustering-in-r
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(data.frame(AIC = D + 2*m*k,
                    BIC = D + log(n)*m*k))
}

cluster_list = list()
IC_df = data.frame(matrix(data=NA, nrow = 49, ncol=1))
colnames(IC_df) = c("AIC")
IC_df$nCluster = 2:50

for(i in 2:50){
  set.seed(240) # setting seed
  kmeans.re <- kmeans(eff_dfs, centers = i, nstart = 50, iter.max = 300)
  cluster_list[[length(cluster_list)+1]] = kmeans.re
  IC = kmeansIC(kmeans.re)
  IC_df[(i-1),1] = IC$AIC
}

# cluster number identification for each observation
kmeans.minAIC = cluster_list[[which(IC_df$AIC==min(IC_df$AIC))]]
nClust.AIC = max(kmeans.minAIC$cluster)

# plot AIC
pdf(file=paste0(res_dir,"/AIC-ncluster_",EXP_pheno,".pdf"), width=5, height = 7)
plot(IC_df$nCluster,IC_df$AIC, xlab="nCluster", ylab="AIC")
abline(v=nClust.AIC,col="red")
dev.off()

# assigning SNPs to clusters
AICclusters_rsid = list()
for(i in 1:nClust.AIC){
  AICclusters_rsid[[i]] = names(kmeans.minAIC$cluster)[which(kmeans.minAIC$cluster==i)]
}
print(paste0("Number of SNPs in each AIC grouped clusters: ", paste0(lengths(AICclusters_rsid), collapse = ", ")))
AICclusters_rsid_df = t(plyr::ldply(AICclusters_rsid, rbind))
write.csv(AICclusters_rsid_df, paste0(res_dir,"/AICclusters_rsid_",EXP_pheno,".csv"), row.names = FALSE)

# MR per cluster - AIC  - first instance of OUTCOME
exp_dat = fread(paste0(res_dir,"/sig-clumped-IVs_",EXP_pheno,".csv")) 
out_dat = fread(paste0(res_dir,"/clumped-IVs_",OUT_pheno,".csv")) 

exp_dat2 = rename(exp_dat,"unstdBeta"="beta", "unstdSE"="se")
out_dat2 = rename(out_dat,"unstdBeta"="beta", "unstdSE"="se")

exp_dat2$beta = exp_dat2$tstat/sqrt(exp_dat2$N)
out_dat2$beta = out_dat2$tstat/sqrt(out_dat2$N)
exp_dat2$se = 1/sqrt(exp_dat2$N)
out_dat2$se = 1/sqrt(out_dat2$N)

exp_dat2 = format_data(exp_dat2, type="exposure", samplesize_col="N", z_col = "tstat", pval_col = "pval.exposure")
out_dat2 = format_data(out_dat2, type="outcome", samplesize_col="N", z_col = "tstat", pval_col = "pval.outcome")

common_SNPs = intersect(exp_dat2$SNP, out_dat2$SNP)
exp_dat2 = exp_dat2[match(common_SNPs, exp_dat2$SNP),]
out_dat2 = out_dat2[match(common_SNPs, out_dat2$SNP),]

if(all(exp_dat2$`effect_allele.exposure` == out_dat2$`effect_allele.outcome`)){
  print("action=1")
  action = 1
} else {
  print("action=2/3")
  aligned = which(exp_dat2$`effect_allele.exposure` == out_dat2$`effect_allele.outcome` &
                    exp_dat2$`other_allele.exposure` == out_dat2$`other_allele.outcome`)
  swapped = which(exp_dat2$`effect_allele.exposure` == out_dat2$`other_allele.outcome` &
                    exp_dat2$`other_allele.exposure` == out_dat2$`effect_allele.outcome`)
  exp_dat2[swapped,'beta.exposure']=exp_dat2[swapped,'beta.exposure']*-1
  exp_dat2 = exp_dat2[c(aligned,swapped),]
  out_dat2 = out_dat2[c(aligned,swapped),]
  action = 1  #made sure all strands are okay
}

# run MR or all SNPS
'%!in%' <- function(x,y)!('%in%'(x,y))
axy_ar = c()
se_ar = c()
MR_output = paste0(res_dir,"/",EXP_pheno,"-",OUT_pheno,"_AICclusters_MR.csv")
write.table("", file = MR_output, append = FALSE)
for(clustN in 0:nClust.AIC){
  # harmonise the exposure and outcome data
  dat <- TwoSampleMR::harmonise_data(
    exposure_dat = exp_dat2,
    outcome_dat = out_dat2, action = action
  )
  dat = dat[which(dat$SNP %!in% sus_SNPs),]
  if(clustN>0){dat = dat[which(dat$SNP %in% AICclusters_rsid[[clustN]]),]}
  # sensitivity - Q-test
  het <- TwoSampleMR::mr_heterogeneity(dat)
  het$I2 = ((het$Q-het$Q_df)/het$Q)*100
  het$Avg_het = het$Q/(het$Q_df-1)
  plei <- TwoSampleMR::mr_pleiotropy_test(dat)
  
  smaller=FALSE
  tryCatch( {res1 <- TwoSampleMR::mr(dat); }, error = function(e) {smaller<<-TRUE}) #print("Bigger MR list")
  if(smaller){
    print("Smaller MR list")
    res1 <- TwoSampleMR::mr(dat, method_list=c("mr_egger_regression","mr_weighted_median","mr_ivw"))
  }else{
    res1 <- TwoSampleMR::mr(dat)
  }
  
  if(nrow(res1)<2 & nrow(res1)>0){
    axy_ar = c(axy_ar,res1[,'b'])
    se_ar = c(se_ar,res1[,'se'])
  }else if(nrow(res1)==0){
    axy_ar = c(axy_ar,0)
    se_ar = c(se_ar,0)
  }else if(nrow(res1)>2){
    axy_ar = c(axy_ar,res1[which(res1$method=="Inverse variance weighted"),'b'])
    se_ar = c(se_ar, res1[which(res1$method=="Inverse variance weighted"),'se'])
  } 
  
  suppressWarnings(write.table(paste0("Cluster", clustN), file=MR_output, sep = ",", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE))
  suppressWarnings(write.table(as.data.frame(res1), file=MR_output, sep = ",", append = TRUE, row.names = FALSE))
  suppressWarnings(write.table(as.data.frame(het), file=MR_output, sep = ",", append = TRUE, row.names = FALSE) )
  suppressWarnings(write.table(as.data.frame(plei), file=MR_output, sep = ",", append = TRUE, row.names = FALSE))
  suppressWarnings(write.table("*", MR_output, sep = ",", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE))
}

# forest plot of IVW of different clusters
cluster_axy = as.data.frame(cbind("axy_ar"=axy_ar,"se_ar"=se_ar))
cluster_axy$clusters = c("All", paste0("Cluster", 1:nClust.AIC))
cluster_axy$clusters = factor(x = cluster_axy$clusters, levels = c(rev(paste0("Cluster", 1:nClust.AIC)),"All"))
cluster_axy$boxsize = c(0.5,lengths(AICclusters_rsid)/sum(lengths(AICclusters_rsid)))

write.csv(cluster_axy, paste0(res_dir,"/",EXP_pheno,"-",OUT_pheno,"_AICclusters_MRsummary.csv"), row.names = FALSE)

by_n <- function(n) { seq(-1000, 1000, by = n) } #https://stackoverflow.com/questions/52579553/scale-x-y-continuous-breaks-by-n-in-r-ggplot2

pdf(file = paste0(res_dir,"/",EXP_pheno,"-",OUT_pheno,"_AICclusters_MR.pdf"), width=8, height = 6)  ###getting corrupt
p <- ggplot(cluster_axy, aes(y=clusters, x=axy_ar, xmin=axy_ar-(1.96*se_ar), xmax=axy_ar+(1.96*se_ar))) +
  geom_errorbarh(height=0, size=0.5, alpha=0.3) +
  geom_point(size=cluster_axy$boxsize*10, shape = 15) + 
  geom_vline(xintercept = cluster_axy[1,1], alpha=0.7, colour = "blue", linetype='dashed') +
  scale_x_continuous(breaks = by_n(0.05)) + #breaks = seq(-0.6,0.4,0.1)
  labs(x='IVW causal effect estimate', y = 'Clusters') +
  geom_vline(xintercept=0, color='black', alpha=.5) +
  ggtitle("IVW causal effect estimates - K-means clustering based on AIC") + 
  theme_classic()
print(p)
dev.off()

# save estimates for Q-test
q_test = q.meta.test(axy_ar,se_ar)
print(paste0("Q-test of heterogeneity between cluster estimates = ", q_test[1],", p-value = ",q_test[2]))
cat(paste0("Q-test of heterogeneity between cluster estimates = ", q_test[1],", p-value = ",q_test[2]),
    file=paste0(res_dir,"/",EXP_pheno,"-",OUT_pheno,"_AICclusters_Qtest.txt"),sep="\n")

## Enrichment analysis: getting traits specific to each cluster
# perSNP heritability in each cluster, get the proportion of that for each trait (enrichment ratio)

trait_h2 = matrix(NA, ncol = ncol(stdBeta_df_noEXP_nosus), nrow = nClust.AIC)
for(i in 1:nClust.AIC){
  c_df =  stdBeta_df_noEXP_nosus[which(rownames(stdBeta_df_noEXP_nosus) %in% AICclusters_rsid[[i]]),]^2
  trait_h2[i,] = colSums(c_df)/nrow(c_df)
}
colnames(trait_h2) = colnames(stdBeta_df_noEXP_nosus)
# normalise perSNP-h2 for each trait: dividing by the col mean
trait_h2_norm = sweep(trait_h2,2,apply(trait_h2, 2,mean), "/")
max_valind = apply(trait_h2_norm, 2, which.max)
max_val = apply(trait_h2_norm, 2, max)
traitPcluster = cbind.data.frame("trait" = colnames(trait_h2_norm), "cluster" = max_valind, "enrichmentRatio" = max_val)
traitPcluster = inner_join(traitPcluster, trait_info_noEXP[,c(1,3)], by=c("trait"="phenotype"))

traitPcluster_top10 = traitPcluster %>%
  group_by(cluster) %>% slice_max(order_by = enrichmentRatio, n = 10)

traitPcluster_top10 %>%
  group_by(cluster) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = cluster, values_from = c(trait, description, enrichmentRatio), names_vary = "slowest") %>%
  select(-row) -> traitPcluster_top10

write.csv(traitPcluster_top10, paste0(res_dir,"/",EXP_pheno,"-",OUT_pheno,"_TraitsTop10_AICclusters_perSNPh2.csv"), row.names = FALSE)
write.csv(trait_h2_norm, paste0(res_dir,"/",EXP_pheno,"-",OUT_pheno,"_Traits_AICclusters_perSNPh2.csv"), row.names = FALSE)
