bi_MR = function(par_exp,par_fpaths, par_out){
  EXP = par_exp
  OUT = par_out
  conf = fread(cmd = paste0("zcat < /data/sgg3/data/neale_files/both_sexes/",par_fpaths))
  colnames(conf) = toupper(colnames(conf))
  
  ## inner functions
  prune_X = function(zX,p_limit=1e-5){
    zX=zX
    z_limit=abs(qnorm(0.5*p_limit))
    ind_keep=which(abs(zX)>z_limit)
    ind_keep=unique(ind_keep)
    ind_keep=list(ind_keep)
    return(ind_keep)
  }
  
  # Taken from Jonathan Sulc to create bins that fit a maximum of 50k SNPs (max threshold for clumping)
  snp_bin  =  function( snp_ranks,
                        chunk_size = 50000 ){
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
  
  # Set threshold values
  pval=2*pnorm(-abs(5.45))
  pval1=2*pnorm(-abs(4))
  reverse_t_threshold  =  qnorm( 5e-2 )
  
  ### BMI
  # Join the exposure and outcome files
  input.df = dplyr::inner_join(conf, OUTdf,
                               by = c("VARIANT"))
  # Remove the HLA region due to highly associated SNPs
  input.df_ind = which(!(input.df$CHR==6 & input.df$POS>=28.5e6 & input.df$POS<=33.5e6))
  input.df = input.df[input.df_ind,]
  input.df$CHR = as.numeric(input.df$CHR)
  # Reorder based on chromosome and position
  input.df %>% dplyr::arrange(CHR, POS) -> input.df
  
  nX = mean(input.df$N_COMPLETE_SAMPLES.x)  #get sample size for trait X
  nY = mean(input.df$N_COMPLETE_SAMPLES.y)  #get sample size for trait Y
  
  bX = input.df$TSTAT.x/sqrt(nX)   #get standardised beta for trait X
  bY = input.df$TSTAT.y/sqrt(nY)   #get standardised beta for trait Y
  
  X = select(input.df, RSID, CHR, ALT, REF, BETA.x, SE.x, PVAL.x, TSTAT.x, N_COMPLETE_SAMPLES.x) %>% rename(unstdb = BETA.x, sderr = SE.x, pval = PVAL.x, TSTAT=TSTAT.x, N = N_COMPLETE_SAMPLES.x, A1 = ALT, A2 = REF)
  Y = select(input.df, RSID, CHR, ALT, REF, BETA.y, SE.y, PVAL.y, TSTAT.y, N_COMPLETE_SAMPLES.y) %>% rename(unstdb = BETA.y, sderr = SE.y, pval = PVAL.y, TSTAT=TSTAT.y, N = N_COMPLETE_SAMPLES.y, A1 = ALT, A2 = REF)
  X$BETA = bX
  Y$BETA = bY
  X$SE = 1/sqrt(nX)
  Y$SE = 1/sqrt(nY)
  
  # Forward MR estimation
  mr_dataX = cbind.data.frame(SNP = X$RSID, beta = X$BETA, se = X$SE, effect_allele = X$A1, other_allele = X$A2, chr=X$CHR, Phenotype=EXP, tstat=X$TSTAT )
  mr_dataY = cbind.data.frame(SNP = Y$RSID, beta = Y$BETA, se = Y$SE, effect_allele = Y$A1, other_allele = Y$A2, chr=Y$CHR, Phenotype=OUT, tstat=Y$TSTAT )
  
  mr_ind=unlist(prune_X(mr_dataX$tstat,pval))
  print(length(mr_ind))
  if(length(mr_ind)==0){
    mr_ind=unlist(prune_X(mr_dataX$tstat,pval1))
    print(length(mr_ind))
  }
  mr_dataX = mr_dataX[mr_ind,]
  mr_dataY = mr_dataY[mr_ind,]
  
  # Remove SNPs that are more strongly associated with the outcome than the exposure
  ind_keep=which((abs(mr_dataX$beta)-abs(mr_dataY$beta))/sqrt(mr_dataX$se^2+mr_dataY$se^2) > reverse_t_threshold)
  print(length(ind_keep))
  mr_dataX = mr_dataX[ind_keep,]
  mr_dataY = mr_dataY[ind_keep,]
  
  exp_dat <- TwoSampleMR::format_data(mr_dataX, type="exposure")  #same rows as mr_dataX
  clump_bin = snp_bin(mr_dataX,50000)
  
  exp_data = c()
  for (x in 1:length(clump_bin)) {
    temp = exp_dat[exp_dat$SNP %in% clump_bin[[x]]$SNP,]
    temp1 = TwoSampleMR::clump_data(temp)
    exp_data=rbind(exp_data,temp1)
  }
  
  dups=which(duplicated(exp_data$SNP)==TRUE)
  if(length(dups)>0){
    exp_dat2 = exp_data[-dups,]
  }else{
    exp_dat2 = exp_data
  }
  
  out_dat <- TwoSampleMR::format_data(mr_dataY, type="outcome")
  out_dat2=out_dat[out_dat$SNP %in% exp_dat2$SNP,]
  
  exp_dat2=exp_dat2[order(exp_dat2$SNP),]
  out_dat2=out_dat2[order(out_dat2$SNP),]
  
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
  # Harmonise the exposure and outcome data
  dat <- TwoSampleMR::harmonise_data(
    exposure_dat = exp_dat2,
    outcome_dat = out_dat2, action = action
  )
  
  # Sensitivity - Q-test
  het <- TwoSampleMR::mr_heterogeneity(dat)
  het$I2 = ((het$Q-het$Q_df)/het$Q)*100
  plei <- TwoSampleMR::mr_pleiotropy_test(dat)
  
  smaller=FALSE
  tryCatch( {res1 <- TwoSampleMR::mr(dat); print("Bigger MR list") }, error = function(e) {smaller<<-TRUE})
  if(smaller){
    print("Smaller MR list")
    res1 <- TwoSampleMR::mr(dat, method_list=c("mr_egger_regression","mr_weighted_median","mr_ivw"))
  }else{
    res1 <- TwoSampleMR::mr(dat)
  }
  ### join res1.BMI with proper exp and out if not already
  
  # Reverse MR estimation
  rm(mr_dataX,mr_dataY,exp_dat,exp_data,exp_dat2,out_dat,out_dat2,dups,ind_keep,mr_ind,clump_bin, action, temp, temp1,dat,het,plei,smaller)
  
  # Reverse the exposure and outcome to Y - X, nothing else besides this needs to change
  mr_dataX = cbind.data.frame(SNP = Y$RSID, beta = Y$BETA, se = Y$SE, effect_allele = Y$A1, other_allele = Y$A2, chr=Y$CHR, Phenotype=OUT, tstat = Y$TSTAT )
  mr_dataY = cbind.data.frame(SNP = X$RSID, beta = X$BETA, se = X$SE, effect_allele = X$A1, other_allele = X$A2, chr=X$CHR, Phenotype=EXP, tstat = X$TSTAT )
  
  mr_ind=unlist(prune_X(mr_dataX$tstat,pval))
  print(length(mr_ind))
  if(length(mr_ind)==0){
    mr_ind=unlist(prune_X(mr_dataX$tstat,pval1))
    print(length(mr_ind))
  }
  mr_dataX = mr_dataX[mr_ind,]
  mr_dataY = mr_dataY[mr_ind,]
  
  ind_keep=which((abs(mr_dataX$beta)-abs(mr_dataY$beta))/sqrt(mr_dataX$se^2+mr_dataY$se^2) > reverse_t_threshold)
  print(length(ind_keep))
  mr_dataX = mr_dataX[ind_keep,]
  mr_dataY = mr_dataY[ind_keep,]
  
  exp_dat <- TwoSampleMR::format_data(mr_dataX, type="exposure")  #same rows as mr_dataX
  clump_bin = snp_bin(mr_dataX,50000)
  
  exp_data = c()
  for (x in 1:length(clump_bin)) {
    temp = exp_dat[exp_dat$SNP %in% clump_bin[[x]]$SNP,]
    temp1 = TwoSampleMR::clump_data(temp)
    exp_data=rbind(exp_data,temp1)
  }
  
  dups=which(duplicated(exp_data$SNP)==TRUE)
  if(length(dups)>0){
    exp_dat2 = exp_data[-dups,]
  }else{
    exp_dat2 = exp_data
  }
  
  out_dat <- TwoSampleMR::format_data(mr_dataY, type="outcome")
  out_dat2=out_dat[out_dat$SNP %in% exp_dat2$SNP,]
  
  exp_dat2=exp_dat2[order(exp_dat2$SNP),]
  out_dat2=out_dat2[order(out_dat2$SNP),]
  
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
  
  # Harmonise the exposure and outcome data
  dat <- TwoSampleMR::harmonise_data(
    exposure_dat = exp_dat2,
    outcome_dat = out_dat2, action = action
  )
  
  # Sensitivity - Q-test
  het <- TwoSampleMR::mr_heterogeneity(dat)
  het$I2 = ((het$Q-het$Q_df)/het$Q)*100
  plei <- TwoSampleMR::mr_pleiotropy_test(dat)
  
  smaller=FALSE
  tryCatch( {res2 <- TwoSampleMR::mr(dat); print("Bigger MR list") }, error = function(e) {smaller<<-TRUE})
  if(smaller){
    print("Smaller MR list")
    res2 <- TwoSampleMR::mr(dat, method_list=c("mr_egger_regression","mr_weighted_median","mr_ivw"))
  }else{
    res2 <- TwoSampleMR::mr(dat)
  }
  
  res = rbind(res1,res2)
  return(res)
}