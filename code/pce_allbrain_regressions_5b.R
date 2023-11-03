#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
args=as.numeric(args)
args=args[1]

library(data.table)
library(dplyr)
library(psych)
library(lmerTest)
library(lme4)
library(tidyr)
library(stringr)

file_location = "C:/Users/David/Documents/BRAINlab/ABCD/Package_1195434/ABCD_5/core/"

cbcl_long = fread(paste(file_location,'allimaging_pce_orig.csv',sep=""),data.table = F)

##############
winsorize = function(x,q=3){
  
  mean_x = mean(x,na.rm = T)
  sd_x = sd(x,na.rm = T)
  top_q = mean_x + q*sd_x
  bottom_q = mean_x - q*sd_x
  
  x[x>top_q] = top_q
  x[x<bottom_q] = bottom_q
  
  return(x)
  
}


any_non_zero = function(dat){
  val_out = apply(dat,
                  MARGIN = 1,function(X){
                    (sum(X,na.rm = T)>0)*1
                  })
  return(val_out)
}


inormal <- function(x){qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))}

##############################

cbcl_long$none = 1

y.vars =  colnames(cbcl_long)[c(
  
  grep(colnames(cbcl_long),pattern = "rsfmri_c_ngd_.*_ngd_.*"),
  grep(colnames(cbcl_long),pattern = "rsfmri_cor"),
  
  grep(colnames(cbcl_long),pattern = "^smri"),
  grep(colnames(cbcl_long),pattern = "^dmri_dti|^dmri_rsi")
)]

y.vars = y.vars[-grep(x = y.vars,pattern = "_n$|_n_|none|total|mean|4thvt|3rdvt|lesion|inflatvent|lvtlh|lvtrh|qc|_all|ventricle|whole|intracranial|hint|suprateialv|subcorticalgv|wmatter|scs_cc")]
#y.vars = y.vars[-grep(x = y.vars,pattern = "dmri_dti.*gm")]

cbcl_long = cbcl_long %>% dplyr::filter(imgincl_t1w_include == 1)

#   "ventr|vent|csf|lesion|wmatter|hint|_cc|whole|intracranialv|suprateialv|subcorticalgv|_all_|hemisphere|cerebral_white_matter")]

#y.vars[grep(x = y.vars,pattern = "_n$|_n_|none")]

global.vars = c(
  "smri_thick_cdk_mean",
  "smri_sulc_cdk_mean",
  "smri_area_cdk_total",
  "smri_vol_scs_intracranialv",
  
  "dmri_dtifa_fiberat_allfibers",
  "dmri_dtimd_fiberat_allfibers",
  "dmri_dtitd_fiberat_allfibers",
  "dmri_dtild_fiberat_allfibers",
  "dmri_dtivol_fiberat_allfibers",
  
  "dmri_dtifawm_cdsn_mean",
  "dmri_dtimdwm_cdsn_mean",
  "dmri_dtitdwm_cdsn_mean",
  "dmri_dtildwm_cdsn_mean",
  
  "dmri_dtifagm_cdsn_mean",
  "dmri_dtimdgm_cdsn_mean",
  "dmri_dtitdgm_cdsn_mean",
  "dmri_dtildgm_cdsn_mean",
  
  # need to make
  "dmri_dtifa_scts_mean",
  "dmri_dtimd_scts_mean",
  "dmri_dtitd_scts_mean",
  "dmri_dtild_scts_mean",
  
  "dmri_rsifni_scs_mean",
  "dmri_rsifni_fib_allfib",
  "dmri_rsifniwm_cdk_mean",
  "dmri_rsifnigm_cdk_mean",
  
  "dmri_rsihnd_scs_mean",
  "dmri_rsihnd_fib_allfib",
  "dmri_rsihndwm_cdk_mean",
  "dmri_rsihndgm_cdk_mean",
  
  "dmri_rsihni_scs_mean",
  "dmri_rsihni_fib_allfib",
  "dmri_rsihniwm_cdk_mean",
  "dmri_rsihnigm_cdk_mean",
  
  "dmri_rsihnt_scs_mean",
  "dmri_rsihnt_fib_allfib",
  "dmri_rsihntwm_cdk_mean",
  "dmri_rsihntgm_cdk_mean",
  
  "dmri_rsirnd_scs_mean",
  "dmri_rsirnd_fib_allfib",
  "dmri_rsirndwm_cdk_mean",
  "dmri_rsirndgm_cdk_mean",
  
  "dmri_rsirni_scs_mean",
  "dmri_rsirni_fib_allfib",
  "dmri_rsirniwm_cdk_mean",
  "dmri_rsirnigm_cdk_mean",
  
  "dmri_rsirnt_scs_mean",
  "dmri_rsirnt_fib_allfib",
  "dmri_rsirntwm_cdk_mean",
  "dmri_rsirntgm_cdk_mean"
  
  
  
)

colnames(cbcl_long)[which(colnames(cbcl_long) =="dmri_dtimdgm_cortdesikan_mean" )] = "dmri_dtimdgm_cdsn_mean"

# for (j in 1: length(global.vars)){
# print(paste(global.vars[j],which(colnames(cbcl_long) ==global.vars[j] )))
#   
# }


y.vars = c(y.vars,global.vars)

#

analysis.framework=expand.grid(stringsAsFactors = F,
                               
                               y.var = y.vars
                               
                               ,
                               global.var = "none",
                               include.var=NA,
                               movement.var = NA
                               
)


analysis.framework$include.var[grep(analysis.framework$y.var,pattern = "smri")] = "imgincl_t1w_include"
analysis.framework$include.var[grep(analysis.framework$y.var,pattern = "rsfmri")] = "imgincl_rsfmri_include"
analysis.framework$include.var[grep(analysis.framework$y.var,pattern = "dmri")] = "imgincl_dmri_include"


analysis.framework$movement.var[grep(analysis.framework$y.var,pattern = "smri")] = "none"
analysis.framework$movement.var[grep(analysis.framework$y.var,pattern = "rsfmri")] = "rsfmri_meanmotion"
analysis.framework$movement.var[grep(analysis.framework$y.var,pattern = "dmri")] = "dmri_meanmotion"


analysis.framework$global.var[grep(analysis.framework$y.var,pattern = "rsfmri")] = "none"


for(i in 1:length(global.vars)){
  v.temp = global.vars[i]
  new.temp = str_split(string = v.temp,pattern = "_",simplify =T)[-4] %>% paste(collapse = "_")
  analysis.framework$global.var[grep(analysis.framework$y.var,pattern = new.temp)] = v.temp
}


analysis.framework$global.var[grep(analysis.framework$y.var,pattern = "smri_vol")] = "smri_vol_scs_intracranialv"
analysis.framework$global.var[analysis.framework$global.var == analysis.framework$y.var] = "none"


to_make=c(
  "dmri_dtifa_scts_mean",
  "dmri_dtimd_scts_mean",
  "dmri_dtitd_scts_mean",
  "dmri_dtild_scts_mean",
  "dmri_rsifni_scs_mean",
  "dmri_rsihnd_scs_mean", 
  "dmri_rsihni_scs_mean", 
  "dmri_rsihnt_scs_mean",
  "dmri_rsirnd_scs_mean",
  "dmri_rsirni_scs_mean",  
  "dmri_rsirnt_scs_mean")

new.global.var = matrix(data = NA,nrow = dim(cbcl_long)[1],ncol = length(to_make)) %>% as.data.frame()
colnames(new.global.var) = to_make

for(i in 1:length(to_make)){
  v.temp = to_make[i]
  new.temp = str_split(string = v.temp,pattern = "_",simplify =T)[-4] %>% paste(collapse = "_")
  var.to.mean = analysis.framework$y.var[grep(analysis.framework$y.var,pattern = new.temp)]
  var.to.mean = var.to.mean[-grep(x = var.to.mean,pattern = "mean")]
  
  new.global.var[,i] = apply(cbcl_long[,match(var.to.mean,colnames(cbcl_long))],1,  FUN = function(X){mean(na.omit(X))})
  new.global.var[,i][is.nan( new.global.var[,i])] = NA
  
}

cbcl_long = cbind(cbcl_long,new.global.var)

colnames(cbcl_long) = str_replace(string =  colnames(cbcl_long),pattern = "cortdesikan",replacement = "cdsn")
################

cbcl_long$y.var = cbcl_long[,which(colnames(cbcl_long) == analysis.framework$y.var[args])]
cbcl_long$global.var = cbcl_long[,which(colnames(cbcl_long) == analysis.framework$global.var[args])]
cbcl_long$include.var = cbcl_long[,which(colnames(cbcl_long) == analysis.framework$include.var[args])]
cbcl_long$movement.var = cbcl_long[,which(colnames(cbcl_long) == analysis.framework$movement.var[args])]

##################
#main model

cbcl_long$rel_relationship_1 = ifelse(test = cbcl_long$rel_relationship == 1,yes = 1,no = 0)
cbcl_long$rel_relationship_2 = ifelse(test = cbcl_long$rel_relationship == 2,yes = 1,no = 0)
cbcl_long$rel_relationship_3 = ifelse(test = cbcl_long$rel_relationship == 3,yes = 1,no = 0)


cbcl_dat = cbcl_long %>% dplyr::select(all_of(
  c(
    "src_subject_id","eventname","site_id_l","rel_family_id","MJ_exposure",
    "mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
    "ed","pds","famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
    "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
    "dgs_after",#"try_dg",
    "income","interview_age",
    "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
    "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
    "income_1","income_2","income_3","income_4",
    "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
    "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
    "Prisma_fit","Prisma","DISCOVERY","Achieva",
    "y.var","global.var","include.var","movement.var"
    #,"mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth"
    # excluding in some analyses -  mat_age_knew_preg, total_birth_weight, devhx_10, devhx_6_p, income      
  ))) %>% dplyr::filter(include.var == 1) %>% unique() %>% na.omit()



var.names=c( "mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
             "ed","pds","famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
             "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
             "dgs_after",#"try_dg",
             "income","interview_age",
             "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
             "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
             "income_1","income_2","income_3","income_4",
             "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
             "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
             "Prisma_fit","Prisma","DISCOVERY","Achieva",
             "y.var")

var.names2=c( #"mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
  "ed","pds",
  #"famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
  #"demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
  #"dgs_after",
  "income","interview_age",
  #"famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
  #"famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
  #"income_1","income_2","income_3","income_4",
  #"famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
  #"devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
  #"Prisma_fit","Prisma","DISCOVERY","Achieva",
  "y.var")


cbcl_dat = cbcl_dat %>% 
  mutate_at(.vars = var.names2,.funs = winsorize) %>% 
  mutate_at(.vars = var.names2,.funs = function(X){
    if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
  }) %>% 
  mutate_at(.vars = var.names,.funs = scale,center=T,scale=T)

if(analysis.framework$global.var[args] != "none"){
  
  cbcl_dat = cbcl_dat %>% 
    mutate_at(.vars = "global.var",.funs = winsorize) %>% 
    mutate_at(.vars = "global.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "global.var",.funs = scale,center=T,scale=T)
  
}


if(analysis.framework$movement.var[args] != "none"){
  
  cbcl_dat = cbcl_dat %>% 
    mutate_at(.vars = "movement.var",.funs = winsorize) %>% 
    mutate_at(.vars = "movement.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "movement.var",.funs = scale,center=T,scale=T)
  
}

cbcl_dat$Age1 = poly(cbcl_dat$interview_age,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat$Age2 = poly(cbcl_dat$interview_age,2)[,2] %>% scale(center = T,scale = T)


run_lmer = function(){lmer(y.var ~ mj_pre+mj_pre_post + 
                             Age1+Age2+
                             rel_relationship_1+rel_relationship_2+rel_relationship_3+
                             ed+income_1+income_2+income_3+income_4+pds+sex2 +
                             famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                             famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                             demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                             dgs_before+dgs_after+
                             #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                             
                             devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                             Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                             (1|site_id_l) + (1|rel_family_id) +(1|src_subject_id)
                           ,data =cbcl_dat  )}


run_lmer_base = function(){lmer(y.var ~ #mj_pre+mj_pre_post + 
                                  Age1+Age2+
                                  rel_relationship_1+rel_relationship_2+rel_relationship_3+
                                  ed+income_1+income_2+income_3+income_4+pds+sex2 +
                                  famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                                  famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+
                                  famhx_ss_firstdeg_nrv_p +
                                  demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                                  dgs_before+dgs_after+
                                  #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                                  
                                  devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                                  Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                                  (1|site_id_l) + (1|rel_family_id) +(1|src_subject_id)
                                ,data =cbcl_dat  )}

f1 = tryCatch(expr = run_lmer(),error = run_lmer()) # if error (b/c of memory), try again
f2 = tryCatch(expr = run_lmer_base(),error = run_lmer_base()) # if error (b/c of memory), try again

a1 = tryCatch(expr = anova(f2,f1),error = anova(f2,f1)) # if error (b/c of memory), try again

l1 = summary(f1)

if(a1$`Pr(>Chisq)`[2] < 0.05){
  c2 = confint(f1,parm = c("mj_pre","mj_pre_post"))
}else{
  c2 = matrix(NA,nrow = 2,ncol = 2)  
  rownames(c2) = c("mj_pre","mj_pre_post")
}

##########


##########

fixed_effs=data.frame(Y = analysis.framework$y.var[args])

fixed_effs2 =  c(
  l1$coefficients[which(rownames(l1$coefficients) == "mj_pre"),],
  l1$coefficients[which(rownames(l1$coefficients) == "mj_pre_post"),],
  
  a1$Chisq[2],
  a1$`Pr(>Chisq)`[2],
  c2[which(rownames(c2) == "mj_pre"),],
  c2[which(rownames(c2) == "mj_pre_post"),]
  
) %>% t() %>% as.data.frame() 

fixed_effs = cbind(fixed_effs,fixed_effs2)
## rm(a4)

colnames(fixed_effs) = c("Y",
                         paste("mj_pre_main",colnames(l1$coefficients),sep = "_"),
                         paste("mj_pre_post_main",colnames(l1$coefficients),sep = "_"),
                         "Chisq_main","Chisq_P_main",
                         "L_mj_pre","U_mj_pre",
                         "L_mj_pre_post","U_mj_pre_post"
                        
)  



fixed_effs$N = length(unique(cbcl_dat$src_subject_id))
fixed_effs$Nobs = length((cbcl_dat$src_subject_id))
fixed_effs$N_1 = table(table(cbcl_dat$src_subject_id))[1] %>% unname()
fixed_effs$N_2 = table(table(cbcl_dat$src_subject_id))[2] %>% unname()
fixed_effs$N_pre = cbcl_dat %>% dplyr::select(c("src_subject_id","mj_pre")) %>% unique() %>%  
  dplyr::select(c("mj_pre")) %>% summarise(n=table(mj_pre)[2]) %>% unname()
fixed_effs$N_pre_post = cbcl_dat %>% dplyr::select(c("src_subject_id","mj_pre_post")) %>% unique() %>%  dplyr::select(c("mj_pre_post")) %>% summarise(n=table(mj_pre_post)[2]) %>% unname()

##########################################################
# extra covariates

cbcl_dat = cbcl_long %>% dplyr::select(all_of(
  c(
    "src_subject_id","eventname","site_id_l","rel_family_id",
    "mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
    "ed","pds","famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
    "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
    "dgs_after","income","interview_age",
    "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
    "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
    "income_1","income_2","income_3","income_4",
    "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
    "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
    "Prisma_fit","Prisma","DISCOVERY","Achieva",
    "y.var","global.var","include.var","movement.var"
    ,"mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth"
    # excluding in some analyses -  mat_age_knew_preg, total_birth_weight, devhx_10, devhx_6_p, income      
  ))) %>% dplyr::filter(include.var == 1) %>% unique()
cbcl_dat = na.omit(cbcl_dat)


var.names=c( "mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
             "ed","pds","famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
             "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
             "dgs_after","income","interview_age",
             "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
             "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
             "income_1","income_2","income_3","income_4",
             "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
             "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
             "Prisma_fit","Prisma","DISCOVERY","Achieva",
             "mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth",
             "y.var")

var.names2=c( #"mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
  "ed","pds",
  #"famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
  #"demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
  #"dgs_after",
  "income","interview_age",
  #"famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
  #"famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
  #"income_1","income_2","income_3","income_4",
  #"famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
  #"devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
  #"Prisma_fit","Prisma","DISCOVERY","Achieva",
  "mat_age_knew_preg", "total_birth_weight","mat_age_birth",
  "y.var")


cbcl_dat = cbcl_dat %>% 
  mutate_at(.vars = var.names2,.funs = winsorize) %>% 
  mutate_at(.vars = var.names2,.funs = function(X){
    if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
  }) %>% 
  mutate_at(.vars = var.names,.funs = scale,center=T,scale=T)

if(analysis.framework$global.var[args] != "none"){
  
  cbcl_dat = cbcl_dat %>% 
    mutate_at(.vars = "global.var",.funs = winsorize) %>% 
    mutate_at(.vars = "global.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "global.var",.funs = scale,center=T,scale=T)
  
}



if(analysis.framework$movement.var[args] != "none"){
  
  cbcl_dat = cbcl_dat %>% 
    mutate_at(.vars = "movement.var",.funs = winsorize) %>% 
    mutate_at(.vars = "movement.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "movement.var",.funs = scale,center=T,scale=T)
  
}

cbcl_dat$Age1 = poly(cbcl_dat$interview_age,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat$Age2 = poly(cbcl_dat$interview_age,2)[,2] %>% scale(center = T,scale = T)


run_lmer = function(){lmer(y.var ~ mj_pre+mj_pre_post + 
                             Age1+Age2+
                             rel_relationship_1+rel_relationship_2+rel_relationship_3+
                             ed+income_1+income_2+income_3+income_4+pds+sex2 +
                             famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                             famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                             demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                             dgs_before+dgs_after+
                             mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                             
                             devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                             Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                             (1|site_id_l) + (1|rel_family_id) +(1|src_subject_id)
                           ,data =cbcl_dat  )}


run_lmer_base = function(){lmer(y.var ~ #mj_pre+mj_pre_post + 
                                  Age1+Age2+
                                  rel_relationship_1+rel_relationship_2+rel_relationship_3+
                                  ed+income_1+income_2+income_3+income_4+pds+sex2 +
                                  famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                                  famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                                  demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                                  dgs_before+dgs_after+
                                  mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                                  
                                  devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                                  Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                                  (1|site_id_l) + (1|rel_family_id) +(1|src_subject_id)
                                ,data =cbcl_dat  )}

f1 = tryCatch(expr = run_lmer(),error = run_lmer()) # if error (b/c of memory), try again
f2 = tryCatch(expr = run_lmer_base(),error = run_lmer_base()) # if error (b/c of memory), try again

a1 = tryCatch(expr = anova(f2,f1),error = anova(f2,f1)) # if error (b/c of memory), try again

l1 = summary(f1)
#########################################

##############################################################

fixed_effsb =  c(
  l1$coefficients[which(rownames(l1$coefficients) == "mj_pre"),],
  l1$coefficients[which(rownames(l1$coefficients) == "mj_pre_post"),],
  
  a1$Chisq[2],
  a1$`Pr(>Chisq)`[2]
  
) %>% t() %>% as.data.frame() 

## rm(a4)

colnames(fixed_effsb) = c(
  paste("mj_pre_extra",colnames(l1$coefficients),sep = "_"),
  paste("mj_pre_post_extra",colnames(l1$coefficients),sep = "_"),
  "Chisq_extra","Chisq_P_extra"
)  

fixed_effsb$N_extra = length(unique(cbcl_dat$src_subject_id))
fixed_effsb$Nobs_extra = length((cbcl_dat$src_subject_id))

fixed_effs = cbind(fixed_effs,fixed_effsb)

###########################################
# genetic control

cbcl_dat = cbcl_long %>% dplyr::select(all_of(
  c(
    "src_subject_id","eventname","site_id_l","rel_family_id",
    "mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
    "ed","pds","famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
    # "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h",
    "dgs_before",
    "dgs_after","income","interview_age",
    "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
    "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
    "income_1","income_2","income_3","income_4",
    "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
    "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
    "Prisma_fit","Prisma","DISCOVERY","Achieva",
    "y.var","global.var","include.var","movement.var",
    colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PRSCS")],
    colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PC")]
    #,"mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth"
    # excluding in some analyses -  mat_age_knew_preg, total_birth_weight, devhx_10, devhx_6_p, income      
  ))) %>% dplyr::filter(include.var == 1) %>% unique() %>% na.omit()


var.names=c( "mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
             "ed","pds","famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
             #"demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h",
             "dgs_before",
             "dgs_after","income","interview_age",
             "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
             "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
             "income_1","income_2","income_3","income_4",
             "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
             "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
             "Prisma_fit","Prisma","DISCOVERY","Achieva",
             "y.var",colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PRSCS")],
             colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PC")])

var.names2=c( #"mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
  "ed","pds",
  #"famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
  #"demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
  #"dgs_after",
  "income","interview_age",
  #"famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
  #"famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
  #"income_1","income_2","income_3","income_4",
  #"famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
  #"devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
  #"Prisma_fit","Prisma","DISCOVERY","Achieva",
  "y.var","CUD_PRSCS",
  colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PRSCS")],
  colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PC")])


cbcl_dat = cbcl_dat %>% 
  mutate_at(.vars = var.names2,.funs = winsorize) %>% 
  mutate_at(.vars = var.names2,.funs = function(X){
    if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
  }) %>% 
  mutate_at(.vars = var.names,.funs = scale,center=T,scale=T)

if(analysis.framework$global.var[args] != "none"){
  
  cbcl_dat = cbcl_dat %>% 
    mutate_at(.vars = "global.var",.funs = winsorize) %>% 
    mutate_at(.vars = "global.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "global.var",.funs = scale,center=T,scale=T)
  
}

if(analysis.framework$movement.var[args] != "none"){
  
  cbcl_dat = cbcl_dat %>% 
    mutate_at(.vars = "movement.var",.funs = winsorize) %>% 
    mutate_at(.vars = "movement.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "movement.var",.funs = scale,center=T,scale=T)
  
}

cbcl_dat$Age1 = poly(cbcl_dat$interview_age,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat$Age2 = poly(cbcl_dat$interview_age,2)[,2] %>% scale(center = T,scale = T)


run_lmer = function(){lmer(y.var ~ mj_pre+mj_pre_post + 
                             Age1+Age2+
                             rel_relationship_1+rel_relationship_2+rel_relationship_3+
                             ed+income_1+income_2+income_3+income_4+pds+sex2 +
                             famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                             famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+
                             famhx_ss_firstdeg_nrv_p +
                             #  demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                             dgs_before+dgs_after+
                             #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                             
                             devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                             Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                             CUD_PRSCS + 
                             PC_1+PC_2+PC_3+PC_4+PC_5+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10+
                             (1|site_id_l) + (1|rel_family_id) +(1|src_subject_id)
                           ,data =cbcl_dat  )}


run_lmer_base = function(){lmer(y.var ~ #mj_pre+mj_pre_post + 
                                  Age1+Age2+
                                  rel_relationship_1+rel_relationship_2+rel_relationship_3+
                                  ed+income_1+income_2+income_3+income_4+pds+sex2 +
                                  famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                                  famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                                  # demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                                  dgs_before+dgs_after+
                                  #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                                  
                                  devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                                  Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                                  CUD_PRSCS + 
                                  PC_1+PC_2+PC_3+PC_4+PC_5+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10+
                                  (1|site_id_l) + (1|rel_family_id) +(1|src_subject_id)
                                ,data =cbcl_dat  )}

f1 = tryCatch(expr = run_lmer(),error = run_lmer()) # if error (b/c of memory), try again
f2 = tryCatch(expr = run_lmer_base(),error = run_lmer_base()) # if error (b/c of memory), try again

a1 = tryCatch(expr = anova(f2,f1),error = anova(f2,f1)) # if error (b/c of memory), try again

l1 = summary(f1)


#############

fixed_effsc =  c(
  l1$coefficients[which(rownames(l1$coefficients) == "mj_pre"),],
  l1$coefficients[which(rownames(l1$coefficients) == "mj_pre_post"),],
  
  a1$Chisq[2],
  a1$`Pr(>Chisq)`[2]
  
  
) %>% t() %>% as.data.frame() 

## rm(a4)

colnames(fixed_effsc) = c(
  paste("mj_pre_genetics",colnames(l1$coefficients),sep = "_"),
  paste("mj_pre_post_genetics",colnames(l1$coefficients),sep = "_"),
  "Chisq_genetics","Chisq_P_genetics"
  
)  

fixed_effsc$N_genetics = length(unique(cbcl_dat$src_subject_id))
fixed_effsc$Nobs_genetics = length((cbcl_dat$src_subject_id))

fixed_effs = cbind(fixed_effs,fixed_effsc)


#################################
#basic model

cbcl_dat = cbcl_long %>% dplyr::select(all_of(
  c(
    "src_subject_id","eventname","site_id_l","rel_family_id",
    "mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
    "sex2","pds","interview_age",
    "Prisma_fit","Prisma","DISCOVERY","Achieva",
    "y.var","global.var","include.var","movement.var"
    
  ))) %>% dplyr::filter(include.var == 1) %>% unique()
cbcl_dat = na.omit(cbcl_dat)


var.names=c( "mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
             "sex2","pds","interview_age",
             "Prisma_fit","Prisma","DISCOVERY","Achieva",
             "y.var")

var.names2=c( #"mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
  #"ed",
  #"famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
  #"demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
  #"dgs_after",
  #"income",
  "interview_age","pds",
  #"famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
  #"famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
  #"income_1","income_2","income_3","income_4",
  #"famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
  #"devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
  #"Prisma_fit","Prisma","DISCOVERY","Achieva",
  "y.var")


cbcl_dat = cbcl_dat %>% 
  mutate_at(.vars = var.names2,.funs = winsorize) %>% 
  mutate_at(.vars = var.names2,.funs = function(X){
    if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
  }) %>% 
  mutate_at(.vars = var.names,.funs = scale,center=T,scale=T)

if(analysis.framework$global.var[args] != "none"){
  
  cbcl_dat = cbcl_dat %>% 
    mutate_at(.vars = "global.var",.funs = winsorize) %>% 
    mutate_at(.vars = "global.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "global.var",.funs = scale,center=T,scale=T)
  
}



if(analysis.framework$movement.var[args] != "none"){
  
  cbcl_dat = cbcl_dat %>% 
    mutate_at(.vars = "movement.var",.funs = winsorize) %>% 
    mutate_at(.vars = "movement.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "movement.var",.funs = scale,center=T,scale=T)
  
}

cbcl_dat$Age1 = poly(cbcl_dat$interview_age,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat$Age2 = poly(cbcl_dat$interview_age,2)[,2] %>% scale(center = T,scale = T)


run_lmer = function(){lmer(y.var ~ mj_pre+mj_pre_post + 
                             Age1+Age2+
                             rel_relationship_1+rel_relationship_2+rel_relationship_3+
                             #ed+income_1+income_2+income_3+income_4+
                             pds+sex2 +
                             #famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                             #famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                             #demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                             #dgs_before+dgs_after+
                             #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                             
                             #devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                             Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                             (1|site_id_l) + (1|rel_family_id) +(1|src_subject_id)
                           ,data =cbcl_dat  )}


run_lmer_base = function(){lmer(y.var ~ #mj_pre+mj_pre_post + 
                                  Age1+Age2+
                                  rel_relationship_1+rel_relationship_2+rel_relationship_3+
                                  pds+sex2 +
                                  Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                                  (1|site_id_l) + (1|rel_family_id) +(1|src_subject_id)
                                ,data =cbcl_dat  )}

f1 = tryCatch(expr = run_lmer(),error = run_lmer()) # if error (b/c of memory), try again
f2 = tryCatch(expr = run_lmer_base(),error = run_lmer_base()) # if error (b/c of memory), try again

a1 = tryCatch(expr = anova(f2,f1),error = anova(f2,f1)) # if error (b/c of memory), try again

l1 = summary(f1)
################

fixed_effsd =  c(
  l1$coefficients[which(rownames(l1$coefficients) == "mj_pre"),],
  l1$coefficients[which(rownames(l1$coefficients) == "mj_pre_post"),],
  
  a1$Chisq[2],
  a1$`Pr(>Chisq)`[2]
  
  
  
) %>% t() %>% as.data.frame() 

## rm(a4)

colnames(fixed_effsd) = c(
  paste("mj_pre_basic",colnames(l1$coefficients),sep = "_"),
  paste("mj_pre_post_basic",colnames(l1$coefficients),sep = "_"),
  "Chisq_basic","Chisq_P_basic"
  
)  



fixed_effsd$N_basic = length(unique(cbcl_dat$src_subject_id))
fixed_effsd$Nobs_basic = length((cbcl_dat$src_subject_id))

fixed_effs = cbind(fixed_effs,fixed_effsd)



######################


cbcl_dat = cbcl_long %>% dplyr::select(all_of(
  c(
    "src_subject_id","eventname","site_id_l","rel_family_id",
    "mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
    "ed","pds","famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
    "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
    "dgs_after","income","interview_age",
    "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
    "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
    "income_1","income_2","income_3","income_4",
    "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
    "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
    "Prisma_fit","Prisma","DISCOVERY","Achieva",
    "y.var","global.var","include.var","movement.var"
    #,"mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth"
    # excluding in some analyses -  mat_age_knew_preg, total_birth_weight, devhx_10, devhx_6_p, income      
  ))) %>%
  dplyr::filter(include.var == 1,mj_pre == 1 | mj_pre_post == 1 ) %>% 
  na.omit() %>% unique()



var.names=c( "mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
             "ed","pds","famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
             "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
             "dgs_after","income","interview_age",
             "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
             "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
             "income_1","income_2","income_3","income_4",
             "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
             "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
             "Prisma_fit","Prisma","DISCOVERY","Achieva",
             "y.var")

var.names2=c( #"mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
  "ed","pds",
  #"famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
  #"demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
  #"dgs_after",
  "income","interview_age",
  #"famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
  #"famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
  #"income_1","income_2","income_3","income_4",
  #"famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
  #"devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
  #"Prisma_fit","Prisma","DISCOVERY","Achieva",
  "y.var")


cbcl_dat = cbcl_dat %>% 
  mutate_at(.vars = var.names2,.funs = winsorize) %>% 
  mutate_at(.vars = var.names2,.funs = function(X){
    if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
  }) %>% 
  mutate_at(.vars = var.names,.funs = scale,center=T,scale=T)

if(analysis.framework$global.var[args] != "none"){
  
  cbcl_dat = cbcl_dat %>% 
    mutate_at(.vars = "global.var",.funs = winsorize) %>% 
    mutate_at(.vars = "global.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "global.var",.funs = scale,center=T,scale=T)
  
}


if(analysis.framework$movement.var[args] != "none"){
  
  cbcl_dat = cbcl_dat %>% 
    mutate_at(.vars = "movement.var",.funs = winsorize) %>% 
    mutate_at(.vars = "movement.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "movement.var",.funs = scale,center=T,scale=T)
  
}

cbcl_dat$Age1 = poly(cbcl_dat$interview_age,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat$Age2 = poly(cbcl_dat$interview_age,2)[,2] %>% scale(center = T,scale = T)


run_lmer = function(){lmer(y.var ~ mj_pre+ 
                             Age1+Age2+
                             rel_relationship_1+rel_relationship_2+
                             ed+income_1+income_2+income_3+income_4+pds+sex2 +
                             famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                             famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                             demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                             dgs_before+dgs_after+
                             #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                             
                             devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                             Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                             (1|site_id_l) + (1|rel_family_id) +(1|src_subject_id)
                           ,data =cbcl_dat  )}


f1 = tryCatch(expr = run_lmer(),error = run_lmer()) # if error (b/c of memory), try again

l1 = summary(f1)

###############

fixed_effsg =  c(
  l1$coefficients[which(rownames(l1$coefficients) == "mj_pre"),]
  
) %>% t() %>% as.data.frame() 


colnames(fixed_effsg) = c(
  paste("mj_pre_pce_only",colnames(l1$coefficients),sep = "_")
  
  
)  



fixed_effsg$N_pce_only = length(unique(cbcl_dat$src_subject_id))
fixed_effsg$Nobs_pce_only = length((cbcl_dat$src_subject_id))

fixed_effs = cbind(fixed_effs,fixed_effsg)
#########################

##################
#baseline model

cbcl_dat = cbcl_long %>% dplyr::select(all_of(
  c(
    "src_subject_id","eventname","site_id_l","rel_family_id","MJ_exposure",
    "mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
    "ed","pds","famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
    "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
    "dgs_after",#"try_dg",
    "income","interview_age",
    "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
    "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
    "income_1","income_2","income_3","income_4",
    "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
    "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
    "Prisma_fit","Prisma","DISCOVERY","Achieva",
    "y.var","global.var","include.var","movement.var"
    #,"mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth"
    # excluding in some analyses -  mat_age_knew_preg, total_birth_weight, devhx_10, devhx_6_p, income      
  ))) %>% dplyr::filter(include.var == 1,eventname == "baseline_year_1_arm_1") %>% unique() %>% na.omit()



var.names=c( "mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
             "ed","pds","famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
             "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
             "dgs_after",#"try_dg",
             "income","interview_age",
             "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
             "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
             "income_1","income_2","income_3","income_4",
             "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
             "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
             "Prisma_fit","Prisma","DISCOVERY","Achieva",
             "y.var")

var.names2=c( #"mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
  "ed","pds",
  #"famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
  #"demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
  #"dgs_after",
  "income","interview_age",
  #"famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
  #"famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
  #"income_1","income_2","income_3","income_4",
  #"famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
  #"devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
  #"Prisma_fit","Prisma","DISCOVERY","Achieva",
  "y.var")


cbcl_dat = cbcl_dat %>% 
  mutate_at(.vars = var.names2,.funs = winsorize) %>% 
  mutate_at(.vars = var.names2,.funs = function(X){
    if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
  }) %>% 
  mutate_at(.vars = var.names,.funs = scale,center=T,scale=T)

if(analysis.framework$global.var[args] != "none"){
  
  cbcl_dat = cbcl_dat %>% 
    mutate_at(.vars = "global.var",.funs = winsorize) %>% 
    mutate_at(.vars = "global.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "global.var",.funs = scale,center=T,scale=T)
  
}


if(analysis.framework$movement.var[args] != "none"){
  
  cbcl_dat = cbcl_dat %>% 
    mutate_at(.vars = "movement.var",.funs = winsorize) %>% 
    mutate_at(.vars = "movement.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "movement.var",.funs = scale,center=T,scale=T)
  
}

cbcl_dat$Age1 = poly(cbcl_dat$interview_age,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat$Age2 = poly(cbcl_dat$interview_age,2)[,2] %>% scale(center = T,scale = T)


run_lmer = function(){lmer(y.var ~ mj_pre+mj_pre_post + 
                             Age1+Age2+
                             rel_relationship_1+rel_relationship_2+rel_relationship_3+
                             ed+income_1+income_2+income_3+income_4+pds+sex2 +
                             famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                             famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                             demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                             dgs_before+dgs_after+
                             #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                             
                             devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                             Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                             (1|site_id_l) + (1|rel_family_id) 
                           ,data =cbcl_dat  )}


run_lmer_base = function(){lmer(y.var ~ #mj_pre+mj_pre_post + 
                                  Age1+Age2+
                                  rel_relationship_1+rel_relationship_2+rel_relationship_3+
                                  ed+income_1+income_2+income_3+income_4+pds+sex2 +
                                  famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                                  famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                                  demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                                  dgs_before+dgs_after+
                                  #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                                  
                                  devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                                  Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                                  (1|site_id_l) + (1|rel_family_id) 
                                ,data =cbcl_dat  )}

f1 = tryCatch(expr = run_lmer(),error = run_lmer()) # if error (b/c of memory), try again
f2 = tryCatch(expr = run_lmer_base(),error = run_lmer_base()) # if error (b/c of memory), try again

a1 = tryCatch(expr = anova(f2,f1),error = anova(f2,f1)) # if error (b/c of memory), try again

l1 = summary(f1)


fixed_effsh =  c(
  l1$coefficients[which(rownames(l1$coefficients) == "mj_pre"),],
  l1$coefficients[which(rownames(l1$coefficients) == "mj_pre_post"),],
  
  a1$Chisq[2],a1$`Pr(>Chisq)`[2]
  
) %>% t() %>% as.data.frame() 

## rm(a4)

colnames(fixed_effsh) = c(
                         paste("mj_pre_baseline",colnames(l1$coefficients),sep = "_"),
                         paste("mj_pre_post_baseline",colnames(l1$coefficients),sep = "_"),
                         "Chisq_baseline","Chisq_P_baseline"
                         
)  



fixed_effsh$N_baseline = length(unique(cbcl_dat$src_subject_id))

fixed_effs = cbind(fixed_effs,fixed_effsh)

##################################

write.table(x = fixed_effs,file = paste(file_location,"PCE_regression.",args[1],".csv",sep=""),quote = F,row.names = F,col.names = F,sep = ",")

if(args[1] == 1){
  
  write.table(x = t(colnames(fixed_effs)),file = paste(file_location,"PCE_regression.0.csv",sep=""),quote = F,row.names = F,col.names = F,sep = ",")
  
}


