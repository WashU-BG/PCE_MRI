#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
args=as.numeric(args)
args=args[1]

library(data.table)
library(dplyr)
library(psych)
#library(lmerTest)
library(lme4)
library(tidyr)
library(stringr)
library(mediation)
library(fastDummies)

file_location = "C:/Users/David/Documents/BRAINlab/ABCD/Package_1195434/ABCD_5/core/"

cbcl_long = fread(paste(file_location,'allimaging_pce_mediation.csv',sep=""),data.table = F)


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


m.vars = c("dmri_rsirnt_fib_fmin",
           "dmri_rsirnd_fib_fmin",
           "rsfmri_cor_ngd_au_scs_ptlh",
           "dmri_rsirnt_fib_fmin",
           "dmri_rsirnd_fib_fmin",
           "rsfmri_cor_ngd_au_scs_ptlh",
           "dmri_rsirnt_fib_fmin",
           "dmri_dtitdgm_cdsn_ptisrh",
           "dmri_dtitdgm_cdsn_ptisrh",
           "dmri_rsirnd_fib_pslflh",
           "dmri_dtimdgm_cdsn_parstangrh",
           "dmri_dtimdgm_cdsn_parstangrh",
           "dmri_dtimdgm_cdsn_parstangrh",
           "dmri_rsirnd_fib_fmin",
           "dmri_rsirnd_fib_fmin",
           "dmri_dtitdgm_cdsn_ptisrh",
           "dmri_dtimdgm_cdsn_parstangrh",
           "dmri_dtimdgm_cdsn_parstangrh"
           
)

y.vars = c("cbcl_scr_dsm5_adhd_r",
           "cbcl_scr_dsm5_adhd_r",
           "pps_y_ss_number",
           "cbcl_scr_syn_attention_r",
           "cbcl_scr_syn_attention_r",
           "cbcl_scr_dsm5_conduct_r",
           "cbcl_scr_syn_social_r",
           "cbcl_scr_syn_attention_r",
           "cbcl_scr_dsm5_adhd_r",
           "pps_y_ss_number",
           "cbcl_scr_syn_attention_r",
           "cbcl_scr_dsm5_conduct_r",
           "cbcl_scr_dsm5_adhd_r",
           "cbcl_scr_dsm5_conduct_r",
           "cbcl_scr_syn_social_r",
           "cbcl_scr_dsm5_conduct_r",
           "pps_y_ss_number",
           "cbcl_scr_syn_rulebreak_r")


analysis.framework=data.frame(stringsAsFactors = F,
                               
                              m.var = m.vars,
                               y.var = y.vars,
                               
                               global.var = "none",
                               include.var=NA,
                               movement.var = NA
                               
)



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


analysis.framework$include.var[grep(analysis.framework$m.var,pattern = "smri")] = "imgincl_t1w_include"
analysis.framework$include.var[grep(analysis.framework$m.var,pattern = "rsfmri")] = "imgincl_rsfmri_include"
analysis.framework$include.var[grep(analysis.framework$m.var,pattern = "dmri")] = "imgincl_dmri_include"


analysis.framework$movement.var[grep(analysis.framework$m.var,pattern = "smri")] = "none"
analysis.framework$movement.var[grep(analysis.framework$m.var,pattern = "rsfmri")] = "rsfmri_meanmotion"
analysis.framework$movement.var[grep(analysis.framework$m.var,pattern = "dmri")] = "dmri_meanmotion"


analysis.framework$global.var[grep(analysis.framework$m.var,pattern = "rsfmri")] = "none"


for(i in 1:length(global.vars)){
  v.temp = global.vars[i]
  new.temp = str_split(string = v.temp,pattern = "_",simplify =T)[-4] %>% paste(collapse = "_")
  analysis.framework$global.var[grep(analysis.framework$m.var,pattern = new.temp)] = v.temp
}


analysis.framework$global.var[grep(analysis.framework$m.var,pattern = "smri_vol")] = "smri_vol_scs_intracranialv"
analysis.framework$global.var[analysis.framework$global.var == analysis.framework$m.var] = "none"



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


cbcl_long$rel_relationship_1 = ifelse(test = cbcl_long$rel_relationship == 1,yes = 1,no = 0)
cbcl_long$rel_relationship_2 = ifelse(test = cbcl_long$rel_relationship == 2,yes = 1,no = 0)
cbcl_long$rel_relationship_3 = ifelse(test = cbcl_long$rel_relationship == 3,yes = 1,no = 0)



################
################
cbcl_long$m.var = cbcl_long[,which(colnames(cbcl_long) == analysis.framework$m.var[args])]
cbcl_long$y.var = cbcl_long[,which(colnames(cbcl_long) == analysis.framework$y.var[args])]
cbcl_long$global.var = cbcl_long[,which(colnames(cbcl_long) == analysis.framework$global.var[args])]
cbcl_long$include.var = cbcl_long[,which(colnames(cbcl_long) == analysis.framework$include.var[args])]
cbcl_long$movement.var = cbcl_long[,which(colnames(cbcl_long) == analysis.framework$movement.var[args])]

######################################################################################

cbcl_long$mj_pre_any = cbcl_long$mj_pre + cbcl_long$mj_pre_post
cbcl_long$mj_group = 0
cbcl_long$mj_group[cbcl_long$mj_pre == 1] = 1
cbcl_long$mj_group[cbcl_long$mj_pre_post == 1] = 2
cbcl_long$mj_group = factor(cbcl_long$mj_group)


site.dummies = fastDummies::dummy_columns(cbcl_long$site_id_l,remove_most_frequent_dummy = F,ignore_na = T)
colnames(site.dummies)[-1]=paste("site",c(1:21),sep="_")
site.dummies = site.dummies[,-1]

cbcl_long = cbind(cbcl_long,site.dummies)
######################################################################################
# Main

print("main")

cbcl_dat_base = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id","site_id_l","rel_family_id","eventname",
      paste("site",c(1:21),sep="_"),
      "mj_pre","mj_pre_post","mj_pre_any","mj_group",
      "rel_relationship_1","rel_relationship_2","rel_relationship_3",
      "sex2",
      "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h",
      
      "dgs_before", "dgs_after",
      
      "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
      "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
      
      "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p", "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p",
      "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco"#,
      #"income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      #"Prisma_fit","Prisma","DISCOVERY","Achieva",
      #"y.var","global.var","include.var","movement.var","x.var"
      #,"mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth"
      # excluding in some analyses -  mat_age_knew_preg, total_birth_weight, devhx_10, devhx_6_p, income      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1") %>% 
  na.omit() %>% unique() 
cbcl_dat_base = cbcl_dat_base[,-which(colnames(cbcl_dat_base) == "eventname")]


cbcl_dat_mediator = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      "Prisma_fit","Prisma","DISCOVERY","Achieva",
      "global.var","include.var","movement.var","m.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1" | eventname == "2_year_follow_up_y_arm_1", include.var == 1) %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      
      "y.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "1_year_follow_up_y_arm_1" | eventname == "3_year_follow_up_y_arm_1") %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome$eventname[cbcl_dat_outcome$eventname == "1_year_follow_up_y_arm_1"] = "baseline_year_1_arm_1"
cbcl_dat_outcome$eventname[cbcl_dat_outcome$eventname == "3_year_follow_up_y_arm_1"] = "2_year_follow_up_y_arm_1"

cbcl_dat_outcome$combo_id = paste(cbcl_dat_outcome$src_subject_id,cbcl_dat_outcome$eventname,sep = "_")
cbcl_dat_mediator$combo_id = paste(cbcl_dat_mediator$src_subject_id,cbcl_dat_mediator$eventname,sep = "_")


cbcl_dat_all = merge(cbcl_dat_base,cbcl_dat_mediator)
cbcl_dat_all = merge(cbcl_dat_all,cbcl_dat_outcome,by = "combo_id")

sum(cbcl_dat_all$src_subject_id.x == cbcl_dat_all$src_subject_id.y)
sum(cbcl_dat_all$eventname.x == cbcl_dat_all$eventname.y)

#cbcl_dat_all = cbcl_dat_all %>% filter(eventname.x == "baseline_year_1_arm_1")

length(unique(cbcl_dat_all$src_subject_id.x))

temp.dat1 = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) %>% unique()
temp.dat = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) 

famids = table(temp.dat1$rel_family_id) %>% as.data.frame() 

# remove relatives, keeping any PCE participants

for(f in 1 :dim(famids)[1]){
  
  matching.fam = which(temp.dat$rel_family_id == famids$Var1[f])
  
  t2 = temp.dat[matching.fam,]
  
  t3 = table(t2$src_subject_id.x) %>% as.data.frame()
  colnames(t3) = c("src_subject_id.x","Freq")
  t2 = merge(t2,t3) %>% unique()
  
  
  if(sum(t2$mj_pre_any)>0){ t2 = t2[t2$mj_pre_any==1,]  }
  
  if(max(t2$Freq == 2)){ t2 = t2[t2$Freq==2,] }
  
  keep = t2[sample(x = c(1:dim(t2)[1]),size = 1),]
  

  if(f == 1){keep.ids = keep$src_subject_id.x}else{keep.ids = c(keep.ids,keep$src_subject_id.x)}

}

cbcl_dat_all = cbcl_dat_all [!is.na(match(cbcl_dat_all$src_subject_id.x,keep.ids)),]

length(unique(cbcl_dat_all$src_subject_id.x))
table(cbcl_dat_all$mj_group)

var.names=c( "mj_pre","mj_pre_post","mj_pre_any","rel_relationship_1","rel_relationship_2","rel_relationship_3",
             "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
             "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
             "dgs_after",#"income","interview_age",
             "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
             "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
             
             "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
             "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
             "Prisma_fit","Prisma","DISCOVERY","Achieva","m.var",
             "y.var", "interview_age.x","interview_age.y",
             paste("site",c(2:21),sep="_"),
             "income_1.x","income_2.x","income_3.x","income_4.x","ed.x","pds.x",
             "income_1.y","income_2.y","income_3.y","income_4.y","ed.y","pds.y")

var.names2=c( #"mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
  "ed.x","pds.x","ed.y","pds.y",
  #"famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
  #"demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
  #"dgs_after",
  "interview_age.x","interview_age.y",
  #"famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
  #"famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
  #"income_1","income_2","income_3","income_4",
  #"famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
  #"devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
  #"Prisma_fit","Prisma","DISCOVERY","Achieva",
  "y.var","m.var")


cbcl_dat_all = cbcl_dat_all %>% 
  mutate_at(.vars = var.names2,.funs = winsorize) %>% 
  mutate_at(.vars = var.names2,.funs = function(X){
    if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
  }) %>% 
  mutate_at(.vars = var.names,.funs = scale,center=T,scale=T)

if(analysis.framework$global.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "global.var",.funs = winsorize) %>% 
    mutate_at(.vars = "global.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "global.var",.funs = scale,center=T,scale=T)
  
}


if(analysis.framework$movement.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "movement.var",.funs = winsorize) %>% 
    mutate_at(.vars = "movement.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "movement.var",.funs = scale,center=T,scale=T)
  
}

cbcl_dat_all$Age1.x = poly(cbcl_dat_all$interview_age.x,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.x = poly(cbcl_dat_all$interview_age.x,2)[,2] %>% scale(center = T,scale = T)
cbcl_dat_all$Age1.y = poly(cbcl_dat_all$interview_age.y,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.y = poly(cbcl_dat_all$interview_age.y,2)[,2] %>% scale(center = T,scale = T)


run_med = function(){lmer(m.var ~ mj_pre_any + 
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            ed.y+income_1.y+income_2.y+income_3.y+income_4.y+
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            site_2+site_3+site_4+site_5+site_6+site_7+site_8+site_9+site_10+site_11+
                            site_12+site_13+site_14+site_15+site_16+site_17+site_18+site_19+site_20+site_21+
                            #             (1|site_id_l) #+ 
                            #(1|rel_family_id)# 
                            (1|src_subject_id.x)
                          ,data =cbcl_dat_all  )}


med.fit = tryCatch(expr = run_med(),error = run_med()) # if error (b/c of memory), try again

run_out = function(){lmer(y.var ~ mj_pre_any + m.var+
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            ed.y+income_1.y+income_2.y+income_3.y+income_4.y+
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            site_2+site_3+site_4+site_5+site_6+site_7+site_8+site_9+site_10+site_11+
                            site_12+site_13+site_14+site_15+site_16+site_17+site_18+site_19+site_20+site_21+
                            #  (1|site_id_l)  
                            #(1|rel_family_id) #+
                            (1|src_subject_id.x)
                          ,data =cbcl_dat_all  )}


out.fit = tryCatch(expr = run_out(),error = run_out()) # if error (b/c of memory), try again


med.out <- mediate(model.m = med.fit, 
                   model.y = out.fit, 
                   treat = "mj_pre_any", 
                   mediator = "m.var",
                   sims = 10000)

out.sum.any = c(med.out$d.avg, med.out$d.avg.ci[1],med.out$d.avg.ci[2],med.out$d.avg.p,
                med.out$z.avg, med.out$z.avg.ci[1],med.out$z.avg.ci[2],med.out$z.avg.p,
                med.out$tau.coef,med.out$tau.ci[1],med.out$tau.ci[2],med.out$tau.p,
                med.out$n.avg, med.out$n.avg.ci[1],med.out$n.avg.ci[2],med.out$n.avg.p) %>% t() %>% as.data.frame()
col.st = expand.grid(a = "main","Any",c =c("Est","L","U","p") ,b = c("ACME","ADE","Tot","Prop"))
colnames(out.sum.any) = apply(col.st,1,function(X){paste(X,collapse = "_",sep="")})



out.sum.any$N_sub_main = length(unique(cbcl_dat_all$src_subject_id.x))
out.sum.any$N_obs_main = length((cbcl_dat_all$src_subject_id.x))
out.sum.any$N_1_main = table(table(cbcl_dat_all$src_subject_id.x))[1] %>% unname()
out.sum.any$N_2_main = table(table(cbcl_dat_all$src_subject_id.x))[2] %>% unname()
out.sum.any$N_pre_main =(cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","mj_group")) %>% unique() %>% summarise(n = table(mj_group)))[2,1] %>% unname()
out.sum.any$N_pre_post_main = (cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","mj_group")) %>% unique() %>% summarise(n = table(mj_group)))[3,1] %>% unname()

#############

cbcl_dat_base = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id","site_id_l","rel_family_id","eventname",
      paste("site",c(1:21),sep="_"),
      "mj_pre","mj_pre_post","mj_pre_any","mj_group",
      "rel_relationship_1","rel_relationship_2","rel_relationship_3",
      "sex2",
      "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h",
      
      "dgs_before", "dgs_after",
      
      "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
      "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
      
      "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p", "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p",
      "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco"#,
      #"income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      #"Prisma_fit","Prisma","DISCOVERY","Achieva",
      #"y.var","global.var","include.var","movement.var","x.var"
      #,"mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth"
      # excluding in some analyses -  mat_age_knew_preg, total_birth_weight, devhx_10, devhx_6_p, income      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1", mj_pre_post == 0) %>% 
  na.omit() %>% unique() 
cbcl_dat_base = cbcl_dat_base[,-which(colnames(cbcl_dat_base) == "eventname")]


cbcl_dat_mediator = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      "Prisma_fit","Prisma","DISCOVERY","Achieva",
      "global.var","include.var","movement.var","m.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1" | eventname == "2_year_follow_up_y_arm_1", include.var == 1) %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      
      "y.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "1_year_follow_up_y_arm_1" | eventname == "3_year_follow_up_y_arm_1") %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome$eventname[cbcl_dat_outcome$eventname == "1_year_follow_up_y_arm_1"] = "baseline_year_1_arm_1"
cbcl_dat_outcome$eventname[cbcl_dat_outcome$eventname == "3_year_follow_up_y_arm_1"] = "2_year_follow_up_y_arm_1"

cbcl_dat_outcome$combo_id = paste(cbcl_dat_outcome$src_subject_id,cbcl_dat_outcome$eventname,sep = "_")
cbcl_dat_mediator$combo_id = paste(cbcl_dat_mediator$src_subject_id,cbcl_dat_mediator$eventname,sep = "_")


cbcl_dat_all = merge(cbcl_dat_base,cbcl_dat_mediator)
cbcl_dat_all = merge(cbcl_dat_all,cbcl_dat_outcome,by = "combo_id")

sum(cbcl_dat_all$src_subject_id.x == cbcl_dat_all$src_subject_id.y)
sum(cbcl_dat_all$eventname.x == cbcl_dat_all$eventname.y)

#cbcl_dat_all = cbcl_dat_all %>% filter(eventname.x == "baseline_year_1_arm_1")

length(unique(cbcl_dat_all$src_subject_id.x))

temp.dat1 = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) %>% unique()
temp.dat = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) 

famids = table(temp.dat1$rel_family_id) %>% as.data.frame() 

# remove relatives, keeping any PCE participants

for(f in 1 :dim(famids)[1]){
  
  matching.fam = which(temp.dat$rel_family_id == famids$Var1[f])
  
  t2 = temp.dat[matching.fam,]
  
  t3 = table(t2$src_subject_id.x) %>% as.data.frame()
  colnames(t3) = c("src_subject_id.x","Freq")
  t2 = merge(t2,t3) %>% unique()
  
  
  if(sum(t2$mj_pre_any)>0){ t2 = t2[t2$mj_pre_any==1,]  }
  
  if(max(t2$Freq == 2)){ t2 = t2[t2$Freq==2,] }
  
  keep = t2[sample(x = c(1:dim(t2)[1]),size = 1),]
  
  
  if(f == 1){keep.ids = keep$src_subject_id.x}else{keep.ids = c(keep.ids,keep$src_subject_id.x)}
  
}

cbcl_dat_all = cbcl_dat_all [!is.na(match(cbcl_dat_all$src_subject_id.x,keep.ids)),]
length(unique(cbcl_dat_all$src_subject_id.x))


var.names=c( "mj_pre","mj_pre_post","mj_pre_any","rel_relationship_1","rel_relationship_2","rel_relationship_3",
             "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
             "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
             "dgs_after",#"income","interview_age",
             "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
             "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
             
             "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
             "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
             "Prisma_fit","Prisma","DISCOVERY","Achieva","m.var",
             "y.var", "interview_age.x","interview_age.y",
             paste("site",c(2:21),sep="_"),
             "income_1.x","income_2.x","income_3.x","income_4.x","ed.x","pds.x",
             "income_1.y","income_2.y","income_3.y","income_4.y","ed.y","pds.y")

var.names2=c( #"mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
  "ed.x","pds.x","ed.y","pds.y",
  #"famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
  #"demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
  #"dgs_after",
  "interview_age.x","interview_age.y",
  #"famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
  #"famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
  #"income_1","income_2","income_3","income_4",
  #"famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
  #"devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
  #"Prisma_fit","Prisma","DISCOVERY","Achieva",
  "y.var","m.var")


cbcl_dat_all = cbcl_dat_all %>% 
  mutate_at(.vars = var.names2,.funs = winsorize) %>% 
  mutate_at(.vars = var.names2,.funs = function(X){
    if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
  }) %>% 
  mutate_at(.vars = var.names,.funs = scale,center=T,scale=T)

if(analysis.framework$global.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "global.var",.funs = winsorize) %>% 
    mutate_at(.vars = "global.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "global.var",.funs = scale,center=T,scale=T)
  
}


if(analysis.framework$movement.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "movement.var",.funs = winsorize) %>% 
    mutate_at(.vars = "movement.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "movement.var",.funs = scale,center=T,scale=T)
  
}

cbcl_dat_all$Age1.x = poly(cbcl_dat_all$interview_age.x,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.x = poly(cbcl_dat_all$interview_age.x,2)[,2] %>% scale(center = T,scale = T)
cbcl_dat_all$Age1.y = poly(cbcl_dat_all$interview_age.y,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.y = poly(cbcl_dat_all$interview_age.y,2)[,2] %>% scale(center = T,scale = T)


run_med = function(){lmer(m.var ~ mj_pre_any + 
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            ed.y+income_1.y+income_2.y+income_3.y+income_4.y+
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            site_2+site_3+site_4+site_5+site_6+site_7+site_8+site_9+site_10+site_11+
                            site_12+site_13+site_14+site_15+site_16+site_17+site_18+site_19+site_20+site_21+
                            #             (1|site_id_l) #+ 
                            #(1|rel_family_id)# 
                            (1|src_subject_id.x)
                          ,data =cbcl_dat_all  )}


med.fit = tryCatch(expr = run_med(),error = run_med()) # if error (b/c of memory), try again

run_out = function(){lmer(y.var ~ mj_pre_any + m.var+
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            ed.y+income_1.y+income_2.y+income_3.y+income_4.y+
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            site_2+site_3+site_4+site_5+site_6+site_7+site_8+site_9+site_10+site_11+
                            site_12+site_13+site_14+site_15+site_16+site_17+site_18+site_19+site_20+site_21+
                            #  (1|site_id_l)  
                            #(1|rel_family_id) #+
                            (1|src_subject_id.x)
                          ,data =cbcl_dat_all  )}


out.fit = tryCatch(expr = run_out(),error = run_out()) # if error (b/c of memory), try again


med.out2a <- mediate(model.m = med.fit, 
                     model.y = out.fit, 
                     treat = "mj_pre_any", 
                     mediator = "m.var",
                     sims = 10000)

out.sum.a = c(med.out2a$d.avg, med.out2a$d.avg.ci[1],med.out2a$d.avg.ci[2],med.out2a$d.avg.p,
              med.out2a$z.avg, med.out2a$z.avg.ci[1],med.out2a$z.avg.ci[2],med.out2a$z.avg.p,
              med.out2a$tau.coef,med.out2a$tau.ci[1],med.out2a$tau.ci[2],med.out2a$tau.p,
              med.out2a$n.avg, med.out2a$n.avg.ci[1],med.out2a$n.avg.ci[2],med.out2a$n.avg.p) %>% t() %>% as.data.frame()
col.st = expand.grid(a = "main","pre",c =c("Est","L","U","p") ,b = c("ACME","ADE","Tot","Prop"))
colnames(out.sum.a) = apply(col.st,1,function(X){paste(X,collapse = "_",sep="")})


####

cbcl_dat_base = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id","site_id_l","rel_family_id","eventname",
      paste("site",c(1:21),sep="_"),
      "mj_pre","mj_pre_post","mj_pre_any","mj_group",
      "rel_relationship_1","rel_relationship_2","rel_relationship_3",
      "sex2",
      "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h",
      
      "dgs_before", "dgs_after",
      
      "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
      "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
      
      "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p", "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p",
      "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco"#,
      #"income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      #"Prisma_fit","Prisma","DISCOVERY","Achieva",
      #"y.var","global.var","include.var","movement.var","x.var"
      #,"mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth"
      # excluding in some analyses -  mat_age_knew_preg, total_birth_weight, devhx_10, devhx_6_p, income      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1", mj_pre == 0) %>% 
  na.omit() %>% unique() 
cbcl_dat_base = cbcl_dat_base[,-which(colnames(cbcl_dat_base) == "eventname")]


cbcl_dat_mediator = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      "Prisma_fit","Prisma","DISCOVERY","Achieva",
      "global.var","include.var","movement.var","m.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1" | eventname == "2_year_follow_up_y_arm_1", include.var == 1) %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      
      "y.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "1_year_follow_up_y_arm_1" | eventname == "3_year_follow_up_y_arm_1") %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome$eventname[cbcl_dat_outcome$eventname == "1_year_follow_up_y_arm_1"] = "baseline_year_1_arm_1"
cbcl_dat_outcome$eventname[cbcl_dat_outcome$eventname == "3_year_follow_up_y_arm_1"] = "2_year_follow_up_y_arm_1"

cbcl_dat_outcome$combo_id = paste(cbcl_dat_outcome$src_subject_id,cbcl_dat_outcome$eventname,sep = "_")
cbcl_dat_mediator$combo_id = paste(cbcl_dat_mediator$src_subject_id,cbcl_dat_mediator$eventname,sep = "_")


cbcl_dat_all = merge(cbcl_dat_base,cbcl_dat_mediator)
cbcl_dat_all = merge(cbcl_dat_all,cbcl_dat_outcome,by = "combo_id")

sum(cbcl_dat_all$src_subject_id.x == cbcl_dat_all$src_subject_id.y)
sum(cbcl_dat_all$eventname.x == cbcl_dat_all$eventname.y)

#cbcl_dat_all = cbcl_dat_all %>% filter(eventname.x == "baseline_year_1_arm_1")

length(unique(cbcl_dat_all$src_subject_id.x))

temp.dat1 = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) %>% unique()
temp.dat = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) 

famids = table(temp.dat1$rel_family_id) %>% as.data.frame() 

# remove relatives, keeping any PCE participants

for(f in 1 :dim(famids)[1]){
  
  matching.fam = which(temp.dat$rel_family_id == famids$Var1[f])
  
  t2 = temp.dat[matching.fam,]
  
  t3 = table(t2$src_subject_id.x) %>% as.data.frame()
  colnames(t3) = c("src_subject_id.x","Freq")
  t2 = merge(t2,t3) %>% unique()
  
  
  if(sum(t2$mj_pre_any)>0){ t2 = t2[t2$mj_pre_any==1,]  }
  
  if(max(t2$Freq == 2)){ t2 = t2[t2$Freq==2,] }
  
  keep = t2[sample(x = c(1:dim(t2)[1]),size = 1),]
  
  
  if(f == 1){keep.ids = keep$src_subject_id.x}else{keep.ids = c(keep.ids,keep$src_subject_id.x)}
  
}

cbcl_dat_all = cbcl_dat_all [!is.na(match(cbcl_dat_all$src_subject_id.x,keep.ids)),]

length(unique(cbcl_dat_all$src_subject_id.x))


var.names=c( "mj_pre","mj_pre_post","mj_pre_any","rel_relationship_1","rel_relationship_2","rel_relationship_3",
             "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
             "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
             "dgs_after",#"income","interview_age",
             "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
             "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
             
             "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
             "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
             "Prisma_fit","Prisma","DISCOVERY","Achieva","m.var",
             "y.var", "interview_age.x","interview_age.y",
             paste("site",c(2:21),sep="_"),
             "income_1.x","income_2.x","income_3.x","income_4.x","ed.x","pds.x",
             "income_1.y","income_2.y","income_3.y","income_4.y","ed.y","pds.y")

var.names2=c( #"mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
  "ed.x","pds.x","ed.y","pds.y",
  #"famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
  #"demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
  #"dgs_after",
  "interview_age.x","interview_age.y",
  #"famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
  #"famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
  #"income_1","income_2","income_3","income_4",
  #"famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
  #"devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
  #"Prisma_fit","Prisma","DISCOVERY","Achieva",
  "y.var","m.var")


cbcl_dat_all = cbcl_dat_all %>% 
  mutate_at(.vars = var.names2,.funs = winsorize) %>% 
  mutate_at(.vars = var.names2,.funs = function(X){
    if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
  }) %>% 
  mutate_at(.vars = var.names,.funs = scale,center=T,scale=T)

if(analysis.framework$global.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "global.var",.funs = winsorize) %>% 
    mutate_at(.vars = "global.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "global.var",.funs = scale,center=T,scale=T)
  
}


if(analysis.framework$movement.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "movement.var",.funs = winsorize) %>% 
    mutate_at(.vars = "movement.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "movement.var",.funs = scale,center=T,scale=T)
  
}

cbcl_dat_all$Age1.x = poly(cbcl_dat_all$interview_age.x,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.x = poly(cbcl_dat_all$interview_age.x,2)[,2] %>% scale(center = T,scale = T)
cbcl_dat_all$Age1.y = poly(cbcl_dat_all$interview_age.y,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.y = poly(cbcl_dat_all$interview_age.y,2)[,2] %>% scale(center = T,scale = T)


run_med = function(){lmer(m.var ~ mj_pre_any + 
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            ed.y+income_1.y+income_2.y+income_3.y+income_4.y+
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            site_2+site_3+site_4+site_5+site_6+site_7+site_8+site_9+site_10+site_11+
                            site_12+site_13+site_14+site_15+site_16+site_17+site_18+site_19+site_20+site_21+
                            #             (1|site_id_l) #+ 
                            #(1|rel_family_id)# 
                            (1|src_subject_id.x)
                          ,data =cbcl_dat_all  )}


med.fit = tryCatch(expr = run_med(),error = run_med()) # if error (b/c of memory), try again

run_out = function(){lmer(y.var ~ mj_pre_any + m.var+
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            ed.y+income_1.y+income_2.y+income_3.y+income_4.y+
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            site_2+site_3+site_4+site_5+site_6+site_7+site_8+site_9+site_10+site_11+
                            site_12+site_13+site_14+site_15+site_16+site_17+site_18+site_19+site_20+site_21+
                            #  (1|site_id_l)  
                            #(1|rel_family_id) #+
                            (1|src_subject_id.x)
                          ,data =cbcl_dat_all  )}


out.fit = tryCatch(expr = run_out(),error = run_out()) # if error (b/c of memory), try again


med.out2b <- mediate(model.m = med.fit, 
                     model.y = out.fit, 
                     treat = "mj_pre_any", 
                     mediator = "m.var",
                     sims = 10000)

out.sum.b = c(med.out2b$d.avg, med.out2b$d.avg.ci[1],med.out2b$d.avg.ci[2],med.out2b$d.avg.p,
              med.out2b$z.avg, med.out2b$z.avg.ci[1],med.out2b$z.avg.ci[2],med.out2b$z.avg.p,
              med.out2b$tau.coef,med.out2b$tau.ci[1],med.out2b$tau.ci[2],med.out2b$tau.p,
              med.out2b$n.avg, med.out2b$n.avg.ci[1],med.out2b$n.avg.ci[2],med.out2b$n.avg.p) %>% t() %>% as.data.frame()
col.st = expand.grid(a = "main","prepost",c =c("Est","L","U","p") ,b = c("ACME","ADE","Tot","Prop"))
colnames(out.sum.b) = apply(col.st,1,function(X){paste(X,collapse = "_",sep="")})


out_final_main = cbind(out.sum.any,out.sum.a,out.sum.b)




#####################################
## Extra covariates

print("extra")

cbcl_dat_base = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id","site_id_l","rel_family_id","eventname",
      paste("site",c(1:21),sep="_"),
      "mj_pre","mj_pre_post","mj_pre_any","mj_group",
      "rel_relationship_1","rel_relationship_2","rel_relationship_3",
      "sex2",
      "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h",
      
      "dgs_before", "dgs_after",
      
      "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
      "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
      
      "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p", "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p",
      "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco"#,
      #"income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      ,"mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth"
      
      #"Prisma_fit","Prisma","DISCOVERY","Achieva",
      #"y.var","global.var","include.var","movement.var","x.var"
      #,"mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth"
      # excluding in some analyses -  mat_age_knew_preg, total_birth_weight, devhx_10, devhx_6_p, income      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1") %>% 
  na.omit() %>% unique() 
cbcl_dat_base = cbcl_dat_base[,-which(colnames(cbcl_dat_base) == "eventname")]


cbcl_dat_mediator = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      "Prisma_fit","Prisma","DISCOVERY","Achieva",
      "global.var","include.var","movement.var","m.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1" | eventname == "2_year_follow_up_y_arm_1", include.var == 1) %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      
      "y.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "1_year_follow_up_y_arm_1" | eventname == "3_year_follow_up_y_arm_1") %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome$eventname[cbcl_dat_outcome$eventname == "1_year_follow_up_y_arm_1"] = "baseline_year_1_arm_1"
cbcl_dat_outcome$eventname[cbcl_dat_outcome$eventname == "3_year_follow_up_y_arm_1"] = "2_year_follow_up_y_arm_1"

cbcl_dat_outcome$combo_id = paste(cbcl_dat_outcome$src_subject_id,cbcl_dat_outcome$eventname,sep = "_")
cbcl_dat_mediator$combo_id = paste(cbcl_dat_mediator$src_subject_id,cbcl_dat_mediator$eventname,sep = "_")


cbcl_dat_all = merge(cbcl_dat_base,cbcl_dat_mediator)
cbcl_dat_all = merge(cbcl_dat_all,cbcl_dat_outcome,by = "combo_id")

sum(cbcl_dat_all$src_subject_id.x == cbcl_dat_all$src_subject_id.y)
sum(cbcl_dat_all$eventname.x == cbcl_dat_all$eventname.y)

#cbcl_dat_all = cbcl_dat_all %>% filter(eventname.x == "baseline_year_1_arm_1")

length(unique(cbcl_dat_all$src_subject_id.x))

temp.dat1 = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) %>% unique()
temp.dat = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) 

famids = table(temp.dat1$rel_family_id) %>% as.data.frame() 

# remove relatives, keeping any PCE participants

for(f in 1 :dim(famids)[1]){
  
  matching.fam = which(temp.dat$rel_family_id == famids$Var1[f])
  
  t2 = temp.dat[matching.fam,]
  
  t3 = table(t2$src_subject_id.x) %>% as.data.frame()
  colnames(t3) = c("src_subject_id.x","Freq")
  t2 = merge(t2,t3) %>% unique()
  
  
  if(sum(t2$mj_pre_any)>0){ t2 = t2[t2$mj_pre_any==1,]  }
  
  if(max(t2$Freq == 2)){ t2 = t2[t2$Freq==2,] }
  
  keep = t2[sample(x = c(1:dim(t2)[1]),size = 1),]
  
  
  if(f == 1){keep.ids = keep$src_subject_id.x}else{keep.ids = c(keep.ids,keep$src_subject_id.x)}
  
}

cbcl_dat_all = cbcl_dat_all [!is.na(match(cbcl_dat_all$src_subject_id.x,keep.ids)),]

length(unique(cbcl_dat_all$src_subject_id.x))
table(cbcl_dat_all$mj_group)

var.names=c( "mj_pre","mj_pre_post","mj_pre_any","rel_relationship_1","rel_relationship_2","rel_relationship_3",
             "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
             "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
             "dgs_after",#"income","interview_age",
             "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
             "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
             "mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth",
             
             "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
             "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
             "Prisma_fit","Prisma","DISCOVERY","Achieva","m.var",
             "y.var", "interview_age.x","interview_age.y",
             paste("site",c(2:21),sep="_"),
             "income_1.x","income_2.x","income_3.x","income_4.x","ed.x","pds.x",
             "income_1.y","income_2.y","income_3.y","income_4.y","ed.y","pds.y")

var.names2=c( #"mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
  "ed.x","pds.x","ed.y","pds.y",
  #"famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
  #"demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
  #"dgs_after",
  "interview_age.x","interview_age.y",
  #"famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
  #"famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
  #"income_1","income_2","income_3","income_4",
  #"famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
  #"devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
  #"Prisma_fit","Prisma","DISCOVERY","Achieva",
  "y.var","m.var")


cbcl_dat_all = cbcl_dat_all %>% 
  mutate_at(.vars = var.names2,.funs = winsorize) %>% 
  mutate_at(.vars = var.names2,.funs = function(X){
    if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
  }) %>% 
  mutate_at(.vars = var.names,.funs = scale,center=T,scale=T)

if(analysis.framework$global.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "global.var",.funs = winsorize) %>% 
    mutate_at(.vars = "global.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "global.var",.funs = scale,center=T,scale=T)
  
}


if(analysis.framework$movement.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "movement.var",.funs = winsorize) %>% 
    mutate_at(.vars = "movement.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "movement.var",.funs = scale,center=T,scale=T)
  
}

cbcl_dat_all$Age1.x = poly(cbcl_dat_all$interview_age.x,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.x = poly(cbcl_dat_all$interview_age.x,2)[,2] %>% scale(center = T,scale = T)
cbcl_dat_all$Age1.y = poly(cbcl_dat_all$interview_age.y,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.y = poly(cbcl_dat_all$interview_age.y,2)[,2] %>% scale(center = T,scale = T)


run_med = function(){lmer(m.var ~ mj_pre_any + 
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            ed.y+income_1.y+income_2.y+income_3.y+income_4.y+
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            site_2+site_3+site_4+site_5+site_6+site_7+site_8+site_9+site_10+site_11+
                            site_12+site_13+site_14+site_15+site_16+site_17+site_18+site_19+site_20+site_21+
                            #             (1|site_id_l) #+ 
                            #(1|rel_family_id)# 
                            (1|src_subject_id.x)
                          ,data =cbcl_dat_all  )}


med.fit = tryCatch(expr = run_med(),error = run_med()) # if error (b/c of memory), try again

run_out = function(){lmer(y.var ~ mj_pre_any + m.var+
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            ed.y+income_1.y+income_2.y+income_3.y+income_4.y+
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            site_2+site_3+site_4+site_5+site_6+site_7+site_8+site_9+site_10+site_11+
                            site_12+site_13+site_14+site_15+site_16+site_17+site_18+site_19+site_20+site_21+
                            #  (1|site_id_l)  
                            #(1|rel_family_id) #+
                            (1|src_subject_id.x)
                          ,data =cbcl_dat_all  )}


out.fit = tryCatch(expr = run_out(),error = run_out()) # if error (b/c of memory), try again


med.out <- mediate(model.m = med.fit, 
                   model.y = out.fit, 
                   treat = "mj_pre_any", 
                   mediator = "m.var",
                   sims = 10000)

out.sum.any = c(med.out$d.avg, med.out$d.avg.ci[1],med.out$d.avg.ci[2],med.out$d.avg.p,
                med.out$z.avg, med.out$z.avg.ci[1],med.out$z.avg.ci[2],med.out$z.avg.p,
                med.out$tau.coef,med.out$tau.ci[1],med.out$tau.ci[2],med.out$tau.p,
                med.out$n.avg, med.out$n.avg.ci[1],med.out$n.avg.ci[2],med.out$n.avg.p) %>% t() %>% as.data.frame()
col.st = expand.grid(a = "extra","Any",c =c("Est","L","U","p") ,b = c("ACME","ADE","Tot","Prop"))
colnames(out.sum.any) = apply(col.st,1,function(X){paste(X,collapse = "_",sep="")})



out.sum.any$N_sub_extra = length(unique(cbcl_dat_all$src_subject_id.x))
out.sum.any$N_obs_extra = length((cbcl_dat_all$src_subject_id.x))
out.sum.any$N_1_extra = table(table(cbcl_dat_all$src_subject_id.x))[1] %>% unname()
out.sum.any$N_2_extra = table(table(cbcl_dat_all$src_subject_id.x))[2] %>% unname()
out.sum.any$N_pre_extra =(cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","mj_group")) %>% unique() %>% summarise(n = table(mj_group)))[2,1] %>% unname()
out.sum.any$N_pre_post_extra = (cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","mj_group")) %>% unique() %>% summarise(n = table(mj_group)))[3,1] %>% unname()

#############

cbcl_dat_base = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id","site_id_l","rel_family_id","eventname",
      paste("site",c(1:21),sep="_"),
      "mj_pre","mj_pre_post","mj_pre_any","mj_group",
      "rel_relationship_1","rel_relationship_2","rel_relationship_3",
      "sex2",
      "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h",
      
      "dgs_before", "dgs_after",
      
      "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
      "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
      
      "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p", "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p",
      "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco"#,
      #"income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      ,"mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth"
      
      #"Prisma_fit","Prisma","DISCOVERY","Achieva",
      #"y.var","global.var","include.var","movement.var","x.var"
      #,"mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth"
      # excluding in some analyses -  mat_age_knew_preg, total_birth_weight, devhx_10, devhx_6_p, income      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1", mj_pre_post == 0) %>% 
  na.omit() %>% unique() 
cbcl_dat_base = cbcl_dat_base[,-which(colnames(cbcl_dat_base) == "eventname")]


cbcl_dat_mediator = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      "Prisma_fit","Prisma","DISCOVERY","Achieva",
      "global.var","include.var","movement.var","m.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1" | eventname == "2_year_follow_up_y_arm_1", include.var == 1) %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      
      "y.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "1_year_follow_up_y_arm_1" | eventname == "3_year_follow_up_y_arm_1") %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome$eventname[cbcl_dat_outcome$eventname == "1_year_follow_up_y_arm_1"] = "baseline_year_1_arm_1"
cbcl_dat_outcome$eventname[cbcl_dat_outcome$eventname == "3_year_follow_up_y_arm_1"] = "2_year_follow_up_y_arm_1"

cbcl_dat_outcome$combo_id = paste(cbcl_dat_outcome$src_subject_id,cbcl_dat_outcome$eventname,sep = "_")
cbcl_dat_mediator$combo_id = paste(cbcl_dat_mediator$src_subject_id,cbcl_dat_mediator$eventname,sep = "_")


cbcl_dat_all = merge(cbcl_dat_base,cbcl_dat_mediator)
cbcl_dat_all = merge(cbcl_dat_all,cbcl_dat_outcome,by = "combo_id")

sum(cbcl_dat_all$src_subject_id.x == cbcl_dat_all$src_subject_id.y)
sum(cbcl_dat_all$eventname.x == cbcl_dat_all$eventname.y)

#cbcl_dat_all = cbcl_dat_all %>% filter(eventname.x == "baseline_year_1_arm_1")

length(unique(cbcl_dat_all$src_subject_id.x))

temp.dat1 = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) %>% unique()
temp.dat = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) 

famids = table(temp.dat1$rel_family_id) %>% as.data.frame() 

# remove relatives, keeping any PCE participants

for(f in 1 :dim(famids)[1]){
  
  matching.fam = which(temp.dat$rel_family_id == famids$Var1[f])
  
  t2 = temp.dat[matching.fam,]
  
  t3 = table(t2$src_subject_id.x) %>% as.data.frame()
  colnames(t3) = c("src_subject_id.x","Freq")
  t2 = merge(t2,t3) %>% unique()
  
  
  if(sum(t2$mj_pre_any)>0){ t2 = t2[t2$mj_pre_any==1,]  }
  
  if(max(t2$Freq == 2)){ t2 = t2[t2$Freq==2,] }
  
  keep = t2[sample(x = c(1:dim(t2)[1]),size = 1),]
  
  
  if(f == 1){keep.ids = keep$src_subject_id.x}else{keep.ids = c(keep.ids,keep$src_subject_id.x)}
  
}

cbcl_dat_all = cbcl_dat_all [!is.na(match(cbcl_dat_all$src_subject_id.x,keep.ids)),]

length(unique(cbcl_dat_all$src_subject_id.x))


var.names=c( "mj_pre","mj_pre_post","mj_pre_any","rel_relationship_1","rel_relationship_2","rel_relationship_3",
             "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
             "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
             "dgs_after",#"income","interview_age",
             "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
             "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
             "mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth",
             
             "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
             "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
             "Prisma_fit","Prisma","DISCOVERY","Achieva","m.var",
             "y.var", "interview_age.x","interview_age.y",
             paste("site",c(2:21),sep="_"),
             "income_1.x","income_2.x","income_3.x","income_4.x","ed.x","pds.x",
             "income_1.y","income_2.y","income_3.y","income_4.y","ed.y","pds.y")

var.names2=c( #"mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
  "ed.x","pds.x","ed.y","pds.y",
  #"famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
  #"demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
  #"dgs_after",
  "interview_age.x","interview_age.y",
  #"famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
  #"famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
  #"income_1","income_2","income_3","income_4",
  #"famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
  #"devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
  #"Prisma_fit","Prisma","DISCOVERY","Achieva",
  "y.var","m.var")


cbcl_dat_all = cbcl_dat_all %>% 
  mutate_at(.vars = var.names2,.funs = winsorize) %>% 
  mutate_at(.vars = var.names2,.funs = function(X){
    if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
  }) %>% 
  mutate_at(.vars = var.names,.funs = scale,center=T,scale=T)

if(analysis.framework$global.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "global.var",.funs = winsorize) %>% 
    mutate_at(.vars = "global.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "global.var",.funs = scale,center=T,scale=T)
  
}


if(analysis.framework$movement.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "movement.var",.funs = winsorize) %>% 
    mutate_at(.vars = "movement.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "movement.var",.funs = scale,center=T,scale=T)
  
}

cbcl_dat_all$Age1.x = poly(cbcl_dat_all$interview_age.x,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.x = poly(cbcl_dat_all$interview_age.x,2)[,2] %>% scale(center = T,scale = T)
cbcl_dat_all$Age1.y = poly(cbcl_dat_all$interview_age.y,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.y = poly(cbcl_dat_all$interview_age.y,2)[,2] %>% scale(center = T,scale = T)


run_med = function(){lmer(m.var ~ mj_pre_any + 
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            ed.y+income_1.y+income_2.y+income_3.y+income_4.y+
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            site_2+site_3+site_4+site_5+site_6+site_7+site_8+site_9+site_10+site_11+
                            site_12+site_13+site_14+site_15+site_16+site_17+site_18+site_19+site_20+site_21+
                            #             (1|site_id_l) #+ 
                            #(1|rel_family_id)# 
                            (1|src_subject_id.x)
                          ,data =cbcl_dat_all  )}


med.fit = tryCatch(expr = run_med(),error = run_med()) # if error (b/c of memory), try again

run_out = function(){lmer(y.var ~ mj_pre_any + m.var+
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            ed.y+income_1.y+income_2.y+income_3.y+income_4.y+
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            site_2+site_3+site_4+site_5+site_6+site_7+site_8+site_9+site_10+site_11+
                            site_12+site_13+site_14+site_15+site_16+site_17+site_18+site_19+site_20+site_21+
                            #  (1|site_id_l)  
                            #(1|rel_family_id) #+
                            (1|src_subject_id.x)
                          ,data =cbcl_dat_all  )}


out.fit = tryCatch(expr = run_out(),error = run_out()) # if error (b/c of memory), try again


med.out2a <- mediate(model.m = med.fit, 
                     model.y = out.fit, 
                     treat = "mj_pre_any", 
                     mediator = "m.var",
                     sims = 10000)

out.sum.a = c(med.out2a$d.avg, med.out2a$d.avg.ci[1],med.out2a$d.avg.ci[2],med.out2a$d.avg.p,
              med.out2a$z.avg, med.out2a$z.avg.ci[1],med.out2a$z.avg.ci[2],med.out2a$z.avg.p,
              med.out2a$tau.coef,med.out2a$tau.ci[1],med.out2a$tau.ci[2],med.out2a$tau.p,
              med.out2a$n.avg, med.out2a$n.avg.ci[1],med.out2a$n.avg.ci[2],med.out2a$n.avg.p) %>% t() %>% as.data.frame()
col.st = expand.grid(a = "extra","pre",c =c("Est","L","U","p") ,b = c("ACME","ADE","Tot","Prop"))
colnames(out.sum.a) = apply(col.st,1,function(X){paste(X,collapse = "_",sep="")})


####

cbcl_dat_base = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id","site_id_l","rel_family_id","eventname",
      paste("site",c(1:21),sep="_"),
      "mj_pre","mj_pre_post","mj_pre_any","mj_group",
      "rel_relationship_1","rel_relationship_2","rel_relationship_3",
      "sex2",
      "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h",
      
      "dgs_before", "dgs_after",
      
      "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
      "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
      
      "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p", "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p",
      "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco"#,
      #"income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      ,"mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth"
      
      #"Prisma_fit","Prisma","DISCOVERY","Achieva",
      #"y.var","global.var","include.var","movement.var","x.var"
      #,"mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth"
      # excluding in some analyses -  mat_age_knew_preg, total_birth_weight, devhx_10, devhx_6_p, income      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1", mj_pre == 0) %>% 
  na.omit() %>% unique() 
cbcl_dat_base = cbcl_dat_base[,-which(colnames(cbcl_dat_base) == "eventname")]


cbcl_dat_mediator = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      "Prisma_fit","Prisma","DISCOVERY","Achieva",
      "global.var","include.var","movement.var","m.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1" | eventname == "2_year_follow_up_y_arm_1", include.var == 1) %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      
      "y.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "1_year_follow_up_y_arm_1" | eventname == "3_year_follow_up_y_arm_1") %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome$eventname[cbcl_dat_outcome$eventname == "1_year_follow_up_y_arm_1"] = "baseline_year_1_arm_1"
cbcl_dat_outcome$eventname[cbcl_dat_outcome$eventname == "3_year_follow_up_y_arm_1"] = "2_year_follow_up_y_arm_1"

cbcl_dat_outcome$combo_id = paste(cbcl_dat_outcome$src_subject_id,cbcl_dat_outcome$eventname,sep = "_")
cbcl_dat_mediator$combo_id = paste(cbcl_dat_mediator$src_subject_id,cbcl_dat_mediator$eventname,sep = "_")


cbcl_dat_all = merge(cbcl_dat_base,cbcl_dat_mediator)
cbcl_dat_all = merge(cbcl_dat_all,cbcl_dat_outcome,by = "combo_id")

sum(cbcl_dat_all$src_subject_id.x == cbcl_dat_all$src_subject_id.y)
sum(cbcl_dat_all$eventname.x == cbcl_dat_all$eventname.y)

#cbcl_dat_all = cbcl_dat_all %>% filter(eventname.x == "baseline_year_1_arm_1")

length(unique(cbcl_dat_all$src_subject_id.x))

temp.dat1 = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) %>% unique()
temp.dat = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) 

famids = table(temp.dat1$rel_family_id) %>% as.data.frame() 

# remove relatives, keeping any PCE participants

for(f in 1 :dim(famids)[1]){
  
  matching.fam = which(temp.dat$rel_family_id == famids$Var1[f])
  
  t2 = temp.dat[matching.fam,]
  
  t3 = table(t2$src_subject_id.x) %>% as.data.frame()
  colnames(t3) = c("src_subject_id.x","Freq")
  t2 = merge(t2,t3) %>% unique()
  
  
  if(sum(t2$mj_pre_any)>0){ t2 = t2[t2$mj_pre_any==1,]  }
  
  if(max(t2$Freq == 2)){ t2 = t2[t2$Freq==2,] }
  
  keep = t2[sample(x = c(1:dim(t2)[1]),size = 1),]
  
  
  if(f == 1){keep.ids = keep$src_subject_id.x}else{keep.ids = c(keep.ids,keep$src_subject_id.x)}
  
}

cbcl_dat_all = cbcl_dat_all [!is.na(match(cbcl_dat_all$src_subject_id.x,keep.ids)),]

length(unique(cbcl_dat_all$src_subject_id.x))


var.names=c( "mj_pre","mj_pre_post","mj_pre_any","rel_relationship_1","rel_relationship_2","rel_relationship_3",
             "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
             "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
             "dgs_after",#"income","interview_age",
             "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
             "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
             "mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth",
             
             "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
             "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
             "Prisma_fit","Prisma","DISCOVERY","Achieva","m.var",
             "y.var", "interview_age.x","interview_age.y",
             paste("site",c(2:21),sep="_"),
             "income_1.x","income_2.x","income_3.x","income_4.x","ed.x","pds.x",
             "income_1.y","income_2.y","income_3.y","income_4.y","ed.y","pds.y")

var.names2=c( #"mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
  "ed.x","pds.x","ed.y","pds.y",
  #"famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
  #"demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
  #"dgs_after",
  "interview_age.x","interview_age.y",
  #"famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
  #"famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
  #"income_1","income_2","income_3","income_4",
  #"famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
  #"devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
  #"Prisma_fit","Prisma","DISCOVERY","Achieva",
  "y.var","m.var")


cbcl_dat_all = cbcl_dat_all %>% 
  mutate_at(.vars = var.names2,.funs = winsorize) %>% 
  mutate_at(.vars = var.names2,.funs = function(X){
    if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
  }) %>% 
  mutate_at(.vars = var.names,.funs = scale,center=T,scale=T)

if(analysis.framework$global.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "global.var",.funs = winsorize) %>% 
    mutate_at(.vars = "global.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "global.var",.funs = scale,center=T,scale=T)
  
}


if(analysis.framework$movement.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "movement.var",.funs = winsorize) %>% 
    mutate_at(.vars = "movement.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "movement.var",.funs = scale,center=T,scale=T)
  
}

cbcl_dat_all$Age1.x = poly(cbcl_dat_all$interview_age.x,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.x = poly(cbcl_dat_all$interview_age.x,2)[,2] %>% scale(center = T,scale = T)
cbcl_dat_all$Age1.y = poly(cbcl_dat_all$interview_age.y,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.y = poly(cbcl_dat_all$interview_age.y,2)[,2] %>% scale(center = T,scale = T)


run_med = function(){lmer(m.var ~ mj_pre_any + 
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            ed.y+income_1.y+income_2.y+income_3.y+income_4.y+
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            site_2+site_3+site_4+site_5+site_6+site_7+site_8+site_9+site_10+site_11+
                            site_12+site_13+site_14+site_15+site_16+site_17+site_18+site_19+site_20+site_21+
                            #             (1|site_id_l) #+ 
                            #(1|rel_family_id)# 
                            (1|src_subject_id.x)
                          ,data =cbcl_dat_all  )}


med.fit = tryCatch(expr = run_med(),error = run_med()) # if error (b/c of memory), try again

run_out = function(){lmer(y.var ~ mj_pre_any + m.var+
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            ed.y+income_1.y+income_2.y+income_3.y+income_4.y+
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            site_2+site_3+site_4+site_5+site_6+site_7+site_8+site_9+site_10+site_11+
                            site_12+site_13+site_14+site_15+site_16+site_17+site_18+site_19+site_20+site_21+
                            #  (1|site_id_l)  
                            #(1|rel_family_id) #+
                            (1|src_subject_id.x)
                          ,data =cbcl_dat_all  )}


out.fit = tryCatch(expr = run_out(),error = run_out()) # if error (b/c of memory), try again


med.out2b <- mediate(model.m = med.fit, 
                     model.y = out.fit, 
                     treat = "mj_pre_any", 
                     mediator = "m.var",
                     sims = 10000)

out.sum.b = c(med.out2b$d.avg, med.out2b$d.avg.ci[1],med.out2b$d.avg.ci[2],med.out2b$d.avg.p,
              med.out2b$z.avg, med.out2b$z.avg.ci[1],med.out2b$z.avg.ci[2],med.out2b$z.avg.p,
              med.out2b$tau.coef,med.out2b$tau.ci[1],med.out2b$tau.ci[2],med.out2b$tau.p,
              med.out2b$n.avg, med.out2b$n.avg.ci[1],med.out2b$n.avg.ci[2],med.out2b$n.avg.p) %>% t() %>% as.data.frame()
col.st = expand.grid(a = "extra","prepost",c =c("Est","L","U","p") ,b = c("ACME","ADE","Tot","Prop"))
colnames(out.sum.b) = apply(col.st,1,function(X){paste(X,collapse = "_",sep="")})


out_final_extra = cbind(out.sum.any,out.sum.a,out.sum.b)






#################################################################
## Genetics

print("genetics")

cbcl_dat_base = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id","site_id_l","rel_family_id","eventname",
      paste("site",c(1:21),sep="_"),
      "mj_pre","mj_pre_post","mj_pre_any","mj_group",
      "rel_relationship_1","rel_relationship_2","rel_relationship_3",
      "sex2",
      "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h",
      
      "dgs_before", "dgs_after",
      
      "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
      "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
      
      "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p", "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p",
      "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco"#,
      #"income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      ,  colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PRSCS")],
      colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PC")]
      #"Prisma_fit","Prisma","DISCOVERY","Achieva",
      #"y.var","global.var","include.var","movement.var","x.var"
      #,"mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth"
      # excluding in some analyses -  mat_age_knew_preg, total_birth_weight, devhx_10, devhx_6_p, income      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1") %>% 
  na.omit() %>% unique() 
cbcl_dat_base = cbcl_dat_base[,-which(colnames(cbcl_dat_base) == "eventname")]


cbcl_dat_mediator = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      "Prisma_fit","Prisma","DISCOVERY","Achieva",
      "global.var","include.var","movement.var","m.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1" | eventname == "2_year_follow_up_y_arm_1", include.var == 1) %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      
      "y.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "1_year_follow_up_y_arm_1" | eventname == "3_year_follow_up_y_arm_1") %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome$eventname[cbcl_dat_outcome$eventname == "1_year_follow_up_y_arm_1"] = "baseline_year_1_arm_1"
cbcl_dat_outcome$eventname[cbcl_dat_outcome$eventname == "3_year_follow_up_y_arm_1"] = "2_year_follow_up_y_arm_1"

cbcl_dat_outcome$combo_id = paste(cbcl_dat_outcome$src_subject_id,cbcl_dat_outcome$eventname,sep = "_")
cbcl_dat_mediator$combo_id = paste(cbcl_dat_mediator$src_subject_id,cbcl_dat_mediator$eventname,sep = "_")


cbcl_dat_all = merge(cbcl_dat_base,cbcl_dat_mediator)
cbcl_dat_all = merge(cbcl_dat_all,cbcl_dat_outcome,by = "combo_id")

sum(cbcl_dat_all$src_subject_id.x == cbcl_dat_all$src_subject_id.y)
sum(cbcl_dat_all$eventname.x == cbcl_dat_all$eventname.y)

#cbcl_dat_all = cbcl_dat_all %>% filter(eventname.x == "baseline_year_1_arm_1")

length(unique(cbcl_dat_all$src_subject_id.x))

temp.dat1 = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) %>% unique()
temp.dat = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) 

famids = table(temp.dat1$rel_family_id) %>% as.data.frame() 

# remove relatives, keeping any PCE participants

for(f in 1 :dim(famids)[1]){
  
  matching.fam = which(temp.dat$rel_family_id == famids$Var1[f])
  
  t2 = temp.dat[matching.fam,]
  
  t3 = table(t2$src_subject_id.x) %>% as.data.frame()
  colnames(t3) = c("src_subject_id.x","Freq")
  t2 = merge(t2,t3) %>% unique()
  
  
  if(sum(t2$mj_pre_any)>0){ t2 = t2[t2$mj_pre_any==1,]  }
  
  if(max(t2$Freq == 2)){ t2 = t2[t2$Freq==2,] }
  
  keep = t2[sample(x = c(1:dim(t2)[1]),size = 1),]
  
  
  if(f == 1){keep.ids = keep$src_subject_id.x}else{keep.ids = c(keep.ids,keep$src_subject_id.x)}
  
}

cbcl_dat_all = cbcl_dat_all [!is.na(match(cbcl_dat_all$src_subject_id.x,keep.ids)),]

length(unique(cbcl_dat_all$src_subject_id.x))
table(cbcl_dat_all$mj_group)

var.names=c( "mj_pre","mj_pre_post","mj_pre_any","rel_relationship_1","rel_relationship_2","rel_relationship_3",
             "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
             "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
             "dgs_after",#"income","interview_age",
             "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
             "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
             colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PRSCS")],
             colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PC")],
             "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
             "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
             "Prisma_fit","Prisma","DISCOVERY","Achieva","m.var",
             "y.var", "interview_age.x","interview_age.y",
             paste("site",c(2:21),sep="_"),
             "income_1.x","income_2.x","income_3.x","income_4.x","ed.x","pds.x",
             "income_1.y","income_2.y","income_3.y","income_4.y","ed.y","pds.y")

var.names2=c( #"mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
  "ed.x","pds.x","ed.y","pds.y",
  #"famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
  #"demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
  #"dgs_after",
  "interview_age.x","interview_age.y",
  colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PRSCS")],
  colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PC")],
  #"famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
  #"famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
  #"income_1","income_2","income_3","income_4",
  #"famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
  #"devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
  #"Prisma_fit","Prisma","DISCOVERY","Achieva",
  "y.var","m.var")


cbcl_dat_all = cbcl_dat_all %>% 
  mutate_at(.vars = var.names2,.funs = winsorize) %>% 
  mutate_at(.vars = var.names2,.funs = function(X){
    if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
  }) %>% 
  mutate_at(.vars = var.names,.funs = scale,center=T,scale=T)

if(analysis.framework$global.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "global.var",.funs = winsorize) %>% 
    mutate_at(.vars = "global.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "global.var",.funs = scale,center=T,scale=T)
  
}


if(analysis.framework$movement.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "movement.var",.funs = winsorize) %>% 
    mutate_at(.vars = "movement.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "movement.var",.funs = scale,center=T,scale=T)
  
}

cbcl_dat_all$Age1.x = poly(cbcl_dat_all$interview_age.x,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.x = poly(cbcl_dat_all$interview_age.x,2)[,2] %>% scale(center = T,scale = T)
cbcl_dat_all$Age1.y = poly(cbcl_dat_all$interview_age.y,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.y = poly(cbcl_dat_all$interview_age.y,2)[,2] %>% scale(center = T,scale = T)


run_med = function(){lmer(m.var ~ mj_pre_any + 
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            #ed.y+
                            income_1.y+income_2.y+income_3.y+income_4.y+
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            #demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            CUD_PRSCS + 
                            PC_1+PC_2+PC_3+PC_4+PC_5+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10+
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            site_2+site_3+site_4+site_5+site_6+site_7+site_8+site_9+site_10+site_11+
                            site_12+site_13+site_14+site_15+site_16+site_17+site_18+site_19+site_20+site_21+
                            #             (1|site_id_l) #+ 
                            #(1|rel_family_id)# 
                            (1|src_subject_id.x)
                          ,data =cbcl_dat_all  )}


med.fit = tryCatch(expr = run_med(),error = run_med()) # if error (b/c of memory), try again

run_out = function(){lmer(y.var ~ mj_pre_any + m.var+
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+
                            pds.x+
                            #ed.y+
                            income_1.y+income_2.y+income_3.y+income_4.y+
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            #demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            CUD_PRSCS + 
                            PC_1+PC_2+PC_3+PC_4+PC_5+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10+
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            site_2+site_3+site_4+site_5+site_6+site_7+site_8+site_9+site_10+site_11+
                            site_12+site_13+site_14+site_15+site_16+site_17+site_18+site_19+site_20+site_21+
                            #  (1|site_id_l)  
                            #(1|rel_family_id) #+
                            (1|src_subject_id.x)
                          ,data =cbcl_dat_all  )}


out.fit = tryCatch(expr = run_out(),error = run_out()) # if error (b/c of memory), try again


med.out <- mediate(model.m = med.fit, 
                   model.y = out.fit, 
                   treat = "mj_pre_any", 
                   mediator = "m.var",
                   sims = 10000)

out.sum.any = c(med.out$d.avg, med.out$d.avg.ci[1],med.out$d.avg.ci[2],med.out$d.avg.p,
                med.out$z.avg, med.out$z.avg.ci[1],med.out$z.avg.ci[2],med.out$z.avg.p,
                med.out$tau.coef,med.out$tau.ci[1],med.out$tau.ci[2],med.out$tau.p,
                med.out$n.avg, med.out$n.avg.ci[1],med.out$n.avg.ci[2],med.out$n.avg.p) %>% t() %>% as.data.frame()
col.st = expand.grid(a = "genetics","Any",c =c("Est","L","U","p") ,b = c("ACME","ADE","Tot","Prop"))
colnames(out.sum.any) = apply(col.st,1,function(X){paste(X,collapse = "_",sep="")})



out.sum.any$N_sub_genetics = length(unique(cbcl_dat_all$src_subject_id.x))
out.sum.any$N_obs_genetics = length((cbcl_dat_all$src_subject_id.x))
out.sum.any$N_1_genetics = table(table(cbcl_dat_all$src_subject_id.x))[1] %>% unname()
out.sum.any$N_2_genetics = table(table(cbcl_dat_all$src_subject_id.x))[2] %>% unname()
out.sum.any$N_pre_genetics =(cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","mj_group")) %>% unique() %>% summarise(n = table(mj_group)))[2,1] %>% unname()
out.sum.any$N_pre_post_genetics = (cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","mj_group")) %>% unique() %>% summarise(n = table(mj_group)))[3,1] %>% unname()

#############

cbcl_dat_base = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id","site_id_l","rel_family_id","eventname",
      paste("site",c(1:21),sep="_"),
      "mj_pre","mj_pre_post","mj_pre_any","mj_group",
      "rel_relationship_1","rel_relationship_2","rel_relationship_3",
      "sex2",
      "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h",
      
      "dgs_before", "dgs_after",
      colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PRSCS")],
      colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PC")],
      "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
      "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
      
      "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p", "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p",
      "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco"#,
      #"income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      #"Prisma_fit","Prisma","DISCOVERY","Achieva",
      #"y.var","global.var","include.var","movement.var","x.var"
      #,"mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth"
      # excluding in some analyses -  mat_age_knew_preg, total_birth_weight, devhx_10, devhx_6_p, income      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1", mj_pre_post == 0) %>% 
  na.omit() %>% unique() 
cbcl_dat_base = cbcl_dat_base[,-which(colnames(cbcl_dat_base) == "eventname")]


cbcl_dat_mediator = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      "Prisma_fit","Prisma","DISCOVERY","Achieva",
      "global.var","include.var","movement.var","m.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1" | eventname == "2_year_follow_up_y_arm_1", include.var == 1) %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      
      "y.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "1_year_follow_up_y_arm_1" | eventname == "3_year_follow_up_y_arm_1") %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome$eventname[cbcl_dat_outcome$eventname == "1_year_follow_up_y_arm_1"] = "baseline_year_1_arm_1"
cbcl_dat_outcome$eventname[cbcl_dat_outcome$eventname == "3_year_follow_up_y_arm_1"] = "2_year_follow_up_y_arm_1"

cbcl_dat_outcome$combo_id = paste(cbcl_dat_outcome$src_subject_id,cbcl_dat_outcome$eventname,sep = "_")
cbcl_dat_mediator$combo_id = paste(cbcl_dat_mediator$src_subject_id,cbcl_dat_mediator$eventname,sep = "_")


cbcl_dat_all = merge(cbcl_dat_base,cbcl_dat_mediator)
cbcl_dat_all = merge(cbcl_dat_all,cbcl_dat_outcome,by = "combo_id")

sum(cbcl_dat_all$src_subject_id.x == cbcl_dat_all$src_subject_id.y)
sum(cbcl_dat_all$eventname.x == cbcl_dat_all$eventname.y)

#cbcl_dat_all = cbcl_dat_all %>% filter(eventname.x == "baseline_year_1_arm_1")

length(unique(cbcl_dat_all$src_subject_id.x))

temp.dat1 = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) %>% unique()
temp.dat = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) 

famids = table(temp.dat1$rel_family_id) %>% as.data.frame() 

# remove relatives, keeping any PCE participants

for(f in 1 :dim(famids)[1]){
  
  matching.fam = which(temp.dat$rel_family_id == famids$Var1[f])
  
  t2 = temp.dat[matching.fam,]
  
  t3 = table(t2$src_subject_id.x) %>% as.data.frame()
  colnames(t3) = c("src_subject_id.x","Freq")
  t2 = merge(t2,t3) %>% unique()
  
  
  if(sum(t2$mj_pre_any)>0){ t2 = t2[t2$mj_pre_any==1,]  }
  
  if(max(t2$Freq == 2)){ t2 = t2[t2$Freq==2,] }
  
  keep = t2[sample(x = c(1:dim(t2)[1]),size = 1),]
  
  
  if(f == 1){keep.ids = keep$src_subject_id.x}else{keep.ids = c(keep.ids,keep$src_subject_id.x)}
  
}

cbcl_dat_all = cbcl_dat_all [!is.na(match(cbcl_dat_all$src_subject_id.x,keep.ids)),]

length(unique(cbcl_dat_all$src_subject_id.x))


var.names=c( "mj_pre","mj_pre_post","mj_pre_any","rel_relationship_1","rel_relationship_2","rel_relationship_3",
             "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
             "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
             "dgs_after",#"income","interview_age",
             "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
             "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
             colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PRSCS")],
             colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PC")],
             "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
             "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
             "Prisma_fit","Prisma","DISCOVERY","Achieva","m.var",
             "y.var", "interview_age.x","interview_age.y",
             paste("site",c(2:21),sep="_"),
             "income_1.x","income_2.x","income_3.x","income_4.x","ed.x","pds.x",
             "income_1.y","income_2.y","income_3.y","income_4.y","ed.y","pds.y")

var.names2=c( #"mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
  "ed.x","pds.x","ed.y","pds.y",
  #"famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
  #"demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
  #"dgs_after",
  "interview_age.x","interview_age.y",
  colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PRSCS")],
  colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PC")],
  #"famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
  #"famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
  #"income_1","income_2","income_3","income_4",
  #"famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
  #"devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
  #"Prisma_fit","Prisma","DISCOVERY","Achieva",
  "y.var","m.var")


cbcl_dat_all = cbcl_dat_all %>% 
  mutate_at(.vars = var.names2,.funs = winsorize) %>% 
  mutate_at(.vars = var.names2,.funs = function(X){
    if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
  }) %>% 
  mutate_at(.vars = var.names,.funs = scale,center=T,scale=T)

if(analysis.framework$global.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "global.var",.funs = winsorize) %>% 
    mutate_at(.vars = "global.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "global.var",.funs = scale,center=T,scale=T)
  
}


if(analysis.framework$movement.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "movement.var",.funs = winsorize) %>% 
    mutate_at(.vars = "movement.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "movement.var",.funs = scale,center=T,scale=T)
  
}

cbcl_dat_all$Age1.x = poly(cbcl_dat_all$interview_age.x,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.x = poly(cbcl_dat_all$interview_age.x,2)[,2] %>% scale(center = T,scale = T)
cbcl_dat_all$Age1.y = poly(cbcl_dat_all$interview_age.y,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.y = poly(cbcl_dat_all$interview_age.y,2)[,2] %>% scale(center = T,scale = T)


run_med = function(){lmer(m.var ~ mj_pre_any + 
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            #ed.y+
                            income_1.y+income_2.y+income_3.y+income_4.y+
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            #demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            CUD_PRSCS + 
                            PC_1+PC_2+PC_3+PC_4+PC_5+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10+
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            site_2+site_3+site_4+site_5+site_6+site_7+site_8+site_9+site_10+site_11+
                            site_12+site_13+site_14+site_15+site_16+site_17+site_18+site_19+site_20+site_21+
                            #             (1|site_id_l) #+ 
                            #(1|rel_family_id)# 
                            (1|src_subject_id.x)
                          ,data =cbcl_dat_all  )}


med.fit = tryCatch(expr = run_med(),error = run_med()) # if error (b/c of memory), try again

run_out = function(){lmer(y.var ~ mj_pre_any + m.var+
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            #ed.y+
                            income_1.y+income_2.y+income_3.y+income_4.y+
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            #demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            CUD_PRSCS + 
                            PC_1+PC_2+PC_3+PC_4+PC_5+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10+
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            site_2+site_3+site_4+site_5+site_6+site_7+site_8+site_9+site_10+site_11+
                            site_12+site_13+site_14+site_15+site_16+site_17+site_18+site_19+site_20+site_21+
                            #  (1|site_id_l)  
                            #(1|rel_family_id) #+
                            (1|src_subject_id.x)
                          ,data =cbcl_dat_all  )}


out.fit = tryCatch(expr = run_out(),error = run_out()) # if error (b/c of memory), try again


med.out2a <- mediate(model.m = med.fit, 
                     model.y = out.fit, 
                     treat = "mj_pre_any", 
                     mediator = "m.var",
                     sims = 10000)

out.sum.a = c(med.out2a$d.avg, med.out2a$d.avg.ci[1],med.out2a$d.avg.ci[2],med.out2a$d.avg.p,
              med.out2a$z.avg, med.out2a$z.avg.ci[1],med.out2a$z.avg.ci[2],med.out2a$z.avg.p,
              med.out2a$tau.coef,med.out2a$tau.ci[1],med.out2a$tau.ci[2],med.out2a$tau.p,
              med.out2a$n.avg, med.out2a$n.avg.ci[1],med.out2a$n.avg.ci[2],med.out2a$n.avg.p) %>% t() %>% as.data.frame()
col.st = expand.grid(a = "genetics","pre",c =c("Est","L","U","p") ,b = c("ACME","ADE","Tot","Prop"))
colnames(out.sum.a) = apply(col.st,1,function(X){paste(X,collapse = "_",sep="")})


####

cbcl_dat_base = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id","site_id_l","rel_family_id","eventname",
      paste("site",c(1:21),sep="_"),
      "mj_pre","mj_pre_post","mj_pre_any","mj_group",
      "rel_relationship_1","rel_relationship_2","rel_relationship_3",
      "sex2",
      "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h",
      
      "dgs_before", "dgs_after",
      colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PRSCS")],
      colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PC")],
      "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
      "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
      
      "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p", "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p",
      "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco"#,
      #"income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      #"Prisma_fit","Prisma","DISCOVERY","Achieva",
      #"y.var","global.var","include.var","movement.var","x.var"
      #,"mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth"
      # excluding in some analyses -  mat_age_knew_preg, total_birth_weight, devhx_10, devhx_6_p, income      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1", mj_pre == 0) %>% 
  na.omit() %>% unique() 
cbcl_dat_base = cbcl_dat_base[,-which(colnames(cbcl_dat_base) == "eventname")]


cbcl_dat_mediator = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      "Prisma_fit","Prisma","DISCOVERY","Achieva",
      "global.var","include.var","movement.var","m.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1" | eventname == "2_year_follow_up_y_arm_1", include.var == 1) %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      
      "y.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "1_year_follow_up_y_arm_1" | eventname == "3_year_follow_up_y_arm_1") %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome$eventname[cbcl_dat_outcome$eventname == "1_year_follow_up_y_arm_1"] = "baseline_year_1_arm_1"
cbcl_dat_outcome$eventname[cbcl_dat_outcome$eventname == "3_year_follow_up_y_arm_1"] = "2_year_follow_up_y_arm_1"

cbcl_dat_outcome$combo_id = paste(cbcl_dat_outcome$src_subject_id,cbcl_dat_outcome$eventname,sep = "_")
cbcl_dat_mediator$combo_id = paste(cbcl_dat_mediator$src_subject_id,cbcl_dat_mediator$eventname,sep = "_")


cbcl_dat_all = merge(cbcl_dat_base,cbcl_dat_mediator)
cbcl_dat_all = merge(cbcl_dat_all,cbcl_dat_outcome,by = "combo_id")

sum(cbcl_dat_all$src_subject_id.x == cbcl_dat_all$src_subject_id.y)
sum(cbcl_dat_all$eventname.x == cbcl_dat_all$eventname.y)

#cbcl_dat_all = cbcl_dat_all %>% filter(eventname.x == "baseline_year_1_arm_1")

length(unique(cbcl_dat_all$src_subject_id.x))

temp.dat1 = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) %>% unique()
temp.dat = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) 

famids = table(temp.dat1$rel_family_id) %>% as.data.frame() 

# remove relatives, keeping any PCE participants

for(f in 1 :dim(famids)[1]){
  
  matching.fam = which(temp.dat$rel_family_id == famids$Var1[f])
  
  t2 = temp.dat[matching.fam,]
  
  t3 = table(t2$src_subject_id.x) %>% as.data.frame()
  colnames(t3) = c("src_subject_id.x","Freq")
  t2 = merge(t2,t3) %>% unique()
  
  
  if(sum(t2$mj_pre_any)>0){ t2 = t2[t2$mj_pre_any==1,]  }
  
  if(max(t2$Freq == 2)){ t2 = t2[t2$Freq==2,] }
  
  keep = t2[sample(x = c(1:dim(t2)[1]),size = 1),]
  
  
  if(f == 1){keep.ids = keep$src_subject_id.x}else{keep.ids = c(keep.ids,keep$src_subject_id.x)}
  
}

cbcl_dat_all = cbcl_dat_all [!is.na(match(cbcl_dat_all$src_subject_id.x,keep.ids)),]

length(unique(cbcl_dat_all$src_subject_id.x))


var.names=c( "mj_pre","mj_pre_post","mj_pre_any","rel_relationship_1","rel_relationship_2","rel_relationship_3",
             "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
             "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
             "dgs_after",#"income","interview_age",
             "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
             "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
             colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PRSCS")],
             colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PC")],
             "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
             "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
             "Prisma_fit","Prisma","DISCOVERY","Achieva","m.var",
             "y.var", "interview_age.x","interview_age.y",
             paste("site",c(2:21),sep="_"),
             "income_1.x","income_2.x","income_3.x","income_4.x","ed.x","pds.x",
             "income_1.y","income_2.y","income_3.y","income_4.y","ed.y","pds.y")

var.names2=c( #"mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
  "ed.x","pds.x","ed.y","pds.y",
  #"famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
  #"demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
  #"dgs_after",
  "interview_age.x","interview_age.y",
  colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PRSCS")],
  colnames(cbcl_long)[grep(x = colnames(cbcl_long),pattern = "PC")],
  #"famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
  #"famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
  #"income_1","income_2","income_3","income_4",
  #"famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
  #"devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
  #"Prisma_fit","Prisma","DISCOVERY","Achieva",
  "y.var","m.var")


cbcl_dat_all = cbcl_dat_all %>% 
  mutate_at(.vars = var.names2,.funs = winsorize) %>% 
  mutate_at(.vars = var.names2,.funs = function(X){
    if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
  }) %>% 
  mutate_at(.vars = var.names,.funs = scale,center=T,scale=T)

if(analysis.framework$global.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "global.var",.funs = winsorize) %>% 
    mutate_at(.vars = "global.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "global.var",.funs = scale,center=T,scale=T)
  
}


if(analysis.framework$movement.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "movement.var",.funs = winsorize) %>% 
    mutate_at(.vars = "movement.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "movement.var",.funs = scale,center=T,scale=T)
  
}

cbcl_dat_all$Age1.x = poly(cbcl_dat_all$interview_age.x,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.x = poly(cbcl_dat_all$interview_age.x,2)[,2] %>% scale(center = T,scale = T)
cbcl_dat_all$Age1.y = poly(cbcl_dat_all$interview_age.y,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.y = poly(cbcl_dat_all$interview_age.y,2)[,2] %>% scale(center = T,scale = T)


run_med = function(){lmer(m.var ~ mj_pre_any + 
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            #ed.y+
                            income_1.y+income_2.y+income_3.y+income_4.y+
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            #demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            CUD_PRSCS + 
                            PC_1+PC_2+PC_3+PC_4+PC_5+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10+
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            site_2+site_3+site_4+site_5+site_6+site_7+site_8+site_9+site_10+site_11+
                            site_12+site_13+site_14+site_15+site_16+site_17+site_18+site_19+site_20+site_21+
                            #             (1|site_id_l) #+ 
                            #(1|rel_family_id)# 
                            (1|src_subject_id.x)
                          ,data =cbcl_dat_all  )}


med.fit = tryCatch(expr = run_med(),error = run_med()) # if error (b/c of memory), try again

run_out = function(){lmer(y.var ~ mj_pre_any + m.var+
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            #ed.y+
                            income_1.y+income_2.y+income_3.y+income_4.y+
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            #demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            CUD_PRSCS + 
                            PC_1+PC_2+PC_3+PC_4+PC_5+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10+
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            site_2+site_3+site_4+site_5+site_6+site_7+site_8+site_9+site_10+site_11+
                            site_12+site_13+site_14+site_15+site_16+site_17+site_18+site_19+site_20+site_21+
                            #  (1|site_id_l)  
                            #(1|rel_family_id) #+
                            (1|src_subject_id.x)
                          ,data =cbcl_dat_all  )}


out.fit = tryCatch(expr = run_out(),error = run_out()) # if error (b/c of memory), try again


med.out2b <- mediate(model.m = med.fit, 
                     model.y = out.fit, 
                     treat = "mj_pre_any", 
                     mediator = "m.var",
                     sims = 10000)

out.sum.b = c(med.out2b$d.avg, med.out2b$d.avg.ci[1],med.out2b$d.avg.ci[2],med.out2b$d.avg.p,
              med.out2b$z.avg, med.out2b$z.avg.ci[1],med.out2b$z.avg.ci[2],med.out2b$z.avg.p,
              med.out2b$tau.coef,med.out2b$tau.ci[1],med.out2b$tau.ci[2],med.out2b$tau.p,
              med.out2b$n.avg, med.out2b$n.avg.ci[1],med.out2b$n.avg.ci[2],med.out2b$n.avg.p) %>% t() %>% as.data.frame()
col.st = expand.grid(a = "genetics","prepost",c =c("Est","L","U","p") ,b = c("ACME","ADE","Tot","Prop"))
colnames(out.sum.b) = apply(col.st,1,function(X){paste(X,collapse = "_",sep="")})


out_final_genetics = cbind(out.sum.any,out.sum.a,out.sum.b)

#####################


# baseline

print("baseline")

cbcl_dat_base = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id","site_id_l","rel_family_id","eventname",
      paste("site",c(1:21),sep="_"),
      "mj_pre","mj_pre_post","mj_pre_any","mj_group",
      "rel_relationship_1","rel_relationship_2","rel_relationship_3",
      "sex2",
      "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h",
      
      "dgs_before", "dgs_after",
      
      "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
      "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
      
      "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p", "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p",
      "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco"#,
      #"income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      #"Prisma_fit","Prisma","DISCOVERY","Achieva",
      #"y.var","global.var","include.var","movement.var","x.var"
      #,"mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth"
      # excluding in some analyses -  mat_age_knew_preg, total_birth_weight, devhx_10, devhx_6_p, income      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1") %>% 
  na.omit() %>% unique() 
cbcl_dat_base = cbcl_dat_base[,-which(colnames(cbcl_dat_base) == "eventname")]


cbcl_dat_mediator = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      "Prisma_fit","Prisma","DISCOVERY","Achieva",
      "global.var","include.var","movement.var","m.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1" , include.var == 1) %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      
      "y.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "1_year_follow_up_y_arm_1" ) %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome$eventname[cbcl_dat_outcome$eventname == "1_year_follow_up_y_arm_1"] = "baseline_year_1_arm_1"
#

cbcl_dat_outcome$combo_id = paste(cbcl_dat_outcome$src_subject_id,cbcl_dat_outcome$eventname,sep = "_")
cbcl_dat_mediator$combo_id = paste(cbcl_dat_mediator$src_subject_id,cbcl_dat_mediator$eventname,sep = "_")


cbcl_dat_all = merge(cbcl_dat_base,cbcl_dat_mediator)
cbcl_dat_all = merge(cbcl_dat_all,cbcl_dat_outcome,by = "combo_id")

sum(cbcl_dat_all$src_subject_id.x == cbcl_dat_all$src_subject_id.y)
sum(cbcl_dat_all$eventname.x == cbcl_dat_all$eventname.y)

#cbcl_dat_all = cbcl_dat_all %>% filter(eventname.x == "baseline_year_1_arm_1")

length(unique(cbcl_dat_all$src_subject_id.x))

temp.dat1 = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) %>% unique()
temp.dat = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) 

famids = table(temp.dat1$rel_family_id) %>% as.data.frame() 

# remove relatives, keeping any PCE participants

for(f in 1 :dim(famids)[1]){
  
  matching.fam = which(temp.dat$rel_family_id == famids$Var1[f])
  
  t2 = temp.dat[matching.fam,]
  
  t3 = table(t2$src_subject_id.x) %>% as.data.frame()
  colnames(t3) = c("src_subject_id.x","Freq")
  t2 = merge(t2,t3) %>% unique()
  
  
  if(sum(t2$mj_pre_any)>0){ t2 = t2[t2$mj_pre_any==1,]  }
  
  if(max(t2$Freq == 2)){ t2 = t2[t2$Freq==2,] }
  
  keep = t2[sample(x = c(1:dim(t2)[1]),size = 1),]
  
  
  if(f == 1){keep.ids = keep$src_subject_id.x}else{keep.ids = c(keep.ids,keep$src_subject_id.x)}
  
}

cbcl_dat_all = cbcl_dat_all [!is.na(match(cbcl_dat_all$src_subject_id.x,keep.ids)),]

length(unique(cbcl_dat_all$src_subject_id.x))
table(cbcl_dat_all$mj_group)

var.names=c( "mj_pre","mj_pre_post","mj_pre_any","rel_relationship_1","rel_relationship_2","rel_relationship_3",
             "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
             "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
             "dgs_after",#"income","interview_age",
             "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
             "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
             
             "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
             "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
             "Prisma_fit","Prisma","DISCOVERY","Achieva","m.var",
             "y.var", "interview_age.x","interview_age.y",
             paste("site",c(2:21),sep="_"),
             "income_1.x","income_2.x","income_3.x","income_4.x","ed.x","pds.x",
             "income_1.y","income_2.y","income_3.y","income_4.y","ed.y","pds.y")

var.names2=c( #"mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
  "ed.x","pds.x","ed.y","pds.y",
  #"famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
  #"demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
  #"dgs_after",
  "interview_age.x","interview_age.y",
  #"famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
  #"famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
  #"income_1","income_2","income_3","income_4",
  #"famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
  #"devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
  #"Prisma_fit","Prisma","DISCOVERY","Achieva",
  "y.var","m.var")


cbcl_dat_all = cbcl_dat_all %>% 
  mutate_at(.vars = var.names2,.funs = winsorize) %>% 
  mutate_at(.vars = var.names2,.funs = function(X){
    if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
  }) %>% 
  mutate_at(.vars = var.names,.funs = scale,center=T,scale=T)

if(analysis.framework$global.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "global.var",.funs = winsorize) %>% 
    mutate_at(.vars = "global.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "global.var",.funs = scale,center=T,scale=T)
  
}


if(analysis.framework$movement.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "movement.var",.funs = winsorize) %>% 
    mutate_at(.vars = "movement.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "movement.var",.funs = scale,center=T,scale=T)
  
}

cbcl_dat_all$Age1.x = poly(cbcl_dat_all$interview_age.x,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.x = poly(cbcl_dat_all$interview_age.x,2)[,2] %>% scale(center = T,scale = T)
cbcl_dat_all$Age1.y = poly(cbcl_dat_all$interview_age.y,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.y = poly(cbcl_dat_all$interview_age.y,2)[,2] %>% scale(center = T,scale = T)


run_med = function(){lmer(m.var ~ mj_pre_any + 
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            
                            
                            #             (1|site_id_l) #+ 
                            #(1|rel_family_id)# 
                            (1|site_id_l) 
                          ,data =cbcl_dat_all  )}


med.fit = tryCatch(expr = run_med(),error = run_med()) # if error (b/c of memory), try again

run_out = function(){lmer(y.var ~ mj_pre_any + m.var+
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            
                            
                            #  (1|site_id_l)  
                            #(1|rel_family_id) #+
                            (1|site_id_l) 
                          ,data =cbcl_dat_all  )}


out.fit = tryCatch(expr = run_out(),error = run_out()) # if error (b/c of memory), try again


med.out <- mediate(model.m = med.fit, 
                   model.y = out.fit, 
                   treat = "mj_pre_any", 
                   mediator = "m.var",
                   sims = 1000)

out.sum.any = c(med.out$d.avg, med.out$d.avg.ci[1],med.out$d.avg.ci[2],med.out$d.avg.p,
                med.out$z.avg, med.out$z.avg.ci[1],med.out$z.avg.ci[2],med.out$z.avg.p,
                med.out$tau.coef,med.out$tau.ci[1],med.out$tau.ci[2],med.out$tau.p,
                med.out$n.avg, med.out$n.avg.ci[1],med.out$n.avg.ci[2],med.out$n.avg.p) %>% t() %>% as.data.frame()
col.st = expand.grid(a = "baseline","Any",c =c("Est","L","U","p") ,b = c("ACME","ADE","Tot","Prop"))
colnames(out.sum.any) = apply(col.st,1,function(X){paste(X,collapse = "_",sep="")})



out.sum.any$N_sub_baseline = length(unique(cbcl_dat_all$src_subject_id.x))
out.sum.any$N_obs_baseline = length((cbcl_dat_all$src_subject_id.x))
out.sum.any$N_1_baseline = table(table(cbcl_dat_all$src_subject_id.x))[1] %>% unname()
out.sum.any$N_2_baseline = table(table(cbcl_dat_all$src_subject_id.x))[2] %>% unname()
out.sum.any$N_pre_baseline =(cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","mj_group")) %>% unique() %>% summarise(n = table(mj_group)))[2,1] %>% unname()
out.sum.any$N_pre_post_baseline = (cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","mj_group")) %>% unique() %>% summarise(n = table(mj_group)))[3,1] %>% unname()

#############

cbcl_dat_base = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id","site_id_l","rel_family_id","eventname",
      paste("site",c(1:21),sep="_"),
      "mj_pre","mj_pre_post","mj_pre_any","mj_group",
      "rel_relationship_1","rel_relationship_2","rel_relationship_3",
      "sex2",
      "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h",
      
      "dgs_before", "dgs_after",
      
      "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
      "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
      
      "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p", "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p",
      "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco"#,
      #"income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      #"Prisma_fit","Prisma","DISCOVERY","Achieva",
      #"y.var","global.var","include.var","movement.var","x.var"
      #,"mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth"
      # excluding in some analyses -  mat_age_knew_preg, total_birth_weight, devhx_10, devhx_6_p, income      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1", mj_pre_post == 0) %>% 
  na.omit() %>% unique() 
cbcl_dat_base = cbcl_dat_base[,-which(colnames(cbcl_dat_base) == "eventname")]


cbcl_dat_mediator = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      "Prisma_fit","Prisma","DISCOVERY","Achieva",
      "global.var","include.var","movement.var","m.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1" , include.var == 1) %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      
      "y.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "1_year_follow_up_y_arm_1" ) %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome$eventname[cbcl_dat_outcome$eventname == "1_year_follow_up_y_arm_1"] = "baseline_year_1_arm_1"


cbcl_dat_outcome$combo_id = paste(cbcl_dat_outcome$src_subject_id,cbcl_dat_outcome$eventname,sep = "_")
cbcl_dat_mediator$combo_id = paste(cbcl_dat_mediator$src_subject_id,cbcl_dat_mediator$eventname,sep = "_")


cbcl_dat_all = merge(cbcl_dat_base,cbcl_dat_mediator)
cbcl_dat_all = merge(cbcl_dat_all,cbcl_dat_outcome,by = "combo_id")

sum(cbcl_dat_all$src_subject_id.x == cbcl_dat_all$src_subject_id.y)
sum(cbcl_dat_all$eventname.x == cbcl_dat_all$eventname.y)

#cbcl_dat_all = cbcl_dat_all %>% filter(eventname.x == "baseline_year_1_arm_1")

length(unique(cbcl_dat_all$src_subject_id.x))

temp.dat1 = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) %>% unique()
temp.dat = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) 

famids = table(temp.dat1$rel_family_id) %>% as.data.frame() 

# remove relatives, keeping any PCE participants

for(f in 1 :dim(famids)[1]){
  
  matching.fam = which(temp.dat$rel_family_id == famids$Var1[f])
  
  t2 = temp.dat[matching.fam,]
  
  t3 = table(t2$src_subject_id.x) %>% as.data.frame()
  colnames(t3) = c("src_subject_id.x","Freq")
  t2 = merge(t2,t3) %>% unique()
  
  
  if(sum(t2$mj_pre_any)>0){ t2 = t2[t2$mj_pre_any==1,]  }
  
  if(max(t2$Freq == 2)){ t2 = t2[t2$Freq==2,] }
  
  keep = t2[sample(x = c(1:dim(t2)[1]),size = 1),]
  
  
  if(f == 1){keep.ids = keep$src_subject_id.x}else{keep.ids = c(keep.ids,keep$src_subject_id.x)}
  
}

cbcl_dat_all = cbcl_dat_all [!is.na(match(cbcl_dat_all$src_subject_id.x,keep.ids)),]
length(unique(cbcl_dat_all$src_subject_id.x))


var.names=c( "mj_pre","mj_pre_post","mj_pre_any","rel_relationship_1","rel_relationship_2","rel_relationship_3",
             "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
             "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
             "dgs_after",#"income","interview_age",
             "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
             "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
             
             "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
             "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
             "Prisma_fit","Prisma","DISCOVERY","Achieva","m.var",
             "y.var", "interview_age.x","interview_age.y",
             paste("site",c(2:21),sep="_"),
             "income_1.x","income_2.x","income_3.x","income_4.x","ed.x","pds.x",
             "income_1.y","income_2.y","income_3.y","income_4.y","ed.y","pds.y")

var.names2=c( #"mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
  "ed.x","pds.x","ed.y","pds.y",
  #"famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
  #"demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
  #"dgs_after",
  "interview_age.x","interview_age.y",
  #"famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
  #"famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
  #"income_1","income_2","income_3","income_4",
  #"famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
  #"devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
  #"Prisma_fit","Prisma","DISCOVERY","Achieva",
  "y.var","m.var")


cbcl_dat_all = cbcl_dat_all %>% 
  mutate_at(.vars = var.names2,.funs = winsorize) %>% 
  mutate_at(.vars = var.names2,.funs = function(X){
    if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
  }) %>% 
  mutate_at(.vars = var.names,.funs = scale,center=T,scale=T)

if(analysis.framework$global.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "global.var",.funs = winsorize) %>% 
    mutate_at(.vars = "global.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "global.var",.funs = scale,center=T,scale=T)
  
}


if(analysis.framework$movement.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "movement.var",.funs = winsorize) %>% 
    mutate_at(.vars = "movement.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "movement.var",.funs = scale,center=T,scale=T)
  
}

cbcl_dat_all$Age1.x = poly(cbcl_dat_all$interview_age.x,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.x = poly(cbcl_dat_all$interview_age.x,2)[,2] %>% scale(center = T,scale = T)
cbcl_dat_all$Age1.y = poly(cbcl_dat_all$interview_age.y,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.y = poly(cbcl_dat_all$interview_age.y,2)[,2] %>% scale(center = T,scale = T)


run_med = function(){lmer(m.var ~ mj_pre_any + 
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            
                            
                            #             (1|site_id_l) #+ 
                            #(1|rel_family_id)# 
                            (1|site_id_l) 
                          ,data =cbcl_dat_all  )}


med.fit = tryCatch(expr = run_med(),error = run_med()) # if error (b/c of memory), try again

run_out = function(){lmer(y.var ~ mj_pre_any + m.var+
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            
                            
                            #  (1|site_id_l)  
                            #(1|rel_family_id) #+
                            (1|site_id_l) 
                          ,data =cbcl_dat_all  )}


out.fit = tryCatch(expr = run_out(),error = run_out()) # if error (b/c of memory), try again


med.out2a <- mediate(model.m = med.fit, 
                     model.y = out.fit, 
                     treat = "mj_pre_any", 
                     mediator = "m.var",
                     sims = 10000)

out.sum.a = c(med.out2a$d.avg, med.out2a$d.avg.ci[1],med.out2a$d.avg.ci[2],med.out2a$d.avg.p,
              med.out2a$z.avg, med.out2a$z.avg.ci[1],med.out2a$z.avg.ci[2],med.out2a$z.avg.p,
              med.out2a$tau.coef,med.out2a$tau.ci[1],med.out2a$tau.ci[2],med.out2a$tau.p,
              med.out2a$n.avg, med.out2a$n.avg.ci[1],med.out2a$n.avg.ci[2],med.out2a$n.avg.p) %>% t() %>% as.data.frame()
col.st = expand.grid(a = "baseline","pre",c =c("Est","L","U","p") ,b = c("ACME","ADE","Tot","Prop"))
colnames(out.sum.a) = apply(col.st,1,function(X){paste(X,collapse = "_",sep="")})


####

cbcl_dat_base = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id","site_id_l","rel_family_id","eventname",
      paste("site",c(1:21),sep="_"),
      "mj_pre","mj_pre_post","mj_pre_any","mj_group",
      "rel_relationship_1","rel_relationship_2","rel_relationship_3",
      "sex2",
      "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h",
      
      "dgs_before", "dgs_after",
      
      "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
      "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
      
      "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p", "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p",
      "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco"#,
      #"income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      #"Prisma_fit","Prisma","DISCOVERY","Achieva",
      #"y.var","global.var","include.var","movement.var","x.var"
      #,"mat_age_knew_preg", "total_birth_weight","devhx_10","devhx_6_p","mat_age_birth"
      # excluding in some analyses -  mat_age_knew_preg, total_birth_weight, devhx_10, devhx_6_p, income      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1", mj_pre == 0) %>% 
  na.omit() %>% unique() 
cbcl_dat_base = cbcl_dat_base[,-which(colnames(cbcl_dat_base) == "eventname")]


cbcl_dat_mediator = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      "Prisma_fit","Prisma","DISCOVERY","Achieva",
      "global.var","include.var","movement.var","m.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "baseline_year_1_arm_1" , include.var == 1) %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome = cbcl_long %>%
  dplyr::select(all_of(
    c(
      "src_subject_id",#"site_id_l","rel_family_id",
      "eventname",
      
      "income","interview_age","ed","pds", "income_1","income_2","income_3","income_4",
      
      
      "y.var"
      
    ))) %>% #dplyr::filter(include.var == 1) %>% 
  dplyr::filter(eventname == "1_year_follow_up_y_arm_1" ) %>% 
  na.omit() %>% unique() 

cbcl_dat_outcome$eventname[cbcl_dat_outcome$eventname == "1_year_follow_up_y_arm_1"] = "baseline_year_1_arm_1"


cbcl_dat_outcome$combo_id = paste(cbcl_dat_outcome$src_subject_id,cbcl_dat_outcome$eventname,sep = "_")
cbcl_dat_mediator$combo_id = paste(cbcl_dat_mediator$src_subject_id,cbcl_dat_mediator$eventname,sep = "_")


cbcl_dat_all = merge(cbcl_dat_base,cbcl_dat_mediator)
cbcl_dat_all = merge(cbcl_dat_all,cbcl_dat_outcome,by = "combo_id")

sum(cbcl_dat_all$src_subject_id.x == cbcl_dat_all$src_subject_id.y)
sum(cbcl_dat_all$eventname.x == cbcl_dat_all$eventname.y)

#cbcl_dat_all = cbcl_dat_all %>% filter(eventname.x == "baseline_year_1_arm_1")

length(unique(cbcl_dat_all$src_subject_id.x))

temp.dat1 = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) %>% unique()
temp.dat = cbcl_dat_all %>% dplyr::select(c("src_subject_id.x","rel_family_id","mj_pre_any")) 

famids = table(temp.dat1$rel_family_id) %>% as.data.frame() 

# remove relatives, keeping any PCE participants

for(f in 1 :dim(famids)[1]){
  
  matching.fam = which(temp.dat$rel_family_id == famids$Var1[f])
  
  t2 = temp.dat[matching.fam,]
  
  t3 = table(t2$src_subject_id.x) %>% as.data.frame()
  colnames(t3) = c("src_subject_id.x","Freq")
  t2 = merge(t2,t3) %>% unique()
  
  
  if(sum(t2$mj_pre_any)>0){ t2 = t2[t2$mj_pre_any==1,]  }
  
  if(max(t2$Freq == 2)){ t2 = t2[t2$Freq==2,] }
  
  keep = t2[sample(x = c(1:dim(t2)[1]),size = 1),]
  
  
  if(f == 1){keep.ids = keep$src_subject_id.x}else{keep.ids = c(keep.ids,keep$src_subject_id.x)}
  
}

cbcl_dat_all = cbcl_dat_all [!is.na(match(cbcl_dat_all$src_subject_id.x,keep.ids)),]

length(unique(cbcl_dat_all$src_subject_id.x))


var.names=c( "mj_pre","mj_pre_post","mj_pre_any","rel_relationship_1","rel_relationship_2","rel_relationship_3",
             "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
             "demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
             "dgs_after",#"income","interview_age",
             "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
             "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
             
             "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
             "devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
             "Prisma_fit","Prisma","DISCOVERY","Achieva","m.var",
             "y.var", "interview_age.x","interview_age.y",
             paste("site",c(2:21),sep="_"),
             "income_1.x","income_2.x","income_3.x","income_4.x","ed.x","pds.x",
             "income_1.y","income_2.y","income_3.y","income_4.y","ed.y","pds.y")

var.names2=c( #"mj_pre","mj_pre_post","rel_relationship_1","rel_relationship_2","rel_relationship_3",
  "ed.x","pds.x","ed.y","pds.y",
  #"famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2",
  #"demo_race_w","demo_race_b","demo_race_na","demo_race_pi","demo_race_h","dgs_before",
  #"dgs_after",
  "interview_age.x","interview_age.y",
  #"famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", "famhx_ss_firstdeg_vs_p", 
  #"famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
  #"income_1","income_2","income_3","income_4",
  #"famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p",
  #"devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco",
  #"Prisma_fit","Prisma","DISCOVERY","Achieva",
  "y.var","m.var")


cbcl_dat_all = cbcl_dat_all %>% 
  mutate_at(.vars = var.names2,.funs = winsorize) %>% 
  mutate_at(.vars = var.names2,.funs = function(X){
    if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
  }) %>% 
  mutate_at(.vars = var.names,.funs = scale,center=T,scale=T)

if(analysis.framework$global.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "global.var",.funs = winsorize) %>% 
    mutate_at(.vars = "global.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "global.var",.funs = scale,center=T,scale=T)
  
}


if(analysis.framework$movement.var[args] != "none"){
  
  cbcl_dat_all = cbcl_dat_all %>% 
    mutate_at(.vars = "movement.var",.funs = winsorize) %>% 
    mutate_at(.vars = "movement.var",.funs = function(X){
      if(abs(psych::skew(X)) > 1){inormal(X)}else{X}
    }) %>% 
    mutate_at(.vars = "movement.var",.funs = scale,center=T,scale=T)
  
}

cbcl_dat_all$Age1.x = poly(cbcl_dat_all$interview_age.x,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.x = poly(cbcl_dat_all$interview_age.x,2)[,2] %>% scale(center = T,scale = T)
cbcl_dat_all$Age1.y = poly(cbcl_dat_all$interview_age.y,2)[,1] %>% scale(center = T,scale = T)
cbcl_dat_all$Age2.y = poly(cbcl_dat_all$interview_age.y,2)[,2] %>% scale(center = T,scale = T)


run_med = function(){lmer(m.var ~ mj_pre_any + 
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            
                            
                            #             (1|site_id_l) #+ 
                            #(1|rel_family_id)# 
                            (1|site_id_l) 
                          ,data =cbcl_dat_all  )}


med.fit = tryCatch(expr = run_med(),error = run_med()) # if error (b/c of memory), try again

run_out = function(){lmer(y.var ~ mj_pre_any + m.var+
                            Age1.x+Age2.x+Age1.y+Age2.y+
                            
                            ed.x+income_1.x+income_2.x+income_3.x+income_4.x+pds.x+
                            
                            pds.y+
                            
                            sex2 +
                            famhx_ss_firstdeg_alc_p+famhx_ss_firstdeg_dg_p+ famhx_ss_firstdeg_dprs_p+
                            famhx_ss_firstdeg_ma_p+  famhx_ss_firstdeg_trb_p+famhx_ss_firstdeg_nrv_p +
                            demo_race_w+demo_race_b+demo_race_na+demo_race_pi+demo_race_h+
                            dgs_before+dgs_after+
                            #mat_age_knew_preg+total_birth_weight+devhx_10+devhx_6_p+mat_age_birth+
                            
                            devhx_8_alcohol+ devhx_9_alcohol+devhx_8_tobacco+ devhx_9_tobacco+
                            Prisma_fit+Prisma+DISCOVERY+Achieva+movement.var+global.var+
                            
                            
                            #  (1|site_id_l)  
                            #(1|rel_family_id) #+
                            (1|site_id_l) 
                          ,data =cbcl_dat_all  )}


out.fit = tryCatch(expr = run_out(),error = run_out()) # if error (b/c of memory), try again


med.out2b <- mediate(model.m = med.fit, 
                     model.y = out.fit, 
                     treat = "mj_pre_any", 
                     mediator = "m.var",
                     sims = 10000)

out.sum.b = c(med.out2b$d.avg, med.out2b$d.avg.ci[1],med.out2b$d.avg.ci[2],med.out2b$d.avg.p,
              med.out2b$z.avg, med.out2b$z.avg.ci[1],med.out2b$z.avg.ci[2],med.out2b$z.avg.p,
              med.out2b$tau.coef,med.out2b$tau.ci[1],med.out2b$tau.ci[2],med.out2b$tau.p,
              med.out2b$n.avg, med.out2b$n.avg.ci[1],med.out2b$n.avg.ci[2],med.out2b$n.avg.p) %>% t() %>% as.data.frame()
col.st = expand.grid(a = "baseline","prepost",c =c("Est","L","U","p") ,b = c("ACME","ADE","Tot","Prop"))
colnames(out.sum.b) = apply(col.st,1,function(X){paste(X,collapse = "_",sep="")})


out_final_baseline = cbind(out.sum.any,out.sum.a,out.sum.b)





##############################################################################################
model.info = data.frame(M =  analysis.framework$m.var[args], Y = analysis.framework$y.var[args] )

out = cbind(model.info,out_final_main,out_final_extra,out_final_genetics,out_final_baseline)




write.table(x = out,file = paste(file_location,"PCE_cbcl_brain_mediation.",args[1],".csv",sep=""),quote = F,row.names = F,col.names = F,sep = ",")

if(args[1] == 1){
  
  write.table(x = t(colnames(out)),file = paste(file_location,"PCE_cbcl_brain_mediation.0.csv",sep=""),quote = F,row.names = F,col.names = F,sep = ",")
  
}


