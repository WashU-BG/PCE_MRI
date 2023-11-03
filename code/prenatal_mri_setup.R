
library(data.table)
library(dplyr)
library(openxlsx)
library(foreach)
library(stringr)


library(data.table)
library(dplyr)
library(openxlsx)
library(foreach)
library(psych)
library(lmerTest)
library(lme4)
library(ggplot2)
library(ggpubr)
#library(dhglm)
library(tidyr)

setwd("C:/Users/David/Documents/BRAINlab/ABCD/Package_1195434/ABCD_5/core")

baseline_data = fread("baseline_data.csv",data.table = F)
long_data = fread("long_data.csv",data.table = F)

###################################################################################


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



############### Setup baseline data


famhx = baseline_data[,grep(colnames(baseline_data),pattern = "famhx_ss_firstdeg")]
baseline_data$famhx_ss_firstdeg_alc_dg_p = any_non_zero(famhx[,c(1,2)])
baseline_data$famhx_ss_firstdeg_mh_p = any_non_zero(famhx[,-c(1,2)])

# 
# famhx2 = baseline_data[,grep(colnames(baseline_data),pattern = "recode")]
# baseline_data$famhx_ss_firstdeg_alc_dg_p2 = any_non_zero(famhx2[,grep(colnames(famhx2),
#                                                                       pattern = "famhx4a_p___0_recode|famhx4b_p___0_recode|.*full_sib.*alc.*|fam_history_q5a_drugs___0_recode|fam_history_q5b_drugs___0_recode|.*full_sib.*drugs.*"
# )])
# 
# table(baseline_data$famhx_ss_firstdeg_alc_dg_p,baseline_data$famhx_ss_firstdeg_alc_dg_p2)

baseline_data$sex2 = baseline_data$demo_sex_v2
#baseline_data$sex2[baseline_data$demo_sex_v2 == "F"] = 2
baseline_data$sex2[baseline_data$sex2 == 3] = NA

baseline_data$demo_race_w = baseline_data$demo_race_a_p___10
baseline_data$demo_race_b = baseline_data$demo_race_a_p___11
baseline_data$demo_race_na = any_non_zero(cbind(baseline_data$demo_race_a_p___12, baseline_data$demo_race_a_p___13))
baseline_data$demo_race_pi = any_non_zero(cbind(baseline_data$demo_race_a_p___14 , baseline_data$demo_race_a_p___15 , 
                                                baseline_data$demo_race_a_p___16 , baseline_data$demo_race_a_p___17) )
baseline_data$demo_race_a = any_non_zero(cbind(baseline_data$demo_race_a_p___18 , baseline_data$demo_race_a_p___19 , 
                                               baseline_data$demo_race_a_p___20 , baseline_data$demo_race_a_p___21,
                                               baseline_data$demo_race_a_p___22 , baseline_data$demo_race_a_p___23,
                                               baseline_data$demo_race_a_p___24) )
baseline_data$demo_ethn_v2[baseline_data$demo_ethn_v2 == 777] = NA
baseline_data$demo_race_h = baseline_data$demo_ethn_v2
baseline_data$demo_race_h[baseline_data$demo_race_h == 2] = 0

baseline_data$demo_race_o = baseline_data$demo_race_a_p___25



#prenatalexposuretotobaccooralcoholbeforeor  aftermaternalknowledgeofpregnancy,
# 
# other_drugs8 = cbind(baseline_data$devhx_8_other1_name_2,
#                     baseline_data$devhx_8_other2_name_2,
#                     baseline_data$devhx_8_other3_name_2,
#                     baseline_data$devhx_8_other4_name_2,
#                     baseline_data$devhx_8_other5_name_2)
# other_drugs8[other_drugs8 == 0] = NA
# other_drugs8[other_drugs8 == 3] = NA
# other_drugs8[other_drugs8 == 12] = NA
# other_drugs8[other_drugs8 == 999] = NA
# 
# other_drugs9 = cbind(baseline_data$devhx_9_other1_name_2,
#                      baseline_data$devhx_9_other2_name_2,
#                      baseline_data$devhx_9_other3_name_2,
#                      baseline_data$devhx_9_other4_name_2,
#                      baseline_data$devhx_9_other5_name_2)
# other_drugs9[other_drugs9 == 0] = NA
# other_drugs9[other_drugs9 == 3] = NA
# other_drugs9[other_drugs9 == 12] = NA
# other_drugs9[other_drugs9 == 999] = NA


baseline_data$dgs_before = any_non_zero(cbind(
  baseline_data$devhx_8_coc_crack,
  
  baseline_data$devhx_8_her_morph,
  baseline_data$devhx_8_oxycont))

baseline_data$dgs_after = 
  any_non_zero(cbind(
    baseline_data$devhx_9_coc_crack,
    
    baseline_data$devhx_9_her_morph,
    baseline_data$devhx_9_oxycont))

baseline_data$any_mj = any_non_zero(cbind(baseline_data$devhx_8_marijuana, 
                                          baseline_data$devhx_9_marijuana))



baseline_data$devhx_10[baseline_data$devhx_10 == -1] = NA

baseline_data$birth_weight_oz[is.na(baseline_data$birth_weight_oz)] = 0
baseline_data$total_birth_weight = (baseline_data$birth_weight_lbs * 16) + baseline_data$birth_weight_oz
baseline_data$total_birth_weight = winsorize(baseline_data$total_birth_weight)
baseline_data$mat_age_birth = winsorize(baseline_data$devhx_3_p)
baseline_data$mat_age_knew_preg = winsorize(log(baseline_data$devhx_7_p+1))
# 
# baseline_data$try_dg = any_non_zero(cbind(baseline_data$tlfb_alc_sip, baseline_data$tlfb_tob_puff,
#                                           baseline_data$tlfb_alc_use  ))



baseline_data$mj_pre_post = baseline_data$devhx_9_marijuana
baseline_data$mj_pre      = baseline_data$devhx_8_marijuana
baseline_data$mj_pre[baseline_data$mj_pre_post == 1] = 0



#dat = na.omit(dat)

###################################################
# Long covariates


baseline_data$demo_prnt_ed_v2[baseline_data$demo_prnt_ed_v2 == 777] = NA
baseline_data$demo_prnt_ed_v2[baseline_data$demo_prnt_ed_v2 == 13] = 12
baseline_data$demo_prnt_ed_v2[baseline_data$demo_prnt_ed_v2 == 14] = 12
baseline_data$demo_prnt_ed_v2[baseline_data$demo_prnt_ed_v2 == 22] = 14
baseline_data$demo_prnt_ed_v2[baseline_data$demo_prnt_ed_v2 == 23] = 14
baseline_data$demo_prnt_ed_v2[baseline_data$demo_prnt_ed_v2 == 17] = 14
baseline_data$demo_prnt_ed_v2[baseline_data$demo_prnt_ed_v2 == 15] = 14
baseline_data$demo_prnt_ed_v2[baseline_data$demo_prnt_ed_v2 == 16] = 14
baseline_data$demo_prnt_ed_v2[baseline_data$demo_prnt_ed_v2 == 18] = 16
baseline_data$demo_prnt_ed_v2[baseline_data$demo_prnt_ed_v2 == 19] = 18
baseline_data$demo_prnt_ed_v2[baseline_data$demo_prnt_ed_v2 == 21] = 20
baseline_data$demo_prnt_ed_v2[baseline_data$demo_prnt_ed_v2  <  9] = 8


baseline_data$demo_comb_income_v2[baseline_data$demo_comb_income_v2 == 777] = NA
baseline_data$demo_comb_income_v2[baseline_data$demo_comb_income_v2 <=  6] = 1
baseline_data$demo_comb_income_v2[baseline_data$demo_comb_income_v2 == 7] = 2
baseline_data$demo_comb_income_v2[baseline_data$demo_comb_income_v2 == 8] = 3
baseline_data$demo_comb_income_v2[baseline_data$demo_comb_income_v2 == 9] = 4
baseline_data$demo_comb_income_v2[baseline_data$demo_comb_income_v2 == 10] = 5



long_data$demo_prnt_ed_v2_2yr_l[long_data$demo_prnt_ed_v2_2yr_l == 777] = NA
long_data$demo_prnt_ed_v2_2yr_l[long_data$demo_prnt_ed_v2_2yr_l == 13] = 12
long_data$demo_prnt_ed_v2_2yr_l[long_data$demo_prnt_ed_v2_2yr_l == 14] = 12
long_data$demo_prnt_ed_v2_2yr_l[long_data$demo_prnt_ed_v2_2yr_l == 22] = 14
long_data$demo_prnt_ed_v2_2yr_l[long_data$demo_prnt_ed_v2_2yr_l == 23] = 14
long_data$demo_prnt_ed_v2_2yr_l[long_data$demo_prnt_ed_v2_2yr_l == 17] = 14
long_data$demo_prnt_ed_v2_2yr_l[long_data$demo_prnt_ed_v2_2yr_l == 15] = 14
long_data$demo_prnt_ed_v2_2yr_l[long_data$demo_prnt_ed_v2_2yr_l == 16] = 14
long_data$demo_prnt_ed_v2_2yr_l[long_data$demo_prnt_ed_v2_2yr_l == 18] = 16
long_data$demo_prnt_ed_v2_2yr_l[long_data$demo_prnt_ed_v2_2yr_l == 19] = 18
long_data$demo_prnt_ed_v2_2yr_l[long_data$demo_prnt_ed_v2_2yr_l == 21] = 20
long_data$demo_prnt_ed_v2_2yr_l[long_data$demo_prnt_ed_v2_2yr_l  <  9] = 8

long_data$demo_prnt_ed_v2_l[long_data$demo_prnt_ed_v2_l == 777] = NA
long_data$demo_prnt_ed_v2_l[long_data$demo_prnt_ed_v2_l == 13] = 12
long_data$demo_prnt_ed_v2_l[long_data$demo_prnt_ed_v2_l == 14] = 12
long_data$demo_prnt_ed_v2_l[long_data$demo_prnt_ed_v2_l == 22] = 14
long_data$demo_prnt_ed_v2_l[long_data$demo_prnt_ed_v2_l == 23] = 14
long_data$demo_prnt_ed_v2_l[long_data$demo_prnt_ed_v2_l == 17] = 14
long_data$demo_prnt_ed_v2_l[long_data$demo_prnt_ed_v2_l == 15] = 14
long_data$demo_prnt_ed_v2_l[long_data$demo_prnt_ed_v2_l == 16] = 14
long_data$demo_prnt_ed_v2_l[long_data$demo_prnt_ed_v2_l == 18] = 16
long_data$demo_prnt_ed_v2_l[long_data$demo_prnt_ed_v2_l == 19] = 18
long_data$demo_prnt_ed_v2_l[long_data$demo_prnt_ed_v2_l == 21] = 20
long_data$demo_prnt_ed_v2_l[long_data$demo_prnt_ed_v2_l  <  9] = 8



long_data$demo_comb_income_v2_l[long_data$demo_comb_income_v2_l == 777] = NA
long_data$demo_comb_income_v2_l[long_data$demo_comb_income_v2_l <=  6] = 1
long_data$demo_comb_income_v2_l[long_data$demo_comb_income_v2_l == 7] = 2
long_data$demo_comb_income_v2_l[long_data$demo_comb_income_v2_l == 8] = 3
long_data$demo_comb_income_v2_l[long_data$demo_comb_income_v2_l == 9] = 4
long_data$demo_comb_income_v2_l[long_data$demo_comb_income_v2_l == 10] = 5



var_keep = c("src_subject_id","rel_family_id",#,
             "site_id_l",
             "rel_group_id",
             "rel_ingroup_order","rel_relationship",
             "mj_pre","mj_pre_post","any_mj","devhx_9_marijuana","devhx_8_marijuana",#"mj_total_base",
             "famhx_ss_firstdeg_alc_dg_p","famhx_ss_firstdeg_mh_p","sex2","demo_race_w","demo_race_b","demo_race_na",
             "demo_race_pi","demo_race_h","demo_race_o","demo_race_a","dgs_before","dgs_after","devhx_10","devhx_6_p","total_birth_weight",
             "famhx_ss_firstdeg_alc_p","famhx_ss_firstdeg_dg_p", "famhx_ss_firstdeg_dprs_p","famhx_ss_firstdeg_ma_p", 
             "famhx_ss_firstdeg_vs_p", "famhx_ss_firstdeg_trb_p","famhx_ss_firstdeg_nrv_p" ,
             "mat_age_birth","mat_age_knew_preg","devhx_8_alcohol","devhx_9_alcohol","devhx_8_tobacco","devhx_9_tobacco")


dat = baseline_data %>% dplyr::select(all_of(var_keep))

all.missing = apply(dat,MARGIN = 2,function(X){
  sum(is.na(X)) 
  
})


all.missing

long_data = long_data %>% dplyr::filter(  eventname == "baseline_year_1_arm_1" | 
                                            eventname == "2_year_follow_up_y_arm_1" )

long_temp = long_data %>% dplyr::select(#c(2:6),
                                        "src_subject_id"  ,    "eventname"  ,
                                      "interview_age",
                                        "demo_prnt_ed_v2_2yr_l","demo_prnt_ed_v2_l",
                                        "demo_comb_income_v2_l",
                                        #                                               "tlfb_cal_scr_num_events",
                                        "pds_p_ss_male_category","pds_p_ss_female_category",
                                        "pds_y_ss_male_category","pds_y_ss_female_category")


base_temp = baseline_data %>% dplyr::select(c(1),"demo_prnt_ed_v2","demo_comb_income_v2")


all_dat = merge(long_temp,base_temp,all=T)

all_dat$income = coalesce(all_dat$demo_comb_income_v2,all_dat$demo_comb_income_v2_l)
all_dat$ed = coalesce(all_dat$demo_prnt_ed_v2,all_dat$demo_prnt_ed_v2_l,all_dat$demo_prnt_ed_v2_2yr_l)
all_dat$pds = coalesce(all_dat$pds_p_ss_male_category,all_dat$pds_p_ss_female_category,
                       all_dat$pds_y_ss_male_category,all_dat$pds_y_ss_female_category)

#all_dat = all_dat[apply(all_dat[,c( (dim(all_dat)[2]-2): dim(all_dat)[2])],1,function(X){sum(is.na(X)) < 3}),]



all_dat = all_dat[order(all_dat$interview_age,decreasing = F),]

all_dat = all_dat %>% dplyr::select(c("src_subject_id","eventname","interview_age","income","ed","pds"))

all_dat = all_dat[!is.na(all_dat$interview_age),]


IDs = unique(all_dat$src_subject_id)

out = matrix(NA,nrow = dim(all_dat)[1],ncol = dim(all_dat)[2]) %>% as.data.frame()


#########################################

r = 1
for(i in 1:length(IDs)){
  #if(r > dim(out)[1]){r = 1}
  temp = all_dat %>% dplyr::filter(src_subject_id == IDs[i])
  
  if(dim(temp)[1] > 1){
    temp = zoo::na.locf.default(object = temp,na.rm = F) #locf
    temp = temp[c(dim(temp)[1]:1),]
    temp = zoo::na.locf.default(object = temp,na.rm = F) #locb
    temp = temp[c(dim(temp)[1]:1),]
  }
  
  
  for(z in 1: dim(temp)[1]){
    out[r,] = temp[z,]
    r=r+1
    #if(r == dim(out)[1]){stop()}
  }
}

colnames(out) = colnames(all_dat)

out[,c(3:6)] = apply(out[,c(3:6)] ,MARGIN = 2,as.numeric)

###################################################################
sui = fread("C:/Users/David/Documents/BRAINlab/ABCD/Package_1195434/ABCD_SUI_R5.csv",header = T,data.table = F)
sui = sui[,c(1,8,9)]
sui = sui %>% pivot_longer(cols = starts_with("any"),names_prefix = "any_",names_to = "eventname",values_to = "sui")

#table(cbcl_long$eventname)
sui$eventname = as.factor(sui$eventname)
levels(sui$eventname) = c("baseline_year_1_arm_1","2_year_follow_up_y_arm_1")
sui$eventname = as.character(sui$eventname)
colnames(sui)[1] = "src_subject_id"

##################################################################


long2 = long_data %>% dplyr::select(c("src_subject_id","eventname","visit_type"),
                                    "rsfmri_meanmotion",
                                    "dmri_meanmotion",
                                    "imgincl_t1w_include",
                                    "imgincl_dmri_include",
                                    "imgincl_rsfmri_include",
                                    "mri_info_manufacturersmn",
                                    "mri_info_softwareversion",

                                    grep(colnames(long_data),pattern = "^cbcl.+_r$"),
                                    grep(colnames(long_data),pattern = "pps_y_ss_number"),
                                    
                                    grep(colnames(long_data),pattern = "rsfmri"),
                                    grep(colnames(long_data),pattern = "^smri"),
                                    #grep(colnames(long_data),pattern = "^mri"),
                                    
                                    grep(colnames(long_data),pattern = "^dmri_dti"),
                                    grep(colnames(long_data),pattern = "^dmri_rsi")
                                    
)

duplicates = long2 %>% t() %>% duplicated()
# rsfmri data have duplicates - it's the full square matrix. 

long2 = long2[,!duplicates]
dim(long2)[2]

#long2 = merge(long2,QC.info,all=T)
long2 = merge(out,long2,all=T)

cbcl_long = merge(dat,long2,by = "src_subject_id",all=T)

cbcl_long = cbcl_long %>% dplyr::filter(  eventname == "baseline_year_1_arm_1" | 
                                            # eventname == "1_year_follow_up_y_arm_1" #| 
                                            eventname == "2_year_follow_up_y_arm_1" 
)




########################################################
cbcl_long = merge(cbcl_long,sui,all=T)

cbcl_long = cbcl_long %>% dplyr::filter(sui ==0)

########################################################


cbcl_long$MJ_exposure = "None"
cbcl_long$MJ_exposure[cbcl_long$mj_pre == max(cbcl_long$mj_pre %>% na.omit())] = "Pre-knowledge only"
cbcl_long$MJ_exposure[cbcl_long$mj_pre_post == max(cbcl_long$mj_pre_post %>% na.omit())] = "Pre & Post-knowledge"

cbcl_long$MJ_any = 0
cbcl_long$MJ_any[cbcl_long$mj_pre == max(cbcl_long$mj_pre)] = 1
cbcl_long$MJ_any[cbcl_long$mj_pre_post == max(cbcl_long$mj_pre_post )] = 198


cbcl_long$demo_race_h[is.na(cbcl_long$demo_race_h)] = min(cbcl_long$demo_race_h,na.rm = T)


cbcl_long$rel_family_id[is.na(cbcl_long$rel_family_id)]=0

cbcl_long$income_1 = (cbcl_long$income ==1)*1
cbcl_long$income_2 = (cbcl_long$income ==2)*1
cbcl_long$income_3 = (cbcl_long$income ==3)*1
cbcl_long$income_4 = (cbcl_long$income ==4)*1


################

cbcl_long$Prisma_fit   <- ifelse(cbcl_long$mri_info_manufacturersmn  == 'Prisma_fit', 1, 0)
cbcl_long$Prisma <- ifelse(cbcl_long$mri_info_manufacturersmn  == 'Prisma', 1, 0)
cbcl_long$DISCOVERY <- ifelse(cbcl_long$mri_info_manufacturersmn  == 'DISCOVERY MR750', 1, 0)
cbcl_long$Achieva  <- ifelse(cbcl_long$mri_info_manufacturersmn  == 'Achieva dStream', 1, 0)

#table(cbcl_long$visit_type)

################################


#,"CROSS","SCZ","RISK","MDD","GAD","ADHD"
traits = c("CUD")

for(i in 1:length(traits)){
  gwas_dat = fread(paste("C:/Users/David/Documents/BRAINlab/ABCD/genetics/ABCD_",traits[i],"_PRSCS.sscore",sep=""),header = T,data.table = F)
  gwas_dat = gwas_dat %>% dplyr::select(c("IID","#FID","SCORE1_AVG"))
  colnames(gwas_dat) = c("src_subject_id","FID",paste(traits[i],"_PRSCS",sep=""))
  
  if(i==1){gwas_dat.out = gwas_dat}else{gwas_dat.out = merge(gwas_dat.out,gwas_dat,all=T)}
  
}

PCs =   gwas_dat = fread("C:/Users/David/Documents/BRAINlab/ABCD/genetics/white_no_hisp_pca_noref.menv.trans.mds",header = F,data.table = F)
PCs = PCs[,c(2,5:24)]
colnames(PCs) = c("src_subject_id",paste("PC",c(1:20),sep="_"))

gwas_dat.out = merge(gwas_dat.out,PCs)

cbcl_long = merge(cbcl_long,gwas_dat.out,all=T)




########################
write.csv(x = cbcl_long,file = "allimaging_pce.csv",quote = F,row.names = F)
