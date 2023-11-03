# clean & merge ABCD data v 5.0
# David Baranger - 7/7/2023

library(data.table)
library(dplyr)
library(openxlsx)
library(foreach)
library(stringr)



setwd("C:/Users/David/Documents/BRAINlab/ABCD/Package_1195434/ABCD_5/core")

datadictionary = readxl::read_xlsx('data.dictionary.xlsx')

#########################################

file.name = c("./physical-health/ph_p_dhx.csv")
file_in = fread(input =file.name,data.table = F,header = T,na.strings = "",#sep="\t",
                stringsAsFactors = F)
file_in = file_in %>% dplyr::filter(eventname == "baseline_year_1_arm_1")

file_in[file_in == 999] = NA
file_in[file_in == 777] = NA

file_in[file_in == ""] = NA

baseline_data = file_in

########################################

file.name = c("./abcd-general/abcd_p_demo.csv")
file_in = fread(input =file.name,data.table = F,header = T,na.strings = "",#sep="\t",
                stringsAsFactors = F)
file_in[file_in == 999] = NA
file_in[file_in == 777] = NA

file_in[file_in == ""] = NA


long.dat = file_in[,c(1,2,grep(x = colnames(file_in),pattern = "_ed_|_income_"))]
file_in = file_in %>% dplyr::filter(eventname == "baseline_year_1_arm_1")

baseline_data = merge(baseline_data,file_in)
long_data = long.dat

######################

file.name = "./mental-health/mh_p_fhx.csv"
file_in = fread(input =file.name,data.table = F,header = T,na.strings = "",#sep="\t",
                stringsAsFactors = F)
file_in[file_in == 999] = NA
file_in[file_in == 777] = NA

file_in[file_in == ""] = NA



prob.names = colnames(file_in)[grep(x = colnames(file_in),pattern = "famhx_ss_fath|famhx_ss_moth|famhx_ss_fulsib")] %>% as.data.frame()
prob.names = str_split(string = prob.names$.,pattern =  "_",simplify = T) %>% as.data.frame()
prob.names = unique(prob.names$V5)

probs = file_in[,grep(x = colnames(file_in),pattern = "famhx_ss_fath|famhx_ss_moth|famhx_ss_fulsib")]

prob.out = matrix(data = NA,ncol = length(prob.names),nrow = dim(probs)[1]) %>% as.data.frame()

colnames(prob.out) = paste("famhx_ss_firstdeg_",prob.names,"_p",sep="")

for(i in 1:length(prob.names)){
  
  any_hist = probs[,grep(colnames(probs),pattern = prob.names[i])]
  
  prob.out[,i] = apply(any_hist,1,function(X){
    
    (sum(na.omit(X)) > 1)*1
    
    })
  
}

prob.out$src_subject_id = file_in$src_subject_id

baseline_data = merge(baseline_data,prob.out,all=T)


##############################

file.name = "./abcd-general/abcd_y_lt.csv"
file_in = fread(input =file.name,data.table = F,header = T,na.strings = "",#sep="\t",
                stringsAsFactors = F)
file_in[file_in == 999] = NA
file_in[file_in == 777] = NA

file_in[file_in == ""] = NA

file_in$site_id_l[file_in$site_id_l == "site22"] = "site21"

long.dat = file_in[,c(1,2,3,grep(x = colnames(file_in),pattern = "age|type"))]

file_in = file_in %>% dplyr::filter(eventname == "baseline_year_1_arm_1")
file_in = file_in[,c(1:5)]

long_data = merge(long_data,long.dat,all=T)
baseline_data = merge(file_in,baseline_data,all=T)


##################

file.name = "./genetics/gen_y_pihat.csv"
file_in = fread(input =file.name,data.table = F,header = T,na.strings = "",#sep="\t",
                stringsAsFactors = F)
file_in[file_in == 999] = NA
file_in[file_in == 777] = NA

file_in[file_in == ""] = NA

baseline_data = merge(file_in,baseline_data,all=T)


################################################

merge_in_long=function(path){
  
  file.name = path
  file_in = fread(input =file.name,data.table = F,header = T,na.strings = "",#sep="\t",
                  stringsAsFactors = F)
  file_in[file_in == 999] = NA
  file_in[file_in == 777] = NA
  
  #file_in[file_in == ""] = NA
  
  long_data = merge(long_data,file_in,all=T)
  return(long_data)
}

long_data = merge_in_long(path = "./physical-health/ph_p_pds.csv")
long_data = merge_in_long(path = "./physical-health/ph_y_pds.csv")


long_data = merge_in_long(path = "./mental-health/mh_p_cbcl.csv")
long_data = merge_in_long(path = "./neurocognition/nc_y_nihtb.csv")
long_data = merge_in_long(path = "./neurocognition/nc_y_ravlt.csv")
long_data = merge_in_long(path = "./mental-health/mh_y_pps.csv")

long_data = merge_in_long(path = "./mental-health/mh_y_upps.csv")
long_data = merge_in_long(path = "./mental-health/mh_y_bisbas.csv")
long_data = merge_in_long(path = "./mental-health/mh_y_poa.csv")
long_data = merge_in_long(path = "./neurocognition/nc_y_cct.csv")
long_data = merge_in_long(path = "./neurocognition/nc_y_lmt.csv")
long_data = merge_in_long(path = "./neurocognition/nc_y_wisc.csv")


long_data = merge_in_long(path = "./imaging/mri_y_adm_info.csv")
long_data = merge_in_long(path = "./imaging/mri_y_qc_incl.csv")

long_data = merge_in_long(path = "./imaging/mri_y_qc_auto_post.csv")
long_data = merge_in_long(path = "./imaging/mri_y_qc_man_fsurf.csv")
long_data = merge_in_long(path = "./imaging/mri_y_qc_man_post_dmr.csv")
long_data = merge_in_long(path = "./imaging/mri_y_qc_man_post_fmr.csv")
long_data = merge_in_long(path = "./imaging/mri_y_qc_motion.csv")


long_data = merge_in_long(path = "./imaging/mri_y_rsfmr_cor_gp_aseg.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsfmr_cor_gp_gp.csv")

long_data = merge_in_long(path = "./imaging/mri_y_rsi_fni_aseg.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsi_fni_at.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsi_fni_wm_dsk.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsi_fni_gm_dsk.csv")

long_data = merge_in_long(path = "./imaging/mri_y_rsi_hnd_aseg.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsi_hnd_at.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsi_hnd_wm_dsk.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsi_hnd_gm_dsk.csv")

long_data = merge_in_long(path = "./imaging/mri_y_rsi_hni_aseg.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsi_hni_at.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsi_hni_wm_dsk.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsi_hni_gm_dsk.csv")

long_data = merge_in_long(path = "./imaging/mri_y_rsi_hnt_aseg.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsi_hnt_at.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsi_hnt_wm_dsk.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsi_hnt_gm_dsk.csv")

long_data = merge_in_long(path = "./imaging/mri_y_rsi_rnd_aseg.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsi_rnd_at.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsi_rnd_wm_dsk.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsi_rnd_gm_dsk.csv")

long_data = merge_in_long(path = "./imaging/mri_y_rsi_rni_aseg.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsi_rni_at.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsi_rni_wm_dsk.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsi_rni_gm_dsk.csv")

long_data = merge_in_long(path = "./imaging/mri_y_rsi_rnt_aseg.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsi_rnt_at.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsi_rnt_wm_dsk.csv")
long_data = merge_in_long(path = "./imaging/mri_y_rsi_rnt_gm_dsk.csv")

long_data = merge_in_long(path = "./imaging/mri_y_smr_area_dsk.csv")
long_data = merge_in_long(path = "./imaging/mri_y_smr_sulc_dsk.csv")
long_data = merge_in_long(path = "./imaging/mri_y_smr_thk_dsk.csv")
long_data = merge_in_long(path = "./imaging/mri_y_smr_vol_aseg.csv")
long_data = merge_in_long(path = "./imaging/mri_y_smr_vol_dsk.csv")



##################################################

long_data = merge_in_long(path = "./imaging/mri_y_dti_fa_is_aseg.csv")
long_data = merge_in_long(path = "./imaging/mri_y_dti_fa_is_at.csv")
long_data = merge_in_long(path = "./imaging/mri_y_dti_fa_is_gm_dsk.csv")
long_data = merge_in_long(path = "./imaging/mri_y_dti_fa_is_wm_dsk.csv")
#long_data = merge_in_long(path = "./imaging/mri_y_dti_fa_is_wm_dsk.csv")
long_data = merge_in_long(path = "./imaging/mri_y_dti_ld_is_aseg.csv")
long_data = merge_in_long(path = "./imaging/mri_y_dti_ld_is_at.csv")
long_data = merge_in_long(path = "./imaging/mri_y_dti_ld_is_gm_dsk.csv")
long_data = merge_in_long(path = "./imaging/mri_y_dti_ld_is_wm_dsk.csv")
#long_data = merge_in_long(path = "./imaging/mri_y_dti_ld_is_wm_dsk.csv")
long_data = merge_in_long(path = "./imaging/mri_y_dti_md_is_aseg.csv")
long_data = merge_in_long(path = "./imaging/mri_y_dti_md_is_at.csv")
long_data = merge_in_long(path = "./imaging/mri_y_dti_md_is_gm_dsk.csv")
long_data = merge_in_long(path = "./imaging/mri_y_dti_md_is_wm_dsk.csv")
#long_data = merge_in_long(path = "./imaging/mri_y_dti_md_is_wm_dsk.csv")
long_data = merge_in_long(path = "./imaging/mri_y_dti_td_is_aseg.csv")
long_data = merge_in_long(path = "./imaging/mri_y_dti_td_is_at.csv")
long_data = merge_in_long(path = "./imaging/mri_y_dti_td_is_gm_dsk.csv")
#long_data = merge_in_long(path = "./imaging/mri_y_dti_td_is_wm_dsk.csv")
long_data = merge_in_long(path = "./imaging/mri_y_dti_td_is_wm_dsk.csv")

long_data = merge_in_long(path = "./imaging/mri_y_dti_vol_is_at.csv")

write.csv(x = baseline_data,file = "baseline_data.csv",row.names = F)
write.csv(x = long_data,file = "long_data.csv",row.names = F)


#######
