library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)

setwd('C:/Users/dbara/Documents/WashU/ABCD/MRI')

all_results = fread("all_PCE_regressions_10.csv",data.table = F)


all_results$type = NA
all_results$global = NA



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

all_results$type[grep(all_results$Y,pattern = "rsfmri_c_ngd_.*_ngd_.*")] = "rsmfri_cortical"
all_results$type[grep(all_results$Y,pattern = "rsfmri_cor")] = "rsmfri_subcortical"



for(i in 1:length(global.vars)){
  v.temp = global.vars[i]
  new.temp = str_split(string = v.temp,pattern = "_",simplify =T)[-4] %>% paste(collapse = "_")
  all_results$type[grep(all_results$Y,pattern = new.temp)] = v.temp
}


all_results$type[grep(all_results$Y,pattern = "smri_vol")] = "smri_vol_scs_intracranialv"


all_results$global[match(global.vars,all_results$Y)] = 1


all_results$type = str_replace(string = all_results$type,pattern = "_scs",replacement = "gm_cdk")
all_results$type = str_replace(string = all_results$type,pattern = "_scts",replacement = "gm_cdsn")

####################################################################################################



all_results$P_fdr = NA
all_results$P_bon = NA


types = names(table(all_results$type))

for(i in 1:length(types)){

  all_results$P_fdr[grep(all_results$type,pattern = types[i])] = all_results$Chisq_P_main[grep(all_results$type,pattern = types[i])] %>% p.adjust(method = "fdr")
  all_results$P_bon[grep(all_results$type,pattern = types[i])] = all_results$Chisq_P_main[grep(all_results$type,pattern = types[i])] %>% p.adjust(method = "bonferroni")

}


all_results$P_fdr[grep(all_results$global,pattern = 1)] = all_results$Chisq_P_main[grep(all_results$global,pattern = 1)] %>% p.adjust(method = "fdr")
all_results$P_bon[grep(all_results$global,pattern = 1)] = all_results$Chisq_P_main[grep(all_results$global,pattern = 1)] %>% p.adjust(method = "bonferroni")



all_results = all_results[order(all_results$P_fdr),]

#all_results$Expanded = str_replace_all(string = all_results$Expanded,pattern = " ",replacement = "_")
#all_results$Expanded = str_replace_all(string = all_results$Expanded,pattern = ",",replacement = "_")



all_results$P_fdr_all = all_results$Chisq_P_main %>% p.adjust(method = "fdr")
all_results$P_bon_all = all_results$Chisq_P_main %>% p.adjust(method = "bonferroni")



write.csv(x = all_results,file ="all_PCE_regressions_10_pfdr.csv",row.names = F,quote = F )



################

setwd('C:/Users/dbara/Documents/WashU/ABCD/MRI')

cbcl_results = fread("all_PCE_cbcl_regression_5.csv",data.table = F)

cbcl_results$type = NA
cbcl_results$global = NA



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

cbcl_results$type[grep(cbcl_results$X,pattern = "rsfmri_c_ngd_.*_ngd_.*")] = "rsmfri_cortical"
cbcl_results$type[grep(cbcl_results$X,pattern = "rsfmri_cor")] = "rsmfri_subcortical"



for(i in 1:length(global.vars)){
  v.temp = global.vars[i]
  new.temp = str_split(string = v.temp,pattern = "_",simplify =T)[-4] %>% paste(collapse = "_")
  cbcl_results$type[grep(cbcl_results$X,pattern = new.temp)] = v.temp
}


cbcl_results$type[grep(cbcl_results$X,pattern = "smri_vol")] = "smri_vol_scs_intracranialv"


cbcl_results$global[match(global.vars,cbcl_results$X)] = 1

####################################################################################################



cbcl_results$P_main_fdr = NA
cbcl_results$P_main_bon = NA


cbcl_results$P_basic_fdr = NA
cbcl_results$P_basic_bon = NA

types = names(table(cbcl_results$type))

for(i in 1:length(types)){
  
  cbcl_results$P_main_fdr[grep(cbcl_results$type,pattern = types[i])] = cbcl_results$`main_Pr(>|t|)`[grep(cbcl_results$type,pattern = types[i])] %>% p.adjust(method = "fdr")
  cbcl_results$P_main_bon[grep(cbcl_results$type,pattern = types[i])] = cbcl_results$`main_Pr(>|t|)`[grep(cbcl_results$type,pattern = types[i])] %>% p.adjust(method = "bonferroni")
  
  cbcl_results$P_basic_fdr[grep(cbcl_results$type,pattern = types[i])] = cbcl_results$`basic_Pr(>|t|)`[grep(cbcl_results$type,pattern = types[i])] %>% p.adjust(method = "fdr")
  cbcl_results$P_basic_bon[grep(cbcl_results$type,pattern = types[i])] = cbcl_results$`basic_Pr(>|t|)`[grep(cbcl_results$type,pattern = types[i])] %>% p.adjust(method = "bonferroni")
  
}


#cbcl_results$P_fdr[grep(cbcl_results$global,pattern = 1)] = cbcl_results$Chisq_P_main[grep(cbcl_results$global,pattern = 1)] %>% p.adjust(method = "fdr")
#cbcl_results$P_bon[grep(cbcl_results$global,pattern = 1)] = cbcl_results$Chisq_P_main[grep(cbcl_results$global,pattern = 1)] %>% p.adjust(method = "bonferroni")



cbcl_results = cbcl_results[order(cbcl_results$P_basic_fdr),]



cbcl_results$Main_all_p_fdr = cbcl_results$`main_Pr(>|t|)` %>% p.adjust(method = "fdr")
cbcl_results$Basic_all_p_fdr = cbcl_results$`basic_Pr(>|t|)` %>% p.adjust(method = "fdr")

write.csv(x =cbcl_results,file = "all_PCE_cbcl_regression_5_fdr.csv", row.names = F,quote = F)


##############################################


setwd('C:/Users/dbara/Documents/WashU/ABCD/MRI')

cbcl_results = fread("allimaging_pce_mediation_2.csv",data.table = F)

cbcl_results$Main_Any_p_ACME_fdr = p.adjust(p = cbcl_results$main_Any_p_ACME,method = "fdr")

write.csv(x =cbcl_results,file = "all_PCE_cbcl_mediation_fdr_3.csv", row.names = F,quote = F)
