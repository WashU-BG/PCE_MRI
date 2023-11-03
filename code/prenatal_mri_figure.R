
library(ggplot2)
library(dplyr)


#
# library(ggcorrplot)
#
 library(ggpubr)
#
# #library(ggalluvial)
# library(data.table)
# library(stringr)
# #library(mygene)
# #library(readxl)
# #library(ACAT)
# #library(biomaRt)
# library(epiR)
# library(ggrepel)
# library(ggpubr)
#
# #library(corrplot)
library(MetBrewer)
library(tidyr)
library(stringr)

dat=read.csv("C:/Users/dbara/Documents/WashU/ABCD/MRI/all_PCE_regressions_10_toplot.csv",header = T)
#colnames(dat)[1]=  "Metric"


dat$Group[dat$Group == "mj_pre"] = "Exposed pre-knowledge only"
dat$Group[dat$Group == "mj_pre_post"] = "Exposed pre- and post-knowledge"

dat$Metric2 = factor(dat$Metric,levels =rev(unique(dat$Metric)))


#dat$Metric2 = factor(dat$Metric2,levels = levels(dat$Metric2)[c(10,2,3,5,7,8,6,9,1,4,11)])

dat$Region2 = factor(dat$Region,levels =(unique(dat$Region)))
#dat$Region2 = factor(dat$Region2,levels = levels(dat$Region2)[c(2,4,3,5,8,6,7,9,1)])

fig1 = ggplot(data = dat,aes(x =Metric2 ,y=Estimate,fill=Group,color=Group))+
  
  scale_fill_manual(values = met.brewer(name = "Egypt",n = 2))+
  scale_color_manual(values = met.brewer(name = "Egypt",n = 2))+
  
  
  scale_y_continuous(limits = c(-.061,.061),breaks = seq(-.08,.06,.02))+
  
  ylab(label = "Standardized effect size (95% CI) relative to no-exposure")+
 # xlab(label = "Child Psychopathology Scales")+
  geom_hline(yintercept = 0,linetype="dashed",color="black",size=.75)+
  geom_hline(yintercept = c(-0.06,-0.04,-0.02,0.02,0.04,0.06),linetype="dotted",color="light grey",size=.75)+
  

  
  coord_flip()+
  
  geom_linerange( aes(xmin=Metric2, xmax=Metric2, ymin=L, ymax=U,color=Group),alpha=.5,
                  stat = "identity",position = position_dodge(width = .75),size=1.25,
                  show.legend = F)+
    geom_vline(xintercept = seq(0.5,8.5,1),color="darkgrey",size=.5)+

  geom_point(stat = "identity",color="black",
             position = position_dodge(width = .75),shape=21,size=2.5)+
  
  theme_classic()+
  
  

  # theme(axis.text.x = element_text(angle=90),
  #       legend.position = "top")+
  theme(axis.text = element_text(color="black"
                                 #,size = 10
  ),
  #axis.title.x = element_text(hjust=1),
  #  legend.text = element_text(size = 10),
  legend.title = element_blank(),
  #plot.margin = margin(40,0,0,5),
  #legend.position = c(-.16,1.12),
  axis.title.y = element_blank(),
  legend.position = "top",
  axis.text.x = element_text(angle = 0)#,
  #axis.title = element_text(size = 10)
  )+
  guides(fill = guide_legend(nrow = 2,reverse = T),color = guide_legend(nrow = 2,reverse = T)) + 
facet_wrap(facets = ~ Region2,nrow = 3)
fig1

ggsave(plot = fig1,filename = "C:/Users/dbara/Documents/WashU/ABCD/MRI/pce_mri_results_fig1.jpeg",
       dpi=500,width = 11,height = 6)
ggsave(plot = fig1,filename = "C:/Users/dbara/Documents/WashU/ABCD/MRI/pce_mri_results_fig1.png",
       dpi=500,width = 11,height = 6)
ggsave(plot = fig1,filename = "C:/Users/dbara/Documents/WashU/ABCD/MRI/pce_mri_results_fig1.tiff",
       dpi=500,width = 11,height = 6)
ggsave(plot = fig1,filename = "C:/Users/dbara/Documents/WashU/ABCD/MRI/pce_mri_results_fig1.pdf",
       dpi=500,width = 11,height = 6)

#####################################################

dat=read.csv("C:/Users/dbara/Documents/WashU/ABCD/MRI/all_PCE_cbcl_regression_4_toplot.csv",header = T)

dat$label = paste(dat$Metric,"of",dat$Region,"  -->  ",dat$Yexp)

# dat$first = paste(dat$Metric,"of",dat$Region,"to")
# 
# dat$first = str_pad(dat$first,width = 85,side = "right",pad = "-")
# dat$second = str_pad(dat$Yexp,width = 25,side = "left",pad = "-")
# dat$label = paste(dat$first,dat$second,sep="")


#dat$label = sapply(dat$label,FUN = function(x){ stringr::str_replace(x, ' to ', strrep(' ', 120 - nchar(x)))})

# 
# mx <- as.numeric(max(nchar(dat$label)))
# 
# target = 100
# 
# for(i in 1: length(dat$label)){
# #first =    paste(dat$Metric[i],"of",dat$Region[i])
# second =   dat$Yexp[i]
# extra = 25-nchar(second)
# 
# #orig = nchar(paste(first,second))
# 
# #pad = round((target - orig - 3)/2)
# 
# #dat$label[i] = paste(first,paste(rep(" ",pad),collapse=""),"-->",paste(rep(" ",pad),collapse=""),second,sep="")
# 
# 
# dat$label[i] =stringr::str_replace(dat$label[i] , ' to ', strrep(' ', target - nchar(dat$label[i]) + extra))
# 
# }


#dat$expanded <- str_pad(dat$expanded, mx, side = "left", pad = " ")

dat$label2 = factor(dat$label,levels = rev(unique(dat$label)))



fig2 = ggplot(data = dat,aes(x =label2 ,y=basic_Estimate))+
  
  scale_fill_manual(values = met.brewer(name = "Egypt")[3])+
  scale_color_manual(values = met.brewer(name = "Egypt")[3])+
  
  
  scale_y_continuous(limits = c(-.115,0.115),breaks = c(seq(-.1,.1,.02)))+
  
  ylab(label = "Standardized effect size (95% CI)")+
  #xlab(label = "Child Psychopathology Scales")+
  
  scale_x_discrete(labels = function(label2){str_replace(string = label2,pattern ="-->",replacement = "\u2194" )})+
  
  geom_hline(yintercept = 0,linetype="dashed",color="black",size=.75)+
  geom_hline(yintercept = c(-.1,-.08,-0.06,-0.04,-0.02,0.02,0.04,0.06,.08,.1),linetype="dotted",color="light grey",size=.75)+
  
  
  
  geom_linerange( aes(xmin=label2, xmax=label2, ymin=L_CI_basic, ymax=U_CI_basic,color = met.brewer(name = "Egypt")[3]),alpha=.5,
                  stat = "identity",position = position_dodge(width = .75),size=1.25,
                  show.legend = F)+
  
  ggtitle(label = "Cross-sectional association")+
  geom_point(stat = "identity",color="black",fill = met.brewer(name = "Egypt")[3],
             position = position_dodge(width = .75),shape=21,size=2.5)+
  theme_classic()+
  coord_flip()+
  geom_vline(xintercept = seq(1.5,12.5,1),color="darkgrey",size=.5)+
  # theme(axis.text.x = element_text(angle=90),
  #       legend.position = "top")+
  theme(axis.text = element_text(color="black"
                                 #,size = 10
  ),
  axis.title.y = element_blank(),
  #  legend.text = element_text(size = 10),
  legend.title = element_blank(),
  #plot.margin = margin(40,0,0,5),
  #legend.position = c(0.45,1.1),
  #legend.position = "top",
  axis.text.x = element_text(angle = 0)#,
  #axis.title = element_text(size = 10)
  )

fig2
ggsave(plot = fig2,filename = "C:/Users/dbara/Documents/WashU/ABCD/MRI/pce_mri_results_fig2.jpeg",
       dpi=500,width = 10,height = 2)
ggsave(plot = fig2,filename = "C:/Users/dbara/Documents/WashU/ABCD/MRI/pce_mri_results_fig2.png",
       dpi=500,width = 10,height = 2)
ggsave(plot = fig2,filename = "C:/Users/dbara/Documents/WashU/ABCD/MRI/pce_mri_results_fig2.tiff",
       dpi=500,width = 10,height = 2)
ggsave(plot = fig2,filename = "C:/Users/dbara/Documents/WashU/ABCD/MRI/pce_mri_results_fig2.pdf",
       dpi=500,width = 10,height = 2)

############################################


dat=read.csv("C:/Users/dbara/Documents/WashU/ABCD/MRI/all_PCE_cbcl_mediation_fdr_toplot_3.csv",header = T)

dat$label = paste(dat$Metric,"of",dat$Region,"  -->  ",dat$Yexp)

dat$Group[dat$Group == "any"] = "Any exposure"
dat$Group[dat$Group == "pre"] = "Exposed pre-knowledge only"
dat$Group[dat$Group == "pre_post"] = "Exposed pre- and post-knowledge"
dat$Group2 = as.factor(dat$Group)
#dat$Metric2 = factor(dat$Metric,levels = unique(dat$Metric))


#dat$Metric2 = factor(dat$Metric2,levels = levels(dat$Metric2)[c(10,2,3,5,7,8,6,9,1,4,11)])


# dat$first = paste(dat$Metric,"of",dat$Region,"to")
# 
# dat$first = str_pad(dat$first,width = 85,side = "right",pad = "-")
# dat$second = str_pad(dat$Yexp,width = 25,side = "left",pad = "-")
# dat$label = paste(dat$first,dat$second,sep="")


#dat$label = sapply(dat$label,FUN = function(x){ stringr::str_replace(x, ' to ', strrep(' ', 120 - nchar(x)))})

# 
# mx <- as.numeric(max(nchar(dat$label)))
# 
# target = 100
# 
# for(i in 1: length(dat$label)){
# #first =    paste(dat$Metric[i],"of",dat$Region[i])
# second =   dat$Yexp[i]
# extra = 25-nchar(second)
# 
# #orig = nchar(paste(first,second))
# 
# #pad = round((target - orig - 3)/2)
# 
# #dat$label[i] = paste(first,paste(rep(" ",pad),collapse=""),"-->",paste(rep(" ",pad),collapse=""),second,sep="")
# 
# 
# dat$label[i] =stringr::str_replace(dat$label[i] , ' to ', strrep(' ', target - nchar(dat$label[i]) + extra))
# 
# }


#dat$expanded <- str_pad(dat$expanded, mx, side = "left", pad = " ")

dat$label2 = factor(dat$label,levels = rev(unique(dat$label)))



fig3 = ggplot(data = dat,aes(x =label2 ,y=Est,fill=Group2,color=Group2))+
  
  scale_fill_manual(values = met.brewer(name = "Egypt",n = 4)[c(4,1,2)])+
  scale_color_manual(values =met.brewer(name = "Egypt",n = 4)[c(4,1,2)])+
  
  
  scale_y_continuous(limits = c(-0.0021,0.0021),breaks = seq(-0.002,0.002,0.0010))+
  
  ylab(label = "Standardized mediation effect (95% CI)")+
  #xlab(label = "Child Psychopathology Scales")+
  scale_x_discrete(labels = function(label2){str_replace(string = label2,pattern ="-->",replacement = "\u2192" )})+
  
# 
  
  geom_hline(yintercept = 0,linetype="dashed",color="black",size=.75)+
  
 geom_hline(yintercept = c(-0.002,-0.0015,-0.0010,-0.0005,0.0005,0.0010,0.0015,0.002   ),
            linetype="dotted",color="light grey",size=.75)+
  
  
  
  geom_linerange( aes(xmin=label2, xmax=label2, ymin=L, ymax=U),
                      alpha=.5,
                  stat = "identity",position = position_dodge(width = .75),size=1.25,
                  show.legend = F)+
  
  ggtitle(label = "Longitudinal mediation",subtitle = "PCE \U2192 brain metrics \U2192 psychopathology")+
  geom_point(stat = "identity",color="black",
             position = position_dodge(width = .75),shape=21,size=2.5)+
  theme_classic()+
  coord_flip()+
  
  #annotate(geom = "text", x = 0, y = -.5, label = "change in psychopathology attributable to effects of PCE on the brain", size = 5) +
  
  geom_vline(xintercept = seq(1.5,12.5,1),color="darkgrey",size=.5)+
  # theme(axis.text.x = element_text(angle=90),
  #       legend.position = "top")+
  theme(axis.text = element_text(color="black"
                                 #,size = 10
  ),
  axis.title.y = element_blank(),plot.subtitle = element_text(size = 10),
  #  legend.text = element_text(size = 10),
  legend.title = element_blank(),
  #plot.margin = margin(5.5,5.5,5.5,5.5),
  #legend.position = c(0.45,1.1),
  legend.position = "top",legend.margin = margin(0,0,0,0),
  axis.text.x = element_text(angle = 0)#,
  #axis.title = element_text(size = 10)
  )+
  guides(fill = guide_legend(nrow = 2,reverse = T),color = guide_legend(nrow = 2,reverse = T)) 
  

fig3


ggsave(plot = fig3,filename = "C:/Users/dbara/Documents/WashU/ABCD/MRI/pce_mri_results_fig3.jpeg",
       dpi=500,width = 10,height = 3)
ggsave(plot = fig3,filename = "C:/Users/dbara/Documents/WashU/ABCD/MRI/pce_mri_results_fig3.tiff",
       dpi=500,width = 10,height = 3)
ggsave(plot = fig3,filename = "C:/Users/dbara/Documents/WashU/ABCD/MRI/pce_mri_results_fig3.png",
       dpi=500,width = 10,height = 3)

fig4 = ggarrange(fig2,fig3,nrow = 2,labels = c("A)","B)"),heights = c(1,1.65),align = "v")
ggsave(plot = fig4,filename = "C:/Users/dbara/Documents/WashU/ABCD/MRI/pce_mri_results_fig4.jpeg",
       dpi=500,width = 10,height = 6)
ggsave(plot = fig4,filename = "C:/Users/dbara/Documents/WashU/ABCD/MRI/pce_mri_results_fig4.pdf",
       dpi=500,width = 10,height = 6)
ggsave(plot = fig4,filename = "C:/Users/dbara/Documents/WashU/ABCD/MRI/pce_mri_results_fig4.png",
       dpi=500,width = 10,height = 6)






sapply(seq(0.01,0.5,0.01),FUN = function(X){pwr::pwr.r.test(r = X/3,sig.level = 0.05,power = .8)$n/pwr::pwr.r.test(r = X,sig.level = 0.05,power = .8)$n})



