library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)
library(viridisLite)

library(gridExtra)
library(cowplot)
library(grid)
library(ggpubr)



setwd('C:/Users/dbara/Documents/WashU/ABCD/MRI')

all_results = fread("all_PCE_regressions_8_pfdr.csv",data.table = F)

group_order = all_results %>% dplyr::select(c("order","label")) %>% unique()


for(i in 1:length(group_order$label)){
  
  temp = all_results %>% dplyr::filter(label == group_order$label[i])
  
  temp = temp[order(temp$Chisq_P_main),]
  temp$order = c(1:length(temp$Chisq_P_main))
  temp$bon_line = 0.05/dim(temp)[1]
  temp$Y = factor(temp$Y,levels = temp$Y )
  if(min(temp$P_fdr)<0.05){
    
    temp$fdr_line = mean(c(
      -log10(max(temp$Chisq_P_main[temp$P_fdr<0.05])),
      -log10(min(temp$Chisq_P_main[temp$P_fdr>0.05]))
    )
    )
    
    
  }else{
    p_fdr = temp$Chisq_P_main[1]
    
    while(min(p.adjust(c(p_fdr,temp$Chisq_P_main[-1]),method = "fdr"))>0.05){
      p_fdr = p_fdr - 0.00001
      
    }
    temp$fdr_line = -log10(p_fdr + 0.00001)
  }
  
  
  if(i ==1){all_results2 = temp}else{all_results2 = rbind(all_results2,temp)}
  
}

all_results2$label2 = factor(all_results2$label,levels = group_order$label)


p_fdr_all = mean(c(
  -log10(max(all_results2$Chisq_P_main[all_results2$P_fdr_all<0.05])),
  -log10(min(all_results2$Chisq_P_main[all_results2$P_fdr_all>0.05]))
)
)

all_results2$col = as.factor(all_results2$col) %>% as.numeric()


# levels(all_results2$col) =  c("rs-fmri","smri","dti fiber tract","rsi fiber tract","dti cortical wm","rsi cortical wm",
#                               "dti cortical gm","rsi cortical gm")

dat =dplyr::filter(.data = all_results2,label2 == unique(all_results2$label2)[40])

###########################################################

man_plot = function(dat){
  
  col.level = unique(dat$col)
  plot.title = unique(dat$label_sub)
  plot.title = stringr::str_replace(string =plot.title,pattern = "_newline_",replacement = "\n" )
  
  #gsub(pattern = "_newline_",replacement = "n" )
  
  plot1 = ggplot(data = dat,aes(x=order,y= -log10(Chisq_P_main)))+

    geom_point(shape=21,show.legend = T,alpha=.5,fill = turbo(n = 8)[col.level] )+
    
    theme_classic()+
    scale_y_continuous(limits = c(0,5))+
    geom_hline(aes(yintercept = p_fdr_all,linetype="p<0.05 fdr - all"),color="darkgrey")+
    geom_hline(aes(yintercept = fdr_line,linetype="p<0.05 fdr - by measure" ),color="darkgrey")+
    geom_hline(aes(yintercept = -log10(0.05),linetype="p<0.05 uncorrected"),color="darkgrey")+
    ggtitle(label = plot.title)+
      
    scale_linetype_manual(name = "p-value cutoff", 
                          values = c(1,2,3), 
                          guide = guide_legend(override.aes = 
                                                 list(color = c("darkgrey", "darkgrey","darkgrey"),shape=NA)))+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          axis.title = element_blank(),
          title = element_text(size=8),
          
          legend.position = "none")
  
    return(plot1)
    
}

##########################################################

plot.list = list()
for(i in 1: length(unique(all_results2$label2))){
  
  plot.list[[i]] = man_plot(dat =dplyr::filter(.data = all_results2,label2 == unique(all_results2$label2)[i])  )
}
#############################################################

plot2 = ggplot(data =dplyr::filter(.data = all_results2,label2 == unique(all_results2$label2)[1])  ,
               aes(x=order,y= -log10(Chisq_P_main)))+
  
  geom_point(shape=21,show.legend = T,alpha=.5,fill = turbo(n = 8)[1] )+
  
  theme_classic()+
  
  geom_hline(aes(yintercept = p_fdr_all,linetype="p<0.05 fdr - all"),color="darkgrey")+
  geom_hline(aes(yintercept = fdr_line,linetype="p<0.05 fdr - by measure" ),color="darkgrey")+
  geom_hline(aes(yintercept = -log10(0.05),linetype="p<0.05 uncorrected"),color="darkgrey")+
 # ggtitle(label = plot.title)+
  
  scale_linetype_manual(name = "p-value cutoff  ", 
                        values = c(1,2,3), 
                        guide = guide_legend(title.position = "top",keywidth = 4,nrow=3,
                                             override.aes = 
                                               list(color = c("darkgrey", "darkgrey","darkgrey"),size = 1.5,shape=NA)))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11),
        legend.position = "top")



leg <- get_legend(plot2 )

# Convert to a ggplot and print
leg2<-as_ggplot(leg)
leg2<-leg2+theme(plot.margin = margin(0,0,0,0))
leg2

###################################################


blank <- grid.rect(gp=gpar(col="white"))

A1 = grid.arrange(ggdraw() + draw_label(label = expression(paste(-log[10]," p-value")),fontfamily = "sans",angle=90,size=9),
                  plot.list[[1]],plot.list[[2]],
                  blank,blank,blank,blank,blank,ncol = 8,widths = c(1,6,6,6,6,6,6,6))
A2 = grid.arrange(  ggdraw() + draw_label(label = "Resting-state fMRI",fontfamily = "sans",x = .025,hjust = 0,size=10),
                        A1,nrow=2,  heights = c(1,6))


B1 = grid.arrange(ggdraw() + draw_label(label = expression(paste(-log[10]," p-value")),fontfamily = "sans",angle=90,size=9),
                  plot.list[[3]],plot.list[[4]],plot.list[[5]],plot.list[[6]],
                  blank,blank,blank,ncol = 8,widths = c(1,6,6,6,6,6,6,6))
B2 = grid.arrange(  ggdraw() + draw_label(label = "Structural MRI",fontfamily = "sans",x = .025,hjust = 0,size=10),
                    B1,nrow=2,  heights = c(1,6))


C1 = grid.arrange(ggdraw() + draw_label(label = expression(paste(-log[10]," p-value")),fontfamily = "sans",angle=90,size=9),
                  plot.list[[7]],plot.list[[8]],plot.list[[9]],plot.list[[10]],plot.list[[11]],
                  blank,blank,ncol = 8,widths = c(1,6,6,6,6,6,6,6))
C2 = grid.arrange(  ggdraw() + draw_label(label = "DTI - white matter tracts",fontfamily = "sans",x = .025,hjust = 0,size=10),
                    C1,nrow=2,  heights = c(1,6))

D1 = grid.arrange(ggdraw() + draw_label(label = expression(paste(-log[10]," p-value")),fontfamily = "sans",angle=90,size=9),
                  plot.list[[12]],plot.list[[13]],plot.list[[14]],plot.list[[15]],plot.list[[16]],plot.list[[17]],plot.list[[18]],
                  ncol = 8,widths = c(1,6,6,6,6,6,6,6))
D2 = grid.arrange(  ggdraw() + draw_label(label = "RSI - white matter tracts",fontfamily = "sans",x = .025,hjust = 0,size=10),
                    D1,nrow=2,  heights = c(1,6))

E1 = grid.arrange(ggdraw() + draw_label(label = expression(paste(-log[10]," p-value")),fontfamily = "sans",angle=90,size=9),
                  plot.list[[19]],plot.list[[20]],plot.list[[21]],plot.list[[22]],
                  blank,blank,blank,ncol = 8,widths = c(1,6,6,6,6,6,6,6))
E2 = grid.arrange(  ggdraw() + draw_label(label = "DTI - cortical white matter",fontfamily = "sans",x = .025,hjust = 0,size=10),
                    E1,nrow=2,  heights = c(1,6))

F1 = grid.arrange(ggdraw() + draw_label(label = expression(paste(-log[10]," p-value")),fontfamily = "sans",angle=90,size=9),
                  plot.list[[23]],plot.list[[24]],plot.list[[25]],plot.list[[26]],plot.list[[27]],plot.list[[28]],plot.list[[29]],
                  ncol = 8,widths = c(1,6,6,6,6,6,6,6))
F2 = grid.arrange(  ggdraw() + draw_label(label = "RSI - cortical white matter",fontfamily = "sans",x = .025,hjust = 0,size=10),
                    F1,nrow=2,  heights = c(1,6))

G1 = grid.arrange(ggdraw() + draw_label(label = expression(paste(-log[10]," p-value")),fontfamily = "sans",angle=90,size=9),
                  plot.list[[30]],plot.list[[31]],plot.list[[32]],plot.list[[33]],
                  blank,blank,blank,ncol = 8,widths = c(1,6,6,6,6,6,6,6))
G2 = grid.arrange(  ggdraw() + draw_label(label = "DTI - cortical grey matter",fontfamily = "sans",x = .025,hjust = 0,size=10),
                    G1,nrow=2,  heights = c(1,6))

H1 = grid.arrange(ggdraw() + draw_label(label = expression(paste(-log[10]," p-value")),fontfamily = "sans",angle=90,size=9),
                  plot.list[[34]],plot.list[[35]],plot.list[[36]],plot.list[[37]],plot.list[[38]],plot.list[[39]],plot.list[[40]],
                  ncol = 8,widths = c(1,6,6,6,6,6,6,6))
H2 = grid.arrange(  ggdraw() + draw_label(label = "RSI - cortical grey matter",fontfamily = "sans",x = .025,hjust = 0,size=10),
                    H1,nrow=2,  heights = c(1,6))

#####################################################################

all = grid.arrange(leg2,
                   A2,B2,C2,D2,E2,F2,G2,H2,
                   ggdraw() + draw_label(label = "Brain metrics - plotted by measure and region",fontfamily = "sans",size=12),
                   nrow = 10,
                   heights = c(5,6,6,6,6,6,6,6,6,1))

ggsave(plot = all,filename = "C:/Users/dbara/Documents/WashU/ABCD/MRI/pce_manhattan2.pdf",
       dpi=500,width = 11,height = 11)
