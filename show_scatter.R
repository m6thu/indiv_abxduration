########Show scatter of different parameters#############
setwd('/Users/moyin/Documents/git_projects/indiv_abxduration/')
rm(list=ls()) # Clean working environment

library(ggplot2)
library(ggpubr)

sd=get(load('runs/LHSdiff_simple_1100_zero21Sep2019_2140BST.Rdata'))
bd=get(load('runs/LHSdiff_binary_1300_zero_23Sep2019_0028BST.Rdata'))
fd=get(load('runs/LHSdiff_frequency_2100_zero_24Sep2019_0845BST.Rdata'))

scaleFUN <- function(x) sprintf("%.2f", x) #keep x axis decimal 2 places

scatter=function(lhs.data){
  sample=lhs.data$data
  y=lhs.data$res[,3,1]
  
  plotlist=list()
  for (i in 1:ncol(sample)){
    x=sample[,i]
    plotdata=data.frame(x=x, y=y)
    
    if (any(x>6)){
    plotlist[[i]]=ggplot(plotdata, aes(x,y))+
      geom_point(alpha=0.05, size=0.1)+
      geom_smooth(method='lm',formula=y~x, se=F)+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      labs(y='Model output', x= colnames(sample)[i])
    } else {
      plotlist[[i]]=ggplot(plotdata, aes(x,y))+
        geom_point(alpha=0.05, size=0.1)+
        geom_smooth(method='lm',formula=y~x, se=F)+
        scale_x_continuous(labels=scaleFUN)+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
        labs(y='Model output', x= colnames(sample)[i])
    }
  }
  no.col=round((ncol(sample))^0.5)
  
  return(ggarrange(plotlist = plotlist, ncol=no.col, nrow=no.col))
}

scatter(sd)
scatter(bd)
scatter(fd)
