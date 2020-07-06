##############################################################################
###########Effect of antibiotic duration on hospitalised patients#############
##show percent of the time duration matters in the parameter space explored###
##############################################################################

setwd("/Users/moyin/Documents/git_projects/indiv_abxduration/")
rm(list=ls()) # Clean working environment

library(pse)
library(plyr)
library(data.table)
library(ggplot2)
library(ggpubr)

#load data from runs 
dA.simple=get(load('runs/LHSdiff_simple_1100_notzero21Sep2019_1937BST.Rdata')) #abxr>0
dA.binary=get(load('runs/LHSdiff_binary_1300_notzero_23Sep2019_0620BST.Rdata'))
dA.freq=get(load('runs/LHSdiff_frequency_1500_notzero_23Sep2019_1241BST.Rdata'))
dB.simple=get(load('runs/LHSdiff_simple_1100_zero21Sep2019_2140BST.Rdata'))#abxr=0
dB.binary=get(load('runs/LHSdiff_binary_1300_zero_23Sep2019_0028BST.Rdata'))
dB.freq=get(load('runs/LHSdiff_frequency_2100_zero_24Sep2019_0845BST.Rdata'))

df1A= data.frame(combi=c(1:length(dA.simple$res)),variable=rep('Model 1A', length(dA.simple$res)) ,value=dA.simple$res[,3,1])
df1B= data.frame(combi=c(1:length(dB.simple$res)),variable=rep('Model 1B', length(dB.simple$res)) ,value=dB.simple$res[,3,1])
df1= rbind.data.frame(df1A,df1B)
ecd1= ddply(df1, .(variable), summarize,
            value = value,
            ecdf = ecdf(value)(value))

df2A= data.frame(combi=c(1:length(dA.binary$res)),variable=rep('Model 2A', length(dA.binary$res)) ,value=dA.binary$res[,3,1])
df2B= data.frame(combi=c(1:length(dB.binary$res)),variable=rep('Model 2B', length(dB.binary$res)) ,value=dB.binary$res[,3,1])
df2= rbind.data.frame(df2A,df2B)
ecd2= ddply(df2, .(variable), summarize,
            value = value,
            ecdf = ecdf(value)(value))

df3A= data.frame(combi=c(1:length(dA.freq$res)),variable=rep('Model 3A', length(dA.freq$res)) ,value=dA.freq$res[,3,1])
df3B= data.frame(combi=c(1:length(dB.freq$res)),variable=rep('Model 3B', length(dB.freq$res)) ,value=dB.freq$res[,3,1])
df3= rbind.data.frame(df3A,df3B)
ecd3= ddply(df3, .(variable), summarize,
            value = value,
            ecdf = ecdf(value)(value))

model1p=ggplot(ecd1,aes(value, ecdf, color = variable)) + 
  geom_step(size=0.5)+
  geom_vline(xintercept=0, linetype='dashed', color='#C03221', size=0.4)+
  scale_colour_manual(name='Model 1',
                      values=c('#3f7478','#8ec0c3'), 
                      labels=c('Scenario A', 'Scenario B'))+
  #annotate("label", x = -0.1, y = 1.05, size = 6, fill='#9AADBF', colour='white', label='Less R carriers with longer duration of antibiotics')+
  #annotate("label", x = 0.1, y = 1.05, size = 6, fill= '#9AADBF', colour='white', label='More R carriers with longer duration of antibiotics')+
  #geom_segment(aes(x=0, xend=min(ecd1$value), y=1.01, yend=1.01), colour='grey50',
  #           arrow = arrow(length = unit(0.2, "cm")))+
  #geom_segment(aes(x=0, xend=max(ecd1$value), y=1.01, yend=1.01), colour='grey50',
  #         arrow = arrow(length = unit(0.2, "cm")))+
  labs(y='Proportion of paramters explored', 
       x= '')+
  theme_minimal()+
  theme(legend.position = 'bottom')

model2p=ggplot(ecd2,aes(value, ecdf, color = variable)) + 
  geom_step(size=0.5)+
  geom_vline(xintercept=0, linetype='dashed', color='#C03221', size=0.4)+
  scale_colour_manual(name='Model 2',
                      values=c('#e8402c','#f5a89f'), labels=c('Scenario A', 'Scenario B'))+
  #annotate("label", x = -0.1, y = 1.05, size = 6, fill='#9AADBF', colour='white', label='Less R carriers with longer duration of antibiotics')+
  #annotate("label", x = 0.1, y = 1.05, size = 6, fill= '#9AADBF', colour='white', label='More R carriers with longer duration of antibiotics')+
  #geom_segment(aes(x=0, xend=min(ecd2$value), y=1.01, yend=1.01), colour='grey50',
  #    arrow = arrow(length = unit(0.2, "cm")))+
  #geom_segment(aes(x=0, xend=max(ecd2$value), y=1.01, yend=1.01), colour='grey50',
  #   arrow = arrow(length = unit(0.2, "cm")))+
  labs(y='', 
       x= 'Difference in number of R carriers in long duration vs short duration/bed/day')+
  theme_minimal()+
  theme(legend.position = 'bottom')

model3p=ggplot(ecd3,aes(value, ecdf, color = variable)) + 
  geom_step(size=0.5)+
  geom_vline(xintercept=0, linetype='dashed', color='#C03221', size=0.4)+
  scale_colour_manual(name='Model 3',
                      values=c('#ea850d','#f9bc97'), labels=c('Scenario A', 'Scenario B'))+
  #annotate("label", x = -0.1, y = 1.05, size = 6, fill='#9AADBF', colour='white', label='Less R carriers with longer duration of antibiotics')+
  #annotate("label", x = 0.1, y = 1.05, size = 6, fill= '#9AADBF', colour='white', label='More R carriers with longer duration of antibiotics')+
  #geom_segment(aes(x=0, xend=min(ecd3$value), y=1.01, yend=1.01), colour='grey50',
  #    arrow = arrow(length = unit(0.2, "cm")))+
  #geom_segment(aes(x=0, xend=max(ecd3$value), y=1.01, yend=1.01), colour='grey50',
  #    arrow = arrow(length = unit(0.2, "cm")))+
  labs(y='', 
       x= '')+
  theme_minimal()+
  theme(legend.position = 'bottom')

ggarrange(model1p,model2p,model3p, ncol=3, labels=c('A', 'B','C'))
