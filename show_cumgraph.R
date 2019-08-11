##############################################################################
###########Effect of antibiotic duration on hospitalised patients#############
##show percent of the time duration matters in the parameter space explored###
##############################################################################

setwd("/Users/moyin/Desktop/indiv_abxduration")
rm(list=ls()) # Clean working environment

library(pse)
library(plyr)
library(data.table)

#load data from runs 
d.simple=get(load('runs/LHS2_simple_60005Aug2019_2050GMT.Rdata'))
d.binary=get(load('runs/LHS2_binary_60006Aug2019_0414GMT.Rdata'))
d.freq=get(load('runs/LHS_frequency_100006Aug2019_1354GMT.Rdata'))

df1= data.frame(combi=c(1:length(d.simple$res)),variable=rep('Model 1', length(d.simple$res)) ,value=d.simple$res)
df2= data.frame(combi=c(1:length(d.binary$res)),variable=rep('Model 2', length(d.binary$res)) ,value=d.binary$res[,1,])
df3= data.frame(combi=c(1:length(d.freq$res)),variable=rep('Model 3', length(d.freq$res)) ,value=d.freq$res[,2,])
df.melt= rbind.data.frame(df1,df2,df3)
ecd= ddply(df.melt, .(variable), summarize,
           value = value,
           ecdf = ecdf(value)(value))

ggplot(ecd,aes(value, ecdf, color = variable)) + 
    geom_step(size=1)+
    scale_colour_manual(name='Model type',
                        values=c('#5CA4A9','#ED6A5A','#f5a94f'), labels=c('Model 1', 'Model 2','Model 3'))+
    annotate("label", x = -0.1, y = 1.05, size = 6, fill='#9AADBF', colour='white', label='Less R carriers with longer duration of antibiotics')+
    annotate("label", x = 0.1, y = 1.05, size = 6, fill= '#9AADBF', colour='white', label='More R carriers with longer duration of antibiotics')+
    geom_segment(aes(x=0, xend=min(ecd$value), y=1.01, yend=1.01), colour='grey50',
                 arrow = arrow(length = unit(0.2, "cm")))+
    geom_segment(aes(x=0, xend=max(ecd$value), y=1.01, yend=1.01), colour='grey50',
                 arrow = arrow(length = unit(0.2, "cm")))+
    labs(y='Proportion of paramters explored', 
         x= 'Difference in number of R carriers in long duration vs short duration/bed/day')+
    theme_minimal()+
    theme(legend.position = 'bottom')

