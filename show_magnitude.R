########Show magnitude of difference over all parameter space#############
setwd('/Users/moyin/Desktop/indiv_abxduration/')
rm(list=ls()) # Clean working environment

library(ggplot2)
library(ggpubr)
require(parallel) # load parallel processing package to use multiple cores on computer (or cluster)

getplotdata=function(raw){
    
    x=raw$data$meanDur
    y=raw$res
    if (dim(raw$data)[2]==14) {
        type= rep('model1', length(y))
    } else if (dim(raw$data)[2]==21){
        type= rep('model2', length(y))
    } else if (dim(raw$data)[2]==16){
        type= rep('model3', length(y))
    } else {
        print('wrong dimensions indicated in the prcc data')
    }
    
    df=data.frame(x=x, y=y, type=type)
    
    return(df)
    
}

s= get(load('runs/LHS_simple_600_zero09Aug2019_2003GMT.Rdata'))
b= get(load('runs/LHS_binary_700_zero_09Aug2019_2221GMT.Rdata'))
f= get(load('runs/LHS_frequency_1300_zero_09Aug2019_2350GMT.Rdata'))

d=rbind.data.frame(getplotdata(s), getplotdata(b),getplotdata(f))

ggplot(d, aes(x,y, colour=type))+
    geom_smooth(method='lm',formula=y~x, se=F)+
    geom_point(alpha=0.1, size=0.7)+
    scale_colour_manual(name='Model type',
                        values=c('#5CA4A9','#ED6A5A','#f5a94f'), 
                        labels=c('Model 1', 'Model 2','Model 3'))+
    labs(x='Antibiotic treatment duration (days)', 
         y='Number of resistant Enterobacteriaceae carriers/bed/day')+
    theme_minimal()+
    theme(legend.position = 'bottom')
