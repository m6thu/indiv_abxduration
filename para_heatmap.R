setwd('/Users/moyin/Desktop/indiv_abxduration')

####Plot a heatmap of the parameters against models 
#libraries
library(reshape)
library(plyr)
library(scales)
library(ggplot2)
library(epiR)

#get parameters
source('model_simple.R')
source('model_binary.R')
source('model_frequency.R')

parameters=as.factor(unique(c(parameters_simple, parameters_binary, parameters_freq)))
levels(parameters) =  rev(c("n.bed", "mean.max.los", #ward 
                            "prob_StartBact_R", "prop_S_nonR","prop_Sr_inR","prop_sr_inR","r_prop", "total_prop", "K", #patient characteristics 
                            "repop.s1", "repop.s2", "s_growth","repop.r1", "repop.r2", "r_growth", #within host dynamics 
                            "r_thres", "r_trans","pi_ssr","bif", #transmission dynamics
                            "mu1","mu2","mu_r",#decolonisation 
                            "abx.s", "abx.r","p.infect", "p.r.day1", "cum.r.1", "short_dur", "long_dur"))
labs.df=data.frame(parameters=parameters, values=NA)
models=c('Scenario A\n Model 1', 'Scenario A\n Model 2', 'Scenario A\n Model 3', 
         'Scenario B\n Model 1', 'Scenario B\n Model 2', 'Scenario B\n Model 3')

#get data

getp<-function(d, para.list, N, labs=labs.df){
    prcc.wo.pvalues=d$prcc[[1]]$PRCC
    dat = cbind.data.frame(d[['data']],unlist(d$res)[1:N])
    prcc= cbind.data.frame(para.list, pvalues= as.numeric(epi.prcc(dat = dat, sided.test = 2)$p.value))
    colnames(prcc)=c('parameters','p')
    out=merge(labs.df,prcc, by='parameters', all.x=TRUE)[,c(1,3)]
    return(out)
}

Amodel1=get(load('runs/LHS_simple_50014Jul2019_0747GMT.Rdata'))
Amodel2=get(load('runs/LHS_binary_50017Jul2019_0706GMT.Rdata'))
Amodel3=get(load('runs/LHS_frequency_80017Jul2019_1520GMT.Rdata'))
Bmodel1=get(load('runs/LHS_simple_50013Jul2019_0057GMT.Rdata'))
Bmodel2=get(load('runs/LHS_binary_50017Jul2019_0042GMT.Rdata')) 
Bmodel3=get(load('runs/LHS_frequency_50017Jul2019_1721GMT.Rdata'))

Amodel1.p=getp(Amodel1, para.list=parameters_simple, 500)
Amodel2.p=getp(Amodel2, para.list=parameters_binary, 500)
Amodel3.p=getp(Amodel3, para.list=parameters_freq, 800)
Bmodel1.p=getp(Bmodel1, para.list=parameters_simple, 500)
Bmodel2.p=getp(Bmodel2, para.list=parameters_binary, 500)
Bmodel3.p=getp(Bmodel3, para.list=parameters_freq, 500)

##prepare data
forplot=rbind.data.frame(Amodel1.p, Amodel2.p, Amodel3.p, 
                         Bmodel1.p, Bmodel2.p, Bmodel3.p)
forplot$model= rep(models, each=length(parameters))
forplot$p[grep('abx.r',forplot$parameters)[4:6]]=NA
forplot$`p value`=rescale(forplot$p, to=c(0,1))

base_size=9
ggplot(forplot, aes(model, parameters)) + 
    geom_tile(aes(fill = `p value`), colour = "white") + 
    scale_fill_gradientn(colours=c("#7A306C","#068D9D","#9BC1BC"),
                         values=c(0,0.001,0.05, 1), na.value = "white")+
    theme_grey(base_size = base_size) + 
    labs(x = "", y="")+ 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + 
    theme(legend.position = "right", 
          axis.ticks = element_blank(), 
          axis.text.x = element_text(size = base_size, angle = 330, hjust = 0, colour = "grey50"))+
        geom_vline(xintercept = 3.5,  color = "black", size=0.5)
