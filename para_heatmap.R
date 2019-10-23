setwd('/Users/moyin/Documents/nBox/git_projects/indiv_abxduration/')

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

parameters=as.factor(unique(c(parameters_prevalence_simple, parameters_prevalence_binary, parameters_prevalence_freq)))
levels(parameters) =  rev(c("n.bed", "max.los", #ward level
                            "prop_S", "prop_Sr","prop_r","prop_R", "total_prop", "K", #Carriage status on day one of admission to the ward  
                            "repop.s","s_growth", #Regrowth of susceptible Enterobacteriaceae 
                            "repop.r", "r_growth",  #Regrowth of resistant Enterobacteriaceae 
                            "r_thres", "r_trans","pi_ssr","bif", #Transmission of resistant Enterobacteriaceae 
                            "mu",#decolonisation 
                            "abx.s", "abx.r", #antibiotics killing
                            "p.infect", "p.r.day1", "cum.r.1", #Number of antibiotic prescriptions 
                            "meanDur")) #antibiotics duration 
labs.df=data.frame(parameters=parameters, values=NA)
models=c('Scenario A\n Model 1', 'Scenario A\n Model 2', 'Scenario A\n Model 3', 
         'Scenario B\n Model 1', 'Scenario B\n Model 2', 'Scenario B\n Model 3', 
         '')

#get data

# getp<-function(d, para.list, N, labs=labs.df){
#     prcc.wo.pvalues=d$prcc[[1]]$PRCC
#     dat = cbind.data.frame(d[['data']],unlist(d$res)[1:N])
#     prcc= cbind.data.frame(para.list, pvalues= as.numeric(epi.prcc(dat = dat, sided.test = 2)$p.value))
#     colnames(prcc)=c('parameters','p')
#     out=merge(labs.df,prcc, by='parameters', all.x=TRUE)[,c(1,3)]
#     return(out)
# }

getposition<-function(data,labs=labs.df){
  prcc= data$prcc[[1]]$PRCC
  #row.names(prcc)[4]='prop_S'
  prcc$ranking=NA
  prcc$ranking[which(prcc$original<0)] = scales::rescale(prcc$original[which(prcc$original<0)], 
                                                                          to=c(0,0.5))
  prcc$ranking[which(prcc$original>0)] = scales::rescale(prcc$original[prcc$original>0], 
                                                                          to=c(0.5,1))
  
  prcc= prcc[order(prcc$ranking),] #arrange in ranking order
  # low = which(prcc$`min. c.i.` < 0 & prcc$`max. c.i.`< 0)
  # high= which(prcc$`min. c.i.` > 0 & prcc$`max. c.i.`> 0)
  # none= setdiff(1:nrow(prcc),c(low,high))
  # prcc$ranking[low]=prcc$`min. c.i.`[low]
  # prcc$ranking[high]=prcc$`max. c.i.`[high]
  # prcc$ranking[none]=0
  
  df=data.frame(parameters=rownames(prcc), ranking=prcc$ranking)
  out=merge(labs.df, df, by='parameters', all.x=TRUE)[,c(1,3)]
  
  return(out)
}

Amodel1=get(load('runs/LHS_simple_900_notzero16Oct2019_1837BST.Rdata')) #abx_r>0
colnames(Amodel1$data)[4]='prop_S'
Amodel2=get(load('runs/LHS_binary_900_notzero_18Oct2019_1627BST.Rdata'))
Amodel3=get(load('runs/LHS_frequency_1300_notzero_17Oct2019_1157BST.Rdata'))
Amodel3=get(load('runs/LHS_frequency_500_notzero_22Oct2019_0735BST.Rdata'))
Bmodel1=get(load('runs/LHS_simple_900_zero16Oct2019_1651BST.Rdata')) #abx_r>0
colnames(Bmodel1$data)[4]='prop_S'
Bmodel2=get(load('runs/LHS_binary_900_zero_18Oct2019_1357BST.Rdata')) 
Bmodel3=get(load('runs/LHS_frequency_1500_zero_21Oct2019_1051BST.Rdata'))
Bmodel3=get(load('runs/LHS_frequency_501_zero_22Oct2019_1222BST.Rdata'))

Amodel1.p= getposition(Amodel1)
Amodel2.p= getposition(Amodel2)
Amodel3.p= getposition(Amodel3)
Bmodel1.p= getposition(Bmodel1)
Bmodel2.p= getposition(Bmodel2)
Bmodel3.p= getposition(Bmodel3)
forlabels= cbind.data.frame(parameters=Bmodel3.p[,1], ranking=NA)

##prepare data
forplot=rbind.data.frame(Amodel1.p, Amodel2.p, Amodel3.p, 
                         Bmodel1.p, Bmodel2.p, Bmodel3.p, 
                         forlabels)
forplot$model= rep(models, each=length(parameters))
forplot$ranking[grep('abx.r',forplot$parameters)[4:6]]=NA

base_size=9
ggplot(forplot, aes(model, parameters)) + 
  geom_tile(aes(fill = ranking), colour = "white") + 
  scale_fill_gradientn(colours=c("#388697",'#fbfae5',"#EB5160"),
                       na.value = "white", 
                       breaks=c(min(forplot$ranking, na.rm = T), 0.5, max(forplot$ranking, na.rm = T)),
                       labels=c('Higher value decreases\nprevalence of resistant carriers',
                                #'Higher value favours long duration',
                                'Does not affect \nprevalence of resistant carriers', 
                                #'Does not affect difference in\nprevalence of resistant carriers', 
                                'Higher value increases\nprevalence of resistant carriers'))+
                                #'Higher value favours short duration'))+
  theme_grey(base_size = base_size) + 
  labs(x = "", y="", fill = "")+ 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0),
                   labels=c('Duration of antibiotics', 
                            '','','Antibiotic prescriptions', 
                            '','Antibiotic killing',
                            'Decolonisation',
                            '','', 'Resistance transmission', 
                            '','', 
                            '','Enterobacteriaceae growth ',
                            '','','','','','','Baseline carriage status', 
                            '','Ward characteristics'))+
  theme(legend.position ='bottom', 
        legend.justification = c(0.6, 1),
        legend.key.width = unit(2,'cm'),
        axis.ticks = element_blank(), 
        axis.text.y = element_text(face='bold.italic', size=7),
        axis.text.x = element_text(size = base_size, angle = 330, hjust = 0, colour = "grey50"))+
  annotate(geom = 'text', x=1, y=c(1:length(parameters)), label=Amodel1.p[,1], colour='grey40', size=3)+
  guides(fill=guide_colorbar(nbin = 200, raster = F))+
  geom_hline(yintercept=c(1.5, 4.5, 6.5, 7.5, 11.5, 15.5, 21.5), color='grey', size=0.5)+
  geom_vline(xintercept = 4.5,  color = "black", size=0.5)

