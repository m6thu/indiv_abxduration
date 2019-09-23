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

parameters=as.factor(unique(c(parameters_diff_prevalence_simple, parameters_diff_prevalence_binary, parameters_diff_prevalence_freq)))
levels(parameters) =  rev(c("n.bed", "max.los", #ward level
                            "prop_R", "prop_S_nonR","prop_Sr_inR","prop_sr_inR", "total_prop", "K", "r_mean", #Carriage status on day one of admission to the ward  
                            "repop.s","s_growth", #Regrowth of susceptible Enterobacteriaceae 
                            "repop.r", "r_growth",  #Regrowth of resistant Enterobacteriaceae 
                            "r_thres", "pi_ssr","bif", #Transmission of resistant Enterobacteriaceae 
                            "mu",#decolonisation 
                            "abx.s", "abx.r", #antibiotics killing
                            "p.infect", "p.r.day1", "cum.r.1", #Number of antibiotic prescriptions 
                            "short_dur", "long_dur")) #antibiotics duration 
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
  prcc= data$prcc[[3]]$PRCC
  prcc$ranking=NA
  prcc= prcc[order(prcc$original),] #arrange in ranking order
  low = which(prcc$`min. c.i.` < 0 & prcc$`max. c.i.`< 0)
  high= which(prcc$`min. c.i.` > 0 & prcc$`max. c.i.`> 0)
  none= setdiff(1:nrow(prcc),c(low,high))
  prcc$ranking[low]=prcc$`min. c.i.`[low]
  prcc$ranking[high]=prcc$`max. c.i.`[high]
  prcc$ranking[none]=0
  
  df=data.frame(parameters=rownames(prcc), ranking=prcc$ranking)
  out=merge(labs.df, df, by='parameters', all.x=TRUE)[,c(1,3)]
  
  return(out)
}

Amodel1=get(load('runs/LHSdiff_simple_1100_notzero21Sep2019_1937BST.Rdata')) #abx_r>0
Amodel2=get(load('runs/LHSdiff_binary_1300_notzero_23Sep2019_0620BST.Rdata'))
Amodel3=get(load('runs/LHSdiff_frequency_1100_notzero_22Sep2019_1423BST.Rdata'))
Bmodel1=get(load('runs/LHSdiff_simple_1100_zero21Sep2019_2140BST.Rdata')) #abx_r>0
Bmodel2=get(load('runs/LHSdiff_binary_1300_zero_23Sep2019_0028BST.Rdata')) 
Bmodel3=get(load('runs/LHSdiff_frequency_1100_zero_22Sep2019_1646BST.Rdata'))

#Amodel1.p=getp(Amodel1, para.list=parameters_simple, 500)
# Amodel2.p=getp(Amodel2, para.list=parameters_binary, 500)
# Amodel3.p=getp(Amodel3, para.list=parameters_freq, 800)
# Bmodel1.p=getp(Bmodel1, para.list=parameters_simple, 500)
# Bmodel2.p=getp(Bmodel2, para.list=parameters_binary, 500)
# Bmodel3.p=getp(Bmodel3, para.list=parameters_freq, 500)

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
forplot$`ranking value`= scales::rescale(forplot$ranking, to=c(0,1))

base_size=9
ggplot(forplot, aes(model, parameters)) + 
  geom_tile(aes(fill = ranking), colour = "white") + 
  scale_fill_gradientn(colours=c("#388697",'#fbfae5',"#EB5160"),
                       na.value = "white", 
                       breaks=c(min(forplot$ranking, na.rm = T),0,max(forplot$ranking, na.rm = T)),
                       labels=c('Higher value decreases difference in\nprevalence of resistant carriers',
                                'Does not affect difference in\nprevalence of resistant carriers', 
                                'Higher value increases difference in\nprevalence of resistant carriers'))+
  theme_grey(base_size = base_size) + 
  labs(x = "", y="", fill = "")+ 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0),
                   labels=c("",'Duration of antibiotics', 
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
  geom_hline(yintercept=c(2.5, 5.5, 7.5, 8.5, 11.5, 15.5, 22.5), color='grey', size=0.5)+
  geom_vline(xintercept = 4.5,  color = "black", size=0.5)

