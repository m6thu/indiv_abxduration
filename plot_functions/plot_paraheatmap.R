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

parameters = as.factor(unique(c(parameters_diff_prevalence_simple, parameters_diff_prevalence_binary, parameters_diff_prevalence_freq)))

levels(parameters) =  rev(c("n.bed", "max.los", #ward level
                            "p.infect", "p.r.day1", "cum.r.1", #Number of antibiotic prescriptions 
                            "r_thres", "r_trans","pi_ssr","bif", #Transmission of resistant Enterobacteriaceae 
                            "prop_S", "prop_Sr","prop_r","prop_R", "total_prop", "K", #Carriage status on day one of admission to the ward 
                            "repop.s","s_growth", #Regrowth of susceptible Enterobacteriaceae 
                            "repop.r", "r_growth",  #Regrowth of resistant Enterobacteriaceae 
                            "mu",#decolonisation 
                            "abx.s", "abx.r", #antibiotics killing
                            "short_dur", "long_dur")) 

labs.df = data.frame(parameters = parameters, values=NA)
models = c('3GCRE\n Simple 3-state\n model', '3GCRE\n Co-carriage 5-state\n model', '3GCRE\n Population growth\n model', 
           'CRE\n Simple 3-state\n model', 'CRE\n Co-carriage 5-state\n model', 'CRE\n Population growth\n model', 
           '')

getposition <- function(data, labs = labs.df){
  
  prcc = data$prcc[[3]]$PRCC
  prcc$ranking = NA
  prcc$ranking[which(prcc$original<0)] = scales::rescale(prcc$original[which(prcc$original<0)], 
                                                         to=c(0,0.5))
  prcc$ranking[which(prcc$original>0)] = scales::rescale(prcc$original[prcc$original>0], 
                                                         to=c(0.5,1))
  prcc$ranking[which(prcc$`min. c.i.` < 0 & prcc$`max. c.i.`> 0)] = 0.5 #those with CI crossing 0 given 0.5
  
  prcc = prcc[order(prcc$ranking),] #arrange in ranking order
  # low = which(prcc$`min. c.i.` < 0 & prcc$`max. c.i.`< 0)
  # high= which(prcc$`min. c.i.` > 0 & prcc$`max. c.i.`> 0)
  # none= setdiff(1:nrow(prcc),c(low,high))
  # prcc$ranking[low]=prcc$`min. c.i.`[low]
  # prcc$ranking[high]=prcc$`max. c.i.`[high]
  # prcc$ranking[none]=0
  
  df = data.frame(parameters = rownames(prcc), ranking = prcc$ranking)
  out = merge(labs.df, df, by='parameters', all.x=TRUE)[,c(1,3)]
  
  return(out)
}

Amodel1.p = getposition(Amodel1)
Amodel2.p = getposition(Amodel2)
Amodel3.p = getposition(Amodel3)
Bmodel1.p = getposition(Bmodel1)
Bmodel2.p = getposition(Bmodel2)
Bmodel3.p = getposition(Bmodel3)
forlabels = cbind.data.frame(parameters = Bmodel3.p[,1], ranking=NA)

##prepare data
forplot = rbind.data.frame(Amodel1.p, Amodel2.p, Amodel3.p, 
                           Bmodel1.p, Bmodel2.p, Bmodel3.p, 
                           forlabels)
forplot$model = rep(models, each = length(parameters))
forplot$model = factor(forplot$model, levels = c('', 
                                                 '3GCRE\n Simple 3-state\n model', '3GCRE\n Co-carriage 5-state\n model', '3GCRE\n Population growth\n model', 
                                                 'CRE\n Simple 3-state\n model', 'CRE\n Co-carriage 5-state\n model', 'CRE\n Population growth\n model'))
forplot$ranking[grep('abx.r', forplot$parameters)[4:6]] = NA
forplot$ranking[grep('p.r.day1', forplot$parameters)[4:6]] = NA

#remove parameters that need not be shown 
forplot=forplot[!forplot$parameters=='short_dur',]
forplot=forplot[!forplot$parameters=='long_dur',]
forplot=forplot[!forplot$parameters=='prop_S',]
forplot=forplot[!forplot$parameters=='prop_Sr',]
forplot=forplot[!forplot$parameters=='prop_r',]

parameters=parameters[!parameters=='short_dur']
parameters=parameters[!parameters=='long_dur']
parameters=parameters[!parameters=='prop_S']
parameters=parameters[!parameters=='prop_Sr']
parameters=parameters[!parameters=='prop_r']
#parameters=as.character(parameters)


plot_paraheatmap <- function (forplot) { 
  
  base_size=9
  
  p = ggplot(forplot, aes(model, parameters)) + 
    geom_tile(aes(fill = ranking), colour = "white") + 
    scale_fill_gradientn(colours=c("#388697",'#fbfae5',"#EB5160"),
                         na.value = "white", 
                         breaks=c(min(forplot$ranking, na.rm = T), 0.5, max(forplot$ranking, na.rm = T)),
                         labels=c('Higher value decreases resistant\ncarriers with longer duration',
                                  'Does not affect difference in\nprevalence of resistant carriers', 
                                  'Higher value increases resistant\ncarriers with longer duration'))+
    theme_grey(base_size = base_size) + 
    labs(x = "", y="", fill = "")+ 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0),
                     labels=c('','Antibiotic killing', 
                              'Decolonisation', 
                              '','','','Enterobacteriaceae growth',
                              '','','Baseline carriage status', 
                              '','','','Resistance transmission',
                              '','','Antibiotic prescriptions',
                              '','Ward characteristics'))+
    theme(legend.position ='bottom', 
          legend.justification = c(0.6, 1),
          legend.key.width = unit(2,'cm'),
          axis.ticks = element_blank(), 
          axis.text.y = element_text(face='bold.italic', size=7),
          axis.text.x = element_text(size = base_size, angle = 330, hjust = 0, colour = "grey50"))+
    annotate(geom = 'text', x=1, y=c(1:length(parameters)), 
             label= rev(c("n.bed", "max.los", #ward level
                          "p.infect", "p.r.day1", "cum.r.1", 
                          "r_thres", "r_trans","pi_ssr","bif", #Transmission of resistant Enterobacteriaceae 
                          "prop_R", "total_prop", "K", #Carriage status on day one of admission to the ward 
                          "repop.s","s_growth", #Regrowth of susceptible Enterobacteriaceae 
                          "repop.r", "r_growth",  #Regrowth of resistant Enterobacteriaceae 
                          "mu",#decolonisation 
                          "abx.s", "abx.r")),#antibiotics killing
             colour='grey40', size=3)+
    guides(fill=guide_colorbar(nbin = 200, raster = F))+
    geom_hline(yintercept=c(2.5, 3.5, 7.5, 10.5, 14.5, 17.5), color='grey', size=0.5)+
    geom_vline(xintercept = 4.5,  color = "black", size=0.5) + 
    theme(plot.margin = unit(c(0,1.1,0,0), "cm"))
  
  return(p)
}






