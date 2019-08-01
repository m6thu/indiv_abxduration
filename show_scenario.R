#########################################################################
#######Effect of antibiotic duration on hospitalised patients############
#########################################################################
setwd('/Users/moyin/Desktop/indiv_abxduration')
rm(list=ls()) # Clean working environment

library(ggplot2)
library(ggpubr)
library(reshape)

# model can be "simple", "binary", or "frequency"
model <- "simple"

source("default_params.R")
source('los_abx_matrix.R')
source(paste0("model_", model,".R"))

dataformosaic<-function(data,label,n.bed=n.bed, n.day=n.day, timestep=timestep){
    
    data.t=t(data)
    rownames(data.t)=as.factor(1:n.bed)
    colnames(data.t)=as.factor(1:(n.day*timestep))
    data.melt=melt(data.t)
    colnames(data.melt)=c('Bed number','Time (days)',label)
    data.melt[,3]=as.factor(data.melt[,3])
    
    return(data.melt)
    
}

####################6. Visualisation #####################

if(model == "simple"){
    
    #Run model 
    timestep = 10
    n.day=350
    sdDur=1
    
    matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, 
                             p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                             meanDur= short_dur, timestep=timestep)
    patient.matrix.short=matrixes[[1]]
    abx.matrix.short=matrixes[[2]]
    los.array.short = summary.los(patient.matrix=patient.matrix.short)
    colo.matrix.short = colo.table(patient.matrix=patient.matrix.short, los=los.array.short, 
                                   prob_StartBact_R=prob_StartBact_R,prop_S_nonR=prop_S_nonR)
    
    colo_table_filled_short = nextDay(patient.matrix=patient.matrix.short, los.array=los.array.short, 
                                      abx.matrix=abx.matrix.short, colo.matrix=colo.matrix.short, 
                                      bif=bif, pi_ssr=pi_ssr, repop.s1=repop.s1, mu_r=mu_r, abx.s=abx.s, abx.r=abx.r,timestep=timestep)
    
    matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, 
                             p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                             meanDur= long_dur, timestep=timestep)
    patient.matrix.long=matrixes[[1]]
    abx.matrix.long=matrixes[[2]]
    los.array.long = summary.los(patient.matrix=patient.matrix.long)
    colo.matrix.long = colo.table(patient.matrix=patient.matrix.long, los=los.array.long, 
                                  prob_StartBact_R=prob_StartBact_R,prop_S_nonR=prop_S_nonR)
    
    colo_table_filled_long = nextDay(patient.matrix=patient.matrix.long, los.array=los.array.long, 
                                     abx.matrix=abx.matrix.long, colo.matrix=colo.matrix.long, 
                                     bif=bif, pi_ssr=pi_ssr, repop.s1=repop.s1, mu_r=mu_r, abx.s=abx.s, abx.r=abx.r,timestep=timestep)
    
    #Plots 
    ####Abx use plots 
    mosaicdata.abx.short = dataformosaic(data=abx.matrix.short ,label='Antibiotic type',n.bed=n.bed, n.day=n.day, timestep=timestep)
    mosaicdata.abx.long = dataformosaic(data=abx.matrix.long ,label='Antibiotic type',n.bed=n.bed, n.day=n.day, timestep=timestep)
    
    winteralmond=c('white',"#87C2BE","#5E8E7B")
    lgdcol=rgb(0.185, 0.188, 0.154, alpha = .05)
    base_size=9
    
    p.abx.short=ggplot(mosaicdata.abx.short, aes(`Time (days)`, `Bed number`)) + 
        geom_tile(aes(fill = `Antibiotic type`)) + 
        scale_fill_manual(values=winteralmond,labels = c("none", "kill S", "kill R"))+
        theme_bw(base_size = base_size) + 
        labs(title= 'Short antibiotic treatment duration',
             x = "Time (days)", y="Bed number")+ 
        scale_x_continuous(expand = c(0, 0), breaks= c(1000,2000,3000), labels = c('100','200','300')) +
        scale_y_continuous(expand = c(0, 0)) + 
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
              legend.position = "bottom", 
              legend.background = element_rect(fill=lgdcol,size=0.05),
              legend.title = element_text(size = base_size-2), 
              legend.text  = element_text(size = base_size-2),
              axis.text = element_text(size = base_size, colour = "grey50"), 
              axis.title = element_text(size=base_size))
    
    p.abx.long=ggplot(mosaicdata.abx.long, aes(`Time (days)`, `Bed number`)) + 
        geom_tile(aes(fill = `Antibiotic type`)) + 
        scale_fill_manual(values=winteralmond,labels = c("none", "kill S", "kill R"))+
        theme_bw(base_size = base_size) + 
        labs(title= 'Long antibiotic treatment duration',
             x = "Time (days)", y="")+ 
        scale_x_continuous(expand = c(0, 0), breaks= c(1000,2000,3000), labels = c('100','200','300')) +
        scale_y_continuous(expand = c(0, 0)) + 
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
              legend.position = "bottom", 
              legend.background = element_rect(fill=lgdcol,size=0.05),
              legend.title = element_text(size = base_size-2), 
              legend.text  = element_text(size = base_size-2),
              axis.text = element_text(size = base_size, colour = "grey50"), 
              axis.title = element_text(size=base_size))
    
    abx.mosaic=ggarrange(p.abx.short, p.abx.long, ncol=2, common.legend = T, legend = 'bottom')
    
    ##Carriage mosaic 
    mosaicdata.car.short = dataformosaic(data=colo_table_filled_short ,label='Carriage type',n.bed=n.bed, n.day=n.day, timestep=timestep)
    mosaicdata.car.long = dataformosaic(data= colo_table_filled_long ,label='Carriage type',n.bed=n.bed, n.day=n.day, timestep=timestep)
    
    sunflower=c('#F2A359',"#AAC0AF","#d6e1d9")
    
    p.car.short=ggplot(mosaicdata.car.short, aes(`Time (days)`, `Bed number`)) + 
        geom_tile(aes(fill = `Carriage type`)) + 
        scale_fill_manual(values=sunflower)+
        theme_bw(base_size = base_size) + 
        labs(x = "Time (days)", y="Bed number")+ 
        scale_x_continuous(expand = c(0, 0), breaks= c(1000,2000,3000), labels = c('100','200','300')) +
        scale_y_continuous(expand = c(0, 0)) + 
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
              legend.position = "bottom", 
              legend.background = element_rect(fill=lgdcol,size=0.05),
              legend.title = element_text(size = base_size-2), 
              legend.text  = element_text(size = base_size-2),
              axis.text = element_text(size = base_size, colour = "grey50"), 
              axis.title = element_text(size=base_size))
    
    p.car.long=ggplot(mosaicdata.car.long, aes(`Time (days)`, `Bed number`)) + 
        geom_tile(aes(fill = `Carriage type`)) + 
        scale_fill_manual(values=sunflower)+
        theme_bw(base_size = base_size) + 
        labs(x = "Time (days)", y="Bed number")+ 
        scale_x_continuous(expand = c(0, 0), breaks= c(1000,2000,3000), labels = c('100','200','300')) +
        scale_y_continuous(expand = c(0, 0)) + 
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
              legend.position = "bottom", 
              legend.background = element_rect(fill=lgdcol,size=0.05),
              legend.title = element_text(size = base_size-2), 
              legend.text  = element_text(size = base_size-2),
              axis.text = element_text(size = base_size, colour = "grey50"), 
              axis.title = element_text(size=base_size))
    
    car.mosaic=ggarrange(p.car.short, p.car.long, ncol=2, common.legend = T, legend = 'bottom')
    
    ##total R per day 
    totalRperday.short=apply(colo_table_filled_short, 1, function(x) length(which(x=='R')))
    totalRperday.short.avg= rowMeans(matrix(totalRperday.short, ncol=timestep, byrow=T))
    Rperdaydata.short=data.frame(`Time (days)`= 1:n.day, 
                                 `Total R per day`= totalRperday.short.avg)
    Rperdayline.short=ggplot(Rperdaydata.short, aes(x=Time..days., y=Total.R.per.day))+
        geom_point(colour= 'grey50', size=0.1)+ 
        labs(x = "Time (days)", y="Total R per day")+
        ylim(0,n.bed)+
        theme_bw()+
        theme(axis.text = element_text(size = base_size, colour = "grey50"), 
              axis.title = element_text(size=base_size))
    
    totalRperday.long=apply(colo_table_filled_long, 1, function(x) length(which(x=='R')))
    totalRperday.long.avg= rowMeans(matrix(totalRperday.long, ncol=timestep, byrow=T))
    Rperdaydata.long=data.frame(`Time (days)`= 1:n.day, 
                                 `Total R per day`= totalRperday.long.avg)
    Rperdayline.long=ggplot(Rperdaydata.long, aes(x=Time..days., y=Total.R.per.day))+
        geom_point(colour= 'grey50', size=0.1)+ 
        labs(x = "Time (days)", y="")+
        ylim(0,n.bed)+
        theme_bw()+
        theme(axis.text = element_text(size = base_size, colour = "grey50"), 
              axis.title = element_text(size=base_size))
    
    totalRline=ggarrange(Rperdayline.short, Rperdayline.long, ncol=2, common.legend = T, legend = 'bottom')
    
    ##cummulative R 
    cumsum.short=cumsum(totalRperday.short.avg)/n.bed/n.day*30*30*30
    cumsumdata.short=data.frame(`Time (days)`= 1:n.day, 
                                cumsum= cumsum.short)
    cumsumdata.short$dur=rep('Short',nrow(cumsumdata.short))
    cumsum.long=cumsum(totalRperday.long.avg)/n.bed/n.day*30*30
    cumsumdata.long=data.frame(`Time (days)`= 1:n.day, 
                                cumsum= cumsum.long)
    cumsumdata.long$dur=rep('Long',nrow(cumsumdata.long))
    cumsumdata=rbind.data.frame(cumsumdata.short,cumsumdata.long)
    sumplot=ggplot(cumsumdata, aes(x=Time..days., y=cumsum, colour = dur))+
        geom_line()+
        labs(x = "Time (days)", y="Cummulative sum of\nR/30-bed ward/month")+
        scale_color_manual(values =  c('#A0495B', '#0096BC'))+
        theme_bw()+
        theme(axis.text = element_text(size = base_size, colour = "grey50"), 
              axis.title = element_text(size=base_size), 
              legend.position = "bottom", 
              legend.background = element_rect(fill=lgdcol,size=0.05),
              legend.text  = element_text(size = base_size-2),
              legend.title = element_text(size = base_size-2,'Treatment duration'))
    
    (allplots= ggarrange(ggarrange(abx.mosaic, car.mosaic, totalRline, nrow=3), 
                         sumplot, nrow=2, heights = c(3, 1)))
    
}else if(model == "binary"){
    
    #Run model 
    timestep = 10
    n.day=350
    sdDur=1
    
    matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, 
                             p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                             meanDur= short_dur, timestep=timestep)
    patient.matrix.short=matrixes[[1]]
    abx.matrix.short=matrixes[[2]]
    los.array.short = summary.los(patient.matrix=patient.matrix.short)
    colo.matrix.short = colo.table(patient.matrix=patient.matrix.short, los=los.array.short, 
                             prob_StartBact_R=prob_StartBact_R, prop_S_nonR=prop_S_nonR, prop_Sr_inR=prop_Sr_inR, prop_sr_inR=prop_sr_inR)
    
    colo_table_filled_short = nextDay(patient.matrix=patient.matrix.short, abx.matrix=abx.matrix.short, colo.matrix=colo.matrix.short, 
                                      pi_ssr=pi_ssr, bif=bif, mu1=mu1, mu2=mu2, mu_r=mu_r, repop.r1=repop.r1, repop.r2=repop.r2,
                                      repop.s1=repop.s1, repop.s2=repop.s2, abx.r=abx.r, abx.s=abx.s, timestep=timestep)
    
    matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, 
                             p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                             meanDur= long_dur, timestep=timestep)
    patient.matrix.long=matrixes[[1]]
    abx.matrix.long=matrixes[[2]]
    los.array.long = summary.los(patient.matrix=patient.matrix.long)
    colo.matrix.long = colo.table(patient.matrix=patient.matrix.long, los=los.array.long, 
                                   prob_StartBact_R=prob_StartBact_R, prop_S_nonR=prop_S_nonR, prop_Sr_inR=prop_Sr_inR, prop_sr_inR=prop_sr_inR)
    
    colo_table_filled_long = nextDay(patient.matrix=patient.matrix.long, abx.matrix=abx.matrix.long, colo.matrix=colo.matrix.long, 
                                       pi_ssr=pi_ssr, bif=bif, mu1=mu1, mu2=mu2, mu_r=mu_r, repop.r1=repop.r1, repop.r2=repop.r2,
                                       repop.s1=repop.s1, repop.s2=repop.s2, abx.r=abx.r, abx.s=abx.s, timestep=timestep)
    
    #Plots 
    ####Abx use plots 
    mosaicdata.abx.short = dataformosaic(data=abx.matrix.short ,label='Antibiotic type',n.bed=n.bed, n.day=n.day, timestep=timestep)
    mosaicdata.abx.long = dataformosaic(data=abx.matrix.long ,label='Antibiotic type',n.bed=n.bed, n.day=n.day, timestep=timestep)
    
    winteralmond=c('white',"#87C2BE","#5E8E7B")
    lgdcol=rgb(0.185, 0.188, 0.154, alpha = .05)
    base_size=9
    
    p.abx.short=ggplot(mosaicdata.abx.short, aes(`Time (days)`, `Bed number`)) + 
        geom_tile(aes(fill = `Antibiotic type`)) + 
        scale_fill_manual(values=winteralmond,labels = c("none", "kill S", "kill R"))+
        theme_bw(base_size = base_size) + 
        labs(title= 'Short antibiotic treatment duration',
             x = "Time (days)", y="Bed number")+ 
        scale_x_continuous(expand = c(0, 0), breaks= c(1000,2000,3000), labels = c('100','200','300')) +
        scale_y_continuous(expand = c(0, 0)) + 
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
              legend.position = "bottom", 
              legend.background = element_rect(fill=lgdcol,size=0.05),
              legend.title = element_text(size = base_size-2), 
              legend.text  = element_text(size = base_size-2),
              axis.text = element_text(size = base_size, colour = "grey50"), 
              axis.title = element_text(size=base_size))
    
    p.abx.long=ggplot(mosaicdata.abx.long, aes(`Time (days)`, `Bed number`)) + 
        geom_tile(aes(fill = `Antibiotic type`)) + 
        scale_fill_manual(values=winteralmond,labels = c("none", "kill S", "kill R"))+
        theme_bw(base_size = base_size) + 
        labs(title= 'Long antibiotic treatment duration',
             x = "Time (days)", y="")+ 
        scale_x_continuous(expand = c(0, 0), breaks= c(1000,2000,3000), labels = c('100','200','300')) +
        scale_y_continuous(expand = c(0, 0)) + 
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
              legend.position = "bottom", 
              legend.background = element_rect(fill=lgdcol,size=0.05),
              legend.title = element_text(size = base_size-2), 
              legend.text  = element_text(size = base_size-2),
              axis.text = element_text(size = base_size, colour = "grey50"), 
              axis.title = element_text(size=base_size))
    
    abx.mosaic=ggarrange(p.abx.short, p.abx.long, ncol=2, common.legend = T, legend = 'bottom')
    
    ##Carriage mosaic 
    mosaicdata.car.short = dataformosaic(data=colo_table_filled_short,label='Carriage type',n.bed=n.bed, n.day=n.day, timestep=timestep)
    mosaicdata.car.short$`Carriage type`=factor(mosaicdata.car.short$`Carriage type`, levels = c('sR', 'sr','Sr', 'ss','S'))
    mosaicdata.car.long = dataformosaic(data= colo_table_filled_long,label='Carriage type',n.bed=n.bed, n.day=n.day, timestep=timestep)
    mosaicdata.car.long$`Carriage type`=factor(mosaicdata.car.long$`Carriage type`, levels = c('sR', 'sr','Sr', 'ss','S'))
    
    sunflower=c('#ee7a12', "#F2A359", "#f8caa0",  "#d6e1d9", "#AAC0AF")
    
    p.car.short=ggplot(mosaicdata.car.short, aes(`Time (days)`, `Bed number`)) + 
        geom_tile(aes(fill = `Carriage type`)) + 
        scale_fill_manual(values=sunflower)+
        theme_bw(base_size = base_size) + 
        labs(x = "Time (days)", y="Bed number")+ 
        scale_x_continuous(expand = c(0, 0), breaks= c(1000,2000,3000), labels = c('100','200','300')) +
        scale_y_continuous(expand = c(0, 0)) + 
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
              legend.position = "bottom", 
              legend.background = element_rect(fill=lgdcol,size=0.05),
              legend.title = element_text(size = base_size-2), 
              legend.text  = element_text(size = base_size-2),
              axis.text = element_text(size = base_size, colour = "grey50"), 
              axis.title = element_text(size=base_size))
    
    p.car.long=ggplot(mosaicdata.car.long, aes(`Time (days)`, `Bed number`)) + 
        geom_tile(aes(fill = `Carriage type`)) + 
        scale_fill_manual(values=sunflower)+
        theme_bw(base_size = base_size) + 
        labs(x = "Time (days)", y="Bed number")+ 
        scale_x_continuous(expand = c(0, 0), breaks= c(1000,2000,3000), labels = c('100','200','300')) +
        scale_y_continuous(expand = c(0, 0)) + 
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
              legend.position = "bottom", 
              legend.background = element_rect(fill=lgdcol,size=0.05),
              legend.title = element_text(size = base_size-2), 
              legend.text  = element_text(size = base_size-2),
              axis.text = element_text(size = base_size, colour = "grey50"), 
              axis.title = element_text(size=base_size))
    
    car.mosaic=ggarrange(p.car.short, p.car.long, ncol=2, common.legend = T, legend = 'bottom')
    
    ##total R per day 
    totalRperday.short=apply(colo_table_filled_short, 1, function(x) length(which(x=='sR'|x=='sr'|x=='Sr')))
    totalRperday.short.avg= rowMeans(matrix(totalRperday.short, ncol=timestep, byrow=T))
    Rperdaydata.short=data.frame(`Time (days)`= 1:n.day, 
                                 `Total R per day`= totalRperday.short.avg)
    Rperdayline.short=ggplot(Rperdaydata.short, aes(x=Time..days., y=Total.R.per.day))+
        geom_point(colour= 'grey50', size=0.1)+ 
        labs(x = "Time (days)", y="Total R per day")+
        ylim(0,n.bed)+
        theme_bw()+
        theme(axis.text = element_text(size = base_size, colour = "grey50"), 
              axis.title = element_text(size=base_size))
    
    totalRperday.long=apply(colo_table_filled_long, 1, function(x) length(which(x=='sR'|x=='sr'|x=='Sr')))
    totalRperday.long.avg= rowMeans(matrix(totalRperday.long, ncol=timestep, byrow=T))
    Rperdaydata.long=data.frame(`Time (days)`= 1:n.day, 
                                `Total R per day`= totalRperday.long.avg)
    Rperdayline.long=ggplot(Rperdaydata.long, aes(x=Time..days., y=Total.R.per.day))+
        geom_point(colour= 'grey50', size=0.1)+ 
        labs(x = "Time (days)", y="")+
        ylim(0,n.bed)+
        theme_bw()+
        theme(axis.text = element_text(size = base_size, colour = "grey50"), 
              axis.title = element_text(size=base_size))
    
    totalRline=ggarrange(Rperdayline.short, Rperdayline.long, ncol=2, common.legend = T, legend = 'bottom')
    
    ##cummulative R 
    cumsum.short=cumsum(totalRperday.short.avg)/n.bed/n.day*30*30
    cumsumdata.short=data.frame(`Time (days)`= 1:n.day, 
                                cumsum= cumsum.short)
    cumsumdata.short$dur=rep('Short',nrow(cumsumdata.short))
    cumsum.long=cumsum(totalRperday.long.avg)/n.bed/n.day*30*30
    cumsumdata.long=data.frame(`Time (days)`= 1:n.day, 
                               cumsum= cumsum.long)
    cumsumdata.long$dur=rep('Long',nrow(cumsumdata.long))
    cumsumdata=rbind.data.frame(cumsumdata.short,cumsumdata.long)
    sumplot=ggplot(cumsumdata, aes(x=Time..days., y=cumsum, colour = dur))+
        geom_line()+
        labs(x = "Time (days)", y="Cummulative sum of\nR/30-bed ward/month")+
        scale_color_manual(values =  c('#A0495B', '#0096BC'))+
        theme_bw()+
        theme(axis.text = element_text(size = base_size, colour = "grey50"), 
              axis.title = element_text(size=base_size), 
              legend.position = "bottom", 
              legend.background = element_rect(fill=lgdcol,size=0.05),
              legend.text  = element_text(size = base_size-2),
              legend.title = element_text(size = base_size-2,'Treatment duration'))
    
    (allplots= ggarrange(ggarrange(abx.mosaic, car.mosaic, totalRline, nrow=3), 
                         sumplot, nrow=2, heights = c(3, 1)))
    
    
}else if(model == "frequency"){
    
    timestep=1
    
    matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, 
                             p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                             meanDur= short_dur, timestep=timestep)
    patient.matrix.short=matrixes[[1]]
    abx.matrix.short=matrixes[[2]]
    los.array.short = summary.los(patient.matrix=patient.matrix.short)
    colo.matrix.short = colo.table(patient.matrix=patient.matrix.short, los.array=los.array.short, total_prop=total_prop, r_prop=r_prop,K=K)
    colo_table_filled_short = nextDay(patient.matrix=patient.matrix.short, los.array=los.array.short, abx.matrix=abx.matrix.short, colo.matrix=colo.matrix.short, 
                                      pi_ssr=pi_ssr, K=K, r_thres=r_thres, r_growth=r_growth, r_trans=r_trans, s_growth=s_growth,
                                      abx.s=abx.s, abx.r=abx.r, timestep=timestep)[[2]]
    
    matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, 
                             p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                             meanDur= long_dur, timestep=timestep)
    patient.matrix.long=matrixes[[1]]
    abx.matrix.long=matrixes[[2]]
    los.array.long = summary.los(patient.matrix=patient.matrix.long)
    colo.matrix.long = colo.table(patient.matrix=patient.matrix.long, los.array=los.array.long, total_prop=total_prop, r_prop=r_prop,K=K)
    colo_table_filled_long = nextDay(patient.matrix=patient.matrix.long, los.array=los.array.long, abx.matrix=abx.matrix.long, colo.matrix=colo.matrix.long, 
                                       pi_ssr=pi_ssr, K=K, r_thres=r_thres, r_growth=r_growth, r_trans=r_trans, s_growth=s_growth,
                                       abx.s=abx.s, abx.r=abx.r, timestep=timestep)[[2]]
   
    #Plots 
    ####Abx use plots 
    mosaicdata.abx.short = dataformosaic(data=abx.matrix.short ,label='Antibiotic type',n.bed=n.bed, n.day=n.day, timestep=timestep)
    mosaicdata.abx.long = dataformosaic(data=abx.matrix.long ,label='Antibiotic type',n.bed=n.bed, n.day=n.day, timestep=timestep)
    
    winteralmond=c('white',"#87C2BE","#5E8E7B")
    lgdcol=rgb(0.185, 0.188, 0.154, alpha = .05)
    base_size=9
    
    p.abx.short=ggplot(mosaicdata.abx.short, aes(`Time (days)`, `Bed number`)) + 
        geom_tile(aes(fill = `Antibiotic type`)) + 
        scale_fill_manual(values=winteralmond,labels = c("none", "kill S", "kill R"))+
        theme_bw(base_size = base_size) + 
        labs(title= 'Short antibiotic treatment duration',
             x = "Time (days)", y="Bed number")+ 
        scale_x_continuous(expand = c(0, 0), breaks= c(100,200,300)) +
        scale_y_continuous(expand = c(0, 0)) + 
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
              legend.position = "bottom", 
              legend.background = element_rect(fill=lgdcol,size=0.05),
              legend.title = element_text(size = base_size-2), 
              legend.text  = element_text(size = base_size-2),
              axis.text = element_text(size = base_size, colour = "grey50"), 
              axis.title = element_text(size=base_size))
    
    p.abx.long=ggplot(mosaicdata.abx.long, aes(`Time (days)`, `Bed number`)) + 
        geom_tile(aes(fill = `Antibiotic type`)) + 
        scale_fill_manual(values=winteralmond,labels = c("none", "kill S", "kill R"))+
        theme_bw(base_size = base_size) + 
        labs(title= 'Long antibiotic treatment duration',
             x = "Time (days)", y="")+ 
        scale_x_continuous(expand = c(0, 0), breaks= c(100,200,300)) +
        scale_y_continuous(expand = c(0, 0)) + 
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
              legend.position = "bottom", 
              legend.background = element_rect(fill=lgdcol,size=0.05),
              legend.title = element_text(size = base_size-2), 
              legend.text  = element_text(size = base_size-2),
              axis.text = element_text(size = base_size, colour = "grey50"), 
              axis.title = element_text(size=base_size))
    
    abx.mosaic=ggarrange(p.abx.short, p.abx.long, ncol=2, common.legend = T, legend = 'bottom')
    
    ##Carriage mosaic 
    mosaicdata.car.short = dataformosaic(data=colo_table_filled_short,label='Number of R',n.bed=n.bed, n.day=n.day, timestep=timestep)
    mosaicdata.car.short$`Carriage status`=rep('Below R threshold for transmission', nrow(mosaicdata.car.short))
    mosaicdata.car.short$`Carriage status`[which(as.numeric(mosaicdata.car.short$`Number of R`)>=r_thres)]='Above R threshold for transmission'
    mosaicdata.car.long = dataformosaic(data=colo_table_filled_long,label='Number of R',n.bed=n.bed, n.day=n.day, timestep=timestep)
    mosaicdata.car.long$`Carriage status`=rep('Below R threshold for transmission', nrow(mosaicdata.car.long))
    mosaicdata.car.long$`Carriage status`[which(as.numeric(mosaicdata.car.long$`Number of R`)>=r_thres)]='Above R threshold for transmission'
    
    sunflower=c('#ee7a12', "#d6e1d9")
    
    p.car.short=ggplot(mosaicdata.car.short, aes(`Time (days)`, `Bed number`)) + 
        geom_tile(aes(fill = `Carriage status`)) + 
        scale_fill_manual(values=sunflower)+
        theme_bw(base_size = base_size) + 
        labs(x = "Time (days)", y="Bed number")+ 
        scale_x_continuous(expand = c(0, 0), breaks= c(100,200,300)) +
        scale_y_continuous(expand = c(0, 0)) + 
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
              legend.position = "bottom", 
              legend.background = element_rect(fill=lgdcol,size=0.05),
              legend.title = element_text(size = base_size-2), 
              legend.text  = element_text(size = base_size-2),
              axis.text = element_text(size = base_size, colour = "grey50"), 
              axis.title = element_text(size=base_size))
    
    p.car.long=ggplot(mosaicdata.car.long, aes(`Time (days)`, `Bed number`)) + 
        geom_tile(aes(fill = `Carriage status`)) + 
        scale_fill_manual(values=sunflower)+
        theme_bw(base_size = base_size) + 
        labs(x = "Time (days)", y="Bed number")+ 
        scale_x_continuous(expand = c(0, 0), breaks= c(100,200,300)) +
        scale_y_continuous(expand = c(0, 0)) + 
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
              legend.position = "bottom", 
              legend.background = element_rect(fill=lgdcol,size=0.05),
              legend.title = element_text(size = base_size-2), 
              legend.text  = element_text(size = base_size-2),
              axis.text = element_text(size = base_size, colour = "grey50"), 
              axis.title = element_text(size=base_size))
    
    car.mosaic=ggarrange(p.car.short, p.car.long, ncol=2, common.legend = T, legend = 'bottom')
    
    ##total R per day 
    totalRperday.short=apply(colo_table_filled_short, 1, function(x) length(which(x>=r_thres)))
    totalRperday.short.avg= rowMeans(matrix(totalRperday.short, ncol=timestep, byrow=T))
    Rperdaydata.short=data.frame(`Time (days)`= 1:n.day, 
                                 `Total R per day`= totalRperday.short.avg)
    Rperdayline.short=ggplot(Rperdaydata.short, aes(x=Time..days., y=Total.R.per.day))+
        geom_point(colour= 'grey50', size=0.1)+ 
        labs(x = "Time (days)", y="Total no. above R\n threshold per day")+
        ylim(0,n.bed)+
        theme_bw()+
        theme(axis.text = element_text(size = base_size, colour = "grey50"), 
              axis.title = element_text(size=base_size))
    
    totalRperday.long=apply(colo_table_filled_long, 1, function(x) length(which(x>=r_thres)))
    totalRperday.long.avg= rowMeans(matrix(totalRperday.long, ncol=timestep, byrow=T))
    Rperdaydata.long=data.frame(`Time (days)`= 1:n.day, 
                                `Total R per day`= totalRperday.long.avg)
    Rperdayline.long=ggplot(Rperdaydata.long, aes(x=Time..days., y=Total.R.per.day))+
        geom_point(colour= 'grey50', size=0.1)+ 
        labs(x = "Time (days)", y="")+
        ylim(0,n.bed)+
        theme_bw()+
        theme(axis.text = element_text(size = base_size, colour = "grey50"), 
              axis.title = element_text(size=base_size))
    
    totalRline=ggarrange(Rperdayline.short, Rperdayline.long, ncol=2, common.legend = T, legend = 'bottom')
    
    ##cummulative R 
    cumsum.short=cumsum(totalRperday.short.avg)/n.bed/n.day*30*30
    cumsumdata.short=data.frame(`Time (days)`= 1:n.day, 
                                cumsum= cumsum.short)
    cumsumdata.short$dur=rep('Short',nrow(cumsumdata.short))
    cumsum.long=cumsum(totalRperday.long.avg)/n.bed/n.day*30*30
    cumsumdata.long=data.frame(`Time (days)`= 1:n.day, 
                               cumsum= cumsum.long)
    cumsumdata.long$dur=rep('Long',nrow(cumsumdata.long))
    cumsumdata=rbind.data.frame(cumsumdata.short,cumsumdata.long)
    sumplot=ggplot(cumsumdata, aes(x=Time..days., y=cumsum, colour = dur))+
        geom_line()+
        labs(x = "Time (days)", y="Cummulative sum of\nR/30-bed ward/month")+
        scale_color_manual(values =  c('#A0495B', '#0096BC'))+
        theme_bw()+
        theme(axis.text = element_text(size = base_size, colour = "grey50"), 
              axis.title = element_text(size=base_size), 
              legend.position = "bottom", 
              legend.background = element_rect(fill=lgdcol,size=0.05),
              legend.text  = element_text(size = base_size-2),
              legend.title = element_text(size = base_size-2,'Treatment duration'))
    
    (allplots= ggarrange(ggarrange(abx.mosaic, car.mosaic, totalRline, nrow=3), 
                         sumplot, nrow=2, heights = c(3, 1)))

}

allplots
