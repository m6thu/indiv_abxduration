#########################################################################
#######Effect of antibiotic duration on hospitalised patients############
######################Show difference over time##########################
#########################################################################
library(reshape)
library(ggplot2)
library(ggpubr)

setwd("/Users/moyin/Documents/nBox/git_projects/indiv_abxduration/")
rm(list=ls()) # Clean working environment

source("los_abx_matrix.R")

plot.stability <- function (model) { # model can be "simple", "binary", or "frequency"
  
  source("default_params.R")
  source(paste0("model_", model,".R"))
  
  if(model == "simple"){
    
    timestep = 1
    sdDur=1
    iterations=100
    
    iter_totalR = matrix(NA, nrow = n.day, ncol = iterations)
    
    for(iter in 1:iterations){
      
      matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, max.los=max.los, 
                               p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                               meanDur= short_dur, timestep=timestep)
      patient.matrix=matrixes[[1]]
      abx.matrix=matrixes[[2]]
      los.array = summary.los(patient.matrix=patient.matrix)
      colo.matrix = colo.table(patient.matrix=patient.matrix, los=los.array, 
                               prop_R=prop_R,prop_S=prop_S)
      
      colo_table_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, 
                                       abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                       bif=bif, pi_ssr=pi_ssr, repop.s=repop.s, mu=mu, abx.s=abx.s, abx.r=abx.r,timestep=timestep)
      
      #Summary
      df = data.frame(colo_table_filled_iter)
      iter_totalR[, iter] = rowMeans(matrix(rowSums(df == "R")/n.bed, ncol=timestep, byrow = T))
    }
    cumsum_short = apply(iter_totalR, 2, cumsum)
    
    iter_totalR = matrix(NA, nrow = n.day, ncol = iterations)
    
    for(iter in 1:iterations){
      
      matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, max.los=max.los, 
                               p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                               meanDur= long_dur, timestep=timestep)
      patient.matrix=matrixes[[1]]
      abx.matrix=matrixes[[2]]
      los.array = summary.los(patient.matrix=patient.matrix)
      colo.matrix = colo.table(patient.matrix=patient.matrix, los=los.array, 
                               prop_R=prop_R,prop_S=prop_S)
      
      colo_table_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, 
                                       abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                       bif=bif, pi_ssr=pi_ssr, repop.s=repop.s, mu=mu, abx.s=abx.s, abx.r=abx.r,timestep=timestep)
      
      #Summary
      df = data.frame(colo_table_filled_iter)
      iter_totalR[, iter] = rowMeans(matrix(rowSums(df == "R")/n.bed, ncol=timestep, byrow = T))
    }
    cumsum_long = apply(iter_totalR, 2, cumsum)
    
    stabilitydata = data.frame(x=1:(n.day),y=(cumsum_long-cumsum_short)/c(1:n.day))
    colnames(stabilitydata) = c('x', 1:(ncol(stabilitydata)-1))
    stabilitydata.melt=melt(stabilitydata,id.vars='x', variable_name='iter')
    
    (stability.p=ggplot(stabilitydata.melt, aes(x=x, y=value, colour=as.factor(iter)))+
        geom_line(size=0.4, alpha=0.4)+
        scale_color_manual(values=rep('grey50', (ncol(stabilitydata)-1)))+
        labs(y='Cumulative difference in number\nof R carriers in long vs short\n duration/bed/day',
             x='Time (days)')+
        geom_vline(xintercept=150, linetype='dashed', size=0.7, colour='red')+
        theme_bw()+
        theme(legend.position = 'none', 
              axis.text = element_text(size = 15, colour = "grey50"), 
              axis.title = element_text(size=15)))
    
  } else if (model == "binary") {
    
    timestep = 1
    sdDur=1
    iterations=100
    
    iter_totalR = matrix(NA, nrow = n.day, ncol = iterations)
    
    for(iter in 1:iterations){
      
      matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, max.los=max.los, 
                               p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                               meanDur= short_dur, timestep=timestep)
      patient.matrix=matrixes[[1]]
      abx.matrix=matrixes[[2]]
      los.array = summary.los(patient.matrix=patient.matrix)
      colo.matrix = colo.table(patient.matrix=patient.matrix, los=los.array, 
                               prop_R=prop_R, prop_S=prop_S, prop_Sr=prop_Sr, prop_r=prop_r)
      colo_table_filled_iter = nextDay(patient.matrix=patient.matrix, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                       pi_ssr=pi_ssr, bif=bif, mu=mu,
                                       repop.s=repop.s, repop.r=repop.r, abx.r=abx.r, abx.s=abx.s, timestep=timestep)
      
      #Summary
      df = data.frame(colo_table_filled_iter)
      iter_totalR[, iter] = rowMeans(matrix(rowSums(df == "sR")/n.bed, ncol=timestep, byrow = T))
    }
    cumsum_short = apply(iter_totalR, 2, cumsum)
    
    iter_totalR = matrix(NA, nrow = n.day, ncol = iterations)
    
    for(iter in 1:iterations){
      
      matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, max.los=max.los, 
                               p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                               meanDur= long_dur, timestep=timestep)
      patient.matrix=matrixes[[1]]
      abx.matrix=matrixes[[2]]
      los.array = summary.los(patient.matrix=patient.matrix)
      colo.matrix = colo.table(patient.matrix=patient.matrix, los=los.array, 
                               prop_R=prop_R, prop_S=prop_S, prop_r=prop_r, prop_Sr=prop_Sr)
      colo_table_filled_iter = nextDay(patient.matrix=patient.matrix, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                       pi_ssr=pi_ssr, bif=bif, mu=mu,
                                       repop.s=repop.s, repop.r=repop.r, abx.r=abx.r, abx.s=abx.s, timestep=timestep)
      
      #Summary
      df = data.frame(colo_table_filled_iter)
      iter_totalR[, iter] = rowMeans(matrix(rowSums(df == "sR")/n.bed, ncol=timestep, byrow = T)) #averaged over timesteps
    }
    cumsum_long = apply(iter_totalR, 2, cumsum)
    
    stabilitydata = data.frame(x=1:(n.day),y=(cumsum_long-cumsum_short)/c(1:n.day))
    colnames(stabilitydata) = c('x', 1:(ncol(stabilitydata)-1))
    stabilitydata.melt=melt(stabilitydata,id.vars='x', variable_name='iter')
    
    (stability.p=ggplot(stabilitydata.melt, aes(x=x, y=value, colour=as.factor(iter)))+
        geom_line(size=0.4, alpha=0.4)+
        scale_color_manual(values=rep('grey50', (ncol(stabilitydata)-1)))+
        labs(y='Cumulative difference in number\nof R carriers in long vs short\n duration/bed/day',
             x='Time (days)')+
        geom_vline(xintercept=150, linetype='dashed', size=0.7, colour='red')+
        theme_bw()+
        theme(legend.position = 'none', 
              axis.text = element_text(size = 15, colour = "grey50"), 
              axis.title = element_text(size=15)))
    
  } else { #frequency 
    
    iterations = 100
    timestep=1
    
    iter_totalR.thres = matrix(NA, nrow = n.day, ncol = iterations)
    
    for(iter in 1:iterations){
      
      matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, max.los=max.los, 
                               p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                               meanDur= short_dur, timestep=timestep)
      patient.matrix=matrixes[[1]]
      abx.matrix=matrixes[[2]]
      los.array = summary.los(patient.matrix=patient.matrix)
      colo.matrix = colo.table(patient.matrix=patient.matrix, los.array=los.array, total_prop=total_prop, prop_R=prop_R, r_mean=r_mean,r_thres=r_thres,K=K)
      
      colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                        pi_ssr=pi_ssr, total_prop = total_prop, K=K, r_mean=r_mean, r_growth=r_growth, r_thres=r_thres, s_growth=s_growth,
                                        abx.s=abx.s, abx.r=abx.r, timestep=timestep)
      # Summary
      df.R = data.frame(colo.matrix_filled_iter[[2]])
      r_thres_matrix=data.frame(colo.matrix[[4]])
      #for number of people who reached R threshold on a day
      iter_totalR.thres[, iter]= rowMeans(matrix(rowSums(df.R >= r_thres_matrix)/n.bed, ncol=timestep, byrow = T))
    }
    cumsum_short = apply(iter_totalR.thres, 2, cumsum)
    
    iter_totalR.thres = matrix(NA, nrow = n.day, ncol = iterations)
    
    for(iter in 1:iterations){
      
      matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, max.los=max.los, 
                               p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                               meanDur= long_dur, timestep=timestep)
      patient.matrix=matrixes[[1]]
      abx.matrix=matrixes[[2]]
      los.array = summary.los(patient.matrix=patient.matrix)
      colo.matrix = colo.table(patient.matrix=patient.matrix, los.array=los.array, total_prop=total_prop, prop_R=prop_R,r_mean=r_mean,r_thres=r_thres,K=K)
      
      colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                        pi_ssr=pi_ssr, total_prop = total_prop, K=K, r_mean=r_mean, r_growth=r_growth,r_thres=r_thres, s_growth=s_growth,
                                        abx.s=abx.s, abx.r=abx.r, timestep=timestep)
      
      # Summary
      df.R = data.frame(colo.matrix_filled_iter[[2]])
      r_thres_matrix= data.frame(colo.matrix[[4]])
      #for number of people who reached R threshold on a day
      iter_totalR.thres[, iter]= rowMeans(matrix(rowSums(df.R >= r_thres_matrix)/n.bed, ncol=timestep, byrow = T))
      
    }
    cumsum_long = apply(iter_totalR.thres, 2, cumsum)
    
    stabilitydata = data.frame(x=1:(n.day),y=(cumsum_long-cumsum_short)/c(1:n.day))
    colnames(stabilitydata) = c('x', 1:(ncol(stabilitydata)-1))
    stabilitydata.melt=melt(stabilitydata,id.vars='x', variable_name='iter')
    
    stability.p=ggplot(stabilitydata.melt, aes(x=x, y=value, colour=as.factor(iter)))+
        geom_line(size=0.4, alpha=0.4)+
        scale_color_manual(values=rep('grey50', (ncol(stabilitydata)-1)))+
        labs(y='Cumulative difference in number\nof R carriers in long vs short\n duration/bed/day',
             x='Time (days)')+
        geom_vline(xintercept=150, linetype='dashed', size=0.7, colour='red')+
        theme_bw()+
        theme(legend.position = 'none', 
              axis.text = element_text(size = 15, colour = "grey50"), 
              axis.title = element_text(size=15))
  }
  return(stability.p)
}

abx.r = 0.0000001
abx.s = 0.2
simple0=plot.stability(model='simple')
binary0=plot.stability(model='binary')
abx.s=12
freq0 = plot.stability(model='frequency')

abx.s= 0.2
abx.r = 0.2
simplenot0=plot.stability(model='simple')
binarynot0=plot.stability(model='binary')
abx.s = 12
abx.r = 12
freqnot0=plot.stability(model='frequency')

figure=ggarrange(
  ggarrange(simplenot0,simple0, ncol = 2), 
  ggarrange(binarynot0,binary0, ncol = 2), 
  ggarrange(freqnot0,freq0, ncol = 2), 
  nrow=3)
annotate_figure(figure,
                top = text_grob("             Scenario A                                                                                                 Scenario B"),
                left = text_grob("    Model 3                                                                                Model 2                                                                                Model 1   ", rot = 90))

#save as png 1000 by 1200
setwd('/Users/moyin/Desktop/')
png(filename="stability.png", width = 1000, height = 1200)
annotate_figure(figure,
                top = text_grob("             Scenario A                                                                                                                                    Scenario B"),
                left = text_grob("    Model 1                                                                                                                Model 2                                                                                                           Model 3   ", rot = 90))
dev.off() 
