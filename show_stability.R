setwd("/Users/moyin/Desktop/indiv_abxduration")
rm(list=ls()) # Clean working environment

# model can be "simple", "binary", or "frequency"
model <- "frequency"

source("default_params.R")
source("los_abx_matrix.R")
source(paste0("model_", model,".R"))

#########################################################################
######################Show difference over time##########################
#########################################################################

if(model == "simple"){
    
    timestep = 10
    n.day=1000
    sdDur=1
    iterations=100
    
    iter_totalR = matrix(NA, nrow = n.day*timestep, ncol = iterations)
    
    for(iter in 1:iterations){
        
        matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, 
                                 p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                                 meanDur= short_dur, timestep=timestep)
        patient.matrix=matrixes[[1]]
        abx.matrix=matrixes[[2]]
        los.array = summary.los(patient.matrix=patient.matrix)
        colo.matrix = colo.table(patient.matrix=patient.matrix, los=los.array, 
                                 prob_StartBact_R=prob_StartBact_R,prop_S_nonR=prop_S_nonR)
        
        colo_table_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, 
                                         abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                         bif=bif, pi_ssr=pi_ssr, repop.s1=repop.s1, mu_r=mu_r, abx.s=abx.s, abx.r=abx.r,timestep=timestep)
        
        #Summary
        df = data.frame(colo_table_filled_iter)
        iter_totalR[, iter] = rowSums(df == "R")    
    }
    cumsum_short = cumsum(rowMeans(iter_totalR[1:nrow(iter_totalR),])/iterations/n.bed)
    
    iter_totalR = matrix(NA, nrow = n.day*timestep, ncol = iterations)
    
    for(iter in 1:iterations){
        
        matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, 
                                 p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                                 meanDur= long_dur, timestep=timestep)
        patient.matrix=matrixes[[1]]
        abx.matrix=matrixes[[2]]
        los.array = summary.los(patient.matrix=patient.matrix)
        colo.matrix = colo.table(patient.matrix=patient.matrix, los=los.array, 
                                 prob_StartBact_R=prob_StartBact_R,prop_S_nonR=prop_S_nonR)
        
        colo_table_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, 
                                         abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                         bif=bif, pi_ssr=pi_ssr, repop.s1=repop.s1, mu_r=mu_r, abx.s=abx.s, abx.r=abx.r,timestep=timestep)
        
        #Summary
        df = data.frame(colo_table_filled_iter)
        iter_totalR[, iter] = rowSums(df == "R")    
    }
    cumsum_long = cumsum(rowMeans(iter_totalR[1:nrow(iter_totalR),])/iterations/n.bed)
    stabilitydata= data.frame(y=cumsum_long-cumsum_short,
                              x=1:(n.day*timestep))
    (stability.p=ggplot(stabilitydata, aes(x=x, y=y))+
        geom_point())
    
} else if (model == "Binary") {
    
    
    
} else {
    
    
    
}
