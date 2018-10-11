# Adapted from MDP_individualmodel_110918.R

rm(list=ls()) # Clean working environment

# model can be "simple", "binary", or "frequency"
model <- "frequency"

source("default_params.R")
source(paste0("model_", model,".R"))

#########################################################################
#######Effect of antibiotic duration on hospitalised patients############
#########################################################################

mean_short <- 4 #antibiotic short duration
mean_long <- 14 #antibiotic short duration


if(model == "binary"){
    ################################## Multiple runs comparing between both ####################################
    iterations <- 100
    
    abx.short <- list()
    abx.long <- list()
    array_LOS_short <- list()
    array_LOS_long <- list()
    array_StartBact_short <- list()
    array_StartBact_long <- list()
    colo_table_filled_short <- list()
    colo_table_filled_long <- list()
    
    short_totalR <- matrix(NA, nrow = n.days, ncol = iterations)
    long_totalR <- matrix(NA, nrow = n.days, ncol = iterations)
    
    for (i in 1:iterations){
        
        print(i)
        #Generate length of stay and antibiotic duration table
        abx.short[[i]] <- abx.table(n.bed, n.days, mean.max.los, p.s, p.r, meanDur=mean_short)
        abx.long[[i]] <- abx.table(n.bed, n.days, mean.max.los, p.s, p.r, meanDur=mean_long)
        
        #Generate baseline carriage status
        array_LOS_short[[i]] <- array_LOS_func(los_duration=abx.short[[i]][2])
        array_LOS_long[[i]] <- array_LOS_func(los_duration=abx.long[[i]][2])
        
        #Update values for every day
        array_StartBact_short[[i]] <- gen_StartBact(los=array_LOS_short[[i]], prob_StartBact)
        array_StartBact_long[[i]] <- gen_StartBact(los=array_LOS_long[[i]], prob_StartBact)
        
        #output
        colo_table_filled_short[[i]] <- nextDay(bed_table= abx.short[[i]][[1]], array_LOS=array_LOS_short[[i]], 
                                                treat_table=abx.short[[i]][[3]], colo_table=array_StartBact_short[[i]], 
                                                pi_r1=pi_r1, pi_r2= pi_r2, mu1=mu1, mu2=mu2, 
                                                abx.r=abx.r,abx.s=abx.s,
                                                repop.r1 = repop.r1, repop.r2 = repop.r2, repop.r3 = repop.r3, 
                                                repop.s1 = repop.s1, repop.s2 = repop.s2,repop.s3 = repop.s3)
        colo_table_filled_long[[i]] <- nextDay(bed_table= abx.long[[i]][[1]], array_LOS=array_LOS_long[[i]], 
                                               treat_table=abx.long[[i]][[3]], colo_table=array_StartBact_long[[i]], 
                                               pi_r1=pi_r1, pi_r2= pi_r2, mu1=mu1, mu2=mu2, 
                                               abx.r=abx.r,abx.s=abx.s,
                                               repop.r1 = repop.r1, repop.r2 = repop.r2, repop.r3 = repop.r3, 
                                               repop.s1 = repop.s1, repop.s2 = repop.s2,repop.s3 = repop.s3)
    
        #increase overtime
        df.short <- colo_table_filled_short[[i]]
        short_totalR[,i]<- rowSums(df.short == "sR")
        
        df.long <- colo_table_filled_long[[i]]
        long_totalR[,i] <- rowSums(df.long == "sR")
    
    }
    
    # Debug edge cases: 43, 71 from run 1431_28Sept2018.RData
    # Debug cases: 36, 61, 160, 163, 212, 240, 258, 275, 300, 557, 565, 648, 658, 692, 735, 783, 813, 897 from run 1501_28Sept2018.RData
    # A clean run: 1553_28Sept2018.RData, still uses lists in iteration collection instead of matrix that is easier to summarize from
    
    
    # Check between 2 simple scenarios
    par(mfrow=c(1,2))
    
    x <- 1:n.days
    y.short <- short_totalR[,1]/n.bed
    r.plot.short <- plot(x, y.short, type="l", ylim=c(0,1), main=paste("Short mean abx:", mean_short, 
                         "days", "(", round(mean(rowSums(short_totalR)/iterations/n.bed), digits=4), ")"),
                         xlab="Time", ylab="% patients carrying MDRO", col=rgb(0,0,0,alpha=0.1))
    for(i in 2:iterations){
        lines(x, short_totalR[,i]/n.bed, col=rgb(0,0,0,alpha=0.1))
    }
    lines(x, rowSums(short_totalR)/iterations/n.bed, col=2)
    abline(h = mean(rowSums(short_totalR)/iterations/n.bed), col=2, lwd=2)
    
    y.long <- long_totalR[,1]/n.bed
    r.plot.long <- plot(x, y.long, type="l", ylim=c(0,1), main=paste("Long mean abx:", mean_long, 
                        "days", "(", round(mean(rowSums(long_totalR)/iterations/n.bed), digits=4), ")"),
                        xlab="Time", ylab="% patients carrying MDRO", col=rgb(0,0,0,alpha=0.1))
    for(i in 2:iterations){
        lines(x, long_totalR[,i]/n.bed, col=rgb(0,0,0,alpha=0.1))
    }
    lines(x, rowSums(long_totalR)/iterations/n.bed, col=2)
    abline(h = mean(rowSums(long_totalR)/iterations/n.bed), col=2, lwd=2)

}else if(model == "frequency"){
    ################################## Multiple runs comparing between both ####################################
    iterations <- 100
    
    abx.short <- list()
    abx.long <- list()
    array_LOS_short <- list()
    array_LOS_long <- list()
    array_StartBact_short <- list()
    array_StartBact_long <- list()
    colo_table_filled_short <- list()
    colo_table_filled_long <- list()
    
    short_totalR <- matrix(NA, nrow = n.days, ncol = iterations)
    long_totalR <- matrix(NA, nrow = n.days, ncol = iterations)
    
    for (i in 1:iterations){
        
        print(i)
        #Generate length of stay and antibiotic duration table
        abx.short[[i]] <- abx.table(n.bed, n.days, mean.max.los, p.s, p.r, meanDur=mean_short)
        abx.long[[i]] <- abx.table(n.bed, n.days, mean.max.los, p.s, p.r, meanDur=mean_long)
        
        #Generate baseline carriage status
        array_LOS_short[[i]] <- array_LOS_func(los_duration=abx.short[[i]][2])
        array_LOS_long[[i]] <- array_LOS_func(los_duration=abx.long[[i]][2])
        
        #Update values for every day
        array_StartBact_short[[i]] <- gen_StartBact(los=array_LOS_short[[i]], prob_StartBact)
        array_StartBact_long[[i]] <- gen_StartBact(los=array_LOS_long[[i]], prob_StartBact)
        
        #output
        colo_table_filled_short[[i]] <- nextDay(bed_table= abx.short[[i]][[1]], array_LOS=array_LOS_short[[i]], 
                                                treat_table=abx.short[[i]][[3]], colo_table=array_StartBact_short[[i]], 
                                                pi_r1=pi_r1, pi_r2= pi_r2, mu1=mu1, mu2=mu2, 
                                                abx.r=abx.r,abx.s=abx.s,
                                                repop.r1 = repop.r1, repop.r2 = repop.r2, repop.r3 = repop.r3, 
                                                repop.s1 = repop.s1, repop.s2 = repop.s2)
        colo_table_filled_long[[i]] <- nextDay(bed_table= abx.long[[i]][[1]], array_LOS=array_LOS_long[[i]], 
                                               treat_table=abx.long[[i]][[3]], colo_table=array_StartBact_long[[i]], 
                                               pi_r1=pi_r1, pi_r2= pi_r2, mu1=mu1, mu2=mu2, 
                                               abx.r=abx.r,abx.s=abx.s,
                                               repop.r1 = repop.r1, repop.r2 = repop.r2, repop.r3 = repop.r3, 
                                               repop.s1 = repop.s1, repop.s2 = repop.s2)
        
        #increase overtime
        df.short <- data.frame(colo_table_filled_short[[i]][2])
        short_totalR[,i]<- rowSums(df.short)
        
        df.long <- data.frame(colo_table_filled_long[[i]][2])
        long_totalR[,i] <- rowSums(df.long)
        
    }
    
    # Check between 2 simple scenarios
    par(mfrow=c(1,2))
    
    x <- 1:n.days
    y.short <- short_totalR[,1]/n.bed
    r.plot.short <- plot(x, y.short, type="l", ylim=c(0,bact_slots), main=paste("Short mean abx:", mean_short, 
                                                                       "days", "(", round(mean(rowSums(short_totalR)/iterations/n.bed), digits=4), ")"),
                         xlab="Time", ylab="% patients carrying MDRO", col=rgb(0,0,0,alpha=0.1))
    for(i in 2:iterations){
        lines(x, short_totalR[,i]/n.bed, col=rgb(0,0,0,alpha=0.1))
    }
    lines(x, rowSums(short_totalR)/iterations/n.bed, col=2)
    abline(h = mean(rowSums(short_totalR)/iterations/n.bed), col=2, lwd=2)
    
    y.long <- long_totalR[,1]/n.bed
    r.plot.long <- plot(x, y.long, type="l", ylim=c(0,bact_slots), main=paste("Long mean abx:", mean_long, 
                                                                     "days", "(", round(mean(rowSums(long_totalR)/iterations/n.bed), digits=4), ")"),
                        xlab="Time", ylab="% patients carrying MDRO", col=rgb(0,0,0,alpha=0.1))
    for(i in 2:iterations){
        lines(x, long_totalR[,i]/n.bed, col=rgb(0,0,0,alpha=0.1))
    }
    lines(x, rowSums(long_totalR)/iterations/n.bed, col=2)
    abline(h = mean(rowSums(long_totalR)/iterations/n.bed), col=2, lwd=2)
    
}
