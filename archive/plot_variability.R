
rm(list=ls()) # Clean working environment

# model can be "simple", "binary", or "frequency"
model <- "frequency"

source("default_params.R")
source("los_abx_matrix.R")
source(paste0("model_", model,".R"))

iterations <- 100

#########################################################################
#######Effect of antibiotic duration on hospitalised patients############
#########################################################################

if(model == "simple"){
    
    ################################## Multiple runs comparing between both ####################################
    
    patient.matrix<-list()
    los.array<-list()
    abx.short <- list()
    abx.long <- list()
    colo.matrix<- list()
    colo_table_filled_short<-list()
    colo_table_filled_long<-list()
    short_totalR <- matrix(NA, nrow = n.day, ncol = iterations)
    long_totalR <- matrix(NA, nrow = n.day, ncol = iterations)
    
    for (i in 1:iterations){
        
        print(i)
        
        
       
        patient.matrix[[i]]<-patient.table(n.bed = n.bed, n.day=n.day, mean.max.los=mean.max.los, timestep=1)
        los.array[[i]]<-summary.los(patient.matrix=patient.matrix[[i]])
        
        #Generate length of stay and antibiotic duration table
        abx.short[[i]] <- abx.table(patient.matrix = patient.matrix[[i]], los.array=los.array[[i]], p=p, meanDur=short_dur, sdDur=sdDur, timestep=1)
        abx.long[[i]] <- abx.table(patient.matrix = patient.matrix[[i]], los.array=los.array[[i]], p=p, meanDur=long_dur, sdDur=sdDur, timestep=1)
        
        #Generate baseline carriage status
        colo.matrix[[i]] <- colo.table(patient.matrix = patient.matrix[[i]], los.array = los.array[[i]], prob_StartBact_R = prob_StartBact_R, prop_S_nonR = prop_S_nonR)
        
        #output
        colo_table_filled_short[[i]] <- nextDay(patient.matrix=patient.matrix[[i]], los.array=los.array[[i]], abx.matrix=abx.short[[i]], colo.matrix=colo.matrix[[i]], 
                                                bif=bif, pi_ssr=pi_ssr, repop.s1=repop.s1, mu_r=mu_r, abx.clear=abx.clear, timestep=1)
        colo_table_filled_long[[i]] <- nextDay(patient.matrix=patient.matrix[[i]], los.array=los.array[[i]], abx.matrix=abx.long[[i]], colo.matrix=colo.matrix[[i]], 
                                               bif=bif, pi_ssr=pi_ssr, repop.s1=repop.s1, mu_r=mu_r, abx.clear=abx.clear, timestep=1)
        
        #increase overtime
        df.short <- colo_table_filled_short[[i]]
        short_totalR[,i]<- rowSums(df.short == "R")
        
        df.long <- colo_table_filled_long[[i]]
        long_totalR[,i] <- rowSums(df.long == "R")
        
    }
    
    # Check between 2 simple scenarios
    par(mfrow=c(1,2))
    
    x <- 1:n.day
    y.short <- short_totalR[,1]/n.bed
    r.plot.short <- plot(x, y.short, type="l", ylim=c(0,1), main=paste("Short mean abx:", short_dur, 
                                                                       "days", "(", round(mean(rowSums(short_totalR)/iterations/n.bed), digits=4), ")"),
                         xlab="Time", ylab="% patients carrying MDRO", col=rgb(0,0,0,alpha=0.1))
    for(i in 2:iterations){
        lines(x, short_totalR[,i]/n.bed, col=rgb(0,0,0,alpha=0.1))
    }
    lines(x, rowSums(short_totalR)/iterations/n.bed, col=2)
    abline(h = mean(rowSums(short_totalR)/iterations/n.bed), col=2, lwd=2)
    
    y.long <- long_totalR[,1]/n.bed
    r.plot.long <- plot(x, y.long, type="l", ylim=c(0,1), main=paste("Long mean abx:", long_dur, 
                                                                     "days", "(", round(mean(rowSums(long_totalR)/iterations/n.bed), digits=4), ")"),
                        xlab="Time", ylab="% patients carrying MDRO", col=rgb(0,0,0,alpha=0.1))
    for(i in 2:iterations){
        lines(x, long_totalR[,i]/n.bed, col=rgb(0,0,0,alpha=0.1))
    }
    lines(x, rowSums(long_totalR)/iterations/n.bed, col=2)
    abline(h = mean(rowSums(long_totalR)/iterations/n.bed), col=2, lwd=2)


}else if(model == "binary"){
    
    ################################## Multiple runs comparing between both ####################################

    patient.matrix<-list()
    los.array<-list()
    abx.short <- list()
    abx.long <- list()
    colo.matrix<- list()
    colo_table_filled_short<-list()
    colo_table_filled_long<-list()
    short_totalR <- matrix(NA, nrow = n.day, ncol = iterations)
    long_totalR <- matrix(NA, nrow = n.day, ncol = iterations)
    
    for (i in 1:iterations){
        
        print(i)
       
        patient.matrix[[i]]<-patient.table(n.bed = n.bed, n.day=n.day, mean.max.los=mean.max.los, timestep=1)
        los.array[[i]]<-summary.los(patient.matrix=patient.matrix[[i]])
        
        #Generate length of stay and antibiotic duration table
        abx.short[[i]] <- abx.table(patient.matrix = patient.matrix[[i]], los.array=los.array[[i]], p.s=p.s, p.r.day1=p.r.day1, p.r.dayafter=p.r.dayafter,
                                    meanDur.s=short_dur.s, meanDur.r=short_dur.r, sdDur=sdDur, timestep=1)
        abx.long[[i]] <- abx.table(patient.matrix = patient.matrix[[i]], los.array=los.array[[i]], p.s=p.s, p.r.day1=p.r.day1, p.r.dayafter=p.r.dayafter,
                                   meanDur.s=long_dur.s, meanDur.r=long_dur.r, sdDur=sdDur, timestep=1)
        
        #Generate baseline carriage status
        colo.matrix[[i]] <- colo.table(patient.matrix = patient.matrix[[i]], los.array = los.array[[i]], prob_StartBact_R = prob_StartBact_R, prop_S_nonR = prop_S_nonR,prop_Sr_inR=prop_Sr_inR, prop_sr_inR=prop_sr_inR)
        
        
        #output
        colo_table_filled_short[[i]] <- nextDay(patient.matrix = patient.matrix[[i]], abx.matrix=abx.short[[i]], colo.matrix=colo.matrix[[i]], 
                                                pi_r1=pi_r1, bif=bif, mu1=mu1, mu2=mu2, repop.r=repop.r, 
                                                repop.s1=repop.s1, repop.s2=repop.s2,depop.r=depop.r, abx.r=abx.r, abx.s=abx.s, timestep=1)
        colo_table_filled_long[[i]] <- nextDay(patient.matrix = patient.matrix[[i]], abx.matrix=abx.long[[i]], colo.matrix=colo.matrix[[i]], 
                                               pi_r1=pi_r1, bif=bif, mu1=mu1, mu2=mu2, repop.r=repop.r, 
                                               repop.s1=repop.s1, repop.s2=repop.s2,depop.r=depop.r, abx.r=abx.r, abx.s=abx.s, timestep=1)
    
        #increase overtime
        df.short <- colo_table_filled_short[[i]]
        short_totalR[,i]<- rowSums(df.short == "sR")
        
        df.long <- colo_table_filled_long[[i]]
        long_totalR[,i] <- rowSums(df.long == "sR")
    
    }
    
    # Check between 2 simple scenarios
    par(mfrow=c(1,2))
    
    x <- 1:n.day
    y.short <- short_totalR[,1]/n.bed
    r.plot.short <- plot(x, y.short, type="l", ylim=c(0,1), main=paste("Short mean abx:", short_dur, 
                         "days", "(", round(mean(rowSums(short_totalR)/iterations/n.bed), digits=4), ")"),
                         xlab="Time", ylab="% patients carrying MDRO", col=rgb(0,0,0,alpha=0.1))
    for(i in 2:iterations){
        lines(x, short_totalR[,i]/n.bed, col=rgb(0,0,0,alpha=0.1))
    }
    lines(x, rowSums(short_totalR)/iterations/n.bed, col=2)
    abline(h = mean(rowSums(short_totalR)/iterations/n.bed), col=2, lwd=2)
    
    y.long <- long_totalR[,1]/n.bed
    r.plot.long <- plot(x, y.long, type="l", ylim=c(0,1), main=paste("Long mean abx:", long_dur, 
                        "days", "(", round(mean(rowSums(long_totalR)/iterations/n.bed), digits=4), ")"),
                        xlab="Time", ylab="% patients carrying MDRO", col=rgb(0,0,0,alpha=0.1))
    for(i in 2:iterations){
        lines(x, long_totalR[,i]/n.bed, col=rgb(0,0,0,alpha=0.1))
    }
    lines(x, rowSums(long_totalR)/iterations/n.bed, col=2)
    abline(h = mean(rowSums(long_totalR)/iterations/n.bed), col=2, lwd=2)

}else if(model == "frequency"){
    
    ################################## Multiple runs comparing between both ####################################
    
    patient.matrix<-list()
    los.array<-list()
    abx.short <- list()
    abx.long <- list()
    colo.matrix<- list()
    colo_table_filled_short<-list()
    colo_table_filled_long<-list()
    short_totalR <- matrix(NA, nrow = n.day, ncol = iterations)
    long_totalR <- matrix(NA, nrow = n.day, ncol = iterations)
    
    for (i in 1:iterations){
        
        print(i)
        
        patient.matrix[[i]]<-patient.table(n.bed = n.bed, n.day=n.day, mean.max.los=mean.max.los, timestep=1)
        los.array[[i]]<-summary.los(patient.matrix=patient.matrix[[i]])
        
        #Generate length of stay and antibiotic duration table
        abx.short[[i]] <- abx.table(patient.matrix = patient.matrix[[i]], los.array=los.array[[i]], p.s=p.s, p.r.day1=p.r.day1, p.r.dayafter=p.r.dayafter,
                                    meanDur.s=short_dur.s, meanDur.r=short_dur.r, sdDur=sdDur, timestep=1)
        abx.long[[i]] <- abx.table(patient.matrix = patient.matrix[[i]], los.array=los.array[[i]], p.s=p.s, p.r.day1=p.r.day1, p.r.dayafter=p.r.dayafter,
                                   meanDur.s=long_dur.s, meanDur.r=long_dur.r, sdDur=sdDur, timestep=1)
        
        #Generate baseline carriage status
        colo.matrix[[i]] <- colo.table(patient.matrix = patient.matrix[[i]], los.array = los.array[[i]], t_mean=t_mean, t_sd=t_sd, r_mean=r_mean, r_sd=r_sd)
        
        #output
        colo_table_filled_short[[i]] <- nextDay(patient.matrix = patient.matrix[[i]], los.array = los.array[[i]], abx.matrix=abx.short[[i]], colo.matrix=colo.matrix[[i]], 
                                                pi_r=pi_r, K=K, r_thres=r_thres, r_growth=r_growth, r_trans=r_trans, r_thres=r_thres,
                                                abxr_killr=abxr_killr, abxr_kills=abxr_kills, abxs_kills=abxs_kills, timestep=1)
        colo_table_filled_long[[i]] <- nextDay(patient.matrix = patient.matrix[[i]], los.array = los.array[[i]], abx.matrix=abx.long[[i]], colo.matrix=colo.matrix[[i]], 
                                               pi_r=pi_r, K=K, r_thres=r_thres, r_growth=r_growth, r_trans=r_trans, r_thres=r_thres,
                                               abxr_killr=abxr_killr, abxr_kills=abxr_kills, abxs_kills=abxs_kills, timestep=1)
        
        #increase overtime
        df.short <- data.frame(colo_table_filled_short[[i]])
        short_totalR[,i]<- rowSums(df.short)
        
        df.long <- data.frame(colo_table_filled_long[[i]])
        long_totalR[,i] <- rowSums(df.long)
        
    }
    
    # Check between 2 simple scenarios
    par(mfrow=c(1,2))
    
    x <- 1:n.day
    y.short <- short_totalR[,1]/n.bed
    y.long <- long_totalR[,1]/n.bed
    r.plot.short <- plot(x, y.short, type="l", ylim=c(0,max(max(y.short),max(y.long))), main=paste("Short mean abx:", short_dur, 
                                                                       "days", "(", round(mean(rowSums(short_totalR)/iterations/n.bed), digits=4), ")"),
                         xlab="Time", ylab="% patients carrying MDRO", col=rgb(0,0,0,alpha=0.1))
    for(i in 2:iterations){
        lines(x, short_totalR[,i]/n.bed, col=rgb(0,0,0,alpha=0.1))
    }
    lines(x, rowSums(short_totalR)/iterations/n.bed, col=2)
    abline(h = mean(rowSums(short_totalR)/iterations/n.bed), col=2, lwd=2)
    
    r.plot.long <- plot(x, y.long, type="l", ylim=c(0,max(max(y.short),max(y.long))), main=paste("Long mean abx:", long_dur, 
                                                                     "days", "(", round(mean(rowSums(long_totalR)/iterations/n.bed), digits=4), ")"),
                        xlab="Time", ylab="% patients carrying MDRO", col=rgb(0,0,0,alpha=0.1))
    for(i in 2:iterations){
        lines(x, long_totalR[,i]/n.bed, col=rgb(0,0,0,alpha=0.1))
    }
    lines(x, rowSums(long_totalR)/iterations/n.bed, col=2)
    abline(h = mean(rowSums(long_totalR)/iterations/n.bed), col=2, lwd=2)
    
}
