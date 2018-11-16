#########################################################################
#######Effect of antibiotic duration on hospitalised patients############
#########################################################################

rm(list=ls()) # Clean working environment

# model can be "simple", "binary", or "frequency"
model <- "simple"

source("default_params.R")
source(paste0("model_", model,".R"))

####################6. Visualisation #####################

if(model == "simple"){
                
    abx.short <- abx.table(n.bed, n.days, mean.max.los, p, meanDur = 4)
    abx.long <- abx.table(n.bed, n.days, mean.max.los, p, meanDur = 14)
    
    array_LOS_short <- array_LOS_func(los_duration=abx.short[2])
    array_LOS_long <- array_LOS_func(los_duration=abx.long[2])
    
    array_StartBact_short <- gen_StartBact(los=array_LOS_short, prob_StartBact_R, prop_S_nonR, n.bed, n.days)
    array_StartBact_long <- gen_StartBact(los=array_LOS_long, prob_StartBact_R, prop_S_nonR, n.bed, n.days)
    
    colo_table_filled_short <- nextDay(bed_table= abx.short[[1]], array_LOS=array_LOS_short, treat_table=abx.short[[3]], 
                                       colo_table=array_StartBact_short, pi_sr=pi_sr, mu_r=mu_r, pi_s=pi_s, pi_r=pi_r, abx.clear=abx.clear)
    colo_table_filled_long <- nextDay(bed_table= abx.long[[1]], array_LOS=array_LOS_long, treat_table=abx.long[[3]], 
                                      colo_table=array_StartBact_long, pi_sr=pi_sr, mu_r=mu_r, pi_s=pi_s, pi_r=pi_r, abx.clear=abx.clear)
    
    #present in one space
    par(mfrow=c(4,2))
    
    #plots
    #antibiotics pixel
    abx.mat.short<- as.matrix(abx.short[[3]])
    abx.mat.long<- as.matrix(abx.long[[3]])
    cols.a<-c('0'='grey', '1'='orange')
    (abx_img.short<-image(1:nrow(abx.mat.short),1:ncol(abx.mat.short),abx.mat.short,col=cols.a, xlab='Time', ylab='Bed No.', main="Short duration of antibiotics"))
    (abx_img.long<-image(1:nrow(abx.mat.long),1:ncol(abx.mat.long),abx.mat.long,col=cols.a, xlab='Time', ylab='Bed No.', main="Long duration of antibiotics"))
    
    #carriage pixel
    o<- colo_table_filled_short
    dim.saved<-dim(o)
    colo_table_filled_short[colo_table_filled_short=="R"]<-2
    colo_table_filled_short[colo_table_filled_short=="S"]<-3
    colo_table_filled_short[colo_table_filled_short=="N"]<-1
    colo_table_filled_short<-as.numeric(colo_table_filled_short)
    dim(colo_table_filled_short)<-dim.saved
    q<- colo_table_filled_long
    colo_table_filled_long[colo_table_filled_long=="R"]<-2
    colo_table_filled_long[colo_table_filled_long=="S"]<-3
    colo_table_filled_long[colo_table_filled_long=="N"]<-1
    colo_table_filled_long<-as.numeric(colo_table_filled_long)
    dim(colo_table_filled_long)<-dim.saved
    cols<-c('grey', 'red', 'chartreuse3')
    (col_img.short<-image(1:nrow(colo_table_filled_short),1:ncol(colo_table_filled_short), colo_table_filled_short, col=cols, xlab='Time', ylab='Bed No.'))
    (col_img.long<-image(1:nrow(colo_table_filled_long),1:ncol(colo_table_filled_long),colo_table_filled_long,col=cols, xlab='Time', ylab='Bed No.'))
    
    #increase overtime
    df.short<-as.data.frame(o)
    df.short$totalR<-rep(0, nrow(df.short))
    for (i in 1:nrow(df.short)) {
        for (j in 1:ncol(df.short)) {
            if (df.short[i,j]=="R") {
                df.short$totalR[i]<- df.short$totalR[i]+1
            } else {df.short$totalR[i]<- df.short$totalR[i]}
        }
    }
    
    df.long<-as.data.frame(q)
    df.long$totalR<-rep(0, nrow(df.long))
    for (i in 1:nrow(df.long)) {
        for (j in 1:ncol(df.long)) {
            if (df.long[i,j]=="R") {
                df.long$totalR[i]<- df.long$totalR[i]+1
            } else {df.long$totalR[i]<- df.long$totalR[i]}
        }
    }
    
    x<-c(1:n.days)
    y.short<-df.short$totalR/n.bed
    y.long<-df.long$totalR/n.bed
    r.plot.short<-plot(x,y.short, type="l", ylim=c(0,1), xlab="Time", ylab="Proportion of patients carrying resistant organisms")
    r.plot.long<-plot(x, y.long, type="l", ylim=c(0,1), xlab="Time", ylab="Proportion of patients carrying resistant organisms")
    
    #total acquisitions 
    df.short$newR<-rep(0, nrow(df.short))
    df.long$newR<-rep(0, nrow(df.long))
    
    for(i in 2:nrow(df.short)){
        for(j in 1:ncol(df.short)){
            if(df.short[i, j] == "R"){
                if(df.short[i-1, j] == "S"|df.short[i-1, j] == "N"){
                    df.short$newR[i]<-df.short$newR[i]+1 
                } else {df.short$newR[i]<-df.short$newR[i]}
            }
        }
    }
    y.short<- df.short$newR
    
    for(i in 2:nrow(df.long)){
        for(j in 1:ncol(df.long)){
            if(df.long[i, j] == "R"){
                if(df.long[i-1, j] == "S"|df.long[i-1, j] == "N"){
                    df.long$newR[i]<-df.long$newR[i]+1 
                } else {df.long$newR[i]<-df.long$newR[i]}
            }
        }
    }
    y.long<- df.long$newR
    
    acq.short<-plot(x,y.short, type='l', xlab="Time", ylab="Number of patients acquiring resistant organisms")
    acq.long<-plot(x,y.long, type='l', xlab="Time", ylab="Number of patients acquiring resistant organisms")

}else if(model == "binary"){
    abx.short <- abx.table(n.bed, n.days, mean.max.los, p.s, p.r.day1, p.r.dayafter, meanDur=4)
    abx.long <- abx.table(n.bed, n.days, mean.max.los, p.s, p.r.day1, p.r.dayafter, meanDur=14)
    
    array_LOS_short <- array_LOS_func(los_duration=abx.short[2])
    array_LOS_long <- array_LOS_func(los_duration=abx.long[2])
    
    array_StartBact_short <- gen_StartBact(los=array_LOS_short, prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR, n.bed, n.days)
    array_StartBact_long <- gen_StartBact(los=array_LOS_long, prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR, n.bed, n.days)
    
    colo_table_filled_short <- nextDay(bed_table= abx.short[[1]], array_LOS=array_LOS_short, 
                                       treat_table=abx.short[[3]], colo_table=array_StartBact_short, 
                                       pi_r1=pi_r1, mu1=mu1, mu2=mu2, pi_r2=pi_r2, 
                                       repop.r1 = repop.r1, repop.r2 = repop.r2,
                                       repop.s1 = repop.s1, repop.s2 = repop.s2, repop.s3 = repop.s3, abx.r = abx.r, abx.s = abx.s)
    colo_table_filled_long <- nextDay(bed_table= abx.long[[1]], array_LOS=array_LOS_long, 
                                      treat_table=abx.long[[3]], colo_table=array_StartBact_long, 
                                      pi_r1=pi_r1, mu1=mu1, mu2=mu2, pi_r2=pi_r2, 
                                      repop.r1 = repop.r1, repop.r2 = repop.r2,
                                      repop.s1 = repop.s1, repop.s2 = repop.s2, repop.s3 = repop.s3, abx.r = abx.r, abx.s = abx.s)
    #present in one space
    par(mfrow=c(4,2))
    
    #plots
    #antibiotics pixel
    abx.mat.short<- abx.short[[3]]
    abx.mat.long<- abx.long[[3]]
    cols.a<-c('0'='grey', '1'='orange', '2'= 'orange3','3'= 'orange3')
    abx_img.short<-image(1:nrow(abx.mat.short),1:ncol(abx.mat.short),abx.mat.short,col=cols.a, xlab='Time', ylab='Bed No.', main="Short duration of antibiotics")
    abx_img.long<-image(1:nrow(abx.mat.long),1:ncol(abx.mat.long),abx.mat.long,col=cols.a, xlab='Time', ylab='Bed No.', main="Long duration of antibiotics")
    
    #carriage pixel
    orig_short <- colo_table_filled_short
    dim.saved<-dim(orig_short)
    colo_table_filled_short[colo_table_filled_short=="S" | colo_table_filled_short=="s"]<-1
    colo_table_filled_short[colo_table_filled_short=="sr"|colo_table_filled_short=="Sr"|colo_table_filled_short=="sR"]<-2
    colo_table_filled_short<-as.numeric(colo_table_filled_short)
    dim(colo_table_filled_short)<-dim.saved
    
    orig_long <- colo_table_filled_long
    colo_table_filled_long[colo_table_filled_long=="S" | colo_table_filled_long=="s"]<-1
    colo_table_filled_long[colo_table_filled_long=="sr"|colo_table_filled_long=="Sr"|colo_table_filled_long=="sR"]<-2
    colo_table_filled_long<-as.numeric(colo_table_filled_long)
    dim(colo_table_filled_long)<-dim.saved
    
    cols<-c('chartreuse3', 'red')
    col_img.short<-image(1:nrow(colo_table_filled_short),1:ncol(colo_table_filled_short), colo_table_filled_short, col=cols, xlab='Time', ylab='Bed No.')
    col_img.long<-image(1:nrow(colo_table_filled_long),1:ncol(colo_table_filled_long),colo_table_filled_long,col=cols, xlab='Time', ylab='Bed No.')
    
    #increase overtime
    df.short<-as.data.frame(orig_short)
    df.short$totalR<-rep(0, nrow(df.short))
    for (i in 1:nrow(df.short)) {
        for (j in 1:ncol(df.short)) {
            if(df.short[i, j]=="sR"|df.short[i, j]=="sr"|df.short[i, j]=="Sr") {
                df.short$totalR[i] <- df.short$totalR[i] + 1
            }
        }
    }
    
    df.long<-as.data.frame(orig_long)
    df.long$totalR<-rep(0, nrow(df.long))
    for (i in 1:nrow(df.long)) {
        for (j in 1:ncol(df.long)) {
            if (df.long[i,j]=="sR"|df.long[i,j]=="sr"|df.long[i,j]=="Sr") {
                df.long$totalR[i] <- df.long$totalR[i] + 1
            }
        }
    }
    
    x <- c(1:n.days)
    y.short <- df.short$totalR/n.bed
    y.long <- df.long$totalR/n.bed
    r.plot.short<-plot(x,y.short, type="l", ylim=c(0,1), xlab="Time", ylab="% patients carrying MDRO")
    lines(x,y.long, col=2)
    r.plot.long<-plot(x, y.long, type="l", ylim=c(0,1), xlab="Time", ylab="% patients carrying MDRO")
    
    #total acquisitions 
    df.short$newR <- rep(0, nrow(df.short))
    df.long$newR <- rep(0, nrow(df.long))
    
    for(i in 2:nrow(df.short)){
        for(j in 1:ncol(df.short)){
            if((df.short[i, j] == "sR"|df.short[i, j] == "sr"|df.short[i, j] == "Sr") & (df.short[i-1, j] == "S" | df.short[i-1, j] == "s")){
                df.short$newR[i] <- df.short$newR[i] + 1 
            }
        }
    }
    y.short<- df.short$newR
    
    for(i in 2:nrow(df.long)){
        for(j in 1:ncol(df.long)){
            if((df.long[i, j] == "sR"|df.long[i, j] == "sr"|df.long[i, j] == "Sr") & (df.long[i-1, j] == "S" | df.long[i-1, j] == "s")){
                df.long$newR[i] <- df.long$newR[i] + 1 
            }
        }
    }
    y.long<- df.long$newR
    
    acq.short<-plot(x,y.short, type='l', xlab="Time", ylab="Number of patients acquiring MDRO")
    lines(x,y.long, col=2)
    acq.long<-plot(x,y.long, type='l', xlab="Time", ylab="Number of patients acquiring MDRO")

}else if(model == "frequency"){
    
    abx.short<-abx.table(n.bed, n.days, mean.max.los, p.s, p.r.day1, p.r.dayafter, meanDur=4)
    abx.long<-abx.table(n.bed, n.days, mean.max.los, p.s, p.r.day1, p.r.dayafter, meanDur=14)
    
    array_LOS_short<-array_LOS_func(los_duration=abx.short[2])
    array_LOS_long<- array_LOS_func(los_duration=abx.long[2])
    
    array_StartBact_short<-gen_StartBact(los=array_LOS_short, K, t_mean = 4.0826, t_sd = 1.1218, r_mean =1.7031, r_sd = 1.8921, n.beds, n.days)
    array_StartBact_long<-gen_StartBact(los=array_LOS_long, K, t_mean = 4.0826, t_sd = 1.1218, r_mean =1.7031, r_sd = 1.8921, n.beds, n.days)
    
    colo_table_filled_short <- nextDay(bed_table= abx.short[[1]], array_LOS=array_LOS_short, 
                                       treat_table=abx.short[[3]], colo_table=array_StartBact_short, 
                                       pi_r1, pi_r2, mu1, mu2, abx.r, abx.s,
                                       repop.r1, repop.r2, repop.r3, repop.s1, repop.s2)
    colo_table_filled_long <- nextDay(bed_table= abx.long[[1]], array_LOS=array_LOS_long, 
                                      treat_table=abx.long[[3]], colo_table=array_StartBact_long, 
                                      pi_r1, pi_r2, mu1, mu2, abx.r, abx.s,
                                      repop.r1, repop.r2, repop.r3, repop.s1, repop.s2)
    
    ####################6. Visualisation #####################
    
    #present in one space
    #par(mfrow=c(4,2))
    layout(rbind(c(1,1,2,2),
                 c(3,3,4,4),
                 c(5,5,6,6),
                 #c(7,7,8,8),
                 c(0,7,7,0)
    ))
    
    #plots
    #antibiotics pixel
    abx.mat.short<- abx.short[[3]]
    abx.mat.long<- abx.long[[3]]
    cols.a<-c('0'='grey', '1'='orange', '2'= 'orange3','3'= 'orange3')
    abx_img.short<-image(1:nrow(abx.mat.short),1:ncol(abx.mat.short),abx.mat.short,col=cols.a, xlab='Time', ylab='Bed No.', main="Short duration of antibiotics")
    abx_img.long<-image(1:nrow(abx.mat.long),1:ncol(abx.mat.long),abx.mat.long,col=cols.a, xlab='Time', ylab='Bed No.', main="Long duration of antibiotics")
    
    #carriage differentiation pixel
    orig_short <- colo_table_filled_short[[2]]
    colo_table_filled_short <- diff(colo_table_filled_short[[2]])
    dim.saved <- dim(colo_table_filled_short)
    colo_table_filled_short[colo_table_filled_short > 0] <- 1
    colo_table_filled_short[colo_table_filled_short <= 0] <- 2
    colo_table_filled_short <- as.numeric(colo_table_filled_short)
    dim(colo_table_filled_short)<-dim.saved
    
    orig_long <- colo_table_filled_long[[2]]
    colo_table_filled_long <- diff(colo_table_filled_long[[2]])
    colo_table_filled_long[colo_table_filled_long > 0]<-1
    colo_table_filled_long[colo_table_filled_long <= 0]<-2
    colo_table_filled_long<-as.numeric(colo_table_filled_long)
    dim(colo_table_filled_long)<-dim.saved
    
    cols<-c('chartreuse3','red')
    col_img.short<-image(1:nrow(colo_table_filled_short),1:ncol(colo_table_filled_short), colo_table_filled_short, col=cols, xlab='Time', ylab='Bed No.')
    col_img.long<-image(1:nrow(colo_table_filled_long),1:ncol(colo_table_filled_long),colo_table_filled_long,col=cols, xlab='Time', ylab='Bed No.')
    
    #carriage pixel
    # colo_table_filled_short <- orig_short
    # 
    # colo_table_filled_long <- orig_long
    # 
    # colfunc <- colorRampPalette(c("black", "white"))
    # cols <- colfunc(0, )
    # col_img.short<-image(1:nrow(colo_table_filled_short),1:ncol(colo_table_filled_short), colo_table_filled_short, col=cols, xlab='Time', ylab='Bed No.')
    # col_img.long<-image(1:nrow(colo_table_filled_long),1:ncol(colo_table_filled_long),colo_table_filled_long,col=cols, xlab='Time', ylab='Bed No.')
    
    #increase overtime
    df.short<-as.data.frame(orig_short)
    df.short$totalR<-rowSums(df.short)
    
    df.long<-as.data.frame(orig_long)
    df.long$totalR<-rowSums(df.long)
    
    x <- c(1:n.days)
    y.short <- df.short$totalR/n.bed
    y.long <- df.long$totalR/n.bed
    r.plot.short<-plot(x,y.short, type="l", ylim=c(0,bact_slots), xlab="Time", ylab="% patients carrying MDRO")
    #lines(x,y.long, col=2)
    r.plot.long<-plot(x, y.long, type="l", ylim=c(0,bact_slots), xlab="Time", ylab="% patients carrying MDRO")
    
    # difference between both plot
    y.diff <- y.long-y.short
    r.plot.long<-plot(x, y.diff, type="l", ylim=c(-bact_slots,bact_slots), xlab="Time", ylab="% patients carrying MDRO")
    abline(h=0)
}
