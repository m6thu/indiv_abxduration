#########################################################################
#######Effect of antibiotic duration on hospitalised patients############
#########################################################################
setwd('Desktop/indiv_abxduration/')
rm(list=ls()) # Clean working environment

# model can be "simple", "binary", or "frequency"
model <- "frequency"

source("default_params.R")
source(paste0("model_", model,".R"))

####################6. Visualisation #####################

if(model == "simple"){
    
    patient.matrix<- patient.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, timestep=1)
    los.array<-summary.los(patient.matrix)
    abx.short<-abx.table(patient.matrix=patient.matrix, los.array=los.array, p=p, meanDur=short_dur, sdDur=sdDur, timestep=1)
    abx.long <-abx.table(patient.matrix=patient.matrix, los.array=los.array, p=p, meanDur=long_dur, sdDur=sdDur, timestep=1)
    colo.matrix<- colo.table(patient.matrix=patient.matrix, los.array=los.array, prob_StartBact_R=prob_StartBact_R, prop_S_nonR=prop_S_nonR)
    colo_table_filled_short<-nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.short, colo.matrix=colo.matrix, 
                              bif=bif, pi_ssr=pi_ssr, repop.s1=repop.s1, mu_r=mu_r, abx.clear=abx.clear, timestep=1)
    colo_table_filled_long<- nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.long, colo.matrix=colo.matrix, 
                                     bif=bif, pi_ssr=pi_ssr, repop.s1=repop.s1, mu_r=mu_r, abx.clear=abx.clear, timestep=1)
    
    #present in one space
    layout(matrix(c(1:8,9,9), 5, 2, byrow = TRUE))
    par( mar=c(5, 6, 4, 2)+ 0.1)
    
    #plots
    #antibiotics pixel
    abx.mat.short<- as.matrix(abx.short)
    abx.mat.long<- as.matrix(abx.long)
    cols.a<-c('0'='grey', '1'='orange')
    (abx_img.short<-image(1:nrow(abx.mat.short),1:ncol(abx.mat.short),abx.mat.short,col=cols.a, xlab='Time', ylab='Bed No.', main="Short duration of antibiotics"))
    (abx_img.long<-image(1:nrow(abx.mat.long),1:ncol(abx.mat.long),abx.mat.long,col=cols.a, xlab='Time', ylab='Bed No.', main="Long duration of antibiotics"))
    
    #carriage pixel
    o<- colo_table_filled_short
    dim.saved<-dim(o)
    colo_table_filled_short[colo_table_filled_short=="R"]<-2
    colo_table_filled_short[colo_table_filled_short=="S"]<-3
    colo_table_filled_short[colo_table_filled_short=="ss"]<-1
    colo_table_filled_short<-as.numeric(colo_table_filled_short)
    dim(colo_table_filled_short)<-dim.saved
    q<- colo_table_filled_long
    colo_table_filled_long[colo_table_filled_long=="R"]<-2
    colo_table_filled_long[colo_table_filled_long=="S"]<-3
    colo_table_filled_long[colo_table_filled_long=="ss"]<-1
    colo_table_filled_long<-as.numeric(colo_table_filled_long)
    dim(colo_table_filled_long)<-dim.saved
    cols<-c('darkseagreen1', 'red', 'chartreuse3')
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
    
    x<-c(1:n.day)
    y.short<-df.short$totalR/n.bed
    y.long<-df.long$totalR/n.bed
    r.plot.short<-plot(x,y.short, type="l", ylim=c(0,1), xlab="Time", ylab="Proportion of patients \ncarrying resistant organisms")
    r.plot.long<-plot(x, y.long, type="l", ylim=c(0,1), xlab="Time", ylab="Proportion of patients \ncarrying resistant organisms")
    
    #total acquisitions 
    df.short$newR<-rep(0, nrow(df.short))
    df.long$newR<-rep(0, nrow(df.long))
    
    for(i in 2:nrow(df.short)){
        for(j in 1:ncol(df.short)){
            if(df.short[i, j] == "R"){
                if(df.short[i-1, j] == "S"|df.short[i-1, j] == "ss"){
                    df.short$newR[i]<-df.short$newR[i]+1 
                } else {df.short$newR[i]<-df.short$newR[i]}
            }
        }
    }
    y.short<- df.short$newR
    
    for(i in 2:nrow(df.long)){
        for(j in 1:ncol(df.long)){
            if(df.long[i, j] == "R"){
                if(df.long[i-1, j] == "S"|df.long[i-1, j] == "ss"){
                    df.long$newR[i]<-df.long$newR[i]+1 
                } else {df.long$newR[i]<-df.long$newR[i]}
            }
        }
    }
    y.long<- df.long$newR
    
    acq.short<-plot(x,y.short, type='l', xlab="Time", ylab="Number of patients \nacquiring resistant organisms")
    acq.long<-plot(x,y.long, type='l', xlab="Time", ylab="Number of patients \nacquiring resistant organisms")

    #final plot for the difference between short and long durations 
    colo_table_filled_short[colo_table_filled_short==2]<-0
    colo_table_filled_short[colo_table_filled_short==1]<-0
    colo_table_filled_short[colo_table_filled_short==3]<-1
    dailyshort<-rowSums(colo_table_filled_short)
    colo_table_filled_long[colo_table_filled_long==2]<-0
    colo_table_filled_long[colo_table_filled_long==1]<-0
    colo_table_filled_long[colo_table_filled_long==3]<-1
    dailylong<-rowSums(colo_table_filled_long)
    dailydiff<- dailylong-dailyshort
    plot(x,dailydiff, type='l', xlab="Time", ylab="Difference in number of patients \ncarrying resistant organisms \ncomparing long vs short duration")
    
}else if(model == "binary"){
    
    patient.matrix<- patient.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, timestep=1)
    los.array<-summary.los(patient.matrix)
    abx.short<-abx.table(patient.matrix=patient.matrix, los.array=los.array, p.s=p.s, p.r.day1=p.r.day1, p.r.dayafter=p.r.dayafter,
                         meanDur.s=short_dur.s, meanDur.r=short_dur.r,  sdDur=sdDur, timestep=1)
    abx.long <-abx.table(patient.matrix=patient.matrix, los.array=los.array, p.s=p.s, p.r.day1=p.r.day1, p.r.dayafter=p.r.dayafter,
                         meanDur.s=short_dur.s, meanDur.r=long_dur.r,  sdDur=sdDur, timestep=1)
    colo.matrix<- colo.table(patient.matrix=patient.matrix, los.array=los.array, prob_StartBact_R=prob_StartBact_R, prop_S_nonR=prop_S_nonR,prop_Sr_inR=prop_Sr_inR, prop_sr_inR=prop_sr_inR)
    colo_table_filled_short<-nextDay(patient.matrix=patient.matrix, abx.matrix=abx.short, colo.matrix=colo.matrix, 
                                     pi_r1=pi_r1, bif=bif, mu1=mu1, mu2=mu2, repop.r=repop.r,
                                     repop.s1=repop.s1, repop.s2=repop.s2, depop.r=depop.r, abx.r=abx.r, abx.s=abx.s, timestep=1)
    colo_table_filled_long<- nextDay(patient.matrix=patient.matrix, abx.matrix=abx.long, colo.matrix=colo.matrix, 
                                     pi_r1=pi_r1, bif=bif, mu1=mu1, mu2=mu2, repop.r=repop.r,
                                     repop.s1=repop.s1, repop.s2=repop.s2, depop.r=depop.r, abx.r=abx.r, abx.s=abx.s, timestep=1)

    #present in one space
    layout(matrix(c(1:8,9,9), 5, 2, byrow = TRUE))
    par( mar=c(5, 6, 4, 2)+ 0.1)
    
    #plots
    #antibiotics pixel
    abx.mat.short<- abx.short
    abx.mat.long<- abx.long
    cols.a<-c('0'='grey', '1'='orange', '2'= 'orange3','3'= 'orange3')
    abx_img.short<-image(1:nrow(abx.mat.short),1:ncol(abx.mat.short),abx.mat.short,col=cols.a, xlab='Time', ylab='Bed No.', main="Short duration of antibiotics")
    abx_img.long<-image(1:nrow(abx.mat.long),1:ncol(abx.mat.long),abx.mat.long,col=cols.a, xlab='Time', ylab='Bed No.', main="Long duration of antibiotics")
    
    #carriage pixel
    orig_short <- colo_table_filled_short
    dim.saved<-dim(orig_short)
    colo_table_filled_short[colo_table_filled_short=="S" | colo_table_filled_short=="ss"]<-1
    colo_table_filled_short[colo_table_filled_short=="sr"|colo_table_filled_short=="Sr"|colo_table_filled_short=="sR"]<-2
    colo_table_filled_short<-as.numeric(colo_table_filled_short)
    dim(colo_table_filled_short)<-dim.saved
    
    orig_long <- colo_table_filled_long
    colo_table_filled_long[colo_table_filled_long=="S" | colo_table_filled_long=="ss"]<-1
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
    
    x <- c(1:n.day)
    y.short <- df.short$totalR/n.bed
    y.long <- df.long$totalR/n.bed
    r.plot.short<-plot(x,y.short, type="l", ylim=c(0,1), xlab="Time", ylab="% patients carrying MDRO")
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
    acq.long<-plot(x,y.long, type='l', xlab="Time", ylab="Number of patients acquiring MDRO")
    
    #final plot for the difference between short and long durations 
    dim.saved<-dim(orig_short)
    orig_short[orig_short=="sR"]<-1
    orig_short[orig_short=="ss"|orig_short=="sR"|orig_short=="Sr"|orig_short=="S"|orig_short=="sr"]<-0
    orig_s<-as.numeric(orig_short)
    dim(orig_s)<-dim.saved
    dailyshort<-rowSums(orig_s)
    orig_long[orig_long=="sR"]<-1
    orig_long[orig_long=="ss"|orig_long=="sR"|orig_long=="Sr"|orig_long=="S"|orig_long=="sr"]<-0
    orig_l<-as.numeric(orig_long)
    dim(orig_l)<-dim.saved
    dailylong<-rowSums(orig_l)
    dailydiff<- dailylong-dailyshort
    plot(x,dailydiff, type='l', xlab="Time", ylab="Difference in number of patients \ncarrying resistant organisms \ncomparing long vs short duration")
    abline(h=0, lty = 3, col='red')
    
}else if(model == "frequency"){
    
    patient.matrix<- patient.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, timestep=1)
    los.array<-summary.los(patient.matrix)
    abx.short<-abx.table(patient.matrix=patient.matrix, los.array=los.array, p.s=p.s, p.r.day1=p.r.day1, p.r.dayafter=p.r.dayafter,
                         meanDur.s=short_dur.s, meanDur.r=short_dur.r,  sdDur=sdDur, timestep=1)
    abx.long <-abx.table(patient.matrix=patient.matrix, los.array=los.array, p.s=p.s, p.r.day1=p.r.day1, p.r.dayafter=p.r.dayafter,
                         meanDur.s=short_dur.s, meanDur.r=long_dur.r,  sdDur=sdDur, timestep=1)
    colo.matrix<- colo.table(patient.matrix=patient.matrix, los.array=los.array, t_mean=t_mean, t_sd=t_sd, r_mean=r_mean, r_sd=r_sd)
    colo_table_filled_short<-nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.short, colo.matrix=colo.matrix, 
                                      pi_r=pi_r, K=K, r_thres=r_thres, r_growth=r_growth, r_trans=r_trans, 
                                      abxr_killr=abxr_killr, abxr_kills=abxr_kills, abxs_kills=abxs_kills, timestep=1)
    colo_table_filled_long<-nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.long, colo.matrix=colo.matrix, 
                                      pi_r=pi_r, K=K, r_thres=r_thres, r_growth=r_growth, r_trans=r_trans, 
                                      abxr_killr=abxr_killr, abxr_kills=abxr_kills, abxs_kills=abxs_kills, timestep=1)
    
    ####################6. Visualisation #####################
    
    #present in one space
    layout(matrix(c(1:6,7,7), 4, 2, byrow = TRUE))
    par( mar=c(5, 6, 4, 2)+ 0.1)
    # layout(rbind(c(1,1,2,2),
    #              c(3,3,4,4),
    #              c(5,5,6,6),
    #              #c(7,7,8,8),
    #              c(0,7,7,0)
    # ))
    # 
    #plots
    #antibiotics pixel
    abx.mat.short<- abx.short
    abx.mat.long<- abx.long
    cols.a<-c('0'='grey', '1'='orange', '2'= 'orange3','3'= 'orange3')
    abx_img.short<-image(1:nrow(abx.mat.short),1:ncol(abx.mat.short),abx.mat.short,col=cols.a, xlab='Time', ylab='Bed No.', main="Short duration of antibiotics")
    abx_img.long<-image(1:nrow(abx.mat.long),1:ncol(abx.mat.long),abx.mat.long,col=cols.a, xlab='Time', ylab='Bed No.', main="Long duration of antibiotics")
    
    #number of patients with R > threshold for transmission
    thres.s<-rowSums(colo_table_filled_short[[2]]>r_thres)
    thres.l<-rowSums(colo_table_filled_long[[2]]>r_thres)
    x <- c(1:n.day)
    thres.plot.short<-plot(x, thres.s, type="l", ylim=c(0,max(max(thres.s),max(thres.l))), xlab="Time", ylab="Number of patients with \nresistant organisms more than \nthreshold of transmission")
    thres.plot.long<-plot(x, thres.l, type="l", ylim=c(0,max(max(thres.s),max(thres.l))), xlab="Time", ylab="Number of patients with \nresistant organisms more than \nthreshold of transmission")
    
    #total burden of resistance 
    abs.r.s<-exp(colo_table_filled_short[[2]])
    abs.r.s.sum<- rowSums(abs.r.s)
    abs.r.l<-exp(colo_table_filled_long[[2]])
    abs.r.l.sum<- rowSums(abs.r.l)
    
    x <- c(1:n.day)
    r.plot.short<-plot(x, abs.r.s.sum, type="l", ylim=c(0,max(max(abs.r.l.sum),max(abs.r.s.sum))), xlab="Time", ylab="Total burden of resistant organisms in ward")
    r.plot.long<-plot(x, abs.r.l.sum, type="l", ylim=c(0,max(max(abs.r.l.sum),max(abs.r.s.sum))), xlab="Time", ylab="Total burden of resistant organisms in ward")
    
    #difference in burden of resistance between short and long duration 
    diff.abs<-abs.r.l.sum-abs.r.s.sum
    diff.plot<-plot(x, diff.abs, type="l", xlab="Time", ylab="Difference in burden of \nresistant organisms in ward between \nlong and short duration")
    abline(h=0, lty = 3, col='red')
    
    # #carriage differentiation pixel
    # orig_short <- colo_table_filled_short
    # colo_table_filled_short <- diff(colo_table_filled_short)
    # dim.saved <- dim(colo_table_filled_short)
    # colo_table_filled_short[colo_table_filled_short > 0] <- 1
    # colo_table_filled_short[colo_table_filled_short <= 0] <- 2
    # colo_table_filled_short <- as.numeric(colo_table_filled_short)
    # dim(colo_table_filled_short)<-dim.saved
    # 
    # orig_long <- colo_table_filled_long
    # colo_table_filled_long <- diff(colo_table_filled_long)
    # colo_table_filled_long[colo_table_filled_long > 0]<-1
    # colo_table_filled_long[colo_table_filled_long <= 0]<-2
    # colo_table_filled_long<-as.numeric(colo_table_filled_long)
    # dim(colo_table_filled_long)<-dim.saved
    # 
    # cols<-c('chartreuse3','red')
    # col_img.short<-image(1:nrow(colo_table_filled_short),1:ncol(colo_table_filled_short), colo_table_filled_short, col=cols, xlab='Time', ylab='Bed No.')
    # col_img.long<-image(1:nrow(colo_table_filled_long),1:ncol(colo_table_filled_long),colo_table_filled_long,col=cols, xlab='Time', ylab='Bed No.')
    # 
    # #carriage pixel
    # # colo_table_filled_short <- orig_short
    # # 
    # # colo_table_filled_long <- orig_long
    # # 
    # # colfunc <- colorRampPalette(c("black", "white"))
    # # cols <- colfunc(0, )
    # # col_img.short<-image(1:nrow(colo_table_filled_short),1:ncol(colo_table_filled_short), colo_table_filled_short, col=cols, xlab='Time', ylab='Bed No.')
    # # col_img.long<-image(1:nrow(colo_table_filled_long),1:ncol(colo_table_filled_long),colo_table_filled_long,col=cols, xlab='Time', ylab='Bed No.')
    # 
    # #increase overtime
    # df.short<-as.data.frame(orig_short)
    # df.short$totalR<-rowSums(df.short)
    # 
    # df.long<-as.data.frame(orig_long)
    # df.long$totalR<-rowSums(df.long)
    # 
    # x <- c(1:n.day)
    # y.short <- df.short$totalR/n.bed
    # y.long <- df.long$totalR/n.bed
    # r.plot.short<-plot(x,y.short, type="l", ylim=c(0,max(y.short|y.long)), xlab="Time", ylab="% patients carrying MDRO")
    # #lines(x,y.long, col=2)
    # r.plot.long<-plot(x, y.long, type="l", ylim=c(0,max(y.short|y.long)), xlab="Time", ylab="% patients carrying MDRO")
    # 
    # # difference between both plot
    # y.diff <- y.long-y.short
    # r.plot.long<-plot(x, y.diff, type="l", ylim=c(-max(y.short|y.long),max(y.short|y.long)), xlab="Time", ylab="% patients carrying MDRO")
    # abline(h=0)
}

