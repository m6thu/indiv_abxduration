#########################################################################
#######Effect of antibiotic duration on hospitalised patients############
#########################################################################

model='frequency'

source(paste0("model_",model,".R"))
source('default_params.R')

if (model=='simple') {
    
    parameters <- list(
        c("qunif", list(min=3, max=50), "n.bed"),             #"n.bed", number of beds in the ward
        c("qunif", list(min=3, max=30), "mean.max.los"),      #"mean.max.los", mean of length of stay
        c("qunif", list(min=0, max=1), "prob_StartBact_R"),    #"prob_StartBact_R",probability of initial carriage of resistant organisms
        c("qunif", list(min=0, max=1), "prop_S_nonR"),         #"prop_S_nonR", proportion of S in the population of S and ss
        c("qunif", list(min=0, max=1), "bif"),                 #"bif", bacterial interference factor
        c("qunif", list(min=0, max=0.1), "pi_ssr"),            # "pi_ssr" probability of being transmitted r to ss (ss—> ssr)
        c("qunif", list(min=0, max=0.05), "repop.s1"),         # "repop.s1" probability of ss repopulated to S (Palleja, Nature Biology, 2018 on gut recovery ~9 months)
        c("qunif", list(min=0, max=0.05), "mu_r"),             # "mu_r", probability of decolonisation (Haggai Bar-Yoseph, JAC, 2016, decreasing colonization rates from 76.7% (95% CI=69.3%–82.8%) at 1 month to 35.2% (95% CI=28.2%–42.9%) at 12 months of follow-up)
        c("qunif", list(min=0.1, max=0.7), "abx.s"),           # "abx.s", probability of S becoming ss after being on narrow spectrum antibiotics
        c("qunif", list(min=0, max=0.7), "abx.r"),             # "abx.r", probability of R becoming ss after being on broad spectrum antibiotics
        c("qunif", list(min=0.1, max=0.9), "p.infect"),        # "p.infect", probability of being prescribed antibiotics
        c("qunif", list(min=10, max=10000), "cum.r.1"),        # admission day when cummulative prabability of HAI requiring abx.r is 1
        c("qunif", list(min=0.1, max=0.9), "p.r.day1"),        #probability of being prescribed broad spectrum antibiotic on admission 
        c("qunif", list(min=5, max=7), "meanDur")              #  mean duration of antibiotics (normal distribution)
        )
    
    prevalence <- function(n.bed, mean.max.los, 
                           prob_StartBact_R, prop_S_nonR, 
                           bif, pi_ssr, repop.s1, mu_r, abx.s, abx.r,
                           p.infect, cum.r.1, p.r.day1, meanDur){
        
        old = Sys.time() # get start time
        # DEBUG
        print(paste(n.bed, mean.max.los, 
                    prob_StartBact_R, prop_S_nonR, 
                    bif, pi_ssr, repop.s1, mu_r, abx.s, abx.r,
                    p.infect, cum.r.1, p.r.day1, meanDur))
        
        timestep = 10
        iterations = 10
        n.day=350
        sdDur=1
        
        iter_totalR = matrix(NA, nrow = n.day*timestep, ncol = iterations)
        
        for(iter in 1:iterations){
            
            matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, 
                                     p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                                     meanDur=meanDur, timestep=timestep)
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
        # Discard first 1/7 runs as burn-in
        totalR = mean(rowSums(iter_totalR[ceiling(n.day*1/7):nrow(iter_totalR), ,drop=FALSE])/iterations)
        
        # print elapsed time
        new = Sys.time() - old # calculate difference
        print(new) # print in nice format
        
        return(totalR)
    }
    
    run.model <- function (data.df) { #data.df is a dataframe of the parameter values in columns 
        return(mapply(prevalence, 
                      data.df[,1], data.df[,2], data.df[,3], 
                      data.df[,4], data.df[,5], data.df[,6], 
                      data.df[,7], data.df[,8], data.df[,9], 
                      data.df[,10], data.df[,11], data.df[,12], 
                      data.df[,13],  data.df[,14]
        ))
    }
    
} else if (model=='binary') {

    parameters <- list(
        c("qunif", list(min=3, max=50), "n.bed"),              #n.bed; number of beds in the ward
        c("qunif", list(min=3, max=30), "mean.max.los"),       #mean.max.los; mean of length of stay (exponential distribution)
        c("qunif", list(min=0, max=1), "prob_StartBact_R"),    #probability of initial carriage of resistant organisms
        c("qunif", list(min=0, max=1), "prop_S_nonR"),         #proportion of S in (S+s): prob_start_S <- prop_S_nonR*(1-prob_StartBact_R)
        c("qunif", list(min=0, max=1), "prop_Sr_inR"),         #proportion of Sr in (r+R): prob_start_Sr <- prop_Sr_inR*prob_StartBact_R
        c("qunif", list(min=0, max=1), "prop_sr_inR"),         #proportion of sr in (r+r): prob_start_sr <- prop_sr_inR*prob_StartBact_R
        c("qunif", list(min=0, max=1), "bif"),                 #bacterial interference factor (pi_ssr = pi_ssr1 * bif )
        c("qunif", list(min=0, max=0.1), "pi_ssr"),             #probability of being transmitted r to ss (ss—> ssr)
        c("qunif", list(min=0, max=0.05), "repop.s1"),         #probability of regrowth of S  (s—>S)
        c("qunif", list(min=0, max=0.05), "repop.s2"),         #probability of regrowth of S  (sr—>Sr)
        c("qunif", list(min=0, max=0.05), "repop.r1"),         #probability of regrowth of s (sr—> sR)
        c("qunif", list(min=0, max=0.05), "repop.r2"),         #probability of regrowth of s (sr—> sR)
        c("qunif", list(min=0, max=0.05), "mu1"),              #probability of being decolonised to S (Sr—> S) 
        c("qunif", list(min=0, max=0.05), "mu2"),              #probability of being decolonised to S (sr—> s) 
        c("qunif", list(min=0, max=0.05), "mu_r"),             #probability of being decolonised to S (Sr—> S) 
        c("qunif", list(min=0.1, max=0.6), "abx.s"),           #probability of clearing S to become s
        c("qunif", list(min=0, max=0.0001), "abx.r"),        #probability of clearing R to become r
        c("qunif", list(min=0.1, max=0.9), "p.infect"),        #probability of being prescribed narrow spectrum antibiotic
        c("qunif", list(min=10, max=10000), "cum.r.1"),        #admission day when cummulative prabability of HAI requiring abx.r is 1
        c("qunif", list(min=0.1, max=0.9), "p.r.day1"),        #probability of being prescribed broad spectrum antibiotic on day 1 of admission 
        c("qunif", list(min=5, max=5), "meanDur")           #mean  duration of antibiotics (normal distribution) 
        )
    
    prevalence <- function(n.bed, mean.max.los, 
                           prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR,
                           bif, pi_ssr, repop.s1, repop.s2, repop.r1, repop.r2,
                           mu1, mu2, mu_r, abx.s, abx.r, 
                           p.infect, cum.r.1, p.r.day1, meanDur){
        
        old = Sys.time() # get start time
        # DEBUG
        print(paste(n.bed, mean.max.los, 
                    prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR,
                    bif, pi_ssr, repop.s1, repop.s2, repop.r1, repop.r2,
                    mu1, mu2, mu_r, abx.s, abx.r, 
                    p.infect, cum.r.1, p.r.day1, meanDur))
        
        timestep = 1
        n.day = 350
        iterations = 1
        
        iter_totalsR = matrix(NA, nrow = n.day*timestep, ncol = iterations)
        iter_totalr_or_R= matrix(NA, nrow = n.day*timestep, ncol = iterations)
        for(iter in 1:iterations){
            matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, 
                                     p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                                     meanDur= meanDur, timestep=timestep)
            patient.matrix=matrixes[[1]]
            abx.matrix=matrixes[[2]]
            los.array = summary.los(patient.matrix=patient.matrix)
            colo.matrix = colo.table(patient.matrix=patient.matrix, los=los.array, 
                                     prob_StartBact_R=prob_StartBact_R, prop_S_nonR=prop_S_nonR, prop_Sr_inR=prop_Sr_inR, prop_sr_inR=prop_sr_inR)
            
            colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                              pi_ssr=pi_ssr, bif=bif, mu1=mu1, mu2=mu2, mu_r=mu_r, repop.r1=repop.r1, repop.r2=repop.r2,
                                              repop.s1=repop.s1, repop.s2=repop.s2, abx.r=abx.r, abx.s=abx.s, timestep=timestep)
            
            #Summary
            df = data.frame(colo.matrix_filled_iter)
            iter_totalsR[, iter] = rowSums(df == "sR")
            iter_totalr_or_R[, iter] = rowSums(df == "sR" | df == "sr" |df == "Sr")
            #print("end iteration loop")
        }
        totalsR = mean(rowSums(iter_totalsR[ceiling(n.day*1/7):nrow(iter_totalsR),, drop=FALSE])/iterations)
        totalr_or_R = mean(rowSums(iter_totalr_or_R[ceiling(n.day*1/7):nrow(iter_totalr_or_R),, drop=FALSE])/iterations)
    
        # print elapsed time
        new = Sys.time() - old # calculate difference
        print(new) # print in nice format
        
        return(c(totalsR, totalr_or_R))
        
    }
    
    run.model<- function (data.df) { #data.df is a dataframe of the parameter values in columns 
        return(mapply(prevalence, 
                      data.df[,1], data.df[,2], data.df[,3], data.df[,4], data.df[,5], 
                      data.df[,6], data.df[,7], data.df[,8], data.df[,9], 
                      data.df[,10], data.df[,11], data.df[,12], data.df[,13], data.df[,14], data.df[,15], 
                      data.df[,16], data.df[,17], data.df[,18],data.df[,19],
                      data.df[,20], data.df[,21]
        ))
    }
    
} else if (model=='frequency') {
    
    parameters <- list(
        c("qunif", list(min=3, max=50), "n.bed"),              #n.bed; number of beds in the ward
        c("qunif", list(min=3, max=30), "mean.max.los"),       #mean.max.los; mean of length of stay (normal distribution)
        c("qunif", list(min=0.1, max=0.9), "p.infect"),        #probability of being prescribed narrow spectrum antibiotic
        c("qunif", list(min=10, max=10000), "cum.r.1"),        #admission day when cummulative prabability of HAI requiring abx.r is 1
        c("qunif", list(min=0.1, max=0.9), "p.r.day1"),          #probability of being prescribed broad spectrum antibiotic on day 1 of admission 
        c("qunif", list(min=3, max=25), "K"),                  # gut holding capacity, on log scale, largest R number possible is exp(300) - typical colonic bacteria 10^14 number/mL content https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4991899/
        c("qunif", list(min=0.1, max=0.9), "total_prop"),             # mean of total starting amount of e coli on log scale
        c("qunif", list(min=0.1,max=0.9), "r_prop"),              # mean of starting amount of resistant gut bacteria on log scale
        c("qunif", list(min=0,max=0.05), "pi_ssr"),              # pi_ssr = daily probability of transmitting resistant E coli
        c("qunif", list(min=1,max=20), "r_thres"),             # r_thres = R threshold level for tranmissibility
        c("qunif", list(min=0.1,max=5), "r_growth"),           # r_growth = growth constant for logistic growth
        c("qunif", list(min=1,max=10), "r_trans"),             # r_trans = amount transmitted on log scale
        c("qunif", list(min=0.1,max=5), "s_growth"),           # s_growth = growth constant for logistic growth
        c("qunif", list(min=1,max=20), "abx.s"),               # abxr_killr = amount of r killed by broad spectrum abx r
        c("qunif", list(min=1,max=20), "abx.r"),               # abxr_kills = amount of s killed by broad spectrum abx r
        c("qunif", list(min=5, max=5), "meanDur")           #mean  duration of narrow spectrum antibiotics (normal distribution) 
     )
    
    prevalence <- function(n.bed, mean.max.los, p.infect, cum.r.1, p.r.day1,
                                K, total_prop, r_prop, pi_ssr, r_thres, r_growth, r_trans, s_growth,
                                abx.s, abx.r, meanDur){
        
        old = Sys.time() # get start time
        # DEBUG
        print(paste(n.bed, mean.max.los, p.infect, cum.r.1, p.r.day1,
                    K, total_prop, r_prop, pi_ssr, r_thres, r_growth, r_trans, s_growth,
                    abx.s, abx.r, meanDur))
        
        timestep = 1
        n.day = 350
        iterations = 1
        
        iter_totalR.no = matrix(NA, nrow = n.day*timestep, ncol = iterations)
        iter_totalR.thres = matrix(NA, nrow = n.day*timestep, ncol = iterations)
        
        for(iter in 1:iterations){
            
            matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, 
                                     p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                                     meanDur= meanDur, timestep=timestep)
            patient.matrix=matrixes[[1]]
            abx.matrix=matrixes[[2]]
            los.array = summary.los(patient.matrix=patient.matrix)
            colo.matrix = colo.table(patient.matrix=patient.matrix, los.array=los.array, total_prop=total_prop, r_prop=r_prop,K=K)
            
            colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                              pi_ssr=pi_ssr, K=K, r_thres=r_thres, r_growth=r_growth, r_trans=r_trans, s_growth=s_growth,
                                              abx.s=abx.s, abx.r=abx.r, timestep=timestep)
            # Summary
            df.R = data.frame(colo.matrix_filled_iter[[2]])
            #print(df.R)
            iter_totalR.no[, iter] = rowSums(df.R)
            
            #for number of people who reached R threshold on a day
            iter_totalR.thres[, iter]=rowSums(df.R >= r_thres)
            #print("end iteration loop")
        }
        totalR_no = mean(rowSums(iter_totalR.no)/iterations)
        totalR_thres = mean(rowSums(iter_totalR.thres)/iterations)
        
        
        # print elapsed time
        new = Sys.time() - old # calculate difference
        print(new) # print in nice format
        
        return(c(totalR_no, totalR_thres))
    }
    
    run.model <- function (data.df) { #data.df is a dataframe of the parameter values in columns 
        return(mapply(prevalence, 
                      data.df[,1], data.df[,2], data.df[,3], 
                      data.df[,4], data.df[,5], data.df[,6], 
                      data.df[,7], data.df[,8], data.df[,9], 
                      data.df[,10], data.df[,11], data.df[,12], 
                      data.df[,13], data.df[,14], data.df[,15], 
                      data.df[,16]
        ))
    }
    
    
}

n.parameters=length(parameters)
parameters_names= unlist(parameters)[seq(4,length(unlist(parameters)), by=4)]

#PLOT 
N.iter=200 #number of iterations 
samples.matrix= as.data.frame(matrix(NA, nrow = N.iter*3, ncol = n.parameters))
colnames(samples.matrix) =c(parameters_names)
samples.matrix$meanDur=rep(c(short_dur,long_dur,long_dur), each=N.iter)
for (i in 1:(n.parameters-1)){ #fill up matrix
    samples.matrix[,which(colnames(samples.matrix)==parameters[[i]][[4]])]=get(parameters[[i]][[4]])
}

#replace parameter values with min and max values 
para='p.infect' #INSERT PARAMETER OF INTEREST HERE
min.para= as.numeric(unlist(parameters)[which(unlist(parameters)==para)-2])
max.para= as.numeric(unlist(parameters)[which(unlist(parameters)==para)-1])
samples.matrix[,which(colnames(samples.matrix)==para)]=rep(c(min.para,min.para,max.para),each=N.iter)
samples.matrix$abx.r=10
samples.matrix$abx.s=10
out=run.model(data.df=samples.matrix)
dat=cbind.data.frame(samples.matrix, out[1,]) #dataframe with columns of parameter values + outcome values 
colnames(dat) =c(parameters_names,'Y')

####plot 
require(ggplot2)
require(hexbin) 

dat$Time=as.factor(ifelse(dat$meanDur<10, 'Short duration', 'Long duration'))
dat$Time = with(dat, factor(Time, levels = rev(levels(Time))))

# calculating mean
meanY= c(mean(dat[1:N.iter,]$Y), mean(dat[(N.iter+1):(N.iter*2),]$Y), mean(dat[(N.iter*2+1):(N.iter*3),]$Y))
means1=data.frame(Time=c('Short duration', 'Long duration'),
                 meanY=meanY[1:2])
means2=data.frame(Time=c('Short duration', 'Long duration'),
                  meanY=meanY[c(1,3)])

ggplot(dat,aes(x=Time,y=Y)) + 
    stat_binhex() + 
    scale_fill_gradientn(colours=c("#16F4D0","#2F9C95"),name = "Frequency",na.value=NA) + 
    # ploting mean dots
    geom_point (aes(Time, meanY), data = means1, pch = 19, col = "#C44536", cex = 3) +
    geom_point (aes(Time, meanY), data = means2, pch = 19, col = "#C44536", cex = 3) + 
    # connecting with line  
    geom_line (aes(x=means1$Time,y=means1$meanY), data = means1, col = "#C44536", lwd = 0.5, group = 1) +
    geom_line (aes(x=means2$Time,y=means2$meanY), data = means2, col = "#C44536", lwd = 0.5, group = 1) +
    #label 
    annotate(geom="text", x=2.3, y=meanY[2], label=paste(para, '=', min.para))+
    annotate(geom="text", x=2.3, y=meanY[3], label=paste(para, '=', max.para))+
    annotate(geom="text", x=0.6, y=meanY[1], label=paste(para, '=', min.para))+
    ylab("Number of patients carrying R in 30-bed ward/month")+
    theme_minimal() 

