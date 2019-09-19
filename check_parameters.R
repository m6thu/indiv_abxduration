#############check parameter space#############
###functions of cumulative risks of events ###
rm(list = ls(all.names = TRUE))
seq_along_admissiondates=function(x){
  return(1:x)
}
seq2 <- Vectorize(seq_along_admissiondates, vectorize.args = "x")

plotlos<-function(n.patient, min, max){
  
  los=runif(200, min=min, max=max)
  all_los = ceiling(rexp(n.patient, 1/max))
  all_los[all_los > 5*max] = max
  dens=density(all_los)
  plot=plot(dens$x,dens$y,type="l", col='blue', ylim=c(0,0.26),
            xlab="Length of stay",ylab="Proportion of patients")
  
  for (iter in 1:length(los)){
  all_los = ceiling(rexp(n.patient, 1/los[iter]))
  all_los[all_los > 5*los[iter]] = los[iter]
  dens=density(all_los)
  lines(x=dens$x, y=dens$y, type = 'l',col='grey50')
  }
  
  all_los = ceiling(rexp(n.patient, 1/min))
  all_los[all_los > 5*los] = min
  dens=density(all_los)
  lines(x=dens$x, y=dens$y, type = 'l',col='red')
  
  all_los = ceiling(rexp(n.patient, 1/max))
  all_los[all_los > 5*max] = max
  dens=density(all_los)
  lines(x=dens$x, y=dens$y, type = 'l',col='blue')
  
  text(60, 0.1, paste0("Red line highlights max.los=",min,'(min)'), cex = .8)
  text(60, 0.12, paste0("Blue line highlights max.los=",max,'(max)'), cex = .8)
}

plotpara<-function(parameter,max, min, nday=30){
  
  ####General parameters
  iter=200
  timesteps=c(3,1)
  models=c('Models 1 and 2','Model 3')
  
  #max value of the parameter 
  r_num=0:50
  
  if (parameter=='pi_ssr') {
    par(mfrow=c(1,2))

    for (type in 1:length(models)) {
      
      #For max value 
      day.risk=list()
      timestep=timesteps[type]
      for (i in 1:length(r_num)){
        pi_ssr = 1-(1-max)^(1/timestep)# max pi_ssr 
        #per timestep, r_num=10, chance of transmission 
        prop_R = 1-((1-pi_ssr)^r_num[i])
        day.risk[[i]]=cumsum((1-prop_R)^(1:max(nday*timestep)-1)*prop_R)
      }
      day.risk.df=data.frame(matrix(unlist(day.risk), ncol=length(day.risk)))
      
      plot(x=1:max(nday*timestep), y=day.risk.df[,2], ylim=c(0:1), type = 'l', 
           col='grey50',
           xaxt='n',
           ylab='Cumulative daily risk of a patient being transmitted R', 
           xlab=paste0('Number of days (timestep=', timestep,')'),
           main=paste0('pi_ssr (', models[type],')'), 
           sub=paste0('max=',max,'  min=',min))
      axis(1, at = seq(0, max(nday*timestep), length=7), las=2, 
           labels=seq(0, max(nday*timestep), length=7)/timestep)
      for (i in 3:ncol(day.risk.df)){
        lines(x=1:max(nday*timestep), y=day.risk.df[,i], type = 'l',col='grey50')
      }
      #label red for the line with r_num=40
      lines(x=1:max(nday*timestep), y=day.risk.df[,41], type = 'l', col='red')
      
      #For min value 
      day.risk=list()
      timestep=timesteps[type]
      for (i in 1:length(r_num)){
        pi_ssr = 1-(1-min)^(1/timestep)# max pi_ssr 
        #per timestep, r_num=10, chance of transmission 
        prop_R = 1-((1-pi_ssr)^r_num[i])
        day.risk[[i]]=cumsum((1-prop_R)^(1:max(nday*timestep)-1)*prop_R)
      }
      day.risk.df=data.frame(matrix(unlist(day.risk), ncol=length(day.risk)))
      
      lines(x=1:max(nday*timestep), y=day.risk.df[,2], type = 'l', 
           col='grey50')
      for (i in 3:ncol(day.risk.df)){
        lines(x=1:max(nday*timestep), y=day.risk.df[,i], type = 'l',col='grey50')
      }
      #label red for the line with r_num=40
      lines(x=1:max(nday*timestep), y=day.risk.df[,41], type = 'l', col='red')
      
    }
    text(250, 0.2, "Each line represents number of R in the ward (0-50)",cex = .8)
    text(250, 0.15, "Red line highlights R=40", cex = .8)
    
  } else {
    par(mfrow=c(1,1))
    timestep=3
    
    day.risk=list()
    
    for (i in 1:iter){
      para=runif(1, min=min, max=max)
      para = 1-(1- para)^(1/timestep)# max 
      day.risk[[i]]=cumsum((1-para)^(1:max(nday*timestep)-1)*para)
    }
    day.risk.df=data.frame(matrix(unlist(day.risk), ncol=length(day.risk)))
    
    plot(x=1:max(nday*timestep), y=day.risk.df[,2], 
         ylim=c(min(day.risk.df), max(day.risk.df)), type = 'l', col='grey50',
         xaxt='n',
         ylab='Cumulative daily risk', 
         xlab=paste0('Number of days (timestep=', timestep,')'),
         main=paste0(parameter,' (Models 1 and 2)'), 
         sub=paste0('max=',max,'  min=',min))
    axis(1, at = seq(0, max(nday*timestep), length=7), las=2, 
         labels=seq(0, max(nday*timestep), length=7)/timestep)
    for (i in 3:ncol(day.risk.df)){
      lines(x=1:max(nday*timestep), y=day.risk.df[,i], type = 'l',col='grey50')
    }
    text(nday*2, min(day.risk.df)-0.01, "Each line represents a single iteration",
         cex = .8)
  }
}

plot.cum.r.1=function(min, max,nday){
  
  par(mfrow=c(1,1))
  
  ####General parameters
  iter=1000
  all_los=nday
  
  prob.r=list()
  
  for (i in 1:iter){
    cum.r.1=runif(1, min, max)
    all_admission_days = seq2(all_los) #patients' admission dates in sequence
    probs=rnorm(10000) #randomly draw 10000 probabilities from a normal distribution
    probs.normalized = (probs - min(probs))/(max(probs)- min(probs)) #normalised probabilities
    p=ecdf(probs.normalized) #cumulative distribution ##SLOW
    
    prob.r.after= p((1/cum.r.1)*all_admission_days)
    prob.r[[i]]= c(0,0,prob.r.after [-c(1,2)])#abx r only can start after 48h
  }
  
  day.risk.df=data.frame(matrix(unlist(prob.r), ncol=length(prob.r)))
  
  plot(x=1:max(all_los), y=day.risk.df[,2], 
       ylim=c(min(day.risk.df),max(day.risk.df)), type = 'l', col='grey50',
       xaxt='n',
       ylab='Cumulative daily risk', 
       xlab=paste0('Number of days'),
       main=paste0('cum.r.1'), 
       sub=paste0('max=',max,'  min=',min))
  axis(1, at = seq(0, max(all_los), length=7), las=2, 
       labels=seq(0, max(all_los), length=7))
  for (i in 3:ncol(day.risk.df)){
    lines(x=1:max(all_los), y=day.risk.df[,i], type = 'l',col='grey50')
  }
  text(250, 1, "Each line represents a single iteration",
       cex = .8)
}

###########Frequency model ###########
setwd('/Users/moyin/Documents/nBox/git_projects/indiv_abxduration/')
source('model_frequency.R')
source('los_abx_matrix.R')
model='Frequency'
source('default_params.R')

###check growth under no antibiotics, r_growth, s_growth
growth.no.abx=function(min, max, k, total_prop, n.day=n.day){

  total_existing= log(total_prop*exp(k))
  
  growlist=list()
  iter=500
  
  for (g in 1:iter){
    
    growth=runif(1, min, max)
    total=c()
    bact=total_existing
    
    for (i in 1:n.day){

      total[i]=bact
      
      # calculate effect of logistic bacteria growth (abs)
      grow = growth*exp(bact)*(1 - (exp(bact)/exp(k))) #assumes no R 
      
      # apply effects to current number
      bact = log(exp(bact) + grow)
    }
    growlist[[g]]= exp(total)
  }
  
  day.risk.df=data.frame(matrix(unlist(growlist), ncol=iter))
  
  par(mfrow=c(1,1))
  plot(x=1:n.day, y=day.risk.df[,2], 
       ylim=c(min(day.risk.df), max(day.risk.df)), type = 'l', col='grey50',
       xaxt='n',
       ylab='Cumulative growth', 
       xlab=paste0('Number of days'),
       main=paste0('s_growth and r_growth (Model 3)'), 
       sub=paste0('max=',max,'  min=',min))
  axis(1, at = seq(0, max(n.day), length=7), las=2, 
       labels=seq(0, max(n.day), length=7))
  for (i in 3:ncol(day.risk.df)){
    lines(x=1:max(n.day), y=day.risk.df[,i], type = 'l',col='grey50')
  }
  text(200, min(day.risk.df)*1.05, "Each line represents a single iteration",
       cex = .8)
  text(200, min(day.risk.df)*1.01, paste0("Given K=22, total_prop=",total_prop),
       cex = .8)
}

#Shaw ISME 2019- microbiome recovery table 3 and 4 
d=data.frame(antibiotic=c('ciprofloxacin', 'clindamycin','amoxicillin','minocycline'), 
             Dmin=c(5.47,6.23,0.13,1.54), 
             Dmax=c(9.75,9.84,6.56,7.82),
             Amin=c(0.28,0.29,-0.66,-1.09),
             Amax=c(1.34,1.42,0.56,0.23), 
             phi1min=c(-0.69,-0.46,-1.96,-1.44),
             phi1max=c(0.16,0.34,0.31,1.29),
             phi2min=c(0.05,0.23,-1.58,1.01),
             phi2max=c(0.92,1.11,1.83,1.97))
t=1:100
D=1.34  #values in paper
A=-.03
phi1=-1.53 
phi2=0.09
#equation 3, figure 1c
y=((D*exp(phi1)*exp(phi2))/(exp(phi2)-exp(phi1)))*(exp(-exp(phi1*t))-exp(-exp(phi2*t)))+A*(1-exp(-exp(phi1*t)))
plot(t,y) #similar distribution as logistic growth

###check antibiotic killing 
abx.killing=function(growth.min, growth.max, abx.min, abx.max, k, total_prop, n.day=n.day){
  
  total_existing= log(total_prop*exp(K))
  
  growlist=list()
  iter=500
  
  for (g in 1:iter){
    
    growth=runif(1, min=growth.min, max=growth.max)
    abx=runif(1, min=abx.min, max=abx.max)
    total=c()
    bact=total_existing
    
    for (i in 1:n.day){

      total[i]=bact
      
      # calculate effect of logistic bacteria growth (abs)
      grow = growth*exp(bact)*(1 - (exp(bact)/exp(K)))
      
      # add effect of abx death 
      bact_abx = -abx*grow
      
      # apply effects to current number
      if (exp(bact) + grow + bact_abx<0) {
        bact=0
        }else{
      bact = log(exp(bact) + grow + bact_abx)
        }
      
    }
    growlist[[g]]= exp(total)
  }
  
  day.risk.df=data.frame(matrix(unlist(growlist), ncol=iter))
  
  par(mfrow=c(1,1))
  plot(x=1:n.day, y=day.risk.df[,2], 
       ylim=c(min(day.risk.df),max(day.risk.df)), type = 'l', col='grey50',
       xaxt='n',
       ylab='Cumulative growth', 
       xlab=paste0('Number of days'),
       main=paste0('growth and abx killing (Model 3)'), 
       sub=paste0('growth(max)=',growth.max,'  growth(min)=',growth.min, 
                  '  abx(max)=',abx.max,'  abx(min)=',abx.min))
  axis(1, at = seq(0, max(n.day), length=7), las=2, 
       labels=seq(0, max(n.day), length=7))
  for (i in 3:ncol(day.risk.df)){
    lines(x=1:max(n.day), y=day.risk.df[,i], type = 'l',col='grey50')
  }
  text(n.day-5, max(day.risk.df)-100, "Each line represents a single iteration",
       cex = .8)
  text(n.day-5, max(day.risk.df)-500, paste0("Given K=22,", ", total_prop=",total_prop),
       cex = .8)
}


####K####
d=read.csv('/Users/moyin/Documents/nBox/git_projects/indiv_abxduration/gutdata/Ini_CTXm_copies_qPCR.csv') #SATURN gut data
mean(d$ini_16S_log) #total capacity in log, K
sd(d$ini_16S_log)
hist(d$ini_16S_log)

#####r_mean#####
##Assuming Enterobacteriaceae 0.1-5% 
enterobacteriaceae.capacity=d$ini_16S_copies*0.05
summary(d$ini_CTXm_copies/enterobacteriaceae.capacity) #proportion of R in total 
#                                                    number of Enterbacterobacteriaceae, r_prop
#                                                    = exponential distribution 
#median 12% 
#png(filename="r_mean_saturn.png", width = 800, height = 350)
hist(d$ini_CTXm_copies/enterobacteriaceae.capacity,
     main='Histogram of CTXM copies/Enterocbacteriaceae copies',
     xlab='',
     breaks =seq(0,105,0.01), xlim = c(0,1))
text(0.5, 25, 'Enterobacteriaceae capacity calculated by 16S*0.05')
#dev.off()

plotrmean<-function(n.patient, min, max){
  
  r_mean=runif(200, min=min, max=max)
  prop.r.ent= rtnorm(n.patient, mean=max, lower=0)
  prop.r.ent.norm= (prop.r.ent-min(prop.r.ent))/(max(prop.r.ent)-min(prop.r.ent))#proportion of R in total population of Enterobacteriaceae (normalised)
  dens=density(prop.r.ent.norm)
  plot=plot(dens$x, dens$y, type="l", col='blue', ylim=c(0,3.5),
            xlab="Proportion of Enterobacteriaceae that is resistant",ylab="Density")
  
  for (iter in 1:length(r_mean)){
    prop.r.ent = rtnorm(n.patient, mean=r_mean[iter], lower=0)
    prop.r.ent.norm= (prop.r.ent-min(prop.r.ent))/(max(prop.r.ent)-min(prop.r.ent))#proportion of R in total population of Enterobacteriaceae (normalised)
    dens=density(prop.r.ent.norm)
    lines(x=dens$x, y=dens$y, type = 'l',col='grey50')
  }

  prop.r.ent= rtnorm(n.patient, mean=max, lower=0)
  prop.r.ent.norm= (prop.r.ent-min(prop.r.ent))/(max(prop.r.ent)-min(prop.r.ent))#proportion of R in total population of Enterobacteriaceae (normalised)
  dens=density(prop.r.ent.norm)
  lines(dens$x, dens$y, type="l", col='blue')
  
  prop.r.ent= rtnorm(n.patient, mean=min, lower=0)
  prop.r.ent.norm= (prop.r.ent-min(prop.r.ent))/(max(prop.r.ent)-min(prop.r.ent))#proportion of R in total population of Enterobacteriaceae (normalised)
  dens=density(prop.r.ent.norm)
  lines(dens$x, dens$y, type="l", col='red')
  
  text(0.6, 2, paste0("Red line highlights r_mean=",min,'(min)'), cex = .8)
  text(0.6, 2.2, paste0("Blue line highlights r_mean=",max,'(max)'), cex = .8)
  
}

###########################################################################
###############################PLOT GRAPHS#################################
###########################################################################
setwd(dir = '/Users/moyin/Desktop/')

#max.los
png(filename="max.los.png", width = 800, height = 350)
plotlos(n.patient=200, min=3, max=20)
dev.off()

#r_mean
png(filename="r_mean.png", width = 800, height = 350)
plotrmean(n.patient=1000, min=0, max=1)
dev.off()

#pi_ssr
png(filename="pi_ssr.png", width = 800, height = 350)
plotpara(parameter='pi_ssr',max=0.002, min=0.000001, nday=720)
dev.off()

#repop.s
png(filename="repop.s.png", width = 800, height = 350)
plotpara(parameter='repop.s',max=0.02, min=0.002, nday=120)
dev.off()

#repop.r
png(filename="repop.r.png", width = 800, height = 350)
plotpara(parameter='repop.r',max=0.05, min=0.01, nday=12)
dev.off()

#mu
png(filename="mu.png", width = 800, height = 350)
plotpara(parameter='mu', max=0.02, min=0.002, nday=360)
dev.off()

png(filename="abx.png", width = 800, height = 350)
plotpara(parameter='abx', max=0.5, min=0.1)
dev.off()

png(filename="cum.r.1.png", width = 800, height = 350)
plot.cum.r.1(min=10, max=1000, nday=60)
dev.off()

png(filename="growth.png", width = 800, height = 350)
growth.no.abx(min=0.01, max=0.1, k=22, total_prop=0.5, n.day=360)
dev.off()

png(filename="kill.png", width = 800, height = 350)
abx.killing(growth.min=0.01, growth.max=0.1, abx.max=15, abx.min=10, 
            k=22, total_prop=0.5, n.day=30)
dev.off()

######timestep
timestep=1
pi_ssr=0.002
pi_ssr = 1-(1-pi_ssr)^(1/timestep)
prop_R = 1-((1-pi_ssr)^50)
abx.r=0.5
prop_R+abx.r

#####Gibson 2018: doubling time of e coli 15h 
# G (generation time) = (time, in minutes or hours)/n(number of generations), G = t/n
# log10(2)/(15/24) #15 hours - a daily growth rate of 115.2 with doubling time 0.625 days
# log10(2)/(25/24) #25 hours - a daily growth rate of 0.29 
# log10(2)/(5/24)  #5 hours - a daily growth rate of 1.44
# 
# t=1:500
# r_growth=0.29
# K=12
# r=200
# s=1000
# r_growth*r*(1 - (r + s)/exp(K)) #growth in 1 day 
