#############check parameter space#############
###functions of cumulative risks of events ###
rm(list = ls(all.names = TRUE))
seq_along_admissiondates=function(x){
  return(1:x)
}
seq2 <- Vectorize(seq_along_admissiondates, vectorize.args = "x")

plotpara<-function(parameter,max, min){
  
  ####General parameters
  iter=200
  nday=1:30
  timesteps=c(3,1)
  models=c('Models 1 and 2','Model 3')
  
  #max value of the parameter 
  r_num=0:50
  
  if (parameter=='pi_ssr') {
    par(mfrow=c(1,2))
    for (type in 1:length(models)) {
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
      #label red for the line with r_num=10
      lines(x=1:max(nday*timestep), y=day.risk.df[,11], type = 'l', col='red')
    }
    text(15, 0.1, "Each line represents number of R in the ward (0-50)",cex = .8)
    text(15, 0.05, "Red line highlights R=10", cex = .8)
    
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
    text(17, min(day.risk.df)-0.01, "Each line represents a single iteration",
         cex = .8)
  }
}

plot.cum.r.1=function(min, max){
  
  par(mfrow=c(1,1))
  
  ####General parameters
  iter=1000
  all_los=30
  
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
  text(5, 1, "Each line represents a single iteration",
       cex = .8)
}

###########Frequency model ###########
setwd('/Users/moyin/Documents/nBox/git_projects/indiv_abxduration/')
source('model_frequency.R')
source('los_abx_matrix.R')
model='Frequency'
source('default_params.R')

###check growth under no antibiotics, r_growth, s_growth
growth.no.abx=function(min, max, k, total_prop, capacity_prop, n.day=n.day){
  
  total_capacity= log(capacity_prop*exp(K)) 
  total_existing= log(total_prop*exp(total_capacity))
  
  growlist=list()
  iter=500
  
  for (k in 1:iter){
    
    growth=runif(1, min, max)
    total=c()
    bact=total_existing
    
    for (i in 1:n.day){

      total[i]=bact
      
      # calculate effect of logistic bacteria growth (abs)
      grow = growth*exp(bact)*(1 - (exp(bact)/exp(total_capacity)))
      
      # apply effects to current number
      bact = log(exp(bact) + grow)
    }
    growlist[[k]]= exp(total)
  }
  
  day.risk.df=data.frame(matrix(unlist(growlist), ncol=iter))
  
  par(mfrow=c(1,1))
  plot(x=1:n.day, y=day.risk.df[,2], 
       ylim=c(min(day.risk.df),max(day.risk.df)), type = 'l', col='grey50',
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
  text(150, min(day.risk.df)+450, "Each line represents a single iteration",
       cex = .8)
  text(150, min(day.risk.df)+100, paste0("Given K=10, capacity_prop=",capacity_prop, ", total_prop=",total_prop),
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
abx.killing=function(growth.min, growth.max, abx.min, abx.max, k, total_prop, capacity_prop, n.day=n.day){
  
  total_capacity= log(capacity_prop*exp(K)) 
  total_existing= log(total_prop*exp(total_capacity))
  
  growlist=list()
  iter=500
  
  for (k in 1:iter){
    
    growth=runif(1, min=growth.min, max=growth.max)
    abx=runif(1, min=abx.min, max=abx.max)
    total=c()
    bact=total_existing
    
    for (i in 1:n.day){

      total[i]=bact
      
      # calculate effect of logistic bacteria growth (abs)
      grow = growth*exp(bact)*(1 - (exp(bact)/exp(total_capacity)))
      
      # add effect of abx death 
      bact_abx = -abx*grow
      
      # apply effects to current number
      if (exp(bact) + grow + bact_abx<0) {
        bact=0
        }else{
      bact = log(exp(bact) + grow + bact_abx)
        }
      
    }
    growlist[[k]]= exp(total)
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
  text(n.day-5, max(day.risk.df)-500, paste0("Given K=10, capacity_prop=",capacity_prop, ", total_prop=",total_prop),
       cex = .8)
}


####K####
d=read.csv('gutdata/Ini_CTXm_copies_qPCR.csv') #SATURN gut data
mean(d$ini_16S_log) #total capacity in log, K
sd(d$ini_16S_log)
hist(d$ini_16S_log)
##Assuming Enterobacteriaceae 0.001-5% 
enterobacteriaceae.capacity=d$ini_16S_copies*0.05
summary(d$ini_CTXm_copies/enterobacteriaceae.capacity) #proportion of R in total 
#                                                    number of Enterbacterobacteriaceae, r_prop
#                                                    = exponential distribution 
#70% of the ratio is less than 1
#median 12% 
hist(d$ini_CTXm_copies/enterobacteriaceae.capacity,breaks =seq(0,105,0.01), xlim = c(0,1))

######PLOT GRAPHS########
setwd(dir = '/Users/moyin/Desktop/')

#repopr=0.02 for ~ 50 days for cummulative risk to become 0.60 -hilty 2012 cid
png(filename="pi_ssr.png", width = 800, height = 350)
plotpara(parameter='pi_ssr',max=0.03, min=0)
dev.off()

png(filename="repop.s.png", width = 800, height = 350)
plotpara(parameter='repop.s',max=0.02, min=0.005)
dev.off()

#returning traveller 25% colonised with esbl after 1-6 weeks of travel
png(filename="repop.r.png", width = 800, height = 350)
plotpara(parameter='repop.r',max=0.05, min=0.005)
dev.off()

#mu_r- Haggai Bar-Yoseph, JAC, 2016, 
##colonization 100 -> 76.7% (95% CI=69.3%–82.8%) 1 month ->35.2% (95% CI=28.2%–42.9%) 12 months 
##decolonized 17.2-30.7 1 month, 57.1-71.8 12 months 
png(filename="mu.png", width = 800, height = 350)
plotpara(parameter='mu',max=0.02, min=0.002)
dev.off()
png(filename="abx.png", width = 800, height = 350)
plotpara(parameter='abx', max=0.5, min=0.1)
dev.off()
png(filename="cum.r.1.png", width = 800, height = 350)
plot.cum.r.1(min=10, max=1000)
dev.off()

png(filename="growth.png", width = 800, height = 350)
growth.no.abx(min=0.01, max=0.5, k=10, total_prop=0.9,capacity_prop=0.1, n.day=18)
dev.off()

png(filename="kill.png", width = 800, height = 350)
abx.killing(growth.min=0.01, growth.max=0.1, abx.max=100, abx.min=20, 
            k=10, total_prop=0.9, capacity_prop=0.1, n.day=18)
dev.off()

#####Gibson 2018: doubling time of e coli 15h 
#G (generation time) = (time, in minutes or hours)/n(number of generations), G = t/n
# log10(2)/(15/24) #15 hours - a daily growth rate of 115.2 with doubling time 0.625 days
# log10(2)/(25/24) #25 hours - a daily growth rate of 0.29 
# log10(2)/(5/24)  #15 hours - a daily growth rate of 1.44
# 
# t=1:500
# r_growth=0.29
# K=12
# total_capacity=log(0.05*(exp(K)))
# r=200
# s=1000
# r_growth*r*(1 - (r + s)/exp(total_capacity)) #growth in 1 day 
