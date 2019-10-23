# Adapted from checkvariability.R

#########################################################################
#######Effect of antibiotic duration on hospitalised patients############
#########################################################################
setwd('/Users/moyin/Documents/nBox/git_projects/indiv_abxduration')
rm(list=ls()) # Clean working environment

library(fields)
source("default_params.R")
model='binary'
source(paste0('model_',model,".R"))

######################### Heatmap exploration of 2 parameters ###########################

pixels=10 #how many pixels per heatmap

## Change parameters here: START
y_name <- "bif"
bif_seq <- seq(0, 1, length.out = pixels)
y_seq <- bif_seq

x_name <- "repop.s"
repop.s_seq <- seq(0, 1, length.out = pixels)
x_seq <- repop.s_seq

# Parameters not added here are considered fixed as defined above in Fixed Parameters
default.para.df= data.frame(para = parameters_diff_prevalence)
for (col in 2: (pixels+1)){
  default.para.df
}

save_x_runs <- list()
save_y_runs <- list()
outcome <- matrix(NA, nrow = length(y_seq), ncol = length(x_seq))

###run model 
if (model == 'simple'){
  
  for (i in x_name) {
    for (j in y_name) {
    
      save_y_runs[j] = diff_prevalence(n.bed=n.bed, max.los=max.los, prop_R=prop_R, prop_S=prop_S, bif=bif, pi_ssr=pissr, repop.s=repop.s, 
                                     mu=mu, abx.s=abx.s, abx.r=abx.r, p.infect=p.infect, cum.r.1=cum.r.1, p.r.day1=p.r.day1, 
                                     short_dur=short_dur, long_dur=long_dur)
    }
    save_runs[[i]] = save_y_runs
  }
  
} if (model == 'binary'){
  
  diff_prevalence(n.bed, max.los, 
                  prop_R, prop_r, prop_Sr, prop_S,
                  bif, pi_ssr, repop.s, repop.r, mu, abx.s, abx.r, 
                  p.infect, cum.r.1, p.r.day1, short_dur, long_dur)
  
} else (model == 'frequency'){
  
  diff_prevalence(n.bed, max.los, p.infect, cum.r.1, p.r.day1,
                  K, total_prop, prop_R, pi_ssr, 
                  r_trans, r_growth, r_thres, s_growth,
                  abx.s, abx.r, short_dur,long_dur)
}

image.plot(data, xlab=x_name, ylab=y_name, axes=F, col=rev(terrain.colors(100)))







i <- 1 # total count
y_count <- 1

for(y in y_seq){
  # calculate parameter dependencies on y paramter here
  pi_r2 <- pi_r1 * y                  # pi_r2= probability of R transmitting to s to become sr 
  #                                     (pi_r1 < pi_r2 if being colonised with S protects colonisation by R)
  x_count <- 1 
  for(x in x_seq){
    # calculate parameter dependencies on y parameter here
    repop.s1 <- x                        # probability of repopulation of s to become S 
    repop.s2 <- x                        # probability of repopulation of sr to become SR 
    repop.s3 <- x                        # probability of repopulation of sR to become sr
    repop.r1 <- x*10                     # probability of repopulation of Sr to become sR 
    repop.r2 <- x*10                     # probability of repopulation of sr to become sR 
    
    # saves for iterations
    abx_iter <- list()
    array_LOS_iter <- list()
    array_StartBact_iter <- list()
    colo_table_filled_iter <- list()
    
    iter_totalS <- matrix(NA, nrow = n.days, ncol = iterations)
    iter_totalsr <- matrix(NA, nrow = n.days, ncol = iterations)
    iter_totalSr <- matrix(NA, nrow = n.days, ncol = iterations)
    iter_totalsR <- matrix(NA, nrow = n.days, ncol = iterations)
    
    for(iter in 1:iterations){
      
      print(paste("iter:", iter, "y:", y_count, '-', y, "x", x_count, '-', x))
      #Generate length of stay and antibiotic duration table
      abx_iter[[iter]] <- abx.table(n.bed=n.bed, n.days=n.days, mean.max.los=mean.max.los, p.s=p.s, p.r.day1=p.r.day1, p.r.dayafter=p.r.dayafter, meanDur=long_dur)
      #Generate baseline carriage status
      array_LOS_iter[[iter]] <- array_LOS_func(los_duration=abx_iter[[iter]][2])
      #Update values for every day
      array_StartBact_iter[[iter]] <- gen_StartBact(los=array_LOS_iter[[iter]],prob_StartBact_R=prob_StartBact_R, 
                                                    prop_S_nonR=prop_S_nonR, prop_sr_inR=prop_sr_inR,
                                                    prop_Sr_inR=prop_Sr_inR, n.bed=n.bed, n.days=n.days)
      #output
      colo_table_filled_iter[[iter]] <- nextDay(bed_table= abx_iter[[iter]][[1]], array_LOS=array_LOS_iter[[iter]], 
                                                treat_table=abx_iter[[iter]][[3]], colo_table=array_StartBact_iter[[iter]], 
                                                pi_r1=pi_r1, pi_r2= pi_r2, mu1=mu1, mu2=mu2, 
                                                abx.r=abx.r,abx.s=abx.s,
                                                repop.r1 = repop.r1, repop.r2 = repop.r2, 
                                                repop.s1 = repop.s1, repop.s2 = repop.s2,repop.s3 = repop.s3)
      
      #Summary
      df <- colo_table_filled_iter[[iter]]
      iter_totalS[,iter] <- rowSums(df == "S")
      iter_totalsr[,iter] <- rowSums(df == "sr")
      iter_totalSr[,iter] <- rowSums(df == "Sr")
      iter_totalsR[,iter] <- rowSums(df == "sR")
      #print("end iteration loop")
    }
    save_runs[[i]] <- list(abx_iter, array_LOS_iter, array_StartBact_iter, colo_table_filled_iter, 
                           iter_totalS, iter_totalsr, iter_totalSr, iter_totalsR)
    # increment index
    i <- i + 1
    
    totalS[y_count, x_count] <- mean(rowSums(iter_totalS)/iterations/n.bed)
    totalsr[y_count, x_count] <- mean(rowSums(iter_totalsr)/iterations/n.bed)
    totalSr[y_count, x_count] <- mean(rowSums(iter_totalSr)/iterations/n.bed)
    totalsR[y_count, x_count] <- mean(rowSums(iter_totalsR)/iterations/n.bed)
    
    x_count <- x_count + 1
    #print("end p.s loop")
  }
  y_count <- y_count + 1
  #print("end pi_r1 loop")
}

# Clean run saved to 1443_2Oct2018.RData
image_name <- paste0(y_name, "_vs_", x_name, "_", format(Sys.time(), "%d%b%Y_%H%M%Z"))
save.image(paste0(image_name, ".Rdata"))

require(fields)

pdf(image_name)
par(mfrow=c(2,2))
image.plot(t(totalS), xlab=x_name, ylab=y_name, axes=F, col=rev(terrain.colors(100)))
title("totalS")
axis(1, at = (1:ncol(totalS)-1)/10, labels=as.character(x_seq))
axis(2, at = (1:nrow(totalS)-1)/10, labels=as.character(y_seq))
image.plot(t(totalsr), xlab=x_name, ylab=y_name, axes=F, col=rev(terrain.colors(100)))
title("totalsr")
axis(1, at = (1:ncol(totalsr)-1)/10, labels=as.character(x_seq))
axis(2, at = (1:nrow(totalsr)-1)/10, labels=as.character(y_seq))
image.plot(t(totalSr), xlab=x_name, ylab=y_name, axes=F, col=rev(terrain.colors(100)))
title("totalSr")
axis(1, at = (1:ncol(totalSr)-1)/10, labels=as.character(x_seq))
axis(2, at = (1:nrow(totalSr)-1)/10, labels=as.character(y_seq))
image.plot(t(totalsR), xlab=x_name, ylab=y_name, axes=F, col=rev(terrain.colors(100)))
title("totalsR")
axis(1, at = (1:ncol(totalsR)-1)/10, labels=as.character(x_seq))
axis(2, at = (1:nrow(totalsR)-1)/10, labels=as.character(y_seq))
dev.off()
