whole_model <- function(n.bed, n.days, mean.max.los, p.s, p.r, mean_dur,
                        prob_StartBact, pi_r1, pi_r2, mu1, mu2, abx.r, abx.s,
                        repop.r1, repop.r2, repop.r3, repop.s1, repop.s2, repop.s3,
                        iterations=10, short_dur, long_dur){
    
    iter_totalsR <- vector()
    for(iter in 1:iterations){
        
        #print(paste("iter:", iter, "y:", y_count, '-', y, "x", x_count, '-', x))
        #Generate length of stay and antibiotic duration table
        abx_iter <- abx.table(n.bed=n.bed, n.days=n.days, mean.max.los=mean.max.los, p.s=p.s, p.r=p.r, meanDur=short_dur)
        #Generate baseline carriage status
        array_LOS_iter <- array_LOS(los_duration=abx_iter[2])
        #Update values for every day
        array_StartBact_iter <- gen_StartBact(los=array_LOS_iter, prob_StartBact)
        #output
        colo_table_filled_iter <- nextDay(bed_table= abx_iter[[1]], array_LOS=array_LOS_iter, 
                                                  treat_table=abx_iter[[3]], colo_table=array_StartBact_iter, 
                                                  pi_r1=pi_r1, pi_r2= pi_r2, mu1=mu1, mu2=mu2, 
                                                  abx.r=abx.r,abx.s=abx.s,
                                                  repop.r1 = repop.r1, repop.r2 = repop.r2, repop.r3 = repop.r3, 
                                                  repop.s1 = repop.s1, repop.s2 = repop.s2,repop.s3 = repop.s3)
        #Summary
        df <- colo_table_filled_iter
        iter_totalsR <- rowSums(df == "sR")
        #print("end iteration loop")
    }
    totalsR_short <- mean(rowSums(iter_totalsR)/iterations/n.bed)
    
    iter_totalsR <- vector()
    for(iter in 1:iterations){
        
        #print(paste("iter:", iter, "y:", y_count, '-', y, "x", x_count, '-', x))
        #Generate length of stay and antibiotic duration table
        abx_iter <- abx.table(n.bed=n.bed, n.days=n.days, mean.max.los=mean.max.los, p.s=p.s, p.r=p.r, meanDur=long_dur)
        #Generate baseline carriage status
        array_LOS_iter <- array_LOS(los_duration=abx_iter[2])
        #Update values for every day
        array_StartBact_iter <- gen_StartBact(los=array_LOS_iter, prob_StartBact)
        #output
        colo_table_filled_iter <- nextDay(bed_table= abx_iter[[1]], array_LOS=array_LOS_iter, 
                                          treat_table=abx_iter[[3]], colo_table=array_StartBact_iter, 
                                          pi_r1=pi_r1, pi_r2= pi_r2, mu1=mu1, mu2=mu2, 
                                          abx.r=abx.r,abx.s=abx.s,
                                          repop.r1 = repop.r1, repop.r2 = repop.r2, repop.r3 = repop.r3, 
                                          repop.s1 = repop.s1, repop.s2 = repop.s2,repop.s3 = repop.s3)
        #Summary
        df <- colo_table_filled_iter
        iter_totalsR <- rowSums(df == "sR")
        #print("end iteration loop")
    }
    totalsR_long <- mean(rowSums(iter_totalsR)/iterations/n.bed)
    
    return(totalsR_long - totalsR_short)
}