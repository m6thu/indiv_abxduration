#Understanding the exponential function
lambda <- seq(0,1000,0.001)
thresholds <- seq(0,1,0.01)

D <- matrix(NA, nrow=length(lambda), ncol = length(thresholds))
for(i in 1:length(lambda)){
    for(j in 1:length(thresholds)){
        D[i, j] <- pexp(thresholds[j], rate = lambda[i])
    }
}
rownames(D) <- lambda
colnames(D) <- thresholds

explookup <- function(density, thresholds){
    
    if(density > 1 && thresholds == 0){
        stop("Density requested higher than 0 when threshold is 0.")}
    idx <- which.min(abs(density-D[,toString(thresholds)]))
    #debug
    #print(paste("ans: ", density, D[idx,toString(thresholds)]))
    #interpolation between values
    
    lambda <- as.numeric(row.names(D)[idx])
    return(lambda)
}

# Test range of generated densities at all threshold levels
testlookup <- function(dresolution, tresolution){
    
    dtest <- seq(0, 1, by=dresolution)
    ttest <- seq(tresolution, 1, by=tresolution)
    error_mat <- matrix(NA, nrow=length(dtest), ncol = length(ttest))
    
    for(i in dtest){
        for(j in ttest){
            res_lambda <- explookup(i, j)
            if(is.nan(D[toString(res_lambda), toString(j)])){
                stop(paste('an answer in test range is NaN:', 
                           i, j, D[toString(res_lambda), toString(j)]))}
            if(is.na(D[toString(res_lambda), toString(j)])){
                stop(paste('an answer in test range is NA:', 
                           i, j, D[toString(res_lambda), toString(j)]))}
            #print(paste(i, D[toString(res_lambda), toString(j)], abs(i-D[toString(res_lambda), toString(j)])))
            # Get lookup error
            error_mat[which(dtest == i),which(ttest == j)] <- 
                abs(i - D[toString(res_lambda), toString(j)])
            
            if(error_mat[which(dtest == i),which(ttest == j)] > 0.01){ print(paste(i, j)) }
        }
    }
    print(paste("Resolution - density:", dresolution, ",threshold:", tresolution))
    print(paste("For resolution generated function has lookup mean error of:", mean(error_mat)))
    print(paste("lookup max error of:", max(error_mat)))
}
testlookup(0.01, 0.01)
# gen res(0.01, 0.001) ~831.6Mb 25 Oct 2019
# testlookup(0.01, 0.01) > avg error 
#[1] "Resolution - density: 0.01 ,threshold: 0.01"
#[1] "For resolution generated function has mean error of: 6.27222569567461e-05"
#[1] "max error of: 0.000464805549974086"

# Rough test because integrate at low values may have issues
# https://stat.ethz.ch/pipermail/r-help/2008-December/182311.html
lambda <- explookup(0.3, 0.5) # get exp distribution param based on density and threshold
dist <- rexp(round(1000), 1/lambda)
pdf <- density(dist)
plot(pdf)
f <- approxfun(pdf$x, pdf$y, yleft=0, yright=0)
cdf <-integrate(f, 0.5, Inf)  # replace '2' by any other value.
cdf

save(list = c('D', 'explookup'), file = 'exp_density_lookup.Rdata')

#test load D
load('exp_density_lookup.Rdata')
