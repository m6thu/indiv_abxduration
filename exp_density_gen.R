#Understanding the exponential function
lambda <- seq(1,500,0.05)
thresholds <- seq(0,1,0.001)

D <- matrix(NA, nrow=length(lambda), ncol = length(thresholds))
for(i in 1:length(lambda)){
    for(j in 1:length(thresholds)){
        D[i, j] <- pexp(thresholds[j], rate = lambda[i])
    }
}
rownames(D) <- lambda
colnames(D) <- thresholds

save(D, file = 'exp_density_lookup.Rdata')

#test load D
load('exp_density_lookup.Rdata')
