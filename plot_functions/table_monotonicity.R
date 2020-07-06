###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
#### Check monotonicity of the parameters - Hoeffding's D and Spearman's ##########
###################################################################################
rm(list=ls()) # Clean working environment

library(xtable)

monotonicity.tab <- function(output) {
  
  parametersamples = output[['data']]
  output = output$res[,3,]
  hd = hd.p = c()
  spm = spm.p = c()
  for (i in 1:ncol(parametersamples)){
    hd[[i]] = hoeffd(x=parametersamples[,i], y = output)$D[1,2]
    hd.p[[i]] = hoeffd(x=parametersamples[,i], y = output)$P[1,2]
    spm[[i]] = cor.test(x=parametersamples[,i], y = output, method = "spearman")$estimate
    spm.p[[i]] = cor.test(x=parametersamples[,i], y = output, method = "spearman")$p.value
  }
  corr.table = data.frame(Parameters = colnames(parametersamples), 
                          HoeffdingD = format(round(hd, 3), nsmall = 3), 
                          HoeffdingD.p = format(round(hd.p, 3), nsmall = 3), 
                          SpearmanRank = format(round(spm, 3), nsmall = 3), 
                          SpearmanRank.p = format(round(spm.p, 3), nsmall = 3))
  colnames(corr.table) = c('Parameters',
                           'Hoeffding`s D measure', 'Hoeffding`s D p-value', 
                           'Spearman`s rank correlation measure', 'Spearman`s rank correlation p-value')
  return(xtable(corr.table))
}

