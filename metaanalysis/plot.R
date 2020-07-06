# ========================================== #
# Plot dose-response curve from model output 
# ========================================== #

rm(list=ls()) # clear environmental 
library(ggplot2); library(rstan) # load required libraries 

##############
# Get data for plot

# outputs from stan model 
fit = readRDS("model0.1.Rdata")

# parameter names 
source('define_data.R') # Get hand hygiene frequency from data 
params.trial = c("a0", "b0", 
                 "a[1]","a[2]","a[3]","a[4]","a[5]","a[6]",
                 "a[7]","a[8]","a[9]","a[10]","a[11]","a[12]",
                 "a[13]","a[14]","a[15]",
                 "b[1]","b[2]","b[3]","b[4]","b[5]","b[6]", 
                 "b[7]","b[8]","b[9]","b[10]","b[11]","b[12]", 
                 "b[13]","b[14]","b[15]",
                 "OR_abx_dur")
p1 = NULL
for(i in 1:13) p1<-c(p1, paste0("predicted_mean[", i, "]") ) # means of the parameters 
params = c(params.trial, p1)

# posterior values from stan model output 
d = rstan::extract(fit, pars = params)

##############
# Prepare posterior values into dataframe for ggplot 

dur = matrix(c(unlist(sapply(data_esbl, "[", 'day_short')), unlist(sapply(data_esbl, "[", 'day_long'))), ncol = 2, byrow = F)

# posterior values for each arm (for bubbles)
df.trial = data.frame(paste(names(stan_data_esbl)[1], 'short') = exp(d$`a[1]` + d$`b[1]`*dur[1,1]),
                      paste(names(stan_data_esbl)[1], 'long') = exp(d$`a[1]` + d$`b[1]`*dur[1,2]),
                      paste(names(stan_data_esbl)[2], 'short') = exp(d$`a[2]` + d$`b[2]`*dur[2,1]),
                      paste(names(stan_data_esbl)[2], 'long') = exp(d$`a[2]` + d$`b[2]`*dur[2,2]),
                      paste(names(stan_data_esbl)[3], 'short') = exp(d$`a[3]` + d$`b[3]`*dur[3,1]),
                      paste(names(stan_data_esbl)[3], 'long') = exp(d$`a[3]` + d$`b[3]`*dur[3,2]),
                      paste(names(stan_data_esbl)[4], 'short') = exp(d$`a[4]` + d$`b[4]`*dur[4,1]),
                      paste(names(stan_data_esbl)[4], 'long') = exp(d$`a[4]` + d$`b[4]`*dur[4,2]),
                      paste(names(stan_data_esbl)[5], 'short') = exp(d$`a[5]` + d$`b[5]`*dur[5,1]),
                      paste(names(stan_data_esbl)[5], 'long') = exp(d$`a[5]` + d$`b[5]`*dur[5,2]),
                      paste(names(stan_data_esbl)[6], 'short') = exp(d$`a[6]` + d$`b[6]`*dur[6,1]),
                      paste(names(stan_data_esbl)[6], 'long') = exp(d$`a[6]` + d$`b[6]`*dur[6,2]),
                      paste(names(stan_data_esbl)[7], 'short') = exp(d$`a[7]` + d$`b[7]`*dur[7,1]),
                      paste(names(stan_data_esbl)[7], 'long') = exp(d$`a[7]` + d$`b[7]`*dur[7,2]),
                      paste(names(stan_data_esbl)[8], 'short') = exp(d$`a[8]` + d$`b[8]`*dur[8,1]),
                      paste(names(stan_data_esbl)[8], 'long') = exp(d$`a[8]` + d$`b[8]`*dur[8,2]),
                      paste(names(stan_data_esbl)[9], 'short') = exp(d$`a[9]` + d$`b[9]`*dur[9,1]),
                      paste(names(stan_data_esbl)[9], 'long') = exp(d$`a[9]` + d$`b[9]`*dur[9,2]),
                      paste(names(stan_data_esbl)[10], 'short') = exp(d$`a[10]` + d$`b[10]`*dur[10,1]),
                      paste(names(stan_data_esbl)[10], 'long') = exp(d$`a[10]` + d$`b[10]`*dur[10,2]),
                      paste(names(stan_data_esbl)[11], 'short') = exp(d$`a[11]` + d$`b[11]`*dur[11,1]),
                      paste(names(stan_data_esbl)[11], 'long') = exp(d$`a[11]` + d$`b[11]`*dur[11,2]),
                      paste(names(stan_data_esbl)[12], 'short') = exp(d$`a[12]` + d$`b[12]`*dur[12,1]),
                      paste(names(stan_data_esbl)[12], 'long') = exp(d$`a[12]` + d$`b[12]`*dur[12,2]),
                      paste(names(stan_data_esbl)[13], 'short') = exp(d$`a[13]` + d$`b[13]`*dur[13,1]),
                      paste(names(stan_data_esbl)[13], 'long') = exp(d$`a[13]` + d$`b[13]`*dur[13,2]),
                      paste(names(stan_data_esbl)[14], 'short') = exp(d$`a[14]` + d$`b[14]`*dur[14,1]),
                      paste(names(stan_data_esbl)[14], 'long') = exp(d$`a[14]` + d$`b[14]`*dur[14,2]),
                      paste(names(stan_data_esbl)[15], 'short') = exp(d$`a[15]` + d$`b[15]`*dur[15,1]),
                      paste(names(stan_data_esbl)[15], 'long') = exp(d$`a[15]` + d$`b[15]`*dur[15,2]))
df.trial = as.data.frame(t(sapply(df.trial, quantile, probs = c(0.1, 0.5, 0.9)))) #plot 80% credible intervals 
dur.counts = c(unlist(t(dur))) # take the abx duration values for each arm 
df.trial$dur = dur.counts
df.trial$trial.names = rep(paste(names(stan_data_esbl)), each = 2)  # trial names (to determine color of the bubbles)

# posterior values for predicted data (for ribbons)
df.predicted = as.data.frame(matrix(unlist(d[grep('predicted', names(d))]), 
                                    byrow = F, nrow = lengths(d)))
colnames(df.predicted) = names(d)[grep('predicted', names(d))]
quants = c(0.1, 0.25, 0.5, 0.75, 0.9) # plot 50% and 80% credible intervals
d.plot.predict = as.data.frame(t(sapply(df.predicted, quantile, probs = quants)))
colnames(d.plot.predict) = paste0("pred_mean", quants)
d.plot.predict$dur = 0:(length(p1)-1)

##############
# Produce plot 

# colors for plots 
ribcol = c(alpha('grey', alph*0.3),
           alpha('grey', alph*0.6))
# plot
ggplot() + 
  geom_ribbon(aes(x = dur, ymin = pred_mean0.1, ymax = pred_mean0.9), data = d.plot.predict, fill = ribcol[1])+
  geom_ribbon(aes(x = dur, ymin = pred_mean0.25, ymax = pred_mean0.75), data = d.plot.predict, fill = ribcol[2])+
  geom_point(aes(x = dur, y = `50%`, color = trial.names,  size = 1/(`90%` - `10%`)), data = df.trial) +
  geom_line(aes(x = dur, y = pred_mean0.5), data = d.plot.predict, color='grey40') +
  scale_color_manual(values = cols) + scale_y_continuous(limits = c(0, 0.02)) + 
  scale_x_continuous(limits = c(0, 18))+ 
  scale_size(guide = 'none')+
  ylab('Daily probability of becoming a resistance carrier') + 
  xlab('Antibiotic duration') + 
  theme_minimal()+
  theme(legend.position = 'bottom', 
        legend.title = element_blank())+
  guides(color = guide_legend(nrow = 1))

# save plot
ggsave(paste0('../../../../Desktop/meta_plot.jpeg'), units = 'cm', width = 30, height= 15)

