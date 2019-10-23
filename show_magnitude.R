########Present magnitude of change in resistance prevalence with time#####
###########################################################################
setwd('/Users/moyin/Documents/nBox/git_projects/indiv_abxduration')
library(ggplot2)

#load data from prevalence prcc 
d.simple.raw=get(load('runs/LHS_simple_900_zero16Oct2019_1651BST.Rdata'))
d.simple= cbind.data.frame(Duration=d.simple.raw$data$meanDur, res=d.simple.raw$res, type='Model 1')
d.binary.raw=get(load('runs/LHS_binary_900_zero_18Oct2019_1357BST.Rdata'))
d.binary= cbind.data.frame(Duration=d.binary.raw$data$meanDur, res=d.binary.raw$res, type='Model 2')
d.freq.raw=get(load('runs/LHS_frequency_500_zero_22Oct2019_0801BST.Rdata'))
d.freq= cbind.data.frame(Duration=d.freq.raw$data$meanDur, res=d.freq.raw$res, type='Model 3')
plotd=rbind.data.frame(d.simple,d.binary,d.freq)

ggplot(aes(x=Duration, y=res, colour=type), data=plotd)+
  geom_smooth(method='lm',formula=y~x, se=F)+
  geom_point(alpha=0.2, size=0.2)+
  ylab('Number of patients carrying CPE/ bed/ day')+
  ylim(0,1)+
  scale_y_continuous(breaks=seq(0,1,0.05))+
  scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  xlab('Duration of antibiotics (days)')+
  theme_minimal()+
  theme(legend.position = 'bottom',legend.title = element_blank())
