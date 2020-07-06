###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
############ plot magnitude of change in resistance prevalence ####################
###################################################################################

library(ggplot2)
library(ggpubr)

get_plotdata <- function(rawlist, name){
  d.simple = cbind.data.frame(Duration = rawlist[['simple']]$data$long_dur, res = rawlist[['simple']]$res[,,1][,3])
  d.binary = cbind.data.frame(Duration = rawlist[['binary']]$data$long_dur, res = rawlist[['binary']]$res[,,1][,3])
  d.freq = cbind.data.frame(Duration = rawlist[['freq']]$data$long_dur, res = rawlist[['freq']]$res[,,1][,3])
  
  d = rbind.data.frame(d.simple, d.binary, d.freq)
  d$name = name
  
  return(d)
}


plot_magnitude <- function(plotlist.left, plotlist.right, n.bed, n.days) {
  
  plotlist.left = lapply(plotlist.left, function(x) { 
    x$res = x$res*n.bed*n.days
    return (x) 
  })
  d = do.call('rbind.data.frame', plotlist.left)
  
  p.left = ggplot(aes(x = Duration, y = res, color = name), data = d) +
    #geom_point(alpha = 0.1, size = 1) +
    geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), alpha = 0.1) +
    coord_cartesian(ylim = c(-120, 40)) +
    scale_color_manual(values=c('#AA6373', '#F4A5AE', '#6DB1BF', '#0075A2'),
                       name = "Scenarios", 
                       labels = c('High resistance prevalence/transmission\nHigh number of prescriptions', 
                                  "Low resistance prevalence/transmission\nHigh number of prescriptions",
                                  "High resistance prevalence/transmission\nLow number of prescriptions",
                                  'Low resistance prevalence/transmission\nLow number of prescriptions'), 
                       limits = c('3gcrehh', 
                                  "3gcrelh",
                                  "3gcrehl", 
                                  '3gcrell')) +
    ylab(paste0('Reduction in patient-days with resistance carriers\nin a ', n.bed ,'-bed ward in ', n.days, ' days')) +
    scale_y_continuous(breaks = seq(-200, 50, by = 10)) + 
    labs(subtitle = '3rd generation cephalosporin resistant Enterobacteriaceae') +
    xlab('Duration of antibiotics (days)')+
    theme_minimal()+
    theme(legend.position = 'bottom') 
  
  plotlist.right = lapply(plotlist.right, function(x) { 
    x$res = x$res*n.bed*n.days
    return (x) 
  })
  d = do.call('rbind.data.frame', plotlist.right)
  
  p.right = ggplot(aes(x = Duration, y = res, color = name), data = d) +
    #geom_point(alpha = 0.1, size = 1) +
    geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), alpha = 0.1) +
    coord_cartesian(ylim = c(-120, 40)) +
    scale_color_manual(values=c('#AA6373', '#F4A5AE', '#6DB1BF',  '#0075A2'),
                       name = "Scenarios", 
                       labels = c('High resistance prevalence/transmission\nHigh number of prescriptions', 
                                  "Low resistance prevalence/transmission\nHigh number of prescriptions",
                                  "High resistance prevalence/transmission\nLow number of prescriptions",
                                  'Low resistance prevalence/transmission\nLow number of prescriptions'),
                       limits = c('crehh', 
                                  "crelh",
                                  "crehl", 
                                  'crell')) +
    scale_y_continuous(breaks = seq(-200, 50, by = 10)) + 
    ylab('') +
    labs(subtitle = 'Carbapenem resistant Enterobacteriaceae') +
    xlab('Duration of antibiotics (days)')+
    theme_minimal()+
    theme(legend.position = 'bottom') 
  
  
  return(ggarrange(p.left, p.right, common.legend = T, legend = 'bottom'))
  
}

