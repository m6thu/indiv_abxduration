###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
###################### show model exploration outputs #############################
###################################################################################

####Plot violin plot if duration matters in the models 
#libraries
library(ggplot2)

make.df <- function(vector, model, scenario) {
  
  df = data.frame(mod = model, 
                  y = vector, 
                  Scenario = scenario)
  return(df)
}

plot_durationsmatters <- function(output.list){
  
  #name the models 
  nam = factor(c('Simple 3-states model', 'Co-carriage 5-states model', 'Population growth model'), levels = c('Simple 3-states model', 'Co-carriage 5-states model','Population growth model'))
  scenarios = factor(c('3GCRE', 'CRE'), levels = c('CRE','3GCRE'))
  
  #put data into dataframe 
  mod1A = make.df (output.list$Amodel1$res[,,1][,3], nam[1], scenarios[1])
  mod2A = make.df (output.list$Amodel2$res[,,1][,3], nam[2], scenarios[1])
  mod3A = make.df (output.list$Amodel3$res[,,1][,3], nam[3], scenarios[1])
  mod1B = make.df (output.list$Bmodel1$res[,,1][,3], nam[1], scenarios[2])
  mod2B = make.df (output.list$Bmodel2$res[,,1][,3], nam[2], scenarios[2])
  mod3B = make.df (output.list$Bmodel3$res[,,1][,3], nam[3], scenarios[2])
  d = rbind(mod1A, mod2A, mod3A, mod1B, mod2B, mod3B)
  
  #plot graph 
  dur.p = ggplot(d, aes(x = Scenario, y = y, fill = Scenario)) +
    geom_violin(aes(color = Scenario, fill = Scenario), trim=FALSE) + 
    geom_hline(yintercept = 0, linetype="dotted", 
               color = "red", size=0.5)+
    scale_x_discrete(expand = c(0, 0.8))+
    scale_fill_manual(values = c('#8B1E3F', '#E2711D' ))+
    scale_color_manual(values = c('#8B1E3F', '#E2711D' ))+
    geom_segment(aes(x = 2.6, xend = 2.6, y = -.2, yend = .2),
                 arrow=arrow(length=unit(0.2,"cm")),colour = "grey") +
    geom_segment(aes(x = 2.6, xend = 2.6, y = .2, yend = -.2),
                 arrow=arrow(length=unit(0.2,"cm")),colour = "grey") +
    annotate('text', x = 2.6, y = .1, label = 'Longer antibiotic duration leads to \nmore resistance carriers', size = 2)+
    annotate('text', x = 2.6, y = -.1, label = 'Shorter antibiotic duration leads to \nmore resistance carriers', size = 2)+
    facet_wrap(~ mod, ncol = 1)+
    labs(x = '', y = 'Difference in proportion of resistance carriers between patients \nreceiving long and short duration of antibiotics per day')+
    theme_bw() +
    theme(legend.position = 'none') +
    coord_flip()
  
  return(dur.p)
}