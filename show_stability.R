#########################################################################
#######Effect of antibiotic duration on hospitalised patients############
#############Show convergence of models over time########################
#########################################################################
rm(list=ls()) # Clean working environment

source("los_abx_matrix.R")
source('plot_functions/plot_stability.R')

iterations = 100

#plot simple cpe 
abx.r = 0.0000001
abx.s = 0.2
model = 'simple'
simple0 = plot.stability(model='simple', iterations = iterations)

#plot binary cpe 
model = 'binary'
binary0 = plot.stability(model='binary', iterations = iterations)

#plot frequency cpe 
abx.s = 1.2
model='frequency'
freq0 = plot.stability(model='frequency', iterations = iterations)

#plot simple esbl
abx.s= 0.2
abx.r = 0.2
model='simple'
simplenot0=plot.stability(model='simple', iterations = iterations)

#plot binary esbl
model='binary'
binarynot0=plot.stability(model='binary', iterations = iterations)

#plot frequency esbl
abx.s = 1.2
abx.r = 1.2
model='frequency'
freqnot0=plot.stability(model='frequency', iterations = iterations)

#arrange all plots 
figure = ggarrange(
  ggarrange(simplenot0,simple0, ncol = 2), 
  ggarrange(binarynot0,binary0, ncol = 2), 
  ggarrange(freqnot0,freq0, ncol = 2), 
  nrow =3)

#save as jpeg
jpeg(filename="../../../../Desktop/stability.jpeg", width = 1000, height = 1200)
annotate_figure(figure,
                top = text_grob("             Scenario A                                                                                                                                    Scenario B"),
                left = text_grob("    Model 1                                                                                                                Model 2                                                                                                           Model 3   ", rot = 90))
dev.off() 
