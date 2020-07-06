###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
######### Plot heatmap of the parameters against models ###########################
###################################################################################
rm(list=ls()) # Clean working environment

#### Input run data required for plotting parameter heatmap 
Amodel1=get(load('runs/LHSdiff_simple_370_notzero23Apr2020_2041BST.Rdata')) #abx_r>0
Amodel2=get(load('runs/LHSdiff_binary_370_notzero_23Apr2020_2328BST.Rdata'))
Amodel3=get(load('runs/LHSdiff_frequency_370_notzero_23Apr2020_2143BST.Rdata'))
Bmodel1=get(load('runs/LHSdiff_simple_370_zero23Apr2020_2000BST.Rdata')) #abx_r>0
Bmodel2=get(load('runs/LHSdiff_binary_370_zero_24Apr2020_0012BST.Rdata')) 
Bmodel3=get(load('runs/LHSdiff_frequency_370_zero_23Apr2020_2243BST.Rdata'))

#### plot parameter heatmap 
source('plot_functions/plot_paraheatmap.R')

plot_paraheatmap(forplot = forplot)

ggsave('../../../../Desktop/para_heatmap.jpeg', units = 'cm', width = 22, height=15)
