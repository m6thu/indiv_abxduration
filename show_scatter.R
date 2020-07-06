###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
######## Plot scatter plot of different parameters to check monotonicity###########
###################################################################################

source('plot_functions/plot_scatter.R')

sd=get(load('runs/LHSdiff_simple_355_notzero15Apr2020_2214BST.Rdata'))
bd=get(load('runs/LHSdiff_binary_355_notzero_15Apr2020_2340BST.Rdata'))
fd=get(load('runs/LHSdiff_frequency_370_notzero_22Apr2020_0152BST.Rdata'))

scatter(sd)
scatter(bd)
scatter(fd)