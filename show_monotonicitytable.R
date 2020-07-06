###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
#### Check monotonicity of the parameters - Hoeffding's D and Spearman's ##########
###################################################################################

source('plot_functions/table_monotonicity.R')

output = get(load('runs/LHSdiff_simple_355_notzero15Apr2020_2214BST.Rdata'))

monotonicity.tab(output)
