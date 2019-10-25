############################################################################
###########Effect of antibiotic duration on hospitalised patients###########
################# Heatmap exploration of 2 parameters ######################
############################################################################

setwd('/Users/moyin/Documents/nBox/git_projects/indiv_abxduration')
rm(list=ls()) # Clean working environment

require(fields) #heatmap

model='simple'
source(paste0('model_',model,".R"))
source("default_params.R")

pixels=30 #how many pixels per heatmap

###################### Change parameters here: START ##################
y_name <- "prop_R"
y_seq <- seq(0, 1, length.out = pixels)

x_name <- "repop.s"
x_seq <- seq(0, 0.2, length.out = pixels)

###################### Change parameters here: END ####################

# Parameters not added here are considered fixed as defined above in Fixed Parameters
default.para.list = para[which(names(para) %in% parameters_diff_prevalence)]
default.para.df = as.data.frame(matrix(as.vector((rep(unlist(default.para.list), each=pixels**2))), ncol=pixels**2, byrow = T))
rownames(default.para.df) = names(default.para.list)
colnames(default.para.df) = 1:pixels**2
default.para.df = default.para.df[match(parameters_diff_prevalence, rownames(default.para.df)),]

default.para.df[which(rownames(default.para.df)==x_name),] = rep(x_seq, each=pixels)
default.para.df[which(rownames(default.para.df)==y_name),] = rep(y_seq, pixels)

# empty vectors and matrix to save runs 
save_runs <- list()

###run model 
if (model == 'simple'){
  
  for (j in 1:ncol(default.para.df)) {
    
    print(j, 'out of', pixels**2, 'runs')
    
    save_runs[j] = diff_prevalence(n.bed=default.para.df[1,j], max.los=default.para.df[2, j], prop_R=default.para.df[3, j], prop_S=default.para.df[4, j],
                                   bif=default.para.df[5, j], pi_ssr=default.para.df[6, j], repop.s=default.para.df[7, j], 
                                   mu=default.para.df[8, j], abx.s=default.para.df[9, j], abx.r=default.para.df[10, j], p.infect=default.para.df[11, j], 
                                   cum.r.1=default.para.df[12, j], p.r.day1=default.para.df[13, j], 
                                   short_dur=default.para.df[14, j], long_dur=default.para.df[15, j])[[3]]
  }

} else if (model == 'binary'){
  
  for (j in 1:ncol(default.para.df)){
    
    print(j, 'out of', pixels**2, 'runs')
    
    save_runs[j] = diff_prevalence(n.bed=default.para.df[1,j], max.los=default.para.df[2,j], prop_R=default.para.df[3,j], prop_r=default.para.df[4,j], 
                                   prop_Sr=default.para.df[5,j], prop_S=default.para.df[6,j],
                                   bif=default.para.df[7,j], pi_ssr=default.para.df[8,j], repop.s=default.para.df[9,j], repop.r=default.para.df[10,j], 
                                   mu=default.para.df[11,j], abx.s=default.para.df[12,j], abx.r=default.para.df[13,j], 
                                   p.infect=default.para.df[14,j], cum.r.1=default.para.df[15,j], p.r.day1=default.para.df[16,j], 
                                   short_dur=default.para.df[17,j], long_dur=default.para.df[18,j])[[3]]
  }
  
} else if (model == 'frequency'){
  
  for (j in 1:ncol(default.para.df)){
    
    print(j, 'out of', pixels**2, 'runs')
    
    save_runs[j] = diff_prevalence(n.bed=default.para.df[1,j], max.los=default.para.df[2,j], p.infect=default.para.df[3,j], 
                                   cum.r.1=default.para.df[4,j], p.r.day1=default.para.df[5,j],
                                   K=default.para.df[6,j], total_prop=default.para.df[7,j], prop_R=default.para.df[8,j], pi_ssr=default.para.df[9,j], 
                                   r_trans=default.para.df[10,j], r_growth=default.para.df[11,j], r_thres=default.para.df[12,j], s_growth=default.para.df[13,j],
                                   abx.s=default.para.df[14,j], abx.r=default.para.df[15,j], short_dur=default.para.df[16,j],long_dur=default.para.df[17,j])[[3]]
  }
}

if (length(save_runs)!=ncol(default.para.df)){
  print('The runs did not complete! Stop here and check the parameter range')
}

outcome.df= cbind.data.frame(x=rep(x_seq, each=pixels), 
                             y=rep(y_seq, pixels) ,
                             outcome=unlist(save_runs))
ggplot(outcome.df, aes(x, y)) +
  geom_tile(aes(fill = outcome)) +
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       midpoint=0, limits=range(outcome.df$outcome),
                       breaks = c(as.numeric(format(round(min(outcome.df$outcome),3),nsmall=3)),
                                  0, 
                                  as.numeric(format(round(max(outcome.df$outcome),3),nsmall=3))),
                       name = "Difference in resistance carriers \nin wards given long vs short \nantibiotic durations") +
  ylab(y_name)+
  xlab(x_name)+
  theme_minimal()+
  theme(legend.position = 'bottom')

