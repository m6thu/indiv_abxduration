##plot 1 heatmap 

model = 'binary'

source(paste0('model_', model, '.R'))
source('default_params.R')

pixels = 8  #how many pixels per heatmap

data_heatmap <- function(model = model, x_name, y_name, x_seq, y_seq, scenario) {
  
  # Parameters not explored are considered fixed as defined in default_params
  if (scenario == 'esbl') {
    para.feed = c(para[[1]]$esbl, para[[2]]$esbl)
  } else {
    para.feed = c(para[[1]]$cpe, para[[2]]$cpe)
  }
  default.para.list = para.feed[which(names(para.feed) %in% parameters_diff_prevalence)]
  
  # format parameter input into lists to run model 
  default.para.df = as.data.frame(matrix(as.vector((rep(unlist(default.para.list), each = pixels**2))), ncol=pixels**2, byrow = T))
  rownames(default.para.df) = names(default.para.list)
  colnames(default.para.df) = 1:pixels**2
  default.para.df = default.para.df[match(parameters_diff_prevalence, rownames(default.para.df)),]
  default.para.df[grep(x_name,rownames(default.para.df)),] = rep(x_seq, each = pixels)
  default.para.df[grep(y_name,rownames(default.para.df)),] = rep(rep(y_seq, pixels), each = length(grep(y_name, rownames(default.para.df))))
  feed.list = as.list(default.para.df)
  
  ###run model 
  if (model == 'simple'){
    
    foo = function (x) {
      diff_prevalence(n.bed= x[1], max.los=x[2], prop_R=x[3], prop_S=x[4],
                      bif=x[5], pi_ssr=x[6], repop.s=x[7], 
                      mu=x[8], abx.s=x[9], abx.r=x[10], p.infect=x[11], 
                      cum.r.1=x[12], p.r.day1=x[13], 
                      short_dur=x[14], long_dur=x[15])[3]
    }
    
    save_runs = mclapply(feed.list, foo,  mc.cores = 11)
    
  } else if (model == 'binary'){
    
    foo=function (x) { 
      diff_prevalence(n.bed=x[1], max.los=x[2], prop_R=x[3], prop_r=x[4], 
                      prop_Sr=x[5], prop_S=x[6],
                      bif=x[7], pi_ssr=x[8], repop.s=x[9], repop.r=x[10], 
                      mu=x[11], abx.s=x[12], abx.r=x[13], 
                      p.infect=x[14], cum.r.1=x[15], p.r.day1=x[16], 
                      short_dur=x[17], long_dur=x[18])[3]
    }
    
    save_runs = mclapply(feed.list, foo, mc.cores = 11)
    
  } else if (model == 'frequency'){
    
    foo = function (x) { 
      diff_prevalence(n.bed=x[1], max.los=x[2], p.infect=x[3], 
                      cum.r.1=x[4], p.r.day1=x[5],
                      K=x[6], total_prop=x[7], prop_R=x[8], pi_ssr=x[9], 
                      r_trans=x[10], r_growth=x[11], r_thres=x[12], s_growth=x[13],
                      abx.s=x[14], abx.r=x[15], short_dur=x[16],long_dur=x[17])[3]
    }
    
    save_runs = mclapply(feed.list, foo,  mc.cores = 11)
  }
  
  outcome.df = cbind.data.frame(x = rep(x_seq, each=pixels), 
                                y = rep(y_seq, pixels),
                                outcome = unlist(save_runs))
  
  out = list (outcome.df, x_name, y_name)
  
  return(out)
}

plot_heatmap <- function (plots_data) {
  
  font.size = 20
  
  outcome.df = plots_data[[1]]
  y_name = plots_data[[3]]
  x_name = plots_data[[2]]
  
  plot = ggplot(outcome.df, aes(x, y)) +
    geom_tile(aes(fill = outcome)) +
    scale_fill_gradient2(low="#388697", mid="grey95", high="#CC2936", 
                         midpoint = 0, limits=range(outcome.df$outcome),
                         breaks = c(as.numeric(format(round(min(outcome.df$outcome),3),nsmall=3))+0.001,
                                    as.numeric(format(round(0),0), nsmall=0), 
                                    as.numeric(format(round(max(outcome.df$outcome),3),nsmall=3))-0.001),
                         name = "Difference in resistance carriers in wards given\nlong vs short antibiotic durations") +
    xlab(x_name)+
    ylab(y_name)+
    theme_minimal()+
    theme(legend.position = 'none', 
          text = element_text(size=font.size))
  
  return(plot)
}

plots_data = data_heatmap(model = model, x_name = 'repop.r', y_name = 'pi_ssr', x_seq = seq(0, 0.15, length.out = pixels), y_seq= seq(0, 0.3, length.out = pixels), scenario = 'cpe')
plot_heatmap(plots_data)
