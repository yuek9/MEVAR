################
## code for real data analysis
################

###############
## use two group comparison
## question: how to determine lag order? here use K=3 by looking at individual marginal ACF; was thinking about using individual VAR order, but for high dimension it is not easy to do lag order selection
###############

rm(list = ls())

## lag order to fit
K=3


## we select random subsample to have equal number of subjects each group
seed_i = 1 ## change this seed to select a different subset of subjects 
usersize = 80 ## set the number of subjects to use each group
sensi_subset = 7 ## change this starting point to select a different time series data, 1~9

time_length = 300

n_node = 200


#setwd('\\\\biostat-fs2\\users\\yuek\\Desktop\\MEVAR')
setwd('~/Desktop/MEVAR')



read_filepath = paste0('../LMM_testing/real_data_analysis/data/3T_HCP1200_MSMAll_d200_ts2_sess1/')

grouping = read.table(paste0(read_filepath, 'THC.txt'), header=F, sep=',')
grouping = grouping[,3]
table(grouping)


## separate analysis for each node and for each group
inputs = commandArgs(T)
infer_row = as.integer(inputs[1])




library(ggplot2)
library(glmnet)
library(hdi)
library(vars)

source('lib.R')


#### perform single group analysis

## read data example
i=1
tmpdata = read.table(paste0(read_filepath, 'subj', i, '.txt'), header=F, sep=',')
tmpdata = tmpdata[seq(1, 1200, 2),] ## down-sample data to control the serial dependence
acf(tmpdata[,1]) ## seems VRA(3) should be good

## select a subsample in each group
set.seed(seed_i)

group_subsample_index = list('0' = sort(sample((1:length(grouping))[grouping == 0], size = usersize, replace=F)),
                             '1' = sort(sample((1:length(grouping))[grouping == 1], size = usersize, replace=F)))


data = list('0' = list(), '1' = list())
for (g in c('0', '1')){
  ii=0
  for(i in group_subsample_index[[g]]){
    ii=ii+1
    tmpdata = read.table(paste0(read_filepath, 'subj', i, '.txt'), header=F, sep=',')
    tmpdata = tmpdata[seq(1, 1200, 2)[1:time_length],] 
    ## standardize the signal for each subject
    tmpdata = scale(tmpdata)
    # tmpdata = cbind(rep(i, nrow(tmpdata)), tmpdata)
    # colnames(tmpdata)<- c('subj', paste0('node', 1:(n_node)))
    colnames(tmpdata)<- c( paste0('node', 1:(n_node)))
    data[[g]][[ii]] = tmpdata
  }
}

#######################


## estimate all edges and test for all edges. 
## test each group on its own for non-zero edges;

a.seq = exp(seq(log(1e-4), log(50), length.out = 20))

tmp = list('beta.hat' = matrix(0, nrow=n_node, ncol=n_node*K),
           'beta.db' = matrix(0, nrow=n_node, ncol=n_node*K), 
           'beta.sd' = matrix(0, nrow=n_node, ncol=n_node*K))
store_res <- list( '0' = tmp, '1'=tmp)


print(infer_row)

inf.coord = 1:(n_node*K) ## inference for all regressed nodes

## compute for each group
for (g in c('0', '1')){
  
  set.seed(1)
  
  input_X = lapply(data[[g]], function(Y){
    input_X_i = t(sapply(0:(nrow(Y)-K-1), function(ii){
      as.vector(t(Y[ii + K:1,]))
    }))
    input_X_i
  })
  
  input_y = lapply(data[[g]], function(Y){
    Y[(K+1):nrow(Y),infer_row]
  })
  
  grp = do.call(c,sapply(1:length(data[[g]]), function(i) rep(i, nrow(input_X[[i]])), simplify = F))
  X = do.call(rbind, input_X)
  y = do.call(c, input_y)
  Z = X
  
  {
    ## implement the proposed method
    cat('Proposed...\n')
    tryCatch(
      {
        
        
        a.opt<-select_a(X=X, y=y, Z=Z, 
                        grp = grp, 
                        a.seq = a.seq,
                        kfold = 4, time_series_cv = T , 
                        cv.max.iter = 5)
        
        res.opta = Fixed_effects_est_inf(X=X, y=y, Z=Z, grp=grp, a=a.opt, inf.coord=inf.coord,
                                         lm=F, 
                                         cv.max.iter = 5,
                                         time_series_cv = T)

        store_res[[g]]$beta.hat[infer_row, ] <- res.opta$beta.hat
        store_res[[g]]$beta.db[infer_row, ] <- res.opta$beta.db
        store_res[[g]]$beta.sd[infer_row, ] <- res.opta$beta.db.sd

      }, error = function(e)print(e)
    )
    
  }
  
  rm(list=c('X', 'Z', 'input_X', 'input_y', 'y'))

  
  save.image(paste0('real_data_analysis/res_two_group_total_', n_node, ' nodes_', usersize, '_groupsize_subjseed_', seed_i, '_timelength_', time_length, '_headnode_', infer_row, '.RData'))
  
}


##########################################################################
## generate plot for paper

if(F){
  setwd('\\\\biostat-fs2\\users\\yuek\\Desktop\\LMM_testing')
  library(ggplot2)
  library(igraph)
  color_cb=c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                      "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                      
  par(mfrow=c(1,3))
  # usersize = 100
  # seed_i = 1
  # sensi_subset = 7 
  # time_length = 120
  # 
  # n_node = 200
  # 
  
  tmp = list('beta.hat' = diag(0,n_node),'beta.db' = diag(0,n_node), 'beta.sd' = diag(0, n_node), 'vc' = diag(0,n_node))
  store_res_all<- store_res_Li_all  <- store_res_lasso_all <- list( '0' = tmp, '1'=tmp)
  rm('tmp')
  
  ## collect results from each node
  for (headnode in 1:n_node){
    print(headnode)
    
    load(paste0('real_data_analysis/res_two_group_total_', n_node, ' nodes_', usersize, '_groupsize_subjseed_', seed_i,'_timestart_', sensi_subset, '_timelength_', time_length, '_headnode_', headnode, '.RData'))
    
    for (g in c('0', '1')){
      
      store_res_lasso_all[[g]]$beta.hat[headnode, -headnode] <- store_res_lasso[[g]]$beta.hat[headnode, -headnode]
      store_res_lasso_all[[g]]$beta.db[headnode, -headnode] <- store_res_lasso[[g]]$beta.db[headnode, -headnode]
      store_res_lasso_all[[g]]$beta.sd[headnode, -headnode] <- store_res_lasso[[g]]$beta.sd[headnode, -headnode]
      
      store_res_all[[g]]$beta.hat[headnode, -headnode] <- store_res[[g]]$beta.hat[headnode, -headnode]
      store_res_all[[g]]$beta.db[headnode, -headnode] <- store_res[[g]]$beta.db[headnode, -headnode]
      store_res_all[[g]]$beta.sd[headnode, -headnode] <- store_res[[g]]$beta.sd[headnode, -headnode]
      
      store_res_all[[g]]$vc[headnode,] <- store_res[[g]]$vc[headnode,]
      
      
      
      store_res_Li_all[[g]]$beta.hat[headnode, -headnode] <- store_res_Li[[g]]$beta.hat[headnode, -headnode]
      store_res_Li_all[[g]]$beta.db[headnode, -headnode] <- store_res_Li[[g]]$beta.db[headnode, -headnode]
      store_res_Li_all[[g]]$beta.sd[headnode, -headnode] <- store_res_Li[[g]]$beta.sd[headnode, -headnode]
    }
  }
  
  ## save variable
  store_res_all_method <- list('Proposed' = store_res_all,
                               'Li' = store_res_Li_all,
                               'lasso' = store_res_lasso_all)
  
  save.image('tmp.RData')
  
  
  setwd('\\\\biostat-fs2\\users\\yuek\\Desktop\\LMM_testing')
  load('\\\\biostat-fs2\\users\\yuek\\Desktop\\LMM_testing/tmp.RData')
  
  library(igraph)
  
  beta_values <- beta_sd_values <- pval <- pval_adj <- vc <- all_edge <- nets <- list()
  
  
  ## the network for proposed method, lasso and Li
  for (method in c('Proposed', 'Li', 'lasso')){
    for (g in c('0', '1')){
      
      beta_values[[method]][[g]] = (store_res_all_method[[method]][[g]]$beta.db+t(store_res_all_method[[method]][[g]]$beta.db))/2
      beta_sd_values[[method]][[g]] = sqrt((store_res_all_method[[method]][[g]]$beta.sd^2+t(store_res_all_method[[method]][[g]]$beta.sd^2))/2)
      
      {
        pval[[method]][[g]] = 2*(1-pnorm(
          abs( beta_values[[method]][[g]] / beta_sd_values[[method]][[g]]  )
        ))
        pval_adj[[method]][[g]] = matrix(0, n_node, n_node)
        pval_adj[[method]][[g]][lower.tri(pval[[method]][[g]], diag=F)] <- 
          p.adjust(pval[[method]][[g]][lower.tri(pval[[method]][[g]], diag=F)], method = 'holm')
        pval_adj[[method]][[g]] = pval_adj[[method]][[g]] + t(pval_adj[[method]][[g]]) + diag(1, n_node)
        
        
        all_edge[[method]][[g]] = (pval_adj[[method]][[g]]<0.05)*1 
        diag(all_edge[[method]][[g]]) <- 0
        
        sum(all_edge[[method]][[g]][lower.tri(all_edge[[method]][[g]], diag=F)])
        
        nets[[method]][[g]] = graph_from_adjacency_matrix(all_edge[[method]][[g]],  mode = 'undirected', weighted=T)
        
        plot(nets[[method]][[g]], layout = layout_in_circle(nets[[method]][[g]]),
             edge.arrow.size=.2, edge.curved=0.7, 
             vertex.size=3,
             edge.color = E(nets[[method]][[g]])$color)
        
        title(paste(method, g))
      }
    }
  }
  
  ## for proposed method, obtain the variance components
  for (method in c('Proposed')){
    for (g in c('0', '1')){
      
      vc[[method]][[g]] = (store_res_all_method[[method]][[g]]$vc+
                             t(store_res_all_method[[method]][[g]]$vc))/2
      
      {
        
        tmp_all_edge = vc[[method]][[g]] > 0
        diag(tmp_all_edge) <- 0
        
        sum(tmp_all_edge[lower.tri(tmp_all_edge, diag=F)])
        
        tmp_nets = graph_from_adjacency_matrix(tmp_all_edge,  mode = 'undirected', weighted=T)
        
        plot(tmp_nets, layout = layout_in_circle(tmp_nets),
             edge.arrow.size=.2, edge.curved=0.7, 
             vertex.size=3,
             edge.color = E(tmp_nets)$color)
        
        title(paste(method, g))      }
    }
  }
  
  
  
  save(list=c('pval', 'pval_adj', 
              'beta_values', 'beta_sd_values',
              'vc'
              
  ), file = paste0('real_data_analysis/matrices_for_plot_nnode_', n_node, '_nsubj_', usersize, '.RData'))
  
  
  
  n_node = 200
  usersize = 80
  load(paste0('real_data_analysis/matrices_for_plot_nnode_', n_node, '_nsubj_', usersize, '.RData'))
  ## print results of estimated networks
  
  sum(pval_adj$Proposed$`0`<0.05)/2
  sum(pval_adj$Proposed$`0`<0.05)/2 / (200*199/2)
  
  sum(pval_adj$Proposed$`1`<0.05)/2
  sum(pval_adj$Proposed$`1`<0.05)/2 / (200*199/2)
  
  
  sum((pval_adj$Proposed$`1`<0.05) & (pval_adj$Proposed$`0` <0.05))/2
  
  
  sum(pval_adj$lasso$`0`<0.05)/2
  sum(pval_adj$lasso$`0`<0.05)/2 / (200*199/2)
  
  sum(pval_adj$Li$`0`<0.05)/2
  sum(pval_adj$Li$`0`<0.05)/2 / (200*199/2)
}





