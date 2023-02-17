########################
## summarize simulation results of MEVAR(d)
## remove summary for VC
########################

rm(list=ls())

setwd('\\\\biostat-fs2\\users\\yuek\\Desktop\\MEVAR')

library(mltools)
library(ggpubr)
library(ggplot2)

total_iter = 200

## expected simulation error
sqrt(0.05*0.95/total_iter)*1.96


## 
file_pre = '_setting6'

# p =30
K = 2
# n =20
# Time = 50
infer_row = 1
# iters=1


################### 
## obtain summary statistics

beta_MSE_res <-  # include total MSE
vc_sd_res <- # include total MSE for sd of random effects, and MSE for error variance and MCC
power_res <- # include power/type-I for each beta tested
typeI_res <-
CI_res <- # include CI coverage for each beta tested
time_res <- data.frame(NULL) ## include average computation time for each step


for(p in c(30)){
  for (n in c(20, 40, 60, 80)){
    for (Time in c(25, 50, 100, 150)){
      cat('processing', 'p', p, 'n', n, 'T', Time, '\n')
      
      {
        ## store estimates from each iteration
        store_all = replicate(2, 
                              list(beta = NULL, 
                                   beta.db = NULL, 
                                   beta.db.sd = NULL, 
                                   vc.b.sd = NULL, 
                                   vc.e.sd = NULL,
                                   CI = NULL), simplify = F)

        time = replicate(2, list(select_t = NULL, fixed_t =NULL, vc_t = NULL, lasso=NULL), simplify = F)
        
        res_CI_all <- res_rej_all <- replicate(2, NULL, simplify = F)
        
        names(store_all) <- names(time) <- names(res_CI_all) <- names(res_rej_all)<- c('proposed', 'lasso')
        
         
        
        for(iters in 1:total_iter){
          read_filename = paste0('res/VAR_', K, 'p', p, 'n', n, 'Time', Time, '/', file_pre, '_iter', iters, '_row_', infer_row, '.RData')
          
          tryCatch({
            
            load(read_filename)
            
            ## append proposed method results
            names(res.opta$beta.hat)  <- 1:p
            
            store_all$proposed$beta = rbind(store_all$proposed$beta, as.vector(res.opta$beta.hat))
            store_all$proposed$beta.db = rbind(store_all$proposed$beta.db, res.opta$beta.db)
            store_all$proposed$beta.db.sd = rbind(store_all$proposed$beta.db.sd, res.opta$beta.db.sd)
            
            if(F){
              names(res.vc.opta$vc.hat) <- 1:p
            store_all$proposed$vc.b.sd = rbind(store_all$proposed$vc.b.sd, sqrt(res.vc.opta$vc.hat))
            store_all$proposed$vc.e.sd = c(store_all$proposed$vc.e.sd, sqrt(res.vc.opta$sig2.eps.hat))
            }
            
            res_CI_all$proposed = rbind(res_CI_all$proposed, proposed_res$CI_cov) 
            res_rej_all$proposed = rbind(res_rej_all$proposed, proposed_res$rej)
            
            time$proposed$select_t = c(time$proposed$select_t, proposed_res$time['time_select_a.elapsed'])
            time$proposed$fixed_t = c(time$proposed$fixed_t, proposed_res$time['time_est_inf.elapsed'])
            time$proposed$vc_t = c(time$proposed$vc_t, proposed_res$time['time_vc.elapsed'])
            
            ## append lasso method results
            names(lasso_res$beta.hat) <- 1:p
            names(lasso_res$beta.db) <- names( lasso_res$beta.db.sd) <- names(lasso_res$CI_cov) <- names(lasso_res$rej)<-  inf.coord
            
            store_all$lasso$beta = rbind(store_all$lasso$beta, lasso_res$beta.hat)
            store_all$lasso$beta.db = rbind(store_all$lasso$beta.db, lasso_res$beta.db)
            store_all$lasso$beta.db.sd = rbind(store_all$lasso$beta.db.sd, lasso_res$beta.db.sd)
            
            
            res_CI_all$lasso = rbind(res_CI_all$lasso, lasso_res$CI_cov)
            res_rej_all$lasso = rbind(res_rej_all$lasso, lasso_res$rej )
            
            time$lasso$lasso = c(time$lasso$lasso, lasso_res$time)
            
          }, error = function(e) print(paste('p', p, 'n', n, 'Time', Time, 'row', infer_row, 'iter', iters))
          )
          
        }
        
        
        
        true_beta = beta
        # true_vc_sd = B_sd_companion[infer_row, inf.coord]
        true_e_sd = Sig_e[infer_row, infer_row]
        
        
        
        ## CI
        CI_res = rbind(CI_res, 
                       as.data.frame(t(c(K=K, p=p, n=n, Time=Time, row = infer_row, method = 'proposed', colMeans(res_CI_all$proposed)))),
                       as.data.frame(t(c(K=K, p=p, n=n, Time=Time, row = infer_row, method = 'lasso', colMeans(res_CI_all$lasso))))
                       )
        ## Power 
        power_res = rbind(power_res, 
                          as.data.frame(t(c(K=K, p=p, n=n, Time=Time, row = infer_row, method = 'proposed', colMeans(res_rej_all$proposed)[true_beta[inf.coord]!=0]))),
                          as.data.frame(t(c(K=K, p=p, n=n, Time=Time, row = infer_row, method = 'lasso', colMeans(res_rej_all$lasso)[true_beta[inf.coord]!=0])))
                          )
        ## Type I
        typeI_res = rbind(typeI_res, 
                          as.data.frame(t(c(K=K, p=p, n=n, Time=Time, row = infer_row, method = 'proposed', colMeans(res_rej_all$proposed)[true_beta[inf.coord]==0]))),
                          as.data.frame(t(c(K=K, p=p, n=n, Time=Time, row = infer_row, method = 'lasso', colMeans(res_rej_all$lasso)[true_beta[inf.coord]==0])))
                          )
        
        ## computation time
        time_res = rbind(time_res, 
                         as.data.frame(t(c(K=K, p=p, n=n, Time=Time, row = infer_row, method= 'proposed', 'select_a'=mean(time$proposed$select_t), 'est_and_inf'=mean(time$proposed$fixed_t), 'vc'=mean(time$proposed$vc_t)))),
                         as.data.frame(t(c(K=K, p=p, n=n, Time=Time, row = infer_row, method= 'lasso', 'select_a'=0, 'est_and_inf'=mean(time$lasso$lasso), 'vc'=0)))
                         )
        
        ## beta MSE (sum for all beta) (bias^2 + var)
        beta_MSE_res = rbind( beta_MSE_res, 
                              as.data.frame(t(c(K=K, p=p, n=n, Time=Time, row = infer_row, method = 'proposed', 
                                                'beta_total_MSE'= sum((apply(store_all$proposed$beta, 2, mean) - true_beta)^2 + apply(store_all$proposed$beta, 2, var))
                                                  ))),
                              as.data.frame(t(c(K=K, p=p, n=n, Time=Time, row = infer_row, method = 'lasso', 
                                                'beta_total_MSE'= sum((apply(store_all$lasso$beta, 2, mean) - true_beta)^2 + apply(store_all$lasso$beta, 2, var))
                                                  )))
                              
        )
        
        ## vc MSE and MCC
        
        if(F){
        vc_sd_res = rbind(vc_sd_res, 
                          as.data.frame(t(c(K=K, p=p, n=n, Time=Time, row = infer_row, method = 'proposed',
                           'vc_sd_total_MSE' =  sum((colMeans(store_all$proposed$vc.b.sd) - true_vc_sd)^2 + apply(store_all$proposed$vc.b.sd, 2, var)),
                           'vc_sd_MCC' = mean(sapply(1:nrow(store_all$proposed$vc.b.sd), function(i) mcc(preds = store_all$proposed$vc.b.sd[i,]>0, actuals = true_vc_sd>0))),
                           'vc_e_sd_MSE'= (mean(store_all$proposed$vc.e.sd) - true_e_sd)^2 + var(store_all$proposed$vc.e.sd)
                           )))
        )
        
        }
        
        
        
      }
    }
  }
}


true_vc_sd = B_sd_companion[infer_row, ]


#################
## generate plots 

color_cb=c("#000000", "#E69F00", "#56B4E9", "#009E73", 
           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plotdata = data.frame(time_res)
plotdata[, -which(colnames(plotdata)=='method')]=sapply(plotdata[, -which(colnames(plotdata)=='method')], as.numeric)
plotdata['y'] = plotdata$select_a + plotdata$est_and_inf + plotdata$vc
plotdata['group'] = plotdata$Time

time_plt = ggplot(data =plotdata, aes(x=n, y=y, group = interaction(group, method), 
                            linetype = method, color = factor(group),  shape = method))+
  geom_line(size = 1.3, alpha=0.5)+
  geom_point(size=3, alpha =0.5)+
  ylab(paste0('time (second)'))+
  scale_color_manual(name= 'Time', values =color_cb[1:length(unique(plotdata$group))])+
  scale_linetype_manual(values=c(1, 3, 4)) +
  labs(shape = 'method', linetype = 'method')+
  theme_bw()
  




plotdata = data.frame(beta_MSE_res)
plotdata[, -which(colnames(plotdata)=='method')]=sapply(plotdata[, -which(colnames(plotdata)=='method')], as.numeric)
plotdata['y'] = plotdata$beta_total_MSE
plotdata['group'] = plotdata$Time

beta_MSE_plt = ggplot(data =plotdata, aes(x=n, y=y, group = interaction(group, method), 
                                      linetype = method, color = factor(group),  shape = method))+
  geom_line(size = 1.3, alpha=0.5)+
  geom_point(size=3, alpha =0.5)+
  ylab((bquote(phi[1] ~" total MSE ")))+
  scale_color_manual(name= 'Time', values =color_cb[1:length(unique(plotdata$group))])+
  scale_linetype_manual(values=c(1, 3, 4)) +
  labs(shape = 'method', linetype = 'method')+
  theme_bw()


if(F){

plotdata = data.frame(vc_sd_res)
plotdata[, -which(colnames(plotdata)=='method')]=sapply(plotdata[, -which(colnames(plotdata)=='method')], as.numeric)
plotdata['group'] = plotdata$Time

vc_sd_MSE_plt = ggplot(data =plotdata, aes(x=n, y=vc_sd_total_MSE, group = interaction(group, method), 
                                          linetype = method, color = factor(group),  shape = method))+
  geom_line(size = 1.3, alpha=0.5)+
  geom_point(size=3, alpha =0.5)+
  ylab(paste0('random effect sd total MSE'))+
  scale_color_manual(name= 'Time', values =color_cb[1:length(unique(plotdata$group))])+
  scale_linetype_manual(values=c(1, 3, 4)) +
  labs(shape = 'method', linetype = 'method')+
  theme_bw()

vc_MCC_plt = ggplot(data =plotdata, aes(x=n, y=vc_sd_MCC, group = interaction(group, method), 
                           linetype = method, color = factor(group),  shape = method))+
  geom_line(size = 1.3, alpha=0.5)+
  geom_point(size=3, alpha =0.5)+
  ylab(paste0('random effect sd MCC'))+
  scale_color_manual(name= 'Time', values =color_cb[1:length(unique(plotdata$group))])+
  scale_linetype_manual(values=c(1, 3, 4)) +
  labs(shape = 'method', linetype = 'method')+
  theme_bw()


vc_e_sd_MSE_plt = ggplot(data =plotdata, aes(x=n, y=vc_e_sd_MSE, group = interaction(group, method), 
                           linetype = method, color = factor(group),  shape = method))+
  geom_line(size = 1.3, alpha=0.5)+
  geom_point(size=3, alpha =0.5)+
  ylab(paste0('error sd MSE'))+
  scale_color_manual(name= 'Time', values =color_cb[1:length(unique(plotdata$group))])+
  scale_linetype_manual(values=c(1, 3, 4)) +
  labs(shape = 'method', linetype = 'method')+
  theme_bw()
}

#-----------

for (get_plt_name in c('CI', 'power', 'typeI')){

plotdata = data.frame(get(paste0(get_plt_name,'_res')))
plotdata[, -which(colnames(plotdata)=='method')]=sapply(plotdata[, -which(colnames(plotdata)=='method')], as.numeric)
nodes = colnames(plotdata)[-c(1:6)]
plotdata['group'] = plotdata$Time



plt_list = lapply(nodes, function(sj){
  j = as.integer(substring(sj, 2))
  plt = ggplot(data =plotdata, aes(x=n, y=get(sj), group = interaction(group, method), 
                             linetype = method, color = factor(group),  shape = method))+
    geom_line(size = 1.3, alpha=0.5)+
    geom_point(size=3, alpha =0.5)+
    ylab(switch(get_plt_name, typeI = 'type I error', power = 'power', CI = 'CI coverage'))+
    scale_color_manual(name= 'Time', values =color_cb[1:length(unique(plotdata$group))])+
    scale_linetype_manual(values=c(1, 3, 4)) +
    labs(shape = 'method', linetype = 'method')+
    theme_bw()+
    ggtitle(bquote(beta[.(j)] ==.(true_beta[j]) ~", "~ 
                     psi[.(j)] == .(round(true_vc_sd[j], 2))~", "~
                   # Sigma[e~","~.(j)] == .(round(Sig_e[infer_row, j], 2))~","~
                   p == .(p)))
  if (get_plt_name == 'CI'){
    plt = plt + geom_abline(intercept=0.95, slope=0, color = color_cb[7])
  } 
  if (get_plt_name == 'typeI'){
    plt = plt + geom_abline(intercept=0.05, slope=0, color = color_cb[7])
  } 
  
  plt
})

ggsave(ggarrange( plotlist = plt_list, nrow=ceiling(length(plt_list)/2),
                  ncol = 2, common.legend = T), 
       filename = paste0('VAR', K, 'plot/row_', infer_row, '/', file_pre, '_',get_plt_name, '_test.pdf'), device = 'pdf', 
       width=10, height=4*ceiling(length(plt_list)/2), units = 'in', limitsize = F)
}


#-------
plt_list = list(
beta_MSE_plt, 
vc_sd_MSE_plt,
vc_MCC_plt,
vc_e_sd_MSE_plt,
time_plt
)

ggsave(ggarrange( plotlist = plt_list, nrow=ceiling(length(plt_list)/2),
                  ncol = 2, common.legend = T), 
       filename = paste0('plot/row_', infer_row, '/', file_pre, '_MSE_MCC_time_test.pdf'), device = 'pdf', 
       width=10, height=4*ceiling(length(plt_list)/2), units = 'in', limitsize = F)



###################################################
## inspect CI coverage failure

## one failure case is Time = 150 p=30 n=80
store_CI_all = NULL
store_beta_all = NULL
inspect = infer_row # 15


for(p in c(30)){
  for (n in c(20, 40, 60, 80)){
    for (Time in c(25, 50, 100, 150)){

     

store_CI = NULL
beta_est = NULL



for(iters in 1:total_iter){
  read_filename = paste0('res/VAR_', K, 'p', p, 'n', n, 'Time', Time, '/', file_pre, '_iter', iters, '_row_', infer_row, '.RData')
  
  tryCatch({
    
    load(read_filename)
    
    # print(beta[inspect])
    
    ## append proposed method results
    store_CI = rbind(store_CI, proposed_res$CI[,inspect])
    beta_est = rbind(beta_est, data.frame(
      betahat = proposed_res$beta.hat[inspect], 
      betadb = proposed_res$beta.db[inspect], 
      betadbsd = proposed_res$beta.db.sd[inspect],
      lassobetadb = lasso_res$beta.db[inspect],
      lassobetahat = lasso_res$beta.hat[inspect],
      Time=Time, n=n))



  }, error = function(e) print(paste('p', p, 'n', n, 'Time', Time, 'row', infer_row, 'iter', iters))
  )
  
}

store_CI = as.data.frame(store_CI)
store_CI$id = 1:nrow(store_CI)
colnames(store_CI)<- c('l', 'u', 'id')
store_CI$Time = Time
store_CI$n = n
store_CI$p = p

store_CI_all = rbind(store_CI_all, store_CI)

store_beta_all = rbind(store_beta_all, data.frame(beta_est))

plot_list = ggplot()+
  geom_hline(yintercept = beta[inspect])+
  geom_errorbar(data = store_CI, aes(ymin = l, ymax=u, x=id), width=.2)+
  ggtitle(paste('p', p, 'Time', Time, 'n', n))

    }}}


inspect_plot = ggplot(data = store_CI_all)+
  geom_hline(yintercept = beta[inspect])+
  geom_errorbar(aes(ymin = l, ymax=u, x=id), width=.2)+
  facet_wrap(~Time + n,labeller = label_both)
  
ggsave(inspect_plot, 
       filename = paste0('Vardplot/row_', infer_row, '/', file_pre, '_CI_inspection_plot_', infer_row, '_', inspect, '.pdf'), device = 'pdf', 
       width=15, height=15, units = 'in', limitsize = F)


inspect_plot2 = ggplot(data = store_beta_all)+
  geom_hline(yintercept = beta[inspect])+
  geom_boxplot(aes(y=betadb))+
  facet_wrap(~Time + n,labeller = label_both)

ggsave(inspect_plot2, 
       filename = paste0('Vardplot/row_', infer_row, '/', file_pre, '_betadb_inspection_plot_', infer_row, '_', inspect, '.pdf'), device = 'pdf', 
       width=15, height=15, units = 'in', limitsize = F)


inspect_plot2 = ggplot(data = store_beta_all)+
  geom_hline(yintercept = beta[inspect])+
  geom_boxplot(aes(y=betadb/betadbsd))+
  facet_wrap(~Time + n,labeller = label_both)

ggsave(inspect_plot2, 
       filename = paste0('Vardplot/row_', infer_row, '/', file_pre, '_betadb_inspection_plot_scaled_', infer_row, '_', inspect, '.pdf'), device = 'pdf', 
       width=15, height=15, units = 'in', limitsize = F)


ggplot(data = store_beta_all)+
  geom_hline(yintercept = beta[inspect])+
  geom_boxplot(aes(y=betadb/betadbsd))+
  facet_wrap(~Time + n,labeller = label_both, scales = 'free')


## how about the lasso estimator? this should be consistent!
plt = ggplot(data = store_beta_all)+
  geom_hline(yintercept = beta[inspect])+
  geom_boxplot(aes(y=lassobetahat))+
  facet_wrap(~Time + n,labeller = label_both)

ggsave(plt, 
       filename = paste0('Vardplot/row_', infer_row, '/', file_pre, '_lassohat.pdf'), device = 'pdf', 
       width=15, height=15, units = 'in', limitsize = F)


plt = ggplot(data = store_beta_all)+
  geom_hline(yintercept = beta[inspect])+
  geom_boxplot(aes(y=lassobetadb))+
  facet_wrap(~Time + n,labeller = label_both)

ggsave(plt, 
       filename = paste0('Vardplot/row_', infer_row, '/', file_pre, '_lassodb.pdf'), device = 'pdf', 
       width=15, height=15, units = 'in', limitsize = F)


plt = ggplot(data = store_beta_all)+
  geom_hline(yintercept = beta[inspect])+
  geom_boxplot(aes(y=betadb))+
  facet_wrap(~Time + n,labeller = label_both)

ggsave(plt, 
       filename = paste0('Vardplot/row_', infer_row, '/', file_pre, '_propdb.pdf'), device = 'pdf', 
       width=15, height=15, units = 'in', limitsize = F)

plt = ggplot(data = store_beta_all)+
  geom_hline(yintercept = beta[inspect])+
  geom_boxplot(aes(y=betahat))+
  facet_wrap(~Time + n,labeller = label_both)

ggsave(plt, 
       filename = paste0('Vardplot/row_', infer_row, '/', file_pre, '_prophat.pdf'), device = 'pdf', 
       width=15, height=15, units = 'in', limitsize = F)
