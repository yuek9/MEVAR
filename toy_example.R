###########################
## toy example to compare two-stage analysis with MEVAR (low dimension)
#############################

rm(list=ls())

library(MASS)
library(tictoc)
library(hdi)
library(BigVAR)

library(vars)
library(poolr)
library(lme4)

library(pbkrtest)


test_on_box = F

if(test_on_box){
  inputs = c(1, 3, 50, 150, 1) #(iter, p, n, Time, infer_row)
  setwd('\\\\biostat-fs2\\users\\yuek\\Desktop\\MEVAR')
}else{
  inputs = commandArgs(T)
  setwd('~/Desktop/MEVAR')
}

ITER = as.integer(inputs[1])

#for (ITER in 1:200)

  {

cat('iter', ITER)
source('lib.R')

a.seq=round(exp((seq(log(1e-4), log(100), length.out=6))),10)

cv.max.iter = 2
time_series_cv = T

## define parameters of VAR
p= as.integer(inputs[2]) ## dimension of response
n = as.integer(inputs[3]) ## number of subjects
Time = as.integer(inputs[4]) ## time length
K = 1 ## generate VAR(K) model

infer_row = as.integer(inputs[5]) 
inf.coord <- 1:p

#---------------
filename = paste0('res/toy_example_VAR_', K, 'p', p, 'n', n, 'Time', Time)
if(!dir.exists(filename)){
  dir.create(filename)
}


file_pre = '_setting1'

if(file.exists(paste0(filename, '/', file_pre, '_iter', ITER, '_row_', infer_row, '.RData'))){
  stop('finished')
}
  
  
  
## setting1: simple diagonal noise cov (value=0.5), and sparse diagonal random effect cov; set (1,1) of Phi being random
if(file_pre=='_setting1'){

## set Sig_e cov value
set.seed(3)
sig.e=0.5
Sig_e = generate_cov_matrix(type = 'diag', para = rep(sig.e, p)) ## diag setting, constant noise component

## set transition matrix population-level true values (sparse)
Phi = matrix(c(0.4, 0, 0.1, 0, 0.6, -0.4, 0.04, 0, 0.8), 3, 3)

if(max(svd(Phi)$d)>=1){
  print('generated population Phi not stationary, processing...')
  Phi = Phi/2
}


## set random coefficient position (whether a transition coefficient is random is pre-set)
set.seed(1)
#B_sd = matrix(c(0.05, 0, 0, 0.04, 0.05, 0, 0.12, 0.05, 0 ), 3, 3)
B_sd = matrix(c(0.05, 0, 0, 0.04, 0.05, 0, 0.08, 0.05, 0 ), 3, 3)
B_pos = B_sd !=0
}


## this is the true value for beta, with picked node as regressor
beta = Phi[infer_row,]

## this is the random effect sd (for the diagonal)
random_cov_diag = B_sd[infer_row,]

#####################
## work under VAR1 setting
#####################


## first generate data from VAR model
all_subj_data = list()
risky_subj = NULL
record_Bi = list()

for (i in 1:n){
  set.seed(ITER*n*2*p*Time + i)
  
  ## for each subject, generate perturbation of the transition matrix; need to make sure the generated coefficient matrix gives a stationary VAR
  B_i = diag(2, p)
  counter=0
  while(max(svd(B_i +  Phi)$d)>=1 & counter < 10){
    counter = counter +1
    B_i = matrix(rnorm(n=p*p, mean=0, sd = as.vector(B_sd)), p, p) * B_pos
  }
  
  if(counter > 1){
    print(paste('subj', i, 'not stationary'))
  }
  
  if(counter ==10 ){
    print(paste('subj', i, 'error'))
    risky_subj = c(risky_subj, i)
  }
  
  record_Bi[[i]] = B_i

  
  ## simulate for a single subject
  Yi = MultVarSim(k=p, A1=B_i+Phi, p=K, Sigma=Sig_e, T = Time)
  colnames(Yi)<-paste0('y', 1:p)
  
  all_subj_data[[i]] = Yi
}

print(risky_subj)

Reduce('+', record_Bi)/n ## this should be very close to 0

#########
## apply algorithm to infer the population-level transition matrix
#########


#########
## apply VAR separately to each individual, and then infer the (1,1) connection
all_est = list()
all_pval = list()

for(i in 1:n){
  
fit_res = VAR(all_subj_data[[i]], p=1, type='none')
tmp = summary(fit_res)

for (target in list(c(1,1), c(1,2), c(1,3))){

est = tmp$varresult[[(paste0('y', target[1]))]]$coefficients[target[2], 1]
pval = tmp$varresult[[(paste0('y', target[1]))]]$coefficients[target[2], 4]

all_est[[paste0(target, collapse = ',')]] = c(all_est[[paste0(target, collapse = ',')]], est)
all_pval[[paste0(target, collapse = ',')]] = c(all_pval[[paste0(target, collapse = ',')]], pval)
}
}

## test using estimates and t test
two_stage_t_test = c('1'=t.test(all_est$`1,1`)$p.value,
'2'=t.test(all_est$`1,2`)$p.value,
'3'=t.test(all_est$`1,3`)$p.value)

## test using Fisher's method on pvalues
two_stage_fisher = c('1' = fisher(all_pval$`1,1`)$p,
'2' = fisher(all_pval$`1,2`)$p,
'3' = fisher(all_pval$`1,3`)$p)

#########
## infer one row for the transition matrix using MEVAR

## prepare input data for the algorithm
Time_list = sapply(all_subj_data, nrow)

input_y = do.call(c, lapply(1:n, function(i){
  all_subj_data[[i]][(K+1):(Time_list[i]), infer_row]
}
))

input_X = do.call(rbind, lapply(1:n, function(i){
  all_subj_data[[i]][1:(Time_list[i]-K), ]
})
)

## input Z=X
grp = do.call(c, lapply(1:n, function(i) rep(i, Time_list[i]-K)))

### note that in lme, model may be singular and gives error
lme_F_test <- lme_t_test <- rep(NA, 3)
 
tryCatch({
data = data.frame(cbind(input_y, input_X, grp))
colnames(data) <-c('y', 'y1', 'y2', 'y3', 'subj')

# fit_res = lme( y~ -1 + y1+ y2+ y3,
#               random = ~ -1 + y1 + y2 + y3| subj,
#               data=data,
#               method = 'REML')
fit_res = lmer(y ~ -1 + y1 + y2 + y3 + 
                 (0 + y1 || subj) + 
                 (0 + y2 || subj) + 
                 (0 + y3 || subj), data=data, REML=TRUE)

tmp = summary(fit_res)
#lme_t_test = tmp$tTable[,'p-value']

lme_t_test = 2*(1-pnorm(abs(tmp$coefficients[,'t value'])))


## another option is F-test, supposed to be more accurate
fit_res1 = lmer(y ~ -1 + y2 + y3 + (0 + y1 || subj) + (0 + y2 || subj) + (0 + y3 || subj), data=data, REML=TRUE)
fit_res2 = lmer(y ~ -1 + y1 + y3 + (0 + y1 || subj) + (0 + y2 || subj) + (0 + y3 || subj), data=data, REML=TRUE)
fit_res3 = lmer(y ~ -1 + y1 + y2 + (0 + y1 || subj) + (0 + y2 || subj) + (0 + y3 || subj), data=data, REML=TRUE)

lme_F_test = c(KRmodcomp(fit_res, fit_res1)$test[1, 5],
               KRmodcomp(fit_res, fit_res2)$test[1, 5],
               KRmodcomp(fit_res, fit_res3)$test[1, 5])

}, error = function(e) print(e))

################
## our method

a.opt = select_a(X=input_X, y = input_y, Z=input_X, grp=grp, a.seq=a.seq, 
                 kfold = 4, time_series_cv = time_series_cv , 
                 cv.max.iter = cv.max.iter, extra_para = list(K=K))

res.opta = Fixed_effects_est_inf(X=input_X, y=input_y, Z=input_X, grp=grp, a=a.opt, lm=F, inf.coord=inf.coord, 
                                 cv.max.iter = cv.max.iter,
                                 time_series_cv = time_series_cv, extra_para = list(K=K) )

beta.db = res.opta$beta.db
beta.db.sd = res.opta$beta.db.sd
CI = sapply(1:length(inf.coord), function(j) 
  c((beta.db[j]-1.96 * beta.db.sd[j]) , ( beta.db[j]+ 1.96 * beta.db.sd[j]))
)
CI_cov = sapply(1:length(inf.coord), function(j) 
  (beta[inf.coord[j]]<= beta.db[j]+1.96 * beta.db.sd[j]) & (beta[inf.coord[j]]>= beta.db[j]- 1.96 * beta.db.sd[j]) # CI coverage
)
res_rej = abs(res.opta$beta.db/res.opta$beta.db.sd)>=1.96


save_res = list(two_stage_fisher = two_stage_fisher,
                two_stage_t_test = two_stage_t_test,
                lme_t_test = lme_t_test,
                lme_F_test = lme_F_test,
                proposed = res_rej)

save(list = c('Phi', 'B_pos', 'B_sd', 'save_res', 'CI_cov'
    ), 
     file = paste0(filename, '/', file_pre, '_iter', ITER, '_row_', infer_row, '.RData'))


}


########################
###

if(test_on_box){
  
  res_all_iter = list()
  infer_row=1
  K=1
  p=3
  n=50
  Time=150
  file_pre = '_setting1'
  filename = paste0('res/toy_example_VAR_', K, 'p', p, 'n', n, 'Time', Time)
  
  check_CI = NULL
  
  for(ITER in 1:200){
    load(paste0(filename, '/', file_pre, '_iter', ITER, '_row_', infer_row, '.RData'))
    res_all_iter[[ITER]] = save_res
    
    check_CI = rbind(check_CI, CI_cov)
  }
  
  ## look at the power for testing each edge
  ## two stage fisher method
  colMeans(do.call(rbind, sapply(res_all_iter, `[`, 1)) < 0.05)
  
  ## two stage t test method
  colMeans(do.call(rbind, sapply(res_all_iter, `[`, 2)) < 0.05)
  
  ## lme method: t-test
  colMeans(do.call(rbind, sapply(res_all_iter, `[`, 3)) < 0.05, na.rm = T)
  
  ## lme method: F-test
  colMeans(do.call(rbind, sapply(res_all_iter, `[`, 4)) < 0.05, na.rm = T)
  
  ## proposed method
  colMeans(do.call(rbind, sapply(res_all_iter, `[`, 5)))
  
  
  ## check CI coverage
  colMeans(check_CI)
}










