###########################
## test doubly HD LMM algorithm on MEVAR
#############################

rm(list=ls())

library(MASS)
library(tictoc)
library(hdi)
library(BigVAR)

test_on_box = F

if(test_on_box){
  inputs = c(1, 30, 80, 150, 3) #(iter, p, n, Time, infer_row)
  setwd('\\\\biostat-fs2\\users\\yuek\\Desktop\\MEVAR')
}else{
  inputs = commandArgs(T)
  setwd('~/Desktop/MEVAR')
}

ITER = as.integer(inputs[1])

source('lib.R')

a.seq=round(exp((seq(log(1e-4), log(100), length.out=6))),10)

cv.max.iter = 8
time_series_cv = F

## define parameters of VAR
p= as.integer(inputs[2]) ## dimension of response
n = as.integer(inputs[3]) ## number of subjects
Time = as.integer(inputs[4]) ## time length
K = 1 ## generate VAR(K) model

infer_row = as.integer(inputs[5]) 
inf.coord <- 1:p

#---------------
filename = paste0('res/VAR_', K, 'p', p, 'n', n, 'Time', Time)
if(!dir.exists(filename)){
  dir.create(filename)
}


file_pre = '_setting6'
#file_pre = '_setting6.2'


# if(file.exists(paste0(filename, '/', file_pre, '_iter', ITER, '_row_', infer_row, '.RData'))){
#   stop('done')
# }

#---------------



## initial state cov matrix
set.seed(1)
Sig_Y0 = generate_cov_matrix(type='sparse-random', param = c(p, 0.7, 0.2))

## noise term cov matrix
#set.seed(3)
#Sig_e = generate_cov_matrix(type='sparse-random', param = c(p, 0.7, 0.2)) ## original setting
#Sig_e = generate_cov_matrix(type='sparse-random', param = c(p, 0.7, 0)) ## diag setting
#Sig_e = generate_cov_matrix(type='sparse-random', param = c(p, 0.4, 0.35)) ## setting 2


## setting1: simple diagonal noise cov (value=0.5), and sparse diagonal random effect cov; set (1,1) of Phi being random
if(file_pre=='_setting1'){
  

## set Sig_e cov value
set.seed(3)
sig.e=0.5
Sig_e = generate_cov_matrix(type = 'diag', para = rep(sig.e, p)) ## diag setting, constant noise component

## set transition matrix population-level true values (sparse)
set.seed(5)
Phi_sparse_p = 0.8
Phi_mag = 0.2
Phi =  (matrix(runif(p*p), p, p)>Phi_sparse_p)* matrix(rnorm(p*p, mean=0, sd = Phi_mag), p, p)
diag(Phi) <- runif(p, min=0.2, max=0.8)
if(max(svd(Phi)$d)>=1){
  print('generated population Phi not stationary, processing...')
  Phi = Phi/2
}


## set random coefficient position (whether a transition coefficient is random is pre-set)
set.seed(13)
B_sparse_p = 0.9
B_pos = (matrix(runif(p*p), p, p)>B_sparse_p)*1
## random coefficient variances (fix the variance for each coefficient across subjects)
B_sd = matrix(runif(p*p, min=0.05, max=0.15), p, p) * B_pos ## this is the sd of the random effects; indicate sparse diagonal random effect cov

B_sd[1,1] <- 0.1 ## manually set this (1,1) entry to be varying with sd=0.1
B_pos[1,1] <- 1
}

## setting2: simple diagonal noise cov (value=0.5), and sparse diagonal random effect cov; set (1,1) of Phi being fixed
if(file_pre=='_setting2'){
  
  
  ## set Sig_e cov value
  set.seed(3)
  sig.e=0.5
  Sig_e = generate_cov_matrix(type = 'diag', para = rep(sig.e, p)) ## diag setting, constant noise component
  
  ## set transition matrix population-level true values (sparse)
  set.seed(5)
  Phi_sparse_p = 0.8
  Phi_mag = 0.2
  Phi =  (matrix(runif(p*p), p, p)>Phi_sparse_p)* matrix(rnorm(p*p, mean=0, sd = Phi_mag), p, p)
  diag(Phi) <- runif(p, min=0.2, max=0.8)
  if(max(svd(Phi)$d)>=1){
    print('generated population Phi not stationary, processing...')
    Phi = Phi/2
  }
  
  
  ## set random coefficient position (whether a transition coefficient is random is pre-set)
  set.seed(13)
  B_sparse_p = 0.9
  B_pos = (matrix(runif(p*p), p, p)>B_sparse_p)*1
  ## random coefficient variances (fix the variance for each coefficient across subjects)
  B_sd = matrix(runif(p*p, min=0.05, max=0.15), p, p) * B_pos ## this is the sd of the random effects; indicate sparse diagonal random effect cov
  
  B_sd[1,1] <- 0 ## manually set this (1,1) entry to be varying with sd=0.1
  B_pos[1,1] <- 0
}

## setting3: simple diagonal noise cov with value 1, and sparse diagonal random effect cov; set (1,1) of Phi being fixed
if(file_pre=='_setting3'){
  
  
  ## set Sig_e cov value
  set.seed(3)
  sig.e = 1
  Sig_e = generate_cov_matrix(type = 'diag', para = rep(sig.e, p)) ## diag setting, constant noise component
  
  ## set transition matrix population-level true values (sparse)
  set.seed(5)
  Phi_sparse_p = 0.8
  Phi_mag = 0.2
  Phi =  (matrix(runif(p*p), p, p)>Phi_sparse_p)* matrix(rnorm(p*p, mean=0, sd = Phi_mag), p, p)
  diag(Phi) <- runif(p, min=0.2, max=0.8)
  if(max(svd(Phi)$d)>=1){
    print('generated population Phi not stationary, processing...')
    Phi = Phi/2
  }
  
  
  ## set random coefficient position (whether a transition coefficient is random is pre-set)
  set.seed(13)
  B_sparse_p = 0.9
  B_pos = (matrix(runif(p*p), p, p)>B_sparse_p)*1
  ## random coefficient variances (fix the variance for each coefficient across subjects)
  B_sd = matrix(runif(p*p, min=0.05, max=0.15), p, p) * B_pos ## this is the sd of the random effects; indicate sparse diagonal random effect cov
  
  B_sd[1,1] <- 0 ## manually set this (1,1) entry to be varying with sd=0.1
  B_pos[1,1] <- 0
}

## setting4: change a.seq to be from 1 to 10; simple diagonal noise cov with value 1, and sparse diagonal random effect cov; set (1,1) of Phi being fixed
if(file_pre=='_setting4'){
  
  a.seq=round(exp((seq(log(1), log(10), length.out=10))),6)
  
  ## set Sig_e cov value
  set.seed(3)
  sig.e = 1
  Sig_e = generate_cov_matrix(type = 'diag', para = rep(sig.e, p)) ## diag setting, constant noise component
  
  ## set transition matrix population-level true values (sparse)
  set.seed(5)
  Phi_sparse_p = 0.8
  Phi_mag = 0.2
  Phi =  (matrix(runif(p*p), p, p)>Phi_sparse_p)* matrix(rnorm(p*p, mean=0, sd = Phi_mag), p, p)
  diag(Phi) <- runif(p, min=0.2, max=0.8)
  if(max(svd(Phi)$d)>=1){
    print('generated population Phi not stationary, processing...')
    Phi = Phi/2
  }
  
  
  ## set random coefficient position (whether a transition coefficient is random is pre-set)
  set.seed(13)
  B_sparse_p = 0.9
  B_pos = (matrix(runif(p*p), p, p)>B_sparse_p)*1
  ## random coefficient variances (fix the variance for each coefficient across subjects)
  B_sd = matrix(runif(p*p, min=0.05, max=0.15), p, p) * B_pos ## this is the sd of the random effects; indicate sparse diagonal random effect cov
  
  B_sd[1,1] <- 0 ## manually set this (1,1) entry to be varying with sd=0.1
  B_pos[1,1] <- 0
}


## setting5: change to fixed VAR1 model: B_sd = 0
if(file_pre=='_setting5'){
  
  ## set Sig_e cov value
  set.seed(3)
  sig.e = 1
  Sig_e = generate_cov_matrix(type = 'diag', para = rep(sig.e, p)) ## diag setting, constant noise component
  
  ## set transition matrix population-level true values (sparse)
  set.seed(5)
  Phi_sparse_p = 0.8
  Phi_mag = 0.2
  Phi =  (matrix(runif(p*p), p, p)>Phi_sparse_p)* matrix(rnorm(p*p, mean=0, sd = Phi_mag), p, p)
  diag(Phi) <- runif(p, min=0.2, max=0.8)
  if(max(svd(Phi)$d)>=1){
    print('generated population Phi not stationary, processing...')
    Phi = Phi/2
  }
  
  
  ## set random coefficient position (whether a transition coefficient is random is pre-set)
  set.seed(13)
  B_sparse_p = 2
  B_pos = (matrix(runif(p*p), p, p)>B_sparse_p)*1
  ## random coefficient variances (fix the variance for each coefficient across subjects)
  B_sd = matrix(runif(p*p, min=0.05, max=0.15), p, p) * B_pos ## this is the sd of the random effects; indicate sparse diagonal random effect cov
  
  B_sd[1,1] <- 0 ## manually set this (1,1) entry to be varying with sd=0.1
  B_pos[1,1] <- 0
}

## setting6: follow setting1, but restrict tuning modification to cv.max.iter = 2
if(file_pre=='_setting6'){
  
  cv.max.iter = 2
  
  ## set Sig_e cov value
  set.seed(3)
  sig.e=0.5
  Sig_e = generate_cov_matrix(type = 'diag', para = rep(sig.e, p)) ## diag setting, constant noise component
  
  ## set transition matrix population-level true values (sparse)
  set.seed(5)
  Phi_sparse_p = 0.8
  Phi_mag = 0.2
  Phi =  (matrix(runif(p*p), p, p)>Phi_sparse_p)* matrix(rnorm(p*p, mean=0, sd = Phi_mag), p, p)
  diag(Phi) <- runif(p, min=0.2, max=0.8)
  if(max(svd(Phi)$d)>=1){
    print('generated population Phi not stationary, processing...')
    Phi = Phi/2
  }
  
  
  ## set random coefficient position (whether a transition coefficient is random is pre-set)
  set.seed(13)
  B_sparse_p = 0.9
  B_pos = (matrix(runif(p*p), p, p)>B_sparse_p)*1
  ## random coefficient variances (fix the variance for each coefficient across subjects)
  B_sd = matrix(runif(p*p, min=0.05, max=0.15), p, p) * B_pos ## this is the sd of the random effects; indicate sparse diagonal random effect cov
  
  B_sd[1,1] <- 0.1 ## manually set this (1,1) entry to be varying with sd=0.1
  B_pos[1,1] <- 1
}

## setting6: follow setting2, but restrict tuning modification to cv.max.iter = 2
if(file_pre=='_setting6.2'){
  
  cv.max.iter = 2
  
  ## set Sig_e cov value
  set.seed(3)
  sig.e=0.5
  Sig_e = generate_cov_matrix(type = 'diag', para = rep(sig.e, p)) ## diag setting, constant noise component
  
  ## set transition matrix population-level true values (sparse)
  set.seed(5)
  Phi_sparse_p = 0.8
  Phi_mag = 0.2
  Phi =  (matrix(runif(p*p), p, p)>Phi_sparse_p)* matrix(rnorm(p*p, mean=0, sd = Phi_mag), p, p)
  diag(Phi) <- runif(p, min=0.2, max=0.8)
  if(max(svd(Phi)$d)>=1){
    print('generated population Phi not stationary, processing...')
    Phi = Phi/2
  }
  
  
  ## set random coefficient position (whether a transition coefficient is random is pre-set)
  set.seed(13)
  B_sparse_p = 0.9
  B_pos = (matrix(runif(p*p), p, p)>B_sparse_p)*1
  ## random coefficient variances (fix the variance for each coefficient across subjects)
  B_sd = matrix(runif(p*p, min=0.05, max=0.15), p, p) * B_pos ## this is the sd of the random effects; indicate sparse diagonal random effect cov
  
  B_sd[1,1] <- 0 ## manually set this (1,1) entry to 0
  B_pos[1,1] <- 0
}

## setting7: follow setting1, but restrict tuning modification to cv.max.iter = 8; use time series special CV; extend a sequence
if(file_pre=='_setting7'){
  time_series_cv = T
  
  a.seq=round(exp((seq(log(1e-7), log(100), length.out=10))),10)
  cv.max.iter = 8
  
  ## set Sig_e cov value
  set.seed(3)
  sig.e=0.5
  Sig_e = generate_cov_matrix(type = 'diag', para = rep(sig.e, p)) ## diag setting, constant noise component
  
  ## set transition matrix population-level true values (sparse)
  set.seed(5)
  Phi_sparse_p = 0.8
  Phi_mag = 0.2
  Phi =  (matrix(runif(p*p), p, p)>Phi_sparse_p)* matrix(rnorm(p*p, mean=0, sd = Phi_mag), p, p)
  diag(Phi) <- runif(p, min=0.2, max=0.8)
  if(max(svd(Phi)$d)>=1){
    print('generated population Phi not stationary, processing...')
    Phi = Phi/2
  }
  
  
  ## set random coefficient position (whether a transition coefficient is random is pre-set)
  set.seed(13)
  B_sparse_p = 0.9
  B_pos = (matrix(runif(p*p), p, p)>B_sparse_p)*1
  ## random coefficient variances (fix the variance for each coefficient across subjects)
  B_sd = matrix(runif(p*p, min=0.05, max=0.15), p, p) * B_pos ## this is the sd of the random effects; indicate sparse diagonal random effect cov
  
  B_sd[1,1] <- 0.1 ## manually set this (1,1) entry to be varying with sd=0.1
  B_pos[1,1] <- 1
}

## setting8: follow setting7, restrict value of a
if(file_pre=='_setting8'){
  time_series_cv = T
  
  a.seq=round(exp((seq(log(0.1), log(50), length.out=10))),10)
  cv.max.iter = 8
  
  ## set Sig_e cov value
  set.seed(3)
  sig.e=0.5
  Sig_e = generate_cov_matrix(type = 'diag', para = rep(sig.e, p)) ## diag setting, constant noise component
  
  ## set transition matrix population-level true values (sparse)
  set.seed(5)
  Phi_sparse_p = 0.8
  Phi_mag = 0.2
  Phi =  (matrix(runif(p*p), p, p)>Phi_sparse_p)* matrix(rnorm(p*p, mean=0, sd = Phi_mag), p, p)
  diag(Phi) <- runif(p, min=0.2, max=0.8)
  if(max(svd(Phi)$d)>=1){
    print('generated population Phi not stationary, processing...')
    Phi = Phi/2
  }
  
  
  ## set random coefficient position (whether a transition coefficient is random is pre-set)
  set.seed(13)
  B_sparse_p = 0.9
  B_pos = (matrix(runif(p*p), p, p)>B_sparse_p)*1
  ## random coefficient variances (fix the variance for each coefficient across subjects)
  B_sd = matrix(runif(p*p, min=0.05, max=0.15), p, p) * B_pos ## this is the sd of the random effects; indicate sparse diagonal random effect cov
  
  B_sd[1,1] <- 0.1 ## manually set this (1,1) entry to be varying with sd=0.1
  B_pos[1,1] <- 1
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
    risky_subj = c(risky_subj, i)
  }
  
  record_Bi[[i]] = B_i

  
  ## simulate for a single subject
  Yi = MultVarSim(k=p, A1=B_i+Phi, p=K, Sigma=Sig_e, T = Time)
  
  
  all_subj_data[[i]] = Yi
}

par(mfrow=c(3, 3))
for (w in 1:9){
  acf(all_subj_data[[1]][,w])
}


#########
## apply algorithm to infer the population-level transition matrix
#########


#########
## infer one row for the transition matrix

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


###################
cat("\nProposed method...\n")


tic()
a.opt = select_a(X=input_X, y = input_y, Z=input_X, grp=grp, a.seq=a.seq, 
                 kfold = 4, time_series_cv = time_series_cv , 
                 cv.max.iter = cv.max.iter, extra_para = list(K=K))
t1 = toc(quiet = T)

tic()
res.opta = Fixed_effects_est_inf(X=input_X, y=input_y, Z=input_X, grp=grp, a=a.opt, lm=F, inf.coord=inf.coord, 
                                 cv.max.iter = cv.max.iter,
                                 time_series_cv = time_series_cv, extra_para = list(K=K) )
t2 = toc(quiet = T)


tic()
res.vc.opta = VC_estimation(X=input_X, y=input_y, Z=input_X, grp=grp, a=a.opt, beta.hat = res.opta$beta.hat)
t3 = toc(quiet=T)

# tic()
# a.opt<-select_a_full_matrix(y=all_subj_data, a.seq = a.seq)
# t1 = toc(quiet = T)
# 
# 
# tic()
# res.opta = Fixed_effects_est_inf_full_matrix( y=all_subj_data, a=a.opt, inf.coord=inf.coord)
# t2 = toc(quiet = T)


# ## need more modification of the variance component estimation; do it later
# res.vc.opta = NULL
# tic()
# res.vc.opta = VC_estimation_full_matrix(y = all_subj_data, a=a.opt,grp=grp, beta.hat = res.opta$beta.hat)
# t3 = toc(quiet=T)

beta.db = res.opta$beta.db
beta.db.sd = res.opta$beta.db.sd
CI = sapply(1:length(inf.coord), function(j) 
  c((beta.db[j]-1.96 * beta.db.sd[j]) , ( beta.db[j]+ 1.96 * beta.db.sd[j]))
)
CI_cov = sapply(1:length(inf.coord), function(j) 
  (beta[inf.coord[j]]<= beta.db[j]+1.96 * beta.db.sd[j]) & (beta[inf.coord[j]]>= beta.db[j]- 1.96 * beta.db.sd[j]) # CI coverage
)
res_rej = abs(res.opta$beta.db/res.opta$beta.db.sd)>=1.96



proposed_res = list(beta.hat = res.opta$beta.hat,
                 beta.db = beta.db,
                 beta.db.sd = beta.db.sd,
                 CI = CI,
                 CI_cov = CI_cov,
                 rej = res_rej,
                 time = c('time_select_a' = t1$toc-t1$tic,
                          'time_est_inf' = t2$toc-t2$tic,
                          'time_vc' = t3$toc-t3$tic)
)


##############
## comparison: with simple lasso

cat("\nLasso method...\n")
tic()
cv.init<-cv.glmnet(input_X, input_y, lambda= NULL)
## inference based on de-sparsified lasso
lasso_fit_res = lasso.proj(input_X, input_y, multiplecorr.method = "none", verbose=F)
lasso_ci_res = confint(lasso_fit_res, level=0.95)[inf.coord,]
lasso_t1 = toc(quiet=T)

lasso_res = list(beta.hat = coef(cv.init, s=cv.init$lambda.min, exact=T)[-1],
                 beta.db = lasso_fit_res$bhat[inf.coord],
                 beta.db.sd = lasso_fit_res$se[inf.coord],
                 CI = confint(lasso_fit_res, level=0.95)[inf.coord,],
                 CI_cov = sapply(1:length(inf.coord), function(j) 
                   (beta[inf.coord[j]]<= lasso_ci_res[j,2] & beta[inf.coord[j]]>= lasso_ci_res[j,1])),
                 rej = lasso_fit_res$pval[inf.coord] <=0.05,
                 time = lasso_t1$toc - lasso_t1$tic
                 )


######################
## save the results


vc.sig.e.sd <- sqrt(res.vc.opta$sig2.eps.hat)
vc.beta.sd <- sqrt(res.vc.opta$vc.hat)




save(list = c('Phi', 'B_pos', 'B_sd', 'Sig_e', 'beta', 'infer_row', 'inf.coord',
    'a.opt', 'res.opta', 'res.vc.opta', 
    'proposed_res',
    'lasso_res',
    'risky_subj'
    # , 'record_Bi'
    ), 
     file = paste0(filename, '/', file_pre, '_iter', ITER, '_row_', infer_row, '.RData'))

















