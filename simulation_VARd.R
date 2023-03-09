##################
## extension to VAR(d) model
#################

## only difference is the processing of the input data matrices

###########################
## test doubly HD LMM algorithm on MEVAR
#############################

rm(list=ls())

library(MASS)
library(tictoc)
library(hdi)
library(BigVAR)

test_on_box = F
K=2 # try a VAR(3) model or VAR(2) model

if(test_on_box){
  inputs = c(1, 30, 80, 150, 1) #(iter, p, n, Time, infer_row)
  setwd('\\\\biostat-fs2\\users\\yuek\\Desktop\\MEVAR')
}else{
  inputs = commandArgs(T)
  setwd('~/Desktop/MEVAR')
}

ITER = as.integer(inputs[1])

source('lib.R')

a.seq=round(exp((seq(log(1e-4), log(100), length.out=8))),10)

cv.max.iter = 3
time_series_cv = T

## define parameters of VAR
p= as.integer(inputs[2]) ## dimension of response
n = as.integer(inputs[3]) ## number of subjects
Time = as.integer(inputs[4]) ## time length


infer_row = as.integer(inputs[5]) 


#---------------
filename = paste0('res/VAR_', K, 'p', p, 'n', n, 'Time', Time)
if(!dir.exists(filename)){
  dir.create(filename)
}

file_setting = as.integer(inputs[6])
file_pre = paste0('_setting', file_setting)

# file_pre = '_setting3' ## use setting4 a.seq, but cv.iter=1
# file_pre = '_setting4' ## setting7 + cv.iter=12 + even longer a.seq
# file_pre = '_setting5' ## setting7 + cv.iter=10 + longer a.seq
# file_pre = '_setting6' ## setting7+cv.iter=6
# file_pre = '_setting7' ## Ali's suggestion: make the second lag more sparse; second lag diagonal all zero


# if(file.exists(paste0(filename, '/', file_pre, '_iter', ITER, '_row_', infer_row, '.RData'))){
#   stop('done')
# }

#---------------


## setting3: 
if(file_pre=='_setting3'){
  
  cv.max.iter = 1
  a.seq=round(exp((seq(log(1e-7), log(1e4), length.out=8))),40)
  
  ## set Sig_e cov value
  set.seed(3)
  sig.e=0.5
  Sig_e = generate_cov_matrix(type = 'diag', para = rep(sig.e, p)) ## diag setting, constant noise component
  
  ## set transition matrix population-level true values (sparse)
  set.seed(5)
  Phi_sparse_p = c(0.8, 0.9, 0.95)
  Phi_mag = c(0.25, 0.15, 0.1)
  diag_max=c(0.8, 0.4, 0.2)
  
  Phi_list = lapply(1:K, function(w){
    tmp = (matrix(runif(p*p), p, p)>Phi_sparse_p[w])* matrix(rnorm(p*p, mean=0, sd = Phi_mag[w]), p, p)
    diag(tmp) <- runif(p, min=0.1, max=diag_max[w])
    tmp
  }
  )
  
  ## make second lag diagonal all zero
  
  diag(Phi_list[[2]]) <- 0
  
  companion_Phi = VarptoVar1MC(do.call(cbind, Phi_list), p=length(Phi_list), k=dim(Phi_list[[1]])[1])
  
  
  
  while(max(abs(eigen(companion_Phi)$value))>=1){
    print('generated population Phi not stationary, processing...')
    for(w in 1:K){Phi_list[[w]] = Phi_list[[w]]/2}
    companion_Phi = VarptoVar1MC(do.call(cbind, Phi_list), p=length(Phi_list), k=dim(Phi_list[[1]])[1])
  }
  
  
  ## set random coefficient position (whether a transition coefficient is random is pre-set)
  set.seed(13)
  B_sparse_p = c(0.7, 0.8, 0.95)
  B_pos_list = lapply(1:K, function(w){
    (matrix(runif(p*p), p, p)>B_sparse_p[w])*1})
  ## random coefficient variances (fix the variance for each coefficient across subjects)
  B_sd_max = c(0.18, 0.08, 0.05)
  B_sd_list = lapply(1:K, function(w) 
    matrix(runif(p*p, min=0.02, max=B_sd_max[w]), p, p) * B_pos_list[[w]]) ## this is the sd of the random effects; indicate sparse diagonal random effect cov
  
  B_sd_list[[1]][1,1] <- 0.1 ## manually set this (1,1) entry to be varying with sd=0.1
  B_pos_list[[1]][1,1] <- 1
  
  B_sd_companion = rbind(do.call(cbind, B_sd_list), matrix(0, nrow=p*(K-1), ncol=p*K))
  B_pos_companion = B_sd_companion>0
  
}

## setting4: 
if(file_pre=='_setting4'){
  
  cv.max.iter = 12
  a.seq=round(exp((seq(log(1e-7), log(1e4), length.out=8))),40)
  
  ## set Sig_e cov value
  set.seed(3)
  sig.e=0.5
  Sig_e = generate_cov_matrix(type = 'diag', para = rep(sig.e, p)) ## diag setting, constant noise component
  
  ## set transition matrix population-level true values (sparse)
  set.seed(5)
  Phi_sparse_p = c(0.8, 0.9, 0.95)
  Phi_mag = c(0.25, 0.15, 0.1)
  diag_max=c(0.8, 0.4, 0.2)
  
  Phi_list = lapply(1:K, function(w){
    tmp = (matrix(runif(p*p), p, p)>Phi_sparse_p[w])* matrix(rnorm(p*p, mean=0, sd = Phi_mag[w]), p, p)
    diag(tmp) <- runif(p, min=0.1, max=diag_max[w])
    tmp
  }
  )
  
  ## make second lag diagonal all zero
  
  diag(Phi_list[[2]]) <- 0
  
  companion_Phi = VarptoVar1MC(do.call(cbind, Phi_list), p=length(Phi_list), k=dim(Phi_list[[1]])[1])
  
  
  
  while(max(abs(eigen(companion_Phi)$value))>=1){
    print('generated population Phi not stationary, processing...')
    for(w in 1:K){Phi_list[[w]] = Phi_list[[w]]/2}
    companion_Phi = VarptoVar1MC(do.call(cbind, Phi_list), p=length(Phi_list), k=dim(Phi_list[[1]])[1])
  }
  
  
  ## set random coefficient position (whether a transition coefficient is random is pre-set)
  set.seed(13)
  B_sparse_p = c(0.7, 0.8, 0.95)
  B_pos_list = lapply(1:K, function(w){
    (matrix(runif(p*p), p, p)>B_sparse_p[w])*1})
  ## random coefficient variances (fix the variance for each coefficient across subjects)
  B_sd_max = c(0.18, 0.08, 0.05)
  B_sd_list = lapply(1:K, function(w) 
    matrix(runif(p*p, min=0.02, max=B_sd_max[w]), p, p) * B_pos_list[[w]]) ## this is the sd of the random effects; indicate sparse diagonal random effect cov
  
  B_sd_list[[1]][1,1] <- 0.1 ## manually set this (1,1) entry to be varying with sd=0.1
  B_pos_list[[1]][1,1] <- 1
  
  B_sd_companion = rbind(do.call(cbind, B_sd_list), matrix(0, nrow=p*(K-1), ncol=p*K))
  B_pos_companion = B_sd_companion>0
  
}

## setting=5: 
if(file_pre=='_setting5'){
  
  cv.max.iter = 10
  a.seq=round(exp((seq(log(1e-5), log(1e3), length.out=8))),25)
  
  ## set Sig_e cov value
  set.seed(3)
  sig.e=0.5
  Sig_e = generate_cov_matrix(type = 'diag', para = rep(sig.e, p)) ## diag setting, constant noise component
  
  ## set transition matrix population-level true values (sparse)
  set.seed(5)
  Phi_sparse_p = c(0.8, 0.9, 0.95)
  Phi_mag = c(0.25, 0.15, 0.1)
  diag_max=c(0.8, 0.4, 0.2)
  
  Phi_list = lapply(1:K, function(w){
    tmp = (matrix(runif(p*p), p, p)>Phi_sparse_p[w])* matrix(rnorm(p*p, mean=0, sd = Phi_mag[w]), p, p)
    diag(tmp) <- runif(p, min=0.1, max=diag_max[w])
    tmp
  }
  )
  
  ## make second lag diagonal all zero
  
  diag(Phi_list[[2]]) <- 0
  
  companion_Phi = VarptoVar1MC(do.call(cbind, Phi_list), p=length(Phi_list), k=dim(Phi_list[[1]])[1])
  
  
  
  while(max(abs(eigen(companion_Phi)$value))>=1){
    print('generated population Phi not stationary, processing...')
    for(w in 1:K){Phi_list[[w]] = Phi_list[[w]]/2}
    companion_Phi = VarptoVar1MC(do.call(cbind, Phi_list), p=length(Phi_list), k=dim(Phi_list[[1]])[1])
  }
  
  
  ## set random coefficient position (whether a transition coefficient is random is pre-set)
  set.seed(13)
  B_sparse_p = c(0.7, 0.8, 0.95)
  B_pos_list = lapply(1:K, function(w){
    (matrix(runif(p*p), p, p)>B_sparse_p[w])*1})
  ## random coefficient variances (fix the variance for each coefficient across subjects)
  B_sd_max = c(0.18, 0.08, 0.05)
  B_sd_list = lapply(1:K, function(w) 
    matrix(runif(p*p, min=0.02, max=B_sd_max[w]), p, p) * B_pos_list[[w]]) ## this is the sd of the random effects; indicate sparse diagonal random effect cov
  
  B_sd_list[[1]][1,1] <- 0.1 ## manually set this (1,1) entry to be varying with sd=0.1
  B_pos_list[[1]][1,1] <- 1
  
  B_sd_companion = rbind(do.call(cbind, B_sd_list), matrix(0, nrow=p*(K-1), ncol=p*K))
  B_pos_companion = B_sd_companion>0
  
}

## setting=6: 
if(file_pre=='_setting6'){
  
  cv.max.iter =6
  # a.seq=round(exp((seq(log(1e-5), log(100), length.out=8))),10)
  
  ## set Sig_e cov value
  set.seed(3)
  sig.e=0.5
  Sig_e = generate_cov_matrix(type = 'diag', para = rep(sig.e, p)) ## diag setting, constant noise component
  
  ## set transition matrix population-level true values (sparse)
  set.seed(5)
  Phi_sparse_p = c(0.8, 0.9, 0.95)
  Phi_mag = c(0.25, 0.15, 0.1)
  diag_max=c(0.8, 0.4, 0.2)
  
  Phi_list = lapply(1:K, function(w){
    tmp = (matrix(runif(p*p), p, p)>Phi_sparse_p[w])* matrix(rnorm(p*p, mean=0, sd = Phi_mag[w]), p, p)
    diag(tmp) <- runif(p, min=0.1, max=diag_max[w])
    tmp
  }
  )
  
  ## make second lag diagonal all zero
  
  diag(Phi_list[[2]]) <- 0
  
  companion_Phi = VarptoVar1MC(do.call(cbind, Phi_list), p=length(Phi_list), k=dim(Phi_list[[1]])[1])
  
  
  
  while(max(abs(eigen(companion_Phi)$value))>=1){
    print('generated population Phi not stationary, processing...')
    for(w in 1:K){Phi_list[[w]] = Phi_list[[w]]/2}
    companion_Phi = VarptoVar1MC(do.call(cbind, Phi_list), p=length(Phi_list), k=dim(Phi_list[[1]])[1])
  }
  
  
  ## set random coefficient position (whether a transition coefficient is random is pre-set)
  set.seed(13)
  B_sparse_p = c(0.7, 0.8, 0.95)
  B_pos_list = lapply(1:K, function(w){
    (matrix(runif(p*p), p, p)>B_sparse_p[w])*1})
  ## random coefficient variances (fix the variance for each coefficient across subjects)
  B_sd_max = c(0.18, 0.08, 0.05)
  B_sd_list = lapply(1:K, function(w) 
    matrix(runif(p*p, min=0.02, max=B_sd_max[w]), p, p) * B_pos_list[[w]]) ## this is the sd of the random effects; indicate sparse diagonal random effect cov
  
  B_sd_list[[1]][1,1] <- 0.1 ## manually set this (1,1) entry to be varying with sd=0.1
  B_pos_list[[1]][1,1] <- 1
  
  B_sd_companion = rbind(do.call(cbind, B_sd_list), matrix(0, nrow=p*(K-1), ncol=p*K))
  B_pos_companion = B_sd_companion>0
  
}


## setting=7: 
if(file_pre=='_setting7'){
  
  cv.max.iter = 2
  # a.seq=round(exp((seq(log(1e-5), log(100), length.out=8))),10)
  
  ## set Sig_e cov value
  set.seed(3)
  sig.e=0.5
  Sig_e = generate_cov_matrix(type = 'diag', para = rep(sig.e, p)) ## diag setting, constant noise component
  
  ## set transition matrix population-level true values (sparse)
  set.seed(5)
  Phi_sparse_p = c(0.8, 0.9, 0.95)
  Phi_mag = c(0.25, 0.15, 0.1)
  diag_max=c(0.8, 0.4, 0.2)
  
  Phi_list = lapply(1:K, function(w){
    tmp = (matrix(runif(p*p), p, p)>Phi_sparse_p[w])* matrix(rnorm(p*p, mean=0, sd = Phi_mag[w]), p, p)
    diag(tmp) <- runif(p, min=0.1, max=diag_max[w])
    tmp
  }
  )
  
  ## make second lag diagonal all zero
  
  diag(Phi_list[[2]]) <- 0
  
  companion_Phi = VarptoVar1MC(do.call(cbind, Phi_list), p=length(Phi_list), k=dim(Phi_list[[1]])[1])
  
  
  
  while(max(abs(eigen(companion_Phi)$value))>=1){
    print('generated population Phi not stationary, processing...')
    for(w in 1:K){Phi_list[[w]] = Phi_list[[w]]/2}
    companion_Phi = VarptoVar1MC(do.call(cbind, Phi_list), p=length(Phi_list), k=dim(Phi_list[[1]])[1])
  }
  
  
  ## set random coefficient position (whether a transition coefficient is random is pre-set)
  set.seed(13)
  B_sparse_p = c(0.7, 0.8, 0.95)
  B_pos_list = lapply(1:K, function(w){
    (matrix(runif(p*p), p, p)>B_sparse_p[w])*1})
  ## random coefficient variances (fix the variance for each coefficient across subjects)
  B_sd_max = c(0.18, 0.08, 0.05)
  B_sd_list = lapply(1:K, function(w) 
    matrix(runif(p*p, min=0.02, max=B_sd_max[w]), p, p) * B_pos_list[[w]]) ## this is the sd of the random effects; indicate sparse diagonal random effect cov
  
  B_sd_list[[1]][1,1] <- 0.1 ## manually set this (1,1) entry to be varying with sd=0.1
  B_pos_list[[1]][1,1] <- 1
  
  B_sd_companion = rbind(do.call(cbind, B_sd_list), matrix(0, nrow=p*(K-1), ncol=p*K))
  B_pos_companion = B_sd_companion>0
  
}


## this is the true value for beta, with picked node as regressor; stack from Phi1 to PhiK
beta = do.call(c, lapply(Phi_list, function(X)X[infer_row,]))

## this is the random effect sd (for the diagonal)
random_cov_diag = do.call(c, lapply(B_sd_list, function(X) X[infer_row,]))

## the inf.coord input will be a list, k-th component corresponds to the k-th lag coefficient matrix, then the vector represent the (inferr_row,l) entry of the coefficient matrix to infer
inf.coord_raw = replicate(K, 1:p, simplify = F)
names(inf.coord_raw) <- as.character(1:K)

## transform into the inf.coord for the stacked coefficient vector
inf.coord = do.call(c, lapply(names(inf.coord_raw), function(s)inf.coord_raw[[s]] + (as.numeric(s)-1) * p))
  
  


#####################
## work under VAR(K) setting
#####################

  
## first generate data from VAR model
all_subj_data = list()
risky_subj = NULL
record_Bi = list()

for (i in 1:n){
  set.seed(ITER*n*2*p*Time + i)
  
  ## for each subject, generate perturbation of the transition matrix; need to make sure the generated coefficient matrix gives a stationary VAR
  B_i = diag(100, dim(B_sd_companion)[1])
    
    
    
  counter=0
  while(max(abs(eigen(B_i + companion_Phi)$value))>=1 & counter < 10){
    counter = counter +1
    B_i = matrix(
      rnorm(n=(dim(B_sd_companion)[1])^2, mean=0, sd = as.vector(B_sd_companion)), 
      p*K, p*K) 
    * B_pos_companion
  }
  
  if(counter > 1){
    print(paste('subj', i, 'not stationary'))
    risky_subj = c(risky_subj, i)
  }
  
  record_Bi[[i]] = B_i
  
  
  ## simulate for a single subject
  Yi = MultVarSim(k=p, A1= companion_Phi + B_i, 
                  p=K, Sigma=Sig_e, T = Time)
  
  
  all_subj_data[[i]] = Yi
}

## want some higher serial correlation in the data
par(mfrow=c(3, 3))
for (w in 1:9){
  acf(all_subj_data[[2]][,w], main=w)
}

par(mfrow=c(3, 3))
for (w in 1:9){
  pacf(all_subj_data[[2]][,w], main=w)
}

##################
## prepare input for VAR(K) model

input_X = lapply(all_subj_data, function(Y){
  input_X_i = t(sapply(0:(nrow(Y)-K-1), function(ii){
    as.vector(t(Y[ii + K:1,]))
  }))
  input_X_i
})

input_y = lapply(all_subj_data, function(Y){
  Y[(K+1):nrow(Y),infer_row]
})

grp = do.call(c, lapply(1:n, function(i) rep(i, nrow(input_X[[i]]))))
input_X = do.call(rbind, input_X)
input_y = do.call(c, input_y)

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
res.vc.opta = NULL # = VC_estimation(X=input_X, y=input_y, Z=input_X, grp=grp, a=a.opt, beta.hat = res.opta$beta.hat)
t3 = toc(quiet=T)


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


# vc.sig.e.sd <- sqrt(res.vc.opta$sig2.eps.hat)
# vc.beta.sd <- sqrt(res.vc.opta$vc.hat)




save(list = c('Phi_list', 'B_pos_list', 'B_sd_list', 
              'companion_Phi', 'B_sd_companion',
              'Sig_e', 'beta', 'infer_row', 'inf.coord', 'inf.coord_raw', 
              'a.opt', 'res.opta', 'res.vc.opta', 
              'proposed_res',
              'lasso_res',
              'risky_subj'
              # , 'record_Bi'
), 
file = paste0(filename, '/', file_pre, '_iter', ITER, '_row_', infer_row, '.RData'))




