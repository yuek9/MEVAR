################################
## estimate fixed effect and variance components in a MEVAR(1) model

## Kun Yue
## update: 10/11/2022
## brute force of SVD for sigmai is slow; utilize the property of kronecker product SVD
################################
library(caret)
library(glmnet)
library(doParallel)
library(foreach)

#######################
## main function for fixed effect estimation and inference
#######################

###################
## the function that only estimate the inference one row of the coefficient matrix Phi (can directly use the previous LMM function)
###################

#Inputs-- 
## X: design for fixed effects; 
## y: response; 
## Z: design for random effects;
## grp: a vector of length N indicating group membership; has to be from 1 to n without gaps
## a : tuning parameter in $\Sig_a$
## lm: linear model fitting (T or F, if T will ignore correlations among observations); 
## inf.coord: the coordinates of beta for inference.

#Outputs-- 
## beta.hat: the Lasso estimator based on psedo-likelihood; 
## beta.db: debiased Lasso estimators for the fixed effects in inf.coord;
## beta.db.sd: standard deviation of the debiased Lasso estimators for the fixed effects in inf.coord.


Fixed_effects_est_inf <- function(X, y, Z, grp, a=1, lm=F, inf.coord=NULL, 
                                  cv.max.iter = 8, 
                                  time_series_cv = T, ## for time series, we can use this special CV setting
                                  extra_para = list(
                                    K=NULL ## K is number of lags in the VAR model
                                    ),
                                  true_cov ## add true VAR based covariance matrix for input Y
                                  ){
  
  N=length(y)
  n=length(unique(grp))
  q=ncol(Z)
  p=ncol(X)
  grp_values = unique(grp)
  
  ## include extra parameter for time series setting
  m_list = table(grp)  ## note that when called in inner CV steps, we may have unequal number of obs per subj
  add_previous = c(0, cumsum(m_list))  ## cumulative number of obs before subject i
  K = extra_para$K
  
  ### preprocessing
  y <- y - mean(y)
  X <- scale (X)
  X.sd <- attributes(X)$`scaled:scale` ## we need to divide these sd of X from the coefficients
  Z <- scale (Z)
  Z.sd <- attributes(Z)$`scaled:scale`

  lam.a.seq<-c(exp(seq(log(100),log(0.1), length.out = 20)))*sqrt(log(p)/N/a^2/q) ## this corresponds to lambda_a scale under glmnet default objective function
  lam.j.seq <- c(exp(seq(log(100),log(0.1), length.out = 20)))*sqrt(log(p)/N/a^2/q^2) ## this is the lambda_j sequence for inference
  
  if(lm){
    ## linear model fitting, no correlation considered
    ## use cvglmnet to get tuning paramter. it is usually of order sqrt(log(p)/N)
    sig.init<-scalreg(X,y)$hsigma ## scaled lasso
    beta.hat <- as.numeric(glmnet(X, y, lambda = sig.init*sqrt(2*log(p)/N))$beta)
    
    return(list(beta.hat=beta.hat/X.sd))
    
  }else{
    ## mixed effect model fitting
    X.a<-X
    y.a<-y
    tr.siga.inv<-0
    for (i in 1:n){
      cur.mem=which(grp==grp_values[i])
      mi=length(cur.mem)
      Zi = as.matrix(Z[cur.mem,])
      
      ## sigmai = a*Zi%*%t(Zi)+diag(1,mi)
      sigmai = true_cov[[as.character(grp_values[i])]]
      sigmai.svd=svd(sigmai)
      Sig.ai.inv.half <- sigmai.svd$u %*% diag(1/sqrt(sigmai.svd$d)) %*% t(sigmai.svd$u) #compute sqrt of Sigmai
      
      X.a[cur.mem,] =  Sig.ai.inv.half %*% as.matrix(X[cur.mem,]) 
      y.a[cur.mem] =  Sig.ai.inv.half %*% as.matrix(y[cur.mem])
      tr.siga.inv<-tr.siga.inv + sum(1/sigmai.svd$d) #compute the effective sample size, which is the trace of Sigma.a.inv
    }
    
    
    ## use scaled-lasso to compute a tuning parameter
    #sig.init<-scalreg(X.a,y.a)$hsigma #scaled-lasso to get tuning parameter    
    #beta.hat= as.numeric(glmnet(X.a, y.a, lambda =sig.init*sqrt(1.2*log(p)/N))$beta)
    
    
    ## select tuning parameter lambda.a
    
    max.iter = cv.max.iter
    count = 0
    while(count < max.iter){
      print(count)
      count = count + 1
      
      
    if(time_series_cv){
      ## create cross validation set: 
      tau = round(max(m_list) * 0.2)
      test_index = do.call(c,lapply(1:n, function(ii){
        add_previous[ii] + seq(sample(1:5, 1), m_list[ii], by = round(m_list[ii]/tau))
      }
        ))
      
      cv.init = glmnet(x=X.a[-test_index,], y = y.a[-test_index], lambda = lam.a.seq, intercept = F)
      cv.mse = sapply(1:length(lam.a.seq), function(ii){
        mean((y[test_index] - X[test_index, ] %*% cv.init$beta[,ii])^2)
      }) 
      
      select_lambda = lam.a.seq[which(cv.mse == min(cv.mse))[1]]
      
      beta.hat = cv.init$beta[,which(cv.mse == min(cv.mse))[1]]
      
      
    }else{
      ## use cv.glmnet to select lambda value
      cv.init<-cv.glmnet(X.a, y.a, nfolds = 4, lam=lam.a.seq, keep=T) ## note that if we directly use MSE from this cv, it is the MSE based on 'decorrelated' data
      cv.mse = cv.init$cvm
      select_lambda = cv.init$lambda.min
      
      beta.hat<-coef(cv.init, s=cv.init$lambda.min)[-1]
      
    }
  
      
      ## if the lambda too large and all coefficients suppressed to 0, we have to make lambda seq smaller
      if(length(unique(cv.mse))==1){
        select_lambda = min(lam.a.seq)
      }
      
      
      if(which.lambda(select_lambda, lam.a.seq) == length(lam.a.seq)){
        if(count == max.iter){
          warning(paste0('lambda.a selected on boundary: smallest lambda selected' ))
        }
        ## in this case lambda need be smaller
        lam.a.seq = sort(update_tuning(sort(lam.a.seq), 'smaller'), decreasing = T)
      }else if(which.lambda(select_lambda, lam.a.seq)== 1){
        if(count == max.iter){
          warning(paste0('lambda.a selected on boundary: largest lambda selected' ))
        }
        ## in this case lambda need be larger
        lam.a.seq =  sort(update_tuning(sort(lam.a.seq), 'larger'), decreasing = T)
      }else{
        ## in this case a good lambda is selected
        break
      }
    }
  

    
    ########
    ## inference for beta (some modification: directly use full-sample estimated beta.hat in this part)
    ########
    ### compute debiased Lasso for LMM, for the specified coordinates
    beta.db.sd<-rep(NA,length=length(inf.coord))
    beta.db<-rep(NA,length=length(inf.coord))
    if(is.null(inf.coord)){
      return( list(beta.hat = beta.hat/X.sd, beta.db=beta.db/X.sd, beta.db.sd=beta.db.sd/X.sd, tr.siga.inv=tr.siga.inv))
    }
    
    
    ## compute for each coordinate of interest
    infer_one_node = function(j){
      col.j<-inf.coord[j]
      X.b<-X
      y.b<-y
      tr.sigb.inv<-0
      ## compute the correction score (wj, kappa.j)
      for (i in 1:n){
        cur.mem=which(grp==grp_values[i])
        mi=length(cur.mem)
        Zi = as.matrix(Z[cur.mem,])
        if(col.j <= q){
          Zi.tmp = Zi[,-col.j]
        }else{
          Zi.tmp = Zi
        }
        
        #sig.bi = a*Zi.tmp%*%t(Zi.tmp)+diag(1,mi)
        sig.bi = true_cov[[as.character(grp_values[i])]]
        sig.bi.svd=svd(sig.bi)
        Sig.bi.inv.half <- sig.bi.svd$u %*% diag(1/sqrt(sig.bi.svd$d)) %*% t(sig.bi.svd$u) #compute sqrt of Sig.bi
        
        X.b[cur.mem,] =  Sig.bi.inv.half %*% as.matrix(X[cur.mem,]) 
        y.b[cur.mem] =  Sig.bi.inv.half %*% as.matrix(y[cur.mem])
        tr.sigb.inv<-tr.sigb.inv + sum(1/sig.bi.svd$d) 
      }
      

      max.iter = cv.max.iter
      count = 0
      while(count < max.iter){
        print(count)
        count = count + 1
        
        if(time_series_cv){
          ## create cross validation set: 
          tau = round(max(m_list) * 0.2)
          test_index = do.call(c,lapply(1:n, function(ii){
            add_previous[ii] + seq(sample(1:5, 1), m_list[ii], by = round(m_list[ii]/tau))
          }
          ))  
          
          cv.init = glmnet(
            x = X.b[-test_index, -col.j], 
            y = X.b[-test_index, col.j], 
            lambda = lam.j.seq, 
            intercept = F)
          
          cv.mse = sapply(1:length(lam.j.seq), function(ii){
            mean((X[test_index, col.j] - X[test_index, -col.j] %*% cv.init$beta[,ii])^2)
          }) 
          
          select_lambda = lam.j.seq[which(cv.mse == min(cv.mse))[1]]
          
          kappa.hat.j = cv.init$beta[,which(cv.mse == min(cv.mse))[1]]
          
          
        }else{
          
          cv.init.j<-cv.glmnet(x = X.b[, -col.j], y = X.b[, col.j], lam=lam.j.seq, nfold=4)
          cv.mse = cv.init.j$cvm
          select_lambda = cv.init.j$lambda.min
          kappa.hat.j<-coef(cv.init.j, s=cv.init.j$lambda.min)[-1]
          
        }
      
        ## double check if the tuning parameter range is reasonable
        ## if the lambda too large and all coefficients suppressed to 0, we have to make lambda seq smaller
        if(length(unique(cv.mse))==1){
          select_lambda = min(lam.j.seq)
        }
        
        if(which.lambda(select_lambda, lam.j.seq) == length(lam.j.seq)){
          if(count == max.iter){
            warning(paste0('lambda.j selected on boundary: smallest lambda selected' ))
          }
          ## in this case lambda need be smaller
          lam.j.seq = sort(update_tuning(sort(lam.j.seq), 'smaller'), decreasing = T)
        }else if(which.lambda(select_lambda, lam.j.seq)== 1){
          if(count == max.iter){
            warning(paste0('lambda.j selected on boundary: largest lambda selected' ))
          }
          ## in this case lambda need be smaller
          lam.j.seq = sort(update_tuning(sort(lam.j.seq), 'larger'), decreasing = T)
          
        }else{
          ## in this case the suitable lambda is selected
          break
        }
      }

      

      
      wj.hat<- X.b[,col.j]-  X.b[,-col.j]%*%kappa.hat.j
      beta.db[j] = beta.hat[col.j] + sum( wj.hat * (y.b - X.b %*% beta.hat ))/sum(wj.hat * X.b[,col.j])
      
      ## and compute sandwich estimates for variance
      num=0
      for(i in 1:n){
        cur.mem=which(grp==grp_values[i])
        num <- num + sum(wj.hat[cur.mem]*(y.b[cur.mem] - X.b[cur.mem,] %*% beta.hat))^2
      }
      beta.db.sd[j] = sqrt(num)/sum(wj.hat*X.b[,col.j])
      
      return(list(inf_index = j, beta.db = beta.db[j], beta.db.sd = beta.db.sd[j]))
    }
    
    
    registerDoParallel(5)
    
    inf_res = foreach(
      j= 1:length(inf.coord),
      .verbose = T,
      .errorhandling = "pass",
      .packages = c('glmnet'),
      .export = c('lam.j.seq', 'which.lambda', 'beta.db', 'beta.db.sd', 'update_tuning')
    ) %dopar% {
      infer_one_node(j)
    }
    stopImplicitCluster()
    
    
    beta.db[do.call(c, sapply(inf_res, `[`, 'inf_index'))] = do.call(c, sapply(inf_res, `[`, 'beta.db'))
    beta.db.sd[do.call(c, sapply(inf_res, `[`, 'inf_index'))] = do.call(c, sapply(inf_res, `[`, 'beta.db.sd'))
    names(beta.db) <- names(beta.db.sd) <- inf.coord
    
    return( list(beta.hat = beta.hat/X.sd, beta.db=beta.db/X.sd[inf.coord], beta.db.sd=beta.db.sd/X.sd[inf.coord], tr.siga.inv=tr.siga.inv))
  }
}

## function to select best value a as a tuning parameter
## default use two-fold cross validation
select_a<- function(X, y, Z, grp, a.seq=seq(0, 10, 0.5), 
                    kfold = 2, # if run simple cv, we separate test and train based on subject index
                    time_series_cv = T, extra_para = list(K=NULL), # for time series case, we specify the cross validation test set index w.r.t. X matrix
                    cv.max.iter = 8,
                    true_cov){

  m_list = table(grp)  ## note that when called in inner CV steps, we may have unequal number of obs per subj
  add_previous = c(0, cumsum(m_list))  ## cumulative number of obs before subject i
  K = extra_para$K
  
  n<-length(unique(grp))
  
  if(time_series_cv){
    ## if use this cv, complicate to divide the true_cov; use the other one for now
    stop('set time_series_cv = F for now')
    
    
    # ## create cross validation set: 
    # tau = round(max(m_list) * 0.2)
    # test_index = do.call(c,lapply(1:n, function(ii){
    #   add_previous[ii] + seq(sample(1:5, 1), m_list[ii], by = round(m_list[ii]/tau))
    # }
    # ))
    # 
    # ## select a based on testing error
    # validation_error = sapply(1:length(a.seq), function(i){
    #   train_index = (1:length(grp))[-test_index]
    #   
    #   est.re<-Fixed_effects_est_inf(X = X[train_index,], 
    #                                 y = y[train_index], 
    #                                 Z = Z[train_index,], 
    #                                 grp = grp[train_index], 
    #                                 a=a.seq[i], 
    #                                 cv.max.iter = cv.max.iter,
    #                                 time_series_cv =  time_series_cv, extra_para = extra_para, ## inherit the cv settings
    #                                 true_cov = true_cov
    #                                 )
    #   
    #   pred.err<-mean((y[test_index]-X[test_index,]%*%matrix(est.re$beta.hat, ncol=1))^2)
    #   
    #   return(pred.err)
    #   
    # }
    # )
    # 
    
  }else{
    ## if simple cv, divide observations based on subject id
    index_fold_list = createFolds(y=unique(grp), k=kfold)
    
    ## select a based on kfold cross validation
    validation_error = sapply(1:length(a.seq), function(i){
      sapply(1:kfold, function(k){
        test_grp = unique(grp)[index_fold_list[[k]]]
        train_grp = unique(grp[!(grp %in% test_grp)])
        test_index = (1:length(grp))[grp %in% (test_grp)]
        train_index = (1:length(grp))[-test_index]
        
        est.re<-Fixed_effects_est_inf(
          X = X[train_index,], 
          y = y[train_index], 
          Z = Z[train_index,], 
          grp = grp[train_index], 
          a=a.seq[i], cv.max.iter = cv.max.iter, 
          time_series_cv = time_series_cv, extra_para = extra_para,
          true_cov = true_cov[as.character(train_grp)])
        
        pred.err<-mean((y[test_index]-X[test_index,]%*%matrix(est.re$beta.hat, ncol=1))^2)
        
        return(pred.err)
      })
    }
    )
    validation_error = colMeans(validation_error)
    
  }

  
  best.a = a.seq[which(validation_error == min(validation_error))[1]] ## in case there are several best a, we choose the first one
  
  return(best.a=best.a)
}


#######################
## main function for variance component estimation (assume independent random effects: diagonal covariance)
#######################


VC_estimation <- function(X, y, Z, grp, a=1, beta.hat, cv.max.iter = 8){
  
  N=length(y)
  n=length(unique(grp))
  q=ncol(Z)
  p=ncol(X)
  m = max(table(grp))
  grp_values = unique(grp)
  
  ### preprocessing
  y <- y - mean(y)
  
  X <- scale (X)
  X.sd <- attributes(X)$`scaled:scale` ## we need to divide these sd of X from the coefficients
  Z <- scale (Z)
  Z.sd <- attributes(Z)$`scaled:scale`
  
  lam.sig.seq<-c(50,10,exp(seq(log(5),log(0.1), length.out = 20)))*q*sqrt(log(q)/n)/ m 
  
  
  ## compute VC
  res = lapply(1:n, function(i){
    
    cur.mem=which(grp==grp_values[i])
    mi=length(cur.mem)
    Zi = as.matrix(Z[cur.mem,])
    ri = y[cur.mem] - as.matrix(X[cur.mem,]) %*% beta.hat
    index.mat.tmp = lower.tri(diag(1, mi), diag=F) ## extract only off-diagonal entries
    
    A_l = lapply(1:q, function(l) Zi[,l, drop=F] %*% t(Zi[,l, drop=F])) ## these are coefficient matrices for the VC parameters
    
    y_vc = (ri%*%t(ri))[index.mat.tmp]
    X_vc = do.call(cbind, lapply(A_l, function(A)A[index.mat.tmp]))
    return(list(y_vc, X_vc))
  })
  
  y_vc = do.call(c, sapply(res, `[`, 1))  
  X_vc = do.call(rbind, sapply(res, `[`, 2))  
  
  
  max.iter = cv.max.iter
  count = 0
  while(count < max.iter){
    count = count + 1
    vc.cv.init<-cv.glmnet(X_vc, y_vc, nfolds = 4, lam=lam.sig.seq)
    
    ## if the lambda too large and all coefficients suppressed to 0, we have to make lambda seq smaller
    if(length(unique(vc.cv.init$cvm))==1){
      vc.cv.init$lambda.min = min(lam.sig.seq)
    }
    
    if(which.lambda(vc.cv.init$lambda.min, lam.sig.seq) == length(lam.sig.seq)){
      if(count == max.iter){
        warning(paste0('lambda VC selected on boundary: smallest lambda selected' ))
      }
      ## in this case lambda need be smaller
      lam.sig.seq = sort(update_tuning(sort(lam.sig.seq), 'smaller'), decreasing = T)
    }else if(which.lambda(vc.cv.init$lambda.min, lam.sig.seq)== 1){
      if(count == max.iter){
        warning(paste0('lambda VC selected on boundary: largest lambda selected' ))
      }
      ## in this case lambda need be larger
      lam.sig.seq = sort(update_tuning(sort(lam.sig.seq), 'larger'), decreasing = T)
    }else{
      ## in this case suitable lambda is selected
      break
    }
  }

  
  
  vc.hat<-coef(vc.cv.init, s=vc.cv.init$lambda.min)[-1]
  
  sig2.eps.hat = 1/N*
    sum(sapply(1:n, function(i){
      
      cur.mem=which(grp==grp_values[i])
      mi=length(cur.mem)
      Zi = as.matrix(Z[cur.mem,])
      ri = y[cur.mem] - as.matrix(X[cur.mem,]) %*% beta.hat
      
      return(sum(ri^2) - sum(diag(Zi %*% diag(vc.hat) %*% t(Zi))) )
    }))
  
  sig2.eps.hat = ifelse(sig2.eps.hat<0, 0, sig2.eps.hat)
  
  ## make estimates non-negative: truncate at 0
  vc.hat[vc.hat<0]<-0
  
  
  return( list(vc.hat = vc.hat/Z.sd^2 , sig2.eps.hat = sig2.eps.hat))
  
}


##################
## the function to estimate and inference the whole coefficient matrix Phi (can be slow when p*T*n is too large)
##################
#Inputs-- 
## y: list of observations for n subjects; observations for each subject, the data is a Time x p matrix from time 1 to time T; 
## a : tuning parameter in $\Sig_a$
## inf.coord: the coordinates of beta for inference, in the vectorized form (by column)

#Outputs-- 
## beta.hat: the Lasso estimator based on psedo-likelihood; 
## beta.db: debiased Lasso estimators for the fixed effects in inf.coord;
## beta.db.sd: standard deviation of the debiased Lasso estimators for the fixed effects in inf.coord.




#########################
## other utility functions
#########################

## function to quickly generate a covariance matrix
generate_cov_matrix = function(type, param){
  ## type = 'diag': param is the diagonal vector
  ## type = 'sparse-random': param = c(matrix dimension, off-diag unif range, off-diag nonzero prob)
  ## type = 'AR1': param = c(matrix dimension, rho)
  ## type = 'const': param = c(matrix dimension, constant for off-diagonal)
  if(type == 'diag'){ 
    Sig = diag(param)
  }else if(type == 'sparse-random'){
    Sig = diag(0, param[1])
    counter=0
    while(any(eigen(Sig)$values<=0) & counter < 10){
      counter = counter+1
      Sig = matrix(runif(param[1]^2),param[1],param[1])
      Sig = 1*((Sig+t(Sig))/2 < param[3])
      Sig = Sig * matrix(runif(param[1]^2, min = -param[2], max=param[2]), param[1], param[1])
      Sig = (Sig + t(Sig))/2
      diag(Sig)<-1
    }
    
    if(any(eigen(Sig)$values<=0)) stop('Cov matrix not pd')
  }else if(type == 'AR1'){
    Sig = toeplitz(param[2]^(0:(param[1]-1)))
  }else if(type == 'const'){
    Sig = matrix(param[2], param[1], param[1])
    diag(Sig)<-1
  }else{
    stop('specify Cov type among diag, sparse-random, AR1, const')
  }
  
  return(Sig)
}


## function to select best value a as a tuning parameter
## default use two-fold cross validation
select_a_full_matrix<- function(y, a.seq=seq(0, 10, 0.5), kfold = 2){
  
  n<-length(y)
  
  index_fold_list = createFolds(y= 1:n, k=kfold)
  
  
  ## select a based on kfold cross validation
  validation_error = sapply(1:length(a.seq), function(i){
    sapply(1:kfold, function(k){
      test_index = (1:n)[index_fold_list[[k]]]
      train_index = (1:n)[-test_index]
      
      est.re<-Fixed_effects_est_inf_full_matrix(y = y[train_index], a=a.seq[i])
      
      ## compute prediction error for all time points 
      pred.err = 
        sum(
          sapply(y[test_index], function(Y){
            sum(
              sapply(1:(nrow(Y)-1), function(ii) sum((Y[ii+1,] - est.re$beta.hat %*% t(Y[ii,, drop=F]))^2))
            )
          })
        )
      
      return(pred.err)
    })
  }
  )
  
  validation_error = colMeans(validation_error)
  best.a = a.seq[which(validation_error == min(validation_error))[1]] ## in case there are several best a, we choose the first one
  
  return(best.a=best.a)
}

## function to determine which tuning parameter lambda is selected (due to rounding error in cv.glmnet, the selected lambda is different from supplied lambda)
which.lambda = function(lambda, lambda_list){
  tmp = abs(lambda_list - lambda)
  which(tmp == min(tmp))
}


## function to adjust the tuning parameter sequence
update_tuning = function(lambda, adjust_direction){
  if(length(lambda)<2){
    stop('need lambda length at least 2')
  }
  if (adjust_direction == 'larger'){
    new_lambda = seq(tail(lambda,2)[1], tail(lambda,1)*10, length.out = length(lambda))
  }  else if (adjust_direction == 'smaller'){
    new_lambda = seq(lambda[1]/10, lambda[2], length.out = length(lambda))
  }else{
    stop('direction needs to be larger or smaller')
  }
  new_lambda
}


## function to generate a perturbed matrix based on a given matrix
## approach to perturb: randomly select a few edges, and add small normally distributed noise
perturb_Sig = function(Sigma, 
                       para, # takes c(prob, sd), prob for probability of having a perturbed edge, and sd for perturbation sd
                       counter_max = 50){
  p = nrow(Sigma)
  continue = T
  counter=0
  while(continue & counter<counter_max){
    counter = counter+1
    vary_edge = matrix(runif(p*p), p, p)
  vary_edge = ((vary_edge + t(vary_edge))< 2*para[1])*1
  diag(vary_edge) <- 0
  
  vary_edge_amount = matrix(rnorm(p*p, sd=para[2]), p, p)
  vary_edge_amount = (vary_edge_amount + t(vary_edge_amount))/2
  
  Sigma.i = Sigma + vary_edge * vary_edge_amount
  
  if(min(eigen(Sigma.i)$values)>0){
    continue = F
  }
  
  if(counter == counter_max & continue){
    stop('Change perturb parameter, Sigma.i not psd')
  }
  }
  
  return(Sigma.i)

}

