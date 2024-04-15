library(rlang)
library(readr)
library(dplyr)
library(tidyr)
library(parallel)
library(rlist)

# Import true regression coefficients into a matrix called beta
# Obtain file path to split data: files_x and files_y

obtain_x <- function(files_x,ind){
  
  data_x <- read_csv(files_x[ind])
  data_x <- data_x %>% drop_na()
  
  return(data_x)
}

obtain_y <- function(files_y,ind){
  
  data_y <- read_csv(files_y[ind])
  data_y <- data_y %>% drop_na()
  
  return(data_y)
  
}

obtain_x_prime_x <- function(ind,files_x){
  options(digits = 20)
  sxx <- 0
  data_use <- obtain_x(files_x,ind)
  
  x <- as.matrix(data_use[,1:ncol(data_use)],nrow(data_use),ncol(data_use))
  
  for(i in 1:nrow(x)){
    sxx <- sxx + t(matrix(x[i,],nrow=1))%*%(matrix(x[i,],nrow=1))
  }
  return(sxx)
}

obtain_k_x_prime_x <- function(k_seq,obtain_x_prime_x,num_files,k){
  options(digits = 20)
  num_blocks <- num_files/k 
  seq <- 1:num_files
  seq_start <- seq(1,num_files,num_blocks)
  seq_end <- seq(num_blocks,num_files,num_blocks)
  sxx <- 0
  
  for(i in seq[-(seq_start[k_seq]:seq_end[k_seq])]){
    sxx <- sxx + obtain_x_prime_x[[i]]
  }
  
  eigenvalues <- eigen(sxx)[[1]]
  eigenvectors <- eigen(sxx)[[2]]
  eigenvalues_diag <- diag(eigenvalues)
  sxx_for_eigen <- list(sxx,eigenvalues_diag,eigenvectors)
  return(sxx_for_eigen)
}

est_params_xx_xy_yy <- function(k_seq,x,y,sxx_eigen){
  options(digits = 20)
  x <- as.matrix(x,nrow = nrow(x), ncol = ncol(x))
  y <- as.matrix(y,nrow = nrow(y),ncol = 1)
  
  sxx <- 0  
  sxy <- 0
  syy <- 0
  
  
  for(i in 1:nrow(x)){
    sxx <- sxx +  (t(sxx_eigen[[k_seq]][[3]])%*%t(matrix(x[i,],nrow=1)))%*%(matrix(x[i,],nrow=1)%*%sxx_eigen[[k_seq]][[3]])  
    sxy <- sxy +  (t(sxx_eigen[[k_seq]][[3]])%*%t(matrix(x[i,],nrow=1)))*y[i,]  
    syy <- syy + y[i,]*y[i,]  
    
  }
  sxx_xy_yy <- list(sxx,sxy,syy)
  return(sxx_xy_yy)  
  
}  

sufficiency_matrix1 <- function(ind,x,y,sxx_eigen,num_files,k,k_seq){
  options(digits = 20)
  x <- obtain_x(files_x,ind)
  y <- obtain_y(files_y,ind)
  
  sxx_xy_yy <- mclapply(k_seq,est_params_xx_xy_yy,x,y,sxx_eigen,mc.cores = 1)
  
  return(sxx_xy_yy)
}

sufficiency_matrix_final <- function(k_seq,num_files,k,sufficiency_matrix1,num_y){
  options(digits = 20)
  num_blocks <- num_files/k 
  seq <- 1:num_files
  seq_start <- seq(1,num_files,num_blocks)
  seq_end <- seq(num_blocks,num_files,num_blocks)
  
  sxx <- 0
  sxy <- 0
  syy <- 0
  
  for(i in seq[-(seq_start[k_seq]:seq_end[k_seq])]){
    sxx <- sxx + sufficiency_matrix1[[i]][[k_seq]][[1]]
    sxy <- sxy + sufficiency_matrix1[[i]][[k_seq]][[2]]
    syy <- syy + sufficiency_matrix1[[i]][[k_seq]][[3]]
    
  }
  
  suff_final <- list(sxx,sxy,syy)
  return(suff_final)
}

estimation_ols <- function(k_seq,num_files,num_obs,k,sufficiency_matrix_final){
  options(digits = 20)  
  n <- num_files*num_obs*(k-1)/k
  alpha <- solve(as.matrix(sufficiency_matrix_final[[k_seq]][[1]]))%*%as.matrix(sufficiency_matrix_final[[k_seq]][[2]])
  sigma <- as.numeric(((1/(n-length(alpha)))*(sufficiency_matrix_final[[k_seq]][[3]] - t(as.matrix(sufficiency_matrix_final[[k_seq]][[2]]))%*%solve(as.matrix(sufficiency_matrix_final[[k_seq]][[1]]))%*%as.matrix(sufficiency_matrix_final[[k_seq]][[2]]))))
  cov_alpha <- sigma*solve(as.matrix(sufficiency_matrix_final[[k_seq]][[1]]))
  ols <- list(alpha,sigma,cov_alpha)
  return(ols)
}


estimation_ridge <- function(k_seq,num_files,num_obs,k,sufficiency_matrix_final,estimation_ols){
  options(digits = 20)
  n <- num_files*num_obs*(k-1)/k
  optimal_kl <- estimation_ols[[k_seq]][[2]]/(max(estimation_ols[[k_seq]][[1]])^2)
  optimal_k <- rep(as.numeric(optimal_kl),ncol(as.matrix(sufficiency_matrix_final[[k_seq]][[1]])))
  optimal_k[1] <- 0
  gzz <- solve(as.matrix(sufficiency_matrix_final[[k_seq]][[1]]) + diag(as.vector(optimal_k)))
  zz <- as.matrix(sufficiency_matrix_final[[k_seq]][[1]])
  zy <- as.matrix(sufficiency_matrix_final[[k_seq]][[2]])
  alpha <- gzz%*%zy
  cov_alpha <- estimation_ols[[k_seq]][[2]]*(gzz%*%zz%*%gzz)
  sigma <- as.numeric(((1/(n-length(alpha)))*(sufficiency_matrix_final[[k_seq]][[3]] - t(as.matrix(sufficiency_matrix_final[[k_seq]][[2]]))%*%solve(as.matrix(sufficiency_matrix_final[[k_seq]][[1]])+ diag(as.vector(optimal_k)))%*%as.matrix(sufficiency_matrix_final[[k_seq]][[2]]))))
  ridge <- list(alpha,cov_alpha,optimal_kl,sigma)
  return(ridge)
  
}


estimation_modified_ridge <- function(k_seq,num_files,num_obs,k,sufficiency_matrix_final,estimation_ols){
  options(digits = 20)
  n <- num_files*num_obs*(k-1)/k
  initial_d <- min(estimation_ols[[k_seq]][[2]]/(estimation_ols[[k_seq]][[1]]^2))
  optimal_kl <- (length(estimation_ols[[k_seq]][[1]])*estimation_ols[[k_seq]][[2]])/(sum((1+initial_d)*(estimation_ols[[k_seq]][[1]]^2)))
  optimal_d1 <- (estimation_ols[[k_seq]][[2]]/(optimal_kl*(estimation_ols[[k_seq]][[1]]^2))) - 1
  optimal_d2 <- length(estimation_ols[[k_seq]][[1]])/(sum((optimal_d1^(-1))))
  optimal_dl = optimal_d2
  optimal_kd <- rep(as.numeric(optimal_kl)*(1 + as.numeric(optimal_dl)),ncol(as.matrix(sufficiency_matrix_final[[k_seq]][[1]])))
  optimal_kd[1] <- 0
  zz <- as.matrix(sufficiency_matrix_final[[k_seq]][[1]])
  gzz <- solve(zz + diag(as.vector(optimal_kd)))
  zy <- as.matrix(sufficiency_matrix_final[[k_seq]][[2]])
  alpha <- gzz%*%zy
  cov_alpha <- estimation_ols[[k_seq]][[2]]*(gzz%*%zz%*%gzz)
  modified_ridge <- list(alpha,cov_alpha,optimal_kl,optimal_dl)
  return(modified_ridge)
}


estimation_liu <- function(k_seq,num_files,num_obs,k,sufficiency_matrix_final,estimation_ols,sxx_eigen){
  options(digits = 20) 
  n <- num_files*num_obs*(k-1)/k
  eigenvalues <- diag(sxx_eigen[[k_seq]][[2]])
  optimal_dl_num <- sum(1/((eigenvalues)*(eigenvalues + 1)))
  optimal_dl_denom <- sum(estimation_ols[[k_seq]][[1]]/((eigenvalues + 1)^2))    
  optimal_dl <- 1 - estimation_ols[[k_seq]][[2]]*(optimal_dl_num/optimal_dl_denom)
  optimal_d <- rep(as.numeric(optimal_dl),ncol(as.matrix(sufficiency_matrix_final[[k_seq]][[1]])))
  optimal_d[1] <- 0
  zz <- as.matrix(sufficiency_matrix_final[[k_seq]][[1]])
  gzz <- solve(zz + diag(ncol(as.matrix(sufficiency_matrix_final[[k_seq]][[1]]))))%*%(zz + diag(as.vector(optimal_d)))%*%solve(zz)
  zy <- as.matrix(sufficiency_matrix_final[[k_seq]][[2]])
  alpha <- gzz%*%zy
  cov_alpha <- estimation_ols[[k_seq]][[2]]*(gzz%*%zz%*%gzz)
  liu <- list(alpha,cov_alpha,optimal_dl)
  return(liu)
  
}


estimation_modified_liu <- function(k_seq,num_files,num_obs,k,sufficiency_matrix_final,estimation_ols,sxx_eigen){
  options(digits = 20)
  n <- num_files*num_obs*(k-1)/k
  eigenvalues <- diag(sxx_eigen[[k_seq]][[2]])
  optimal_dl_num <- ((eigenvalues)*(estimation_ols[[k_seq]][[2]] + estimation_ols[[k_seq]][[1]]^2))
  optimal_dl_denom <- (estimation_ols[[k_seq]][[2]] + eigenvalues*(estimation_ols[[k_seq]][[1]]^2))
  optimal_dl <- min(optimal_dl_num/optimal_dl_denom)
  optimal_d <- rep(as.numeric(optimal_dl),ncol(as.matrix(sufficiency_matrix_final[[k_seq]][[1]])))
  optimal_d[1] <- 0
  zz <- as.matrix(sufficiency_matrix_final[[k_seq]][[1]])
  gzz <- solve(zz + diag(ncol(zz)))%*%(zz - diag(as.vector(optimal_d)))%*%solve(zz)
  zy <- as.matrix(sufficiency_matrix_final[[k_seq]][[2]])
  alpha <- gzz%*%zy
  cov_alpha <- estimation_ols[[k_seq]][[2]]*(gzz%*%zz%*%gzz)
  modified_liu <- list(alpha,cov_alpha,optimal_dl)
  return(modified_liu)
  
}


estimation_kl <- function(k_seq,num_files,num_obs,k,sufficiency_matrix_final,estimation_ols,sxx_eigen){
  options(digits = 20)
  n <- num_files*num_obs*(k-1)/k
  eigenvalues <- diag(sxx_eigen[[k_seq]][[2]])
  optimal_kl <- (estimation_ols[[k_seq]][[2]])/sum((2*(estimation_ols[[k_seq]][[1]]^2)+(estimation_ols[[k_seq]][[2]]/eigenvalues)))
  optimal_k <- rep(as.numeric(optimal_kl),ncol(as.matrix(sufficiency_matrix_final[[k_seq]][[1]])))
  optimal_k[1] <- 0
  zz <- as.matrix(sufficiency_matrix_final[[k_seq]][[1]])
  gzz <- solve(zz + diag(as.vector(optimal_k)))%*%(zz - diag(as.vector(optimal_k)))%*%solve(zz)
  zy <- as.matrix(sufficiency_matrix_final[[k_seq]][[2]])
  alpha <- gzz%*%zy
  cov_alpha <- estimation_ols[[k_seq]][[2]]*(gzz%*%zz%*%gzz)
  kl <- list(alpha,cov_alpha,optimal_kl)
  return(kl)
}

mse_estimators <- function(k_seq,sxx_eigen,estimation_ols,estimation_ridge,estimation_modified_ridge,estimation_liu,estimation_modified_liu,estimation_kl,beta){
  options(digits = 20)
  ols <- sum((solve(t(method2b[[k_seq]][[3]]))%*%estimation_ols[[k_seq]][[1]] - beta)^2)
  ridge <- sum((solve(t(method2b[[k_seq]][[3]]))%*%estimation_ridge[[k_seq]][[1]] - beta)^2)
  modified_ridge <- sum((solve(t(method2b[[k_seq]][[3]]))%*%estimation_modified_ridge[[k_seq]][[1]] - beta)^2)
  liu <- sum((solve(t(method2b[[k_seq]][[3]]))%*%estimation_liu[[k_seq]][[1]] - beta)^2)
  modified_liu <- sum((solve(t(method2b[[k_seq]][[3]]))%*%estimation_modified_liu[[k_seq]][[1]] - beta)^2)
  kl <- sum((solve(t(method2b[[k_seq]][[3]]))%*%estimation_kl[[k_seq]][[1]] - beta)^2)
  estimators <- list(ols,ridge,modified_ridge,liu,modified_liu,kl)
  return(estimators)
}


mse_final <- function(k,mse_estimators){
  options(digits = 20)
  ols <- 0
  ridge <- 0
  modified_ridge <- 0
  liu <-  0
  modified_liu <- 0
  kl <- 0
  
  for(i in 1:k){
    ols <- ols + mse_estimators[[i]][[1]]
    ridge <- ridge + mse_estimators[[i]][[2]]
    modified_ridge <- modified_ridge + mse_estimators[[i]][[3]]
    liu <- liu + mse_estimators[[i]][[4]]
    modified_liu <- modified_liu + mse_estimators[[i]][[5]]
    kl <- kl + mse_estimators[[i]][[6]]
  }
  
  mse <- list(ols,ridge,modified_ridge,liu,modified_liu,kl)
  return(mse)
}

predict_y <- function(k_seq,num_files,k,sxx_eigen,estimation_ols,estimation_ridge,estimation_modified_ridge,
                      estimation_liu,estimation_modified_liu,estimation_kl){
  options(digits = 20)
  num_blocks <- num_files/k 
  seq <- 1:num_files
  seq_start <- seq(1,num_files,num_blocks)
  seq_end <- seq(num_blocks,num_files,num_blocks)
  
  pred_ols_y <- list()
  mse_ols_part <- list()
  SMAPE_ols_part <- list()
  
  pred_ridge_y <- list()
  mse_ridge_part <- list()
  SMAPE_ridge_part <- list()
  
  pred_modified_ridge_y <- list()
  mse_modified_ridge_part <- list()
  SMAPE_modified_ridge_part <- list()
  
  pred_liu_y <- list()
  mse_liu_part <- list()
  SMAPE_liu_part <- list()
  
  pred_modified_liu_y <- list()
  mse_modified_liu_part <- list()
  SMAPE_modified_liu_part <- list()
  
  pred_kl_y <- list()
  mse_kl_part <- list()
  SMAPE_kl_part <- list()
  
  j <- 1
  
  for(i in seq[seq_start[k_seq]:seq_end[k_seq]]){
    x <- obtain_x(files_x,i)
    y <- obtain_y(files_y,i)
    
    
    pred_ols_y[[j]] <- as.matrix(x,ncol = ncol(x))%*%(solve(t(sxx_eigen[[k_seq]][[3]]))%*%estimation_ols[[k_seq]][[1]])
    mse_ols_part[[j]] <- (pred_ols_y[[j]] - as.matrix(y,ncol = 1))^2
    SMAPE_ols_part[[j]] <- abs(pred_ols_y[[j]] - as.matrix(y,ncol = 1))/(0.5*(abs(as.matrix(y,ncol = 1)) + abs(pred_ols_y[[j]])))
    
    pred_ridge_y[[j]] <- as.matrix(x,ncol = ncol(x))%*%(solve(t(sxx_eigen[[k_seq]][[3]]))%*%estimation_ridge[[k_seq]][[1]])
    mse_ridge_part[[j]] <- (pred_ridge_y[[j]] - as.matrix(y,ncol = 1))^2
    SMAPE_ridge_part[[j]] <- abs(pred_ridge_y[[j]] - as.matrix(y,ncol = 1))/(0.5*(abs(as.matrix(y,ncol = 1)) + abs(pred_ridge_y[[j]])))
    
    pred_modified_ridge_y[[j]] <- as.matrix(x,ncol = ncol(x))%*%(solve(t(sxx_eigen[[k_seq]][[3]]))%*%estimation_modified_ridge[[k_seq]][[1]])
    mse_modified_ridge_part[[j]] <- (pred_modified_ridge_y[[j]] - as.matrix(y,ncol = 1))^2
    SMAPE_modified_ridge_part[[j]] <- abs(pred_modified_ridge_y[[j]] - as.matrix(y,ncol = 1))/(0.5*(abs(as.matrix(y,ncol = 1)) + abs(pred_modified_ridge_y[[j]])))
    
    pred_liu_y[[j]] <- as.matrix(x,ncol = ncol(x))%*%(solve(t(sxx_eigen[[k_seq]][[3]]))%*%estimation_liu[[k_seq]][[1]])
    mse_liu_part[[j]] <- (pred_liu_y[[j]] - as.matrix(y,ncol = 1))^2
    SMAPE_liu_part[[j]] <- abs(pred_liu_y[[j]] - as.matrix(y,ncol = 1))/(0.5*(abs(as.matrix(y,ncol = 1)) + abs(pred_liu_y[[j]])))
    
    pred_modified_liu_y[[j]] <- as.matrix(x,ncol = ncol(x))%*%(solve(t(sxx_eigen[[k_seq]][[3]]))%*%estimation_modified_liu[[k_seq]][[1]])
    mse_modified_liu_part[[j]] <- (pred_modified_liu_y[[j]] - as.matrix(y,ncol = 1))^2
    SMAPE_modified_liu_part[[j]] <- abs(pred_modified_liu_y[[j]] - as.matrix(y,ncol = 1))/(0.5*(abs(as.matrix(y,ncol = 1)) + abs(pred_modified_liu_y[[j]])))
    
    pred_kl_y[[j]] <- as.matrix(x,ncol = ncol(x))%*%(solve(t(sxx_eigen[[k_seq]][[3]]))%*%estimation_kl[[k_seq]][[1]])
    mse_kl_part[[j]] <- (pred_kl_y[[j]] - as.matrix(y,ncol = 1))^2
    SMAPE_kl_part[[j]] <- abs(pred_kl_y[[j]] - as.matrix(y,ncol = 1))/(0.5*(abs(as.matrix(y,ncol = 1)) + abs(pred_kl_y[[j]])))
    
    j <- j + 1
  }
  
  pred_diff_y <- list(pred_ols_y,mse_ols_part,SMAPE_ols_part,pred_ridge_y,mse_ridge_part,SMAPE_ridge_part,
                      pred_modified_ridge_y,mse_modified_ridge_part,SMAPE_modified_ridge_part,pred_liu_y,mse_liu_part,SMAPE_liu_part,
                      pred_modified_liu_y,mse_modified_liu_part,SMAPE_modified_liu_part,pred_kl_y,mse_kl_part,SMAPE_kl_part)
  
  return(pred_diff_y)
  
}

mse_y <- function(k_seq,num_files,k,predict_y){
  options(digits = 20)
  num_blocks <- num_files/k   
  
  ols <- 0
  ridge <- 0
  modified_ridge <- 0
  liu <- 0
  modified_liu <- 0
  kl <- 0
  
  for(j in 1:num_blocks){
    ols <- ols + sum(predict_y[[k_seq]][[2]][[j]])
    ridge <- ridge + sum(predict_y[[k_seq]][[5]][[j]])
    modified_ridge <- modified_ridge  + sum(predict_y[[k_seq]][[8]][[j]])
    liu <- liu + sum(predict_y[[k_seq]][[11]][[j]])
    modified_liu <- modified_liu + sum(predict_y[[k_seq]][[14]][[j]])
    kl <- kl + sum(predict_y[[k_seq]][[17]][[j]])
  }
  
  estimators <- list(ols,ridge,modified_ridge,liu,modified_liu,kl)
  return(estimators)
  
}

mse_y_final <- function(k,mse_y){
  options(digits = 20)  
  ols <- 0
  ridge <- 0
  modified_ridge <- 0
  liu <-  0
  modified_liu <- 0
  kl <- 0
  
  for(j in 1:k){
    ols <- ols + as.numeric(mse_y[[j]][[1]][[1]])
    ridge <- ridge + mse_y[[j]][[2]][[1]]
    modified_ridge <- modified_ridge + mse_y[[j]][[3]][[1]]
    liu <- liu + mse_y[[j]][[4]][[1]]
    modified_liu <- modified_liu + mse_y[[j]][[5]][[1]]
    kl <- kl + mse_y[[j]][[6]][[1]]
  }
  
  mse <- list(ols,ridge,modified_ridge,liu,modified_liu,kl)
  return(mse)
}

smape_y <- function(k_seq,num_files,k,predict_y){
  options(digits = 20)  
  num_blocks <- num_files/k   
  
  ols <- 0
  ridge <- 0
  modified_ridge <- 0
  liu <- 0
  modified_liu <- 0
  kl <- 0
  
  for(j in 1:num_blocks){
    ols <- ols + sum(predict_y[[k_seq]][[3]][[j]])
    ridge <- ridge + sum(predict_y[[k_seq]][[6]][[j]])
    modified_ridge <- modified_ridge  + sum(predict_y[[k_seq]][[9]][[j]])
    liu <- liu + sum(predict_y[[k_seq]][[12]][[j]])
    modified_liu <- modified_liu + sum(predict_y[[k_seq]][[15]][[j]])
    kl <- kl + sum(predict_y[[k_seq]][[18]][[j]])
  }
  
  estimators <- list(ols,ridge,modified_ridge,liu,modified_liu,kl)
  return(estimators)
  
}

smape_y_final <- function(k,smape_y){
  options(digits = 20)  
  ols <- 0
  ridge <- 0
  modified_ridge <- 0
  liu <-  0
  modified_liu <- 0
  kl <- 0
  
  for(j in 1:k){
    ols <- ols + smape_y[[j]][[1]]
    ridge <- ridge + smape_y[[j]][[2]]
    modified_ridge <- modified_ridge + smape_y[[j]][[3]]
    liu <- liu + smape_y[[j]][[4]]
    modified_liu <- modified_liu + smape_y[[j]][[5]]
    kl <- kl + smape_y[[j]][[6]]
  }
  
  smape <- list(ols,ridge,modified_ridge,liu,modified_liu,kl)
  return(smape)
}

num_files <- 100
k <- 5
k_seq <- 1:k
ind <- 1:num_files
num_obs <- 100000

time_estimation <- system.time({method2 <- mclapply(1:num_files,obtain_x_prime_x,files_x,mc.cores = 10)
method2b <- mclapply(k_seq,obtain_k_x_prime_x,method2,num_files,k,mc.cores = k)
sufficiency_matrix1c <- mclapply(ind,sufficiency_matrix1,x,y,method2b,num_files,k,k_seq,mc.cores = 10)
sufficiency_matrix_finalc <- mclapply(k_seq,sufficiency_matrix_final,num_files,k,sufficiency_matrix1c,num_y,mc.cores = k)
ols <- mclapply(k_seq,estimation_ols,num_files,num_obs,k,sufficiency_matrix_finalc,mc.cores = k)
ridge <- mclapply(k_seq,estimation_ridge,num_files,num_obs,k,sufficiency_matrix_finalc,ols,mc.cores = k)
modified_ridge <- mclapply(k_seq,estimation_modified_ridge,num_files,num_obs,k,sufficiency_matrix_finalc,ols,mc.cores = k)
liu <- mclapply(k_seq,estimation_liu,num_files,num_obs,k,sufficiency_matrix_finalc,ols,method2b,mc.cores = k)
modified_liu <- mclapply(k_seq,estimation_modified_liu,num_files,num_obs,k,sufficiency_matrix_finalc,ols,method2b,mc.cores = k)
kl <- mclapply(k_seq,estimation_kl,num_files,num_obs,k,sufficiency_matrix_finalc,ols,method2b,mc.cores = k)})

time_performance <- system.time({mse_estimators1 <- mclapply(k_seq,mse_estimators,method2b,ols,ridge,modified_ridge,liu,modified_liu,kl,beta,mc.cores = k)
mse_final1 <- mse_final(k,mse_estimators1)
predict_y1 <- mclapply(k_seq,predict_y,num_files,k,method2b,ols,ridge,modified_ridge,liu,
                       modified_liu,kl,mc.cores = k)
mse_y1 <- mclapply(k_seq,mse_y,num_files,k,predict_y1,mc.cores = k)
mse_y_final1 <- mse_y_final(k,mse_y1)
smape_y1 <- mclapply(k_seq,smape_y,num_files,k,predict_y1,mc.cores = k)
smape_y_final1 <- smape_y_final(k,smape_y1)
})

