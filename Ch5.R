library(corpcor)
library(MASS)
library(readr)
library(rlang)
library(caret)
library(tidyverse)
library(parallel)
library(car)
library(dplyr)
library(tidyr)

# Import true regression coefficients into a matrix called betas

xy_function <- function(ind,x,beta,sigma,nsmall){
  
  y <- x%*%beta[,ind] + rnorm(nsmall,0,sigma)
  
  return(y)
}

simulate_part1 <- function(num_nsmall,nsmall,p,beta,gamma,sigma,num_beta){
  
  mean <- rep(0,p+1)
  cov <- diag(p+1)
  z <- mvrnorm(nsmall,mean,cov)
  x <- matrix(nrow = nrow(z),ncol = (ncol(z)-1))
  
  for(i in 1:(ncol(z)-1)){
    x[,i] <- sqrt(1-gamma*gamma)*(z[,i]) + gamma*z[,ncol(z)]
  }
  
  rm(list = c("z"))  
  
  xy_function <- mclapply(1:num_beta,xy_function,x,beta,sigma,nsmall,mc.cores = 1)
  
  xy <- list(x,xy_function)
  rm(list = c("x"))  
  
  xprimex <- t(xy[[1]])%*%xy[[1]]
  
  xprimey <- matrix(data = NA,ncol = num_beta,nrow = p)
  yprimey <- matrix(data = NA,ncol = num_beta,nrow = 1)

  for(j in 1:num_beta){
    xprimey[,j] <- t(xy[[1]])%*%xy[[2]][[j]]
    yprimey[,j] <- t(xy[[2]][[j]])%*%xy[[2]][[j]]
  }
  sufficient <- list(xprimex,xprimey,yprimey)
  
  return(sufficient)
}

simulate_part2 <- function(sufficient_list){
  
  xprimex_final <- 0
  xprimey_final <- 0
  yprimey_final <- 0

  for(i in 1:length(sufficient_list)){
    xprimex_final <- xprimex_final + sufficient_list[[i]][[1]]
    xprimey_final <- xprimey_final + sufficient_list[[i]][[2]]
    yprimey_final <- yprimey_final + sufficient_list[[i]][[3]]

  }
  
  xyprimexy <- list(xprimex_final,xprimey_final,yprimey_final)
  
  return(xyprimexy)
}


simulate_part3 <- function(it,num_nsmall,nsmall,p,beta,gamma,sigma,num_beta){
  
  sufficient1 <- mclapply(1:num_nsmall,simulate_part1,nsmall,p2,betas,gamma,sigma,num_beta,mc.cores = 1)
  sufficient2 <- simulate_part2(sufficient1)
  
  eigenv <- eigen(sufficient2[[1]])
  beta_hat <- solve(sufficient2[[1]])%*%sufficient2[[2]]
  sigma_hat <- (sufficient2[[3]] - diag(t(sufficient2[[2]])%*%solve(sufficient2[[1]])%*%sufficient2[[2]]))/(nsmall*num_nsmall - nrow(beta_hat))
  beta_hat_cov <- matrix(rep(diag(solve(sufficient2[[1]])),num_beta),ncol=num_beta)*matrix(rep(sigma_hat,p),ncol = num_beta)

  result <- list(beta_hat,sigma_hat,beta_hat_cov,eigenv)
  
}

it <- 1000
n <- 10000000
p <- 100

nsmall <- 1000
num_nsmall <- n/nsmall
num_beta <- ncol(betas)
gamma <- 0.999
sigma <- 1

###############################################################################

n5p1gamma4 <- mclapply(1:it,simulate_part3,num_nsmall,nsmall,
                       p,betas,gamma,sigma,num_beta,mc.cores=50)



