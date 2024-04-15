library(readr)
library(glmnet)
library(cvTools)
library(schoolmath)
library(caret)
library(pROC)
library(parallel)
library(e1071)

# Import data into two matrices x and y

mysd <- function(y) sqrt(sum((y-mean(y))^2)/length(y))

lambda_path <- function(x,y,alpha){
  sx <- scale(x, scale=apply(x, 2, mysd))
  sx <- as.matrix(sx, ncol=ncol(x), nrow=nrow(x)) 
  max_lambda <- max(abs(colSums(sx*rep(ifelse(y==0, -prop.table(table(y))[1], 
                                              prop.table(table(y))[2]),ncol(x)))))/nrow(x)
  max_lambda <- max_lambda*1000
  epsilon <- 0.05
  K <- 100
  lambdapath <- round(exp(seq(log(max_lambda), log(max_lambda*epsilon), 
                              length.out = K)), digits = 10)
  lambdapath_vec <- matrix(lambdapath,nrow = 1)
  if(alpha != 0){lambdapath_vec <- lambdapath_vec/(1000*alpha)}
  
  return(lambdapath_vec)
}

## Obtain knot values
knot_values <- function(x,y){
  
  nfolds <- 5 
  lasso_estimate <- cv.glmnet(x,y,family = "binomial", alpha = 0, 
                              standardize = FALSE,intercept = TRUE, nfolds = nfolds,
                              maxit = 500000 )
  
  beta_knot <-  coef(lasso_estimate, x , s="lambda.min")
  beta_knot1 <- as.numeric(beta_knot)
  beta_knot1 <- beta_knot1 + 0.000000001
  return(beta_knot1)
}



Q1 <- function(knot_values){

  Q1_a <- 2*diag(as.vector(knot_values/abs(knot_values)))
  Q1_a[1,1] <- 0

  return(Q1_a)
}

Q2 <- function(q2,d,knot_values){
  
  Q2_a <- diag(as.vector(q2*(d*knot_values - knot_values)*(abs(d*knot_values - knot_values)**(q2-3))))
  Q2_a[1,1] <- 0 
  return(Q2_a)
}

lambda_vecs <- function(lambda,alpha,x){
  lambda1 <- lambda*(1-alpha)/2
  lambda2 <- lambda*alpha
  
  lambda_vec1 <- rbind(0,as.matrix(rep(lambda1,(ncol(x))),ncol=1))
  lambda_vec2 <- rbind(0,as.matrix(rep(lambda2,(ncol(x))),ncol=1))
  lambda_vecs <- cbind(lambda_vec1,lambda_vec2)
  return(lambda_vecs)
}

predicted_p <- function(x,beta){
  x <- cbind(1,x)
  pred_p <-  matrix((exp(x%*%beta)/(1 + exp(x%*%beta))),ncol = 1)
  return(pred_p)
}

predicted_y <- function(pred_p,threshold){
  pred_y <- matrix(ifelse(pred_p >= threshold,1,0),ncol = 1) 
  return(pred_y)
}

average_prediction_error <- function(pred_p,y){
  mean_pred_error <- mean(abs(y - pred_p))
  return(mean_pred_error)
}

area_under_the_roc_curve <- function(y,predicted_p){
  roc_obj <- roc(as.vector(y), as.vector(predicted_p))
  auc <- auc(roc_obj)
  return(auc)
}

beta_estimates <- function(x,y,lambda,alpha,q2,it,beta_initial,Q1,Q2,knot_values){
  
  beta_old <- as.matrix(beta_initial,ncol = 1)
  x1 <- cbind(1,x)
  lambda1 <- lambda*(1-alpha)/2
  lambda2 <- lambda*alpha
  
  lambda_vec1 <- lambda_vecs(lambda,alpha,x)[,1]
  lambda_vec2 <- lambda_vecs(lambda,alpha,x)[,2]
  
  
  for(i in 1:it){
    p <- as.vector(exp(x1%*%beta_old)/(1 + exp(x1%*%beta_old)))
    w <- diag(p*(1-p))
    
    beta_new <- beta_old + 
      solve((t(x1)%*%w%*%x1)*(1/nrow(x1)) + lambda1*Q1 + lambda2*Q2)%*%
      ((t(x1)%*%(y - p))*(1/nrow(x1)) - lambda_vec1*(Q1%*%beta_old) - lambda_vec2*(Q2%*%beta_old))
    
    converge <- max(abs(beta_new - beta_old))
    beta_old <- beta_new
    
    if(converge<0.00001){break}
    
  }
  return(matrix(beta_old,ncol = 1))
}

beta_estimates_different_lambda <- function(x,y,alpha,q2,it,d)
{
  beta_estimate_lambda_list <- vector()
  beta_estimate_lambda_list2 <- list()
  beta_initial <- rep(0,ncol(x)+1)
  lambda_vec <- lambda_path(x,y,alpha)
  
  knot_values <- knot_values(x,y)
  Q1 <- Q1(knot_values)
  Q2 <- Q2(q2,d,knot_values)
  
  
  for(lambda_pos in 1:ncol(lambda_vec)){
    lambda <- lambda_vec[lambda_pos]
    beta_estimate_lambda <- beta_estimates(x,y,lambda,alpha,q2,it,beta_initial,Q1,Q2,knot_values)
    beta_estimate_lambda_list <- list(lambda,beta_estimate_lambda)
    beta_initial <- beta_estimate_lambda
    beta_estimate_lambda_list2 <- append(beta_estimate_lambda_list2,
                                         list(beta_estimate_lambda_list))
  }
  return(beta_estimate_lambda_list2)
}

cross_validation <- function(num_folds,x,y,alpha,q2,it,d){
  
  list2 <- list()
  k = num_folds
  folds <- cvFolds(nrow(x), K = k)
  
  for(i in 1:k){
    train_x <- x[folds$subsets[folds$which != i], ] #Set the training set
    test_x <- x[folds$subsets[folds$which == i], ] #Set the validation set
    train_y <- y[folds$subsets[folds$which != i], ] #Set the training set
    test_y <- y[folds$subsets[folds$which == i], ] #Set the validation set
    
    estimate_beta <- beta_estimates_different_lambda(train_x,train_y,
                                                     alpha,q2,it,d)
    
    lambda_vec <- vector()
    beta_vec <- vector()
    auc_vec <- vector()
    average_predict_error_vec <- vector()
    
    for(j in 1:length(estimate_beta)){
      
      lambda <- estimate_beta[[j]][[1]]
      lambda_vec <- cbind(lambda_vec,lambda)
      
      beta_use <- matrix(estimate_beta[[j]][[2]], ncol=1)
      beta_vec <- cbind(beta_vec,beta_use)  
      average_predict_error <- average_prediction_error(predicted_p(test_x,beta_use),
                                                        test_y)
      average_predict_error_vec <- cbind(average_predict_error_vec,average_predict_error)
      auc <- area_under_the_roc_curve(train_y,predicted_p(train_x,beta_use))
      auc_vec <- cbind(auc_vec,auc)
      
      
    }
    list1 <- list(lambda_vec,beta_vec,
                  auc_vec,average_predict_error_vec)
    
    list2 <- append(list2,list(list1))
  }
  
  return(list2) 
}

choose_optimal <- function(num_folds,x,y,alpha,q2,it,auc_or_prediction_error,d){
  
  build_crossvalidation <- cross_validation(num_folds,x,y,alpha,
                                            q2,it,d)
  list2 <- list()
  
  for(fold_k in 1:num_folds){
    list1 <- list()
    for(lambda in 1:ncol(build_crossvalidation[[fold_k]][[1]])){
      
      z1 <- which(build_crossvalidation[[fold_k]][[3]] == max(build_crossvalidation[[fold_k]][[3]]))
      z2 <- which(build_crossvalidation[[fold_k]][[4]] == min(build_crossvalidation[[fold_k]][[4]]))
      
      build_optimal <- list(build_crossvalidation[[fold_k]][[1]][,z1[length(z1)]],
                            max(build_crossvalidation[[fold_k]][[3]]),
                            build_crossvalidation[[fold_k]][[2]][,z1[length(z1)]],
                            build_crossvalidation[[fold_k]][[1]][,z2[length(z2)]],                            
                            min(build_crossvalidation[[fold_k]][[4]]),
                            build_crossvalidation[[fold_k]][[2]][,z2[length(z2)]])
      
      list1 <- append(list1, list(build_optimal))
    }
    
    if(auc_or_prediction_error == 2)
    {
      choose_auc <- vector()
      for(i in 1:length(build_crossvalidation))
      {
        choose_auc <- cbind(choose_auc,list1[[i]][[2]])
      }
      choose_lambda <- list1[[which.max(choose_auc)]][1:3]
      
    }    
    else if(auc_or_prediction_error == 5)
    {
      choose_ape <- vector()
      for(j in 1:length(build_crossvalidation))
      {
        choose_ape <- cbind(choose_ape,list1[[j]][[5]])
      }  
      choose_lambda <- list1[[which.min(choose_ape)]][4:6]
      
    }
    
    
    else{choose_lambda <- "Wrong choice of auc_or_prediction_error"}
    
    list2 <- append(list2,list(choose_lambda))
  }
  return(list2)
}

optimal_optimal <- function(num_folds,x,y,alpha,q2,it,auc_or_prediction_error,d){
  
  choose_optimal_1 <- choose_optimal(num_folds,x,y,alpha,q2,it,auc_or_prediction_error,d) 
  print(choose_optimal_1)
  
  if(auc_or_prediction_error == 2)
  {  
    choose_auc <- vector()
    lambda_auc <- vector()
    for(i in 1:num_folds){
      choose_auc <- cbind(choose_auc,choose_optimal_1[[i]][[2]])
      lambda_auc <- cbind(lambda_auc,choose_optimal_1[[i]][[1]])
      
    }
    z3<-which(choose_auc==max(choose_auc))
    lambda_z3 <- lambda_auc[z3]
    z4 <- which(lambda_auc == min(lambda_z3))
    
    optimal <- choose_optimal_1[[z4]]
    
  }
  
  else if(auc_or_prediction_error == 5)
  {
    choose_ape <- vector()
    for(i in 1:num_folds){
      choose_ape <- cbind(choose_ape,choose_optimal_1[[i]][[2]])
    }
    optimal <- choose_optimal_1[[which.min(choose_ape)]]
    
  }
  
  else{optimal <- "wrong choice of auc_or_prediction_error"}
  
  
  optimal_list <- list(optimal)
  return(optimal_list)
}

optimal_with_thresholds <- function(num_folds,x,y,alpha,q2,it,auc_or_prediction_error,quantile_sequence,d){
  
  optimal <- optimal_optimal(num_folds,x,y,alpha,q2,it,auc_or_prediction_error,d)
  lambda <- optimal[[1]][[1]]
  beta <- matrix(optimal[[1]][[3]],ncol = 1)
  beta_estimate_threshold_vec <- beta
  beta_drop_intercept <- beta[2:nrow(beta),]
  threshold_vec <- 0
  
  
  for(j in quantile_sequence){
    beta_drop_intercept_temp <- matrix(beta_drop_intercept,ncol = 1)
    beta_drop_intercept_temp[abs(beta_drop_intercept) <= 
                               quantile(abs(beta_drop_intercept),j),] <- 0
    beta_threshold <- rbind(beta[1,],beta_drop_intercept_temp)
    beta_estimate_threshold_vec <- cbind(beta_estimate_threshold_vec,
                                         beta_threshold)
    threshold_vec <- cbind(threshold_vec, quantile(abs(beta_drop_intercept),j))
  }
  
  beta_estimate_lambda_threshold_list <- list(lambda,threshold_vec,beta_estimate_threshold_vec)
  return(beta_estimate_lambda_threshold_list)
}

optimal_with_thresholds_performance <- function(num_folds,x,y,alpha,q2,it,
                                                auc_or_prediction_error,quantile_sequence,d){
  
  optimal <- optimal_with_thresholds(num_folds,x,y,alpha,q2,it,
                                     auc_or_prediction_error,quantile_sequence,d)
  average_predict_error_vec <- vector()
  auc_vec <- vector()
  for(w in 1:ncol(optimal[[3]]))
  {
    beta_use <- matrix(optimal[[3]][,w], ncol=1)
    average_predict_error <- average_prediction_error(predicted_p(x,beta_use),
                                                      y)
    average_predict_error_vec <- cbind(average_predict_error_vec,average_predict_error)
    
    auc <- area_under_the_roc_curve(y,predicted_p(x,beta_use))
    auc_vec <- cbind(auc_vec,auc)
    
  }
  list1 <- list(optimal[[1]],optimal[[2]],optimal[[3]],auc_vec,average_predict_error_vec)
  return(list1)
}

simulation2 <- function(num_sim, num_folds,x,y,alpha,q2,it,
                        auc_or_prediction_error,quantile_sequence,d){
  set.seed(num_sim)

  sample_ind <- sample(1:nrow(x),100,replace = FALSE)
  training_data_x <- matrix(x[-sample_ind,], ncol = ncol(x))
  training_data_y <- matrix(y[-sample_ind,], ncol = 1)
  
  testing_data_x <- matrix(x[sample_ind,], ncol = ncol(x))
  testing_data_y <- matrix(y[sample_ind,], ncol = 1)
  
  result <-  optimal_with_thresholds_performance(num_folds,training_data_x,training_data_y,
                                                 alpha,q2,it,auc_or_prediction_error,
                                                 quantile_sequence,d)
  
  average_prediction_error_test_vec <- vector()
  accuracy_vec <- vector()
  
  for(j in 1:ncol(result[[3]])){
    
    beta_use_all <- result[[3]][,j]
    pred_p_test <- predicted_p(testing_data_x,matrix(beta_use_all,ncol = 1))
    pred_y_test <- predicted_y(pred_p_test,0.5)
    average_prediction_error_test <- average_prediction_error(pred_p_test,
                                                              testing_data_y)
    average_prediction_error_test_vec <- cbind(average_prediction_error_test_vec,
                                               average_prediction_error_test)
    accuracy <- confusionMatrix(as.factor(pred_y_test),
                                as.factor(testing_data_y))
    accuracy_test1 <- accuracy$overall[1]
    accuracy_vec <- cbind(accuracy_vec,accuracy_test1)
  }
  
  test_performance_all <- list(average_prediction_error_test_vec,accuracy_vec)
  
  named_results_all <- list(lambda = result[[1]],
                            threshold = result[[2]],
                            beta = result[[3]],
                            auc_train = result[[4]],
                            ape_train = result[[5]],
                            ape_test = test_performance_all[[1]],
                            accuracy_test = test_performance_all[[2]] )
  
  
  return(named_results_all)
}

quantile_sequence <- seq(0.01,0.99,0.01)

empty_vec <- 1:500
num_folds = 5
q2 = 0.75
alpha = 0.25
d = 0.1
auc_or_pred = 2
alpha0.25_q0.75_d0.1 <- mclapply(empty_vec, simulation2, num_folds = num_folds, x = x, y = y,
                                           alpha = alpha, q2 = q2, it = 1000,
                                           auc_or_prediction_error = auc_or_pred,
                                           quantile_sequence = quantile_sequence,d = d, mc.cores = 50)

