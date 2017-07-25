library(leaps)

create_folds <- function(data, k) {
  fold_indeces <- sample(1:nrow(data)%%k+1)
  
  train_folds <- list()
  test_folds <- list()
  for (i in 1:k){
    train = data[fold_indeces != i,]
    test = data[fold_indeces == i,]
    train_folds[[i]] <- train
    test_folds[[i]] <- test
  }
  return (list(train=train_folds, test=test_folds))
}

get_fold_errors <- function(k_folds, response_name) {
  num_features <- dim(k_folds$train[[1]])[2]-1 #subtract the response variable
  num_folds <- length(k_folds[[1]])
  model_info <- vector()
  drops_model <- list()
  
  for (p in 2:num_features) {
    fold_errors <- vector()
    drops_folds <- list()
    for (k in 1:num_folds) {
      #split X and y
      X_train <- k_folds$train[[k]][,-which(colnames(k_folds$train[[k]])==response_name)] #drop response
      y_train <- k_folds$train[[k]][,which(colnames(k_folds$train[[k]])==response_name)] #use only response
      
      #find best subset model
      drop_predictors <- find_best_subset(X_train, y_train, p)
      
      #now use the best model on the test data
      X_test <- k_folds$test[[k]][,-which(colnames(k_folds$test[[k]])==response_name)] #drop response
      y_test <- k_folds$test[[k]][,which(colnames(k_folds$test[[k]])==response_name)] #use only response
      
      if (length(drop_predictors) > 0){
        X_test <- X_test[,!(names(X_test) %in% names(drop_predictors))] #drop bad predictors
      }
      
      X_train_subset <- X_train[,!(names(X_train) %in% names(drop_predictors))]
      if(is.null(dim(X_train_subset)[1])){ #this is a vector, not matrix
        X_train_subset <- as.data.frame(X_train_subset)
        X_test <- as.data.frame(X_test)
      } 
      train_subset_model <- lm(y_train ~ ., data=X_train_subset)
      y_hat_test <- predict(train_subset_model, newdata = X_test)
      
      fold_error <- mean((y_hat_test - y_test)^2) #MSE
      fold_errors <- c(fold_errors, fold_error)
      drops_folds[[k]] <- drop_predictors
    }
    mean_folds <- mean(fold_errors)
    sd_folds <- sqrt(var(fold_errors)/length(fold_errors))
    model_info <- rbind(model_info, c(mean=mean_folds, sd=sd_folds))
    drops_model[[p]] <- drops_folds[[1]]
  }
  
  return(list(model_info=model_info, drops=drops_model))
}

get_best_model_ols <- function(k_folds, response_name) {
  num_folds <- length(k_folds[[1]])
  
  fold_errors <- vector()
  models <- list()
  df_values <- list()
  for (k in 1:num_folds) {
      #split X and y
      X_train <- k_folds$train[[k]][,-which(colnames(k_folds$train[[k]])==response_name)] #drop response
      y_train <- k_folds$train[[k]][,which(colnames(k_folds$train[[k]])==response_name)] #use only response
      
      #run ols
      ols_model <- lm(y_train ~ ., data=X_train)
      
      #now use the model on the test data
      X_test <- k_folds$test[[k]][,-which(colnames(k_folds$test[[k]])==response_name)] #drop response
      y_test <- k_folds$test[[k]][,which(colnames(k_folds$test[[k]])==response_name)] #use only response
      y_hat_test <- predict(ols_model, newdata = X_test)
      
      fold_error <- mean((y_hat_test - y_test)^2) #MSE
      fold_errors <- c(fold_errors, fold_error)
      models[[k]] <- ols_model
      df_values[[k]] <- k_folds$train[[k]]
  }
  best_model <- models[[which.min(fold_errors)]]
  df <- df_values[[which.min(fold_errors)]]
  return(list(best_model=best_model, df=df))
}

find_best_subset <- function(X_train, y_train, p){
  #reconstruct the formula since passing in x and y crashes R
  subset_formula <- as.formula("y_train ~ .")
  subsets_obj <- regsubsets(subset_formula, data=X_train, nvmax=(p))
  subsets_summ <- summary(subsets_obj)
  best_model_index <- which.min(subsets_summ$bic)
  #find out which predictors to drop
  drops <- which(!subsets_summ$which[best_model_index,])
  
  return (drops)
}

get_complexity_for_simplest_model <- function(subset_errors) {
  model_info <- subset_errors$model_info
  num_features <- dim(model_info)[1]
  x=1:num_features
  y=model_info[,1]
  sd=model_info[,2]
  best_model_num <- which.min(model_info[,1])
  plus_one_sd <- model_info[best_model_num,][1] + model_info[best_model_num,][2]
  most_parsimonious_complexity <- which(model_info[,1] <= plus_one_sd)[1]
  
  #plot the graph
  ymin <- min(y-sd)
  ymax <- max(y+sd)
  plot(x, y, xlab="Subset Size p", ylab="Mean Prediction Error", ylim=c(ymin, ymax))
  arrows(x,y-sd,x,y+sd, code=3, length=0.02, angle = 90)
  abline(h=plus_one_sd, lty="dashed")
  
  return (most_parsimonious_complexity)
}

cv_best_model_ols <- function(data_train, num_folds, response_name){
  #create folds
  k_folds <- create_folds(data_train, num_folds)
  
  #get best model
  best_model <- get_best_model_ols(k_folds, response_name)
  
  return(best_model)
}

cv_best_model_subset <- function(data_train, num_folds, response_name){
  #create folds
  k_folds <- create_folds(data_train, num_folds)
  
  #get fold errors
  fold_errors <- get_fold_errors(k_folds, response_name)
  
  return(fold_errors)
}

rmse = function(a,b) { sqrt(mean((a-b)^2)) }

#find predictors with unit variance
get_unit_var_predictors <- function(data, threshold){
  var_data <- var(data)
  unit_var_predictors <- vector()
  for(i in 1:ncol(var_data)){
    if(abs(var_data[i,i])<=threshold){
      unit_var_predictors <- c(unit_var_predictors, colnames(var_data)[i])
    }
  }
  return(unit_var_predictors)
}

drop_combine_correlated_predictors <- function(data, response_name){
  cor_matrix=cor(data)
  #find highly-correlated predictors
  cor_matrix = cor(data[ , !(names(data) %in% response_name)])
  ncol=dim(cor_matrix)[2] 
  collinear = list()
  collinear2=list()
  ind1<-1
  ind2<-1
  for (j in 1:(ncol-1)){
    for (i in (j+1):ncol){
      if (abs(cor_matrix[i,j]) >= 0.99){
        collinear[ind1] = colnames(data[ , !(names(data) %in% response_name)])[j]    ##If cor>0.95 put the name of the two preditors sequentially in collinear matrix
        collinear[ind1+1]=colnames(data[ , !(names(data) %in% response_name)])[i]
        collinear[ind1+2]=cor_matrix[i,j]
        ind1=ind1+3
      } 
      if (abs(cor_matrix[i,j]) > 0.95 && abs(cor_matrix[i,j])<0.99){
        collinear2[ind2] = colnames(data[ , !(names(data) %in% response_name)])[j]    ##If cor>0.95 put the name of the two preditors sequentially in collinear matrix
        collinear2[ind2+1]=colnames(data[ , !(names(data) %in% response_name)])[i]
        collinear2[ind2+2]=cor_matrix[i,j]
        ind2=ind2+3
      }
    }
  }
  
  ##Drop the highly correlated predictors in collinear (cor>0.99).  
  ##The way the list is constructed above, need to grab every 3rd element to remove
  drop_index=seq(1,length(collinear),3)
  drop_names=collinear[drop_index]
  print(drop_names)
  flush.console()
  ##Now drop the predictors in the drop_names list
  data<-data[,!(names(data) %in% drop_names)]
  
  #combine the predictors that are less highly correlated .95<x<.99
  combine_index1=seq(1,length(collinear2),3)
  combine_index2=seq(2,length(collinear2),3)
  combine1=collinear2[combine_index1]
  combine2=collinear2[combine_index2]
  
  ##Combine names and add new Columns to the end of the data
  for (i in 1: length(combine1))
  {
    #first check if we already dropped one of the predictors
    if (!(combine1[i] %in% drop_names | combine2[i] %in% drop_names)){
      rm1=which(names(data) %in% combine1[i])
      rm2=which(names(data) %in% combine2[i])
      new_data=(data[,rm1]+data[,rm2])
      new_data=data.frame(new_data)
      colnames(new_data)=paste0(combine1[i],"_", combine2[i])
      data=cbind(data,new_data)
    }
  }
  
  ##Now remove names in combine1 and 2
  drop_names_2=cbind(combine1,combine2)
  print(drop_names_2)
  flush.console()
  data<-data[,!(names(data) %in% drop_names_2)]
  
  return(data)
}