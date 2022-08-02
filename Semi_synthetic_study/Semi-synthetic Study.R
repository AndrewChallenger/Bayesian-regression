#Semi synthetic simulation using all the priors on Bayesreg. Will test the prediction power with 
#AUC score, and testing variable selection with proportion of true positives, true positives, false 
#positives, false negatives

#Will just use SRRCT data for now
library(plsgenomics)
library(pROC)
library(bayesreg)
library(glmnet)

set.seed(1)
data(SRBCT) 

# 83 data entries with 2308 variables. We'll take 75% of the data to train model and the rest to test.
#Will run the models on lots of different scenarios as in Paul's code. 

#This is Paul's function to generate the semisynthetic data
generateSemisyntheticDataset <- function(X, p = 250, s0 = 5, correlation = "Normal" , muddleIrrelevantColumns = FALSE)
{
  
  #Standardise the data, such that columns have mean 0 and s.d. 1:
  X <- scale(X)
  
  originalDimension <- ncol(X)
  
  if(p > originalDimension)
  {
    stop("Dimension of original dataset is less than p")
  }
  else
  {
    sampledIndices <- sample(1:originalDimension, size = p)
  }
  
  X <- X[,sampledIndices]
  
  # X is now of size n x p
  
  if(s0 > p)
  {
    stop("s0 must be less than (or equal to) p")
  }
  else
  {
    if(correlation == "Normal")
    {
      relevantIndices <- sample(1:p, size = s0)
    }
    else
    {
      relevantIndex   <- sample(1:p, size = 1)
      
      # Now add in variables that have the highest absolute correlation with the selected variable
      dataForCurrentVariable            <- X[, relevantIndex]
      absCorrelationWithCurrentVariable <- abs(cor(dataForCurrentVariable,X))
      relevantIndices                   <- sort.int(absCorrelationWithCurrentVariable, decreasing = T, index.return = T)$ix[1:s0]
    }
  }
  
  
  #Create the vector of coefficients, initialising by setting all to 0:
  betaVec <- vector(mode="numeric", length = p)
  #Randomly set the elements given in relevantIndices to -2 or 2
  betaVec[relevantIndices] <- (4*rbinom(s0, 1,.5)) - 2
  
  if(muddleIrrelevantColumns)
  {
    for(ind in setdiff(1:p, relevantIndices))
    {
      X[,ind] <- X[sample(1:nrow(X)),ind]
    }
    
  }
  
  
  #Calculate the linear predictor:
  linearPredictor <- X %*% betaVec
  
  #Use sigmoid function to transform to probabilities
  probs           <- plogis(linearPredictor)
  
  #Sample the binary outcome according to these probs:
  Y <- rbinom(length(probs), 1, probs)
  
  return(list(X = X, Y = Y, probs = probs, linearPredictor = linearPredictor, betaVec = betaVec)) 
}

#Gives back summary stats for bayesreg models
get_summary_stats_bs <- function(model, test_data, variables, current_data){
  
  #Prediction accuracy
  predicted <- predict(model, test_data, type = 'class')
  AUC <- auc(as.numeric(test_data$label), as.numeric(predicted))
  #Confusion matrix - take 1 to be a positive
  reality <- factor(test_data$label)
  predicted <- factor(predicted)
  levels(predicted) <- c(levels(predicted), c('0','1'))
  levels(reality) <- c(levels(reality), c('0','1'))
  cnfmatrix <- table(reality, predicted)
  TPR <- cnfmatrix[2,2]/(cnfmatrix[2,1]+cnfmatrix[2,2])
  FPR <- cnfmatrix[1,2]/(cnfmatrix[1,1]+cnfmatrix[1,2])
  TNR <- 1 - FPR
  FNR <- 1 - TPR
  
  #Variable selection
  real_sig_vbls <- which(current_data$betaVec != 0)
  #True positives - correctly identifies sig variable
  true_pos <- length(variables[variables %in% real_sig_vbls])
  prop_TP <- true_pos/length(real_sig_vbls)
  #False positives - incorrectly identifies what is sig
  fls_pos <- length(variables[!(variables %in% real_sig_vbls)])
  prop_FP <- fls_pos/(length(current_data$betaVec) - length(real_sig_vbls))
  #True negative rate
  prop_TN <- 1 - prop_FP
  #False negative rate
  prop_FN <- 1 - prop_TP
  output <- c(AUC, TPR, FPR, TNR, FNR, prop_TP, prop_FP, prop_TN, prop_FN)
  
  return(output)
  
}

#Same for the glmnet 
get_summary_stats_glm <- function(model, testX, testY, current_data){
  
  predicted <- predict(model, as.matrix(testX), type = 'class')
  #Assessing prediction
  AUC <- auc(as.numeric(testY), as.numeric(predicted)) 
  
  #Confusion matrix - take 1 to be a positive
  reality <- factor(testY)
  predicted <- factor(predicted)
  levels(predicted) <- c(levels(predicted), c('0','1'))
  levels(reality) <- c(levels(reality), c('0','1'))
  cnfmatrix <- table(reality, predicted)
  TPR <- cnfmatrix[2,2]/(cnfmatrix[2,1]+cnfmatrix[2,2])
  FPR <- cnfmatrix[1,2]/(cnfmatrix[1,1]+cnfmatrix[1,2])
  TNR <- 1 - FPR
  FNR <- 1 - TPR
  
  #Assessing variable selection
  sig_vbls <- tail(which(coef(model) !=0) - 1, -1) #accounting for intercept
  real_sig_vbls <- which(current_data$betaVec != 0)
  #True positives - correctly identifies sig variable
  true_pos <- length(sig_vbls[sig_vbls %in% real_sig_vbls])
  prop_TP <- true_pos/length(real_sig_vbls)
  #False postives - incorrectly identifies what is sig - is this code right?
  fls_pos <- length(sig_vbls[!(sig_vbls %in% real_sig_vbls)])
  prop_FP <- fls_pos/(length(testX) - length(real_sig_vbls))
  #True negative rate
  prop_TN <- 1 - prop_FP
  #False negative rate
  prop_FN <- 1 - prop_TP
  
  return(c(AUC, TPR, FPR, TNR, FNR, prop_TP, prop_FP, prop_TN, prop_FN))
  
}

#Gets new model with updated predictors for 75% CI - uses ridge prior as discussed
#Returns summary_stats
get_new_model75 <- function(model, trainX, trainY, testX, testY, current_data){
  
  model_sum <- summary(model)
  variables <- which(head(model_sum$n.stars, -1) >= 1)  # This gives the indices of the variables with 1+ stars
  
  if (length(variables) > 0){
    #Train new model with just those predictors
    new_trainX <- trainX[,variables]
    new_testX <- testX[,variables]
    vbl_names <- rownames(model_sum$n.stars)[which(model_sum$n.stars >= 1)]
    new_train_data <- data.frame(new_trainX, label = trainY)
    colnames(new_train_data) <- c(vbl_names, 'label')
    new_test_data <- data.frame(new_testX, label = testY)
    colnames(new_test_data) <- c(vbl_names, 'label')
    
    new_model <- bayesreg(label ~., data = new_train_data, model = 'logistic', prior = 'ridge', n.samples = n_sam, burnin = burn, thin = thin)
    
    summary_stats <- get_summary_stats_bs(new_model, new_test_data, variables, current_data)
  }
  else{
    summary_stats <- rep(NA, outputs)
  }
  
  return(summary_stats)
}

#Same for 95% CI
get_new_model95 <- function(model, trainX, trainY, testX, testY, current_data){
  
  model_sum <- summary(model)
  variables <- which(head(model_sum$n.stars, -1) == 2)  # This gives the indices of the variables with 1+ stars
  
  if (length(variables) > 0){ #actually have some predictors
    
    #Train new model with just those predictors
    new_trainX <- trainX[,variables]
    new_testX <- testX[,variables]
    vbl_names <- rownames(model_sum$n.stars)[which(model_sum$n.stars == 2)]
    new_train_data <- data.frame(new_trainX, label = trainY)
    colnames(new_train_data) <- c(vbl_names, 'label')
    new_test_data <- data.frame(new_testX, label = testY)
    colnames(new_test_data) <- c(vbl_names, 'label')
    
    new_model <- bayesreg(label ~., data = new_train_data, model = 'logistic', prior = 'ridge', n.samples = n_sam, burnin = burn, thin = thin)
    summary_stats <- get_summary_stats_bs(new_model, new_test_data, variables, current_data)
  }
  
  else { #no predictors
    summary_stats <- rep(NA, outputs)
  }
  
  return(summary_stats)
}



#What we'll be simulating
pRange           <- c(250, 1000)
s0Range          <- c(5, 20)
correlationRange <- c("Normal", "High")
muddleRange      <- c(TRUE, FALSE)

#MCMC parameters
n_sam <- 3000
burn <- 15000
thin <- 5

#Other parameters
test_proportion <- 0.25 #proportion of data in test set
outputs <- 9 #number of outputs to save for each model
c_names <- c('Ridge_soft', 'Ridge_hard75', 'Ridge_hard95', 'Lasso_soft', 'Lasso_hard75', 'Lasso_hard95', 'Horseshoe_soft', 'Horseshoe_hard75', 'Horseshoe_hard95', 'Glmnet_0.1', 'Glmnet_0.5', 'Glmnet_1')
r_names <- c('AUC', 'TPR', 'FPR', 'TNR', 'FNR', 'TPR (variable)',  'FPR (variable)', 'TNR (variable)', 'FNR (variable)')

currentDataset <- SRBCT

#Replicating for loop
p_inds <- c(rep(1,40), rep(2, 40))
s0_inds <- rep(c(rep(1, 20), rep(2,20)), 2)
correlation_inds <- rep(c(rep(1, 10), rep(2, 10)), 4)
muddle_inds <- rep(c(rep(1, 5), rep(2, 5)), 8)
sims_inds <- rep(1:5,16)



#Function that takes an index, and saves a file with the summary stats that corresponds to that index

analysis <- function(n){
  
  p_ind <- p_inds[n]
  s0_ind <- s0_inds[n]
  correlation_ind <- correlation_inds[n]
  muddle_ind <- muddle_inds[n]
  
  p_ <- pRange[p_ind]
  s0_ <- s0Range[s0_ind]
  correlation_ <- correlationRange[correlation_ind]
  muddle_ <- muddleRange[muddle_ind]
  sim_id <- sims_inds[n]
  
  currentDatasetX <- currentDataset$X
  currentSyntheticDataset <- generateSemisyntheticDataset(currentDatasetX, p = p_, s0 = s0_, 
                                                          correlation = correlation_,
                                                          muddleIrrelevantColumns = muddle_)
  #Setting up train and test data
  N <- length(currentSyntheticDataset$Y)
  test_indices <- sample(1:N, floor(N*test_proportion))
  train_indices <- setdiff(1:N, test_indices)
  
  trainX = as.data.frame(currentSyntheticDataset$X[train_indices,])
  trainY = as.factor(currentSyntheticDataset$Y[train_indices])
  
  testX = as.data.frame(currentSyntheticDataset$X[test_indices,])
  testY = as.factor(currentSyntheticDataset$Y[test_indices])
  
  train_data <- data.frame(trainX, label = trainY)
  test_data <- data.frame(testX, label = testY)
  
  #Setting up dataframe where we'll store the outputs for each n
  results <- as.data.frame(matrix(rep(0, 12*outputs), ncol = 12))
  colnames(results) <- c_names
  rownames(results) <- r_names
  
  
  
  
  
  #Ridge prior - soft
  model <- bayesreg(label ~., data = train_data, model = 'logistic', prior = 'ridge', n.samples = n_sam, burnin = burn, thin = thin)
  variables <- which(model$mu.beta < Inf) #very weird I know
  summary_stats <- get_summary_stats_bs(model, test_data, variables, currentSyntheticDataset)
  results$Ridge_soft <- summary_stats
 
  
  
  #Ridge prior - hard (less strict 75% CI acceptance)
  summary_stats <- get_new_model75(model, trainX = trainX, trainY = trainY, testX = testX, testY = testY, currentSyntheticDataset)
  results$Ridge_hard75 <- summary_stats
  
  
  #Ridge prior - hard (more strict 95% CI acceptance)
  summary_stats <- get_new_model95(model, trainX = trainX, trainY = trainY, testX = testX, testY = testY, currentSyntheticDataset)
  results$Ridge_hard95 <- summary_stats
  
  
  #Lasso prior - soft
  model <- bayesreg(label ~., data = train_data, model = 'logistic', prior = 'lasso', n.samples = n_sam, burnin = burn, thin = thin)
  variables <- which(model$mu.beta < Inf) #very weird I know
  summary_stats <- get_summary_stats_bs(model, test_data, variables, currentSyntheticDataset)
  results$Lasso_soft <- summary_stats
  
  
  #Lasso prior - hard(less strict)
  summary_stats <- get_new_model75(model, trainX = trainX, trainY = trainY, testX = testX, testY = testY, currentSyntheticDataset)
  results$Lasso_hard75 <- summary_stats
  
  
  #Lasso prior - hard(more strict)
  summary_stats <- get_new_model95(model, trainX = trainX, trainY = trainY, testX = testX, testY = testY, currentSyntheticDataset)
  results$Lasso_hard95 <- summary_stats
  
  
  #Horseshoe prior - soft
  model <- bayesreg(label ~., data = train_data, model = 'logistic', prior = 'hs', n.samples = n_sam, burnin = burn, thin = thin)
  variables <- which(model$mu.beta < Inf) #very weird I know
  summary_stats <- get_summary_stats_bs(model, test_data, variables, currentSyntheticDataset)
  results$Horseshoe_soft <- summary_stats
  
  
  #Horseshoe prior - hard (less strict)
  summary_stats <- get_new_model75(model, trainX = trainX, trainY = trainY, testX = testX, testY = testY, currentSyntheticDataset)
  results$Horseshoe_hard75 <- summary_stats
  
  
  #Horseshoe prior - hard (more strict)
  summary_stats <- get_new_model95(model, trainX = trainX, trainY = trainY, testX = testX, testY = testY, currentSyntheticDataset)
  results$Horseshoe_hard95 <- summary_stats
  
  
  
  #Glmnet alpha = 0.1
  model <- cv.glmnet(as.matrix(trainX), trainY, family = 'binomial', alpha = 0.1)
  sum_stats <- get_summary_stats_glm(model, as.matrix(testX), testY, currentSyntheticDataset)
  results$Glmnet_0.1 <- sum_stats
  
  
  
  #Glmnet alpha = 0.5
  model <- cv.glmnet(as.matrix(trainX), trainY, family = 'binomial', alpha = 0.5)
  sum_stats <- get_summary_stats_glm(model, as.matrix(testX), testY, currentSyntheticDataset)
  results$Glmnet_0.5 <- sum_stats
  
  
  #Glmnet alpha = 1
  model <- cv.glmnet(as.matrix(trainX), trainY, family = 'binomial', alpha = 1)
  sum_stats <- get_summary_stats_glm(model, as.matrix(trainX), trainY, currentSyntheticDataset)
  results$Glmnet_1 <- sum_stats
  
  results <- round(results, 3)
  
  setwd("output")
  save(results, file = paste0(p_, "-", s0_, "-", correlation_, "-", muddle_, "-", sim_id, ".RData"))
}

task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id <- as.numeric(task_id_string)
analysis(task_id)
























