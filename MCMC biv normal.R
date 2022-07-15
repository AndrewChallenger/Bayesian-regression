#Create a function to streamline


compare_models <- function(n0, n1, p, mu0 = rep(0,p) , mu1, Sigma0 = diag(1, nrow = p), Sigma1 = diag(1, nrow = p)){
  

library(mvtnorm)

#Draw some data from a multivariate normal for Class 0:
x0 <- rmvnorm(n0, mu0, Sigma0)
x1 <- rmvnorm(n1, mu1, Sigma1)

y0 <- rep(0,n0) # Class 0 labels
y1 <- rep(1,n1) # Class 1 labels

x <- rbind(x0, x1)
y <- as.factor(c(y0, y1))

myData <- data.frame(x, label = y)

#Splitting data
library(caTools)
split <- sample.split(myData$label, 0.25)
trData <- subset(myData, split == TRUE)
tstData <- subset(myData, split == FALSE)

#MLE 
myglm <- glm(label ~ ., data = trData, family = "binomial")
intercept_only <- glm(label ~ 1, data = trData, family = "binomial")
myglm <- step(intercept_only, direction = 'both', scope = formula(myglm), trace = 0)


#Applying to the test data

library(dplyr)
preds <- predict(myglm, select(tstData, -label), type = 'response')
preds <- ifelse(preds > 0.5, 1, 0)
#confusion matrix for MLE
MLE_T <- table(preds, tstData$label)

MLE_accuracy <- (MLE_T[1,1]+MLE_T[2,2])/length(tstData$label)

cat('Accuarcy for MLE is', MLE_accuracy, '\n')

#See if bayesreg can improve this

library(bayesreg)

#First ridge prior
rmodel <- bayesreg(label ~., data = trData, model = 'logistic', prior = 'ridge', n.samples = 2000)

rpreds <- predict(rmodel, select(tstData, -label), type = 'response')
rpreds <- ifelse(rpreds > 0.5, 1, 0)
rid_T <- table(rpreds, tstData$label)

rid_accuracy <- (rid_T[1,1]+rid_T[2,2])/length(tstData$label)

cat('Accuarcy for ridge is', rid_accuracy, '\n')

#Lasso prior

lmodel <- bayesreg(label ~., data = trData, model = 'logistic', prior = 'lasso', n.samples = 2000)

lpreds <- predict(lmodel, select(tstData, -label), type = 'response')
lpreds <- ifelse(lpreds > 0.5, 1, 0)
las_T <- table(lpreds, tstData$label)

las_accuracy <- (las_T[1,1]+las_T[2,2])/length(tstData$label)

cat('Accuarcy for lasso is', las_accuracy, '\n')

#Horseshoe prior

hmodel <- bayesreg(label ~., data = trData, model = 'logistic', prior = 'hs', n.samples = 2000)

hpreds <- predict(hmodel, select(tstData, -label), type = 'response')
hpreds <- ifelse(hpreds > 0.5, 1, 0)
hs_T <- table(hpreds, tstData$label)

hs_accuracy <- (hs_T[1,1]+hs_T[2,2])/length(tstData$label)

cat('Accuarcy for horeshoe is', hs_accuracy, '\n')



#See how a random forest compares

library(randomForest)

rfmodel <- randomForest(label~., data = trData)
rfmodel

rfpreds <- predict(rfmodel, select(tstData, -label))
rf_T <- table(rfpreds, tstData$label)

rf_accuracy <- (rf_T[1,1]+rf_T[2,2])/length(tstData$label)

cat('Accuarcy for random forest is', rf_accuracy, '\n')
}

n0 <- 200
n1 <- 200
p <- 50
mu1 <- c(rep(0,45), rep(0.5,5))

compare_models(n0=n0, n1=n1,p=p, mu1=mu1)
