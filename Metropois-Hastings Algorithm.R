#Creating a MH algorithm for logistic regression data with the ridge prior
library(bayesreg)
library(extraDistr)

#Some simulated data to test - using bivariate normal example first

#assumes mu0 just zeros and no correlation for now
create_MVNdata <- function(n0, n1, mu1){
  #n0, n1 number of draws from each dist, mu1 mean vector for second dist
p <- length(mu1)
mu0 <- rep(0,p)
Sigma0 <- Sigma1 <- diag(1,nrow = p)

library(mvtnorm)

x0 <- rmvnorm(n0, mu0, Sigma0)
x1 <- rmvnorm(n1, mu1, Sigma1)

y0 <- rep(0,n0) # Class 0 labels
y1 <- rep(1,n1) # Class 1 labels

x <- rbind(x0, x1)
y <- c(y0, y1)

myData <- data.frame(x, label = y)
return(myData)
}

data <- create_MVNdata(100, 100, c(2,1,0)) #making it more difficult

# Let's quickly visualise this data
ggplot(data, aes(x=X1, y=X2, group = label) ) +
 stat_density_2d(aes(color = as.factor(label), fill = as.factor(label), alpha = ..level..), geom = "polygon")

# See how MLE logistic regression performs first
mle <- glm(label ~., data = data, family = binomial(link = 'logit')) 
summary(mle)
mle$coefficients #seems good - coef 3 not sig but still quite large as expected with MLE

#See how bayesreg does - all inputs match my algorithm
library(bayesreg)
breg <- bayesreg(as.factor(label) ~.,data = data, model = 'logistic', prior = 'ridge', n.samples = 2000, thin = 10)
beta1 <- breg$beta[1,]
plot(beta1) #looks good
hist(beta1)
summary(breg)
breg$mu.beta
#coefs have shrunk a bit from MLE which is good - but thinks coef 3 is sig - problem 
# with ridge

#Use Paul's t-test check for convergence for beta1
sample1 <- head(beta1, 100)
sample2 <- tail(beta1, 100)
t.test(sample1, sample2) 
#Bayesreg not working great - taking longer to converge and not shrinking coef3 enough
#due to ridge prior



#Set up for my algorithm
library(dplyr)
Y <- data$label
X <- select(data, -label)
burnin <- 1000
n_iters <- 2000
thin <- 10
sd <- c(1.2,1.5,1.5, 1.5,7.5,6) #seemed to be good

#First create function that calculates pointwise log posteriors (proportional)


post <- function(Y, X, beta, beta0, sigma, lambda){ 
  if (sigma > 0 & lambda > 0){
    
  #Calculating log priors - not including beta0 as always vanishes
    pr_sigma <- log(1/sigma) #maybe sigma^2 not sure?
    pr_lambda <- log(dhcauchy(lambda))
    pr_beta <- log(dnorm(beta, 0, sigma/sqrt(lambda)))
    
    #log prior assuming independence
    l_prior <- pr_sigma+pr_lambda+sum(pr_beta)
    logit <- beta0 + as.matrix(X) %*% beta #is this right? 
    prob <- 1/(1+exp(-logit))
    
    l_like <- sum(log(dbinom(Y, size = 1, prob = prob))) #likelihood -is this right?
    
    l_pst <- l_like + l_prior #posterior
  }
  else{
    l_pst <- -Inf
  }
  return(l_pst)
}


#Y is response vector, X is matrix of predictors, number of its, burnin amount, 
#thin, sd is standard deviations for markov chain jumping
MH <- function(Y, X, n_iters, burnin, thin, sd){
  p <- length(X)
  
  #for now choose all of beta starting values as 1, as well as the rest
  initial_beta <- rep(1, p)
  initial_sigma <- initial_lambda <- initial_beta0 <- 1
  
  #Counters for acceptance rate
  count <- rep(0, p+3) #for beta0, lambda, sigma too
  
  total_its <- burnin + n_iters*thin #total iterations
  
  #Now making matrix/vectors to save iterations
  sigma <- lambda <- beta0 <- rep(0, total_its)
  beta <- matrix(rep(0, total_its * p), nrow = p)
  
  beta0[1] <- initial_beta0
  sigma[1] <- initial_sigma
  lambda[1] <- initial_lambda
  beta[,1] <- initial_beta
  
  #Now loop through algorithm
  
  for (n in 2:total_its){
    
    #current values
    cur_beta0 <- beta0[n-1]
    cur_sigma <- sigma[n-1]
    cur_lambda <- lambda[n-1]
    cur_beta <- beta[,n-1]
    
    #generating new value of lambda
    new_lambda <- cur_lambda + rnorm(n=1, mean = 0, sd = sd[p+3])
    #Calculate log(A), where A is acceptance prob
    new_post <- post(Y, X, beta = cur_beta, beta0 = cur_beta0, sigma = cur_sigma, lambda = new_lambda)
    old_post <- post(Y, X, beta = cur_beta, beta0 = cur_beta0, sigma = cur_sigma, lambda = cur_lambda)
    logA <- new_post - old_post
    A <- exp(logA)
    
    if (runif(1) < A){
      lambda[n] <- new_lambda
      count[p+3] <- count[p+3] + 1 #as we've accepted
    }
    else {
      lambda[n] <- cur_lambda
    }
    
    #And repeat for rest
    
    new_sigma <- cur_sigma + rnorm(n=1, mean = 0, sd = sd[p+2])
    new_post <- post(Y, X, beta = cur_beta, beta0 = cur_beta0, sigma = new_sigma, lambda = lambda[n])
    old_post <- post(Y, X, beta = cur_beta, beta0 = cur_beta0, sigma = cur_sigma, lambda = lambda[n])
    logA <- new_post - old_post
    A <- exp(logA)
    
    if (runif(1) < A){
      sigma[n] <- new_sigma
      count[p+2] <- count[p+2] + 1 #as we've accepted
    }
    else {
      sigma[n] <- cur_sigma
    }
    
    new_beta0 <- cur_beta0 + rnorm(n=1, mean = 0, sd = sd[p+1])
    new_post <- post(Y, X, beta = cur_beta, beta0 = new_beta0, sigma = sigma[n], lambda = lambda[n])
    old_post <- post(Y, X, beta = cur_beta, beta0 = cur_beta0, sigma = sigma[n], lambda = lambda[n])
    logA <- new_post - old_post
    A <- exp(logA)
    
    if (runif(1) < A){
      beta0[n] <- new_beta0
      count[p+1] <- count[p+1] + 1 
    }
    else {
      beta0[n] <- cur_beta0
    }
    
    
    for (m in 1:p){ #for each variable in beta
      new_beta <- cur_beta
      new_beta[m] <- cur_beta[m] + rnorm(n=1, mean = 0, sd = sd[m])
      
      new_post <- post(Y, X, beta = new_beta, beta0 = beta0[n], sigma = sigma[n], lambda = lambda[n])
      old_post <- post(Y, X, beta = cur_beta, beta0 = beta0[n], sigma = sigma[n], lambda = lambda[n])
      logA <- new_post - old_post
      A <- exp(logA)
      
      if (runif(1) < A){
        beta[,n] <- new_beta
        count[m] <- count[m] + 1 
      }
      else {
        beta[,n] <- cur_beta
      }
      cur_beta <- beta[,n]
    }
  }
  
  #Getting rid of burnin and thinning
  
  lambda <- lambda[seq(burnin+1, total_its, thin)]
  sigma <- sigma[seq(burnin+1, total_its, thin)]
  beta0 <- beta0[seq(burnin+1, total_its, thin)]
  beta <- beta[,seq(burnin+1, total_its, thin)]
  
  #Finding acceptance rate - over ALL iterations (including burnin?)
  acc_rate <- count/total_its
  
  #Finding means and estimate of credibility intervals
  mu_beta <- apply(beta, 1, mean)
  CI <- matrix(rep(0,2*p), nrow = p)
  for (m in 1:p){
    CI[m,] <- quantile(beta[m,], probs = c(0.025, 0.975))
  }
  
  return(list(beta = beta, beta0 = beta0, sigma = sigma, lambda = lambda, mu_beta = mu_beta, CI = CI, acc_rate = acc_rate))
}

result <- MH(Y,X,n_iters,burnin,thin, sd)
#Look at beta 1 again
betaMH <- result$beta
beta1MH <- betaMH[1,]
beta2MH <- betaMH[2,]
beta3MH <- betaMH[3,]
plot(beta1MH)
plot(beta2MH)
plot(beta3MH)
hist(beta1MH)
hist(beta2MH)
hist(beta3MH)

result$mu_beta
result$CI

#Test for convergence for beta1
sample1 <- head(beta1MH, 100)
sample2 <- tail(beta1MH, 100)
t.test(sample1, sample2) #converged in this case - more thinning needed

#Test for convergence for beta2
sample1 <- head(beta2MH, 100)
sample2 <- tail(beta2MH, 100)
t.test(sample1, sample2) #not sig - seems to work!

#Beta 3
sample1 <- head(beta3MH, 100)
sample2 <- tail(beta3MH, 100)
t.test(sample1, sample2) #not sig - seems to work!


#Comparing parameter estimates again
mle$coefficients #mle
breg$mu.beta #bayesreg - uses mean
result$mu_beta #my algorithm - uses mean
#Mine shrinking more than MLE, and sometimes more than bayesreg - sometimes not

#Comparing distribution of beta1 for bayesreg and my algorithm
library(ggplot2)
library(gridExtra)
pl <- ggplot(data = as.data.frame(beta1)) + geom_histogram(aes(x=beta1),binwidth = 0.2) + coord_cartesian(xlim = c(1.3,3.8))
pl2 <- ggplot(data = as.data.frame(beta1MH)) + geom_histogram(aes(x=beta1MH),binwidth = 0.2) + coord_cartesian(xlim = c(1.3,3.8))
grid.arrange(pl, pl2)


#So needed more thinning in this case but not more its - still working fairly well






