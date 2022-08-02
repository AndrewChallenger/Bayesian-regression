# Example 4: Logistic regression on spambase
data(spambase)

# bayesreg expects binary targets to be factors
spambase$is.spam <- factor(spambase$is.spam)

# First take a subset of the data (1/10th) for training, reserve the rest for testing
spambase.tr  = spambase[seq(1,nrow(spambase),10),]
spambase.tst = spambase[-seq(1,nrow(spambase),10),]

# Fit a model using logistic horseshoe for 2,000 samples
rv <- bayesreg(is.spam ~ ., spambase.tr, model = "logistic", prior = "horseshoe", n.samples = 2e3)

# Summarise
summary(rv)

# Save summary:
rv.s <- summary(rv)

#Figure out the 2-star variables:
which(rv.s$n.stars == 2)  # This gives the indices of the variables with 2 stars
rownames(rv.s$n.stars)[which(rv.s$n.stars == 2)]  # This gives the  names of the variables with 2 stars

#Figure out the 1-or-more-star variables:
which(rv.s$n.stars >= 1)  # This gives the indices of the variables with 1+ stars
rownames(rv.s$n.stars)[which(rv.s$n.stars >= 1)]  # This gives the names of the variables with 1+ stars

