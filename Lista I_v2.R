library(mfx)
library(tidyverse)
library(lmtest)
library(ggplot2)
library(matlib)

#-----------------------------------------------------------------------------------------

# Seed
set.seed(20201)

# Simulate independent variables
X <- data.frame(
  # Normal with 4 different groups of means and std.
  income = rnorm(n = 1000, mean = c(20, 30, 40, 50), sd = c(6, 6, 9, 9)), 
  vote_share = runif(n = 1000, min = 0, max = 1), # U[0, 1]
  age = rpois(n = 1000, lambda = 50), # Poisson
  race = rep(letters[1:4], times = 250) # Factor variable with groups having different incomes
)


# True parameters
true_beta = c(0.05, 0.2, 0.04, 0.06, 0.08, 0.1)
true_delta = 0.0

# We need to turn the groups into dummies. 
X_dummy <- fastDummies::dummy_cols(X)


# Heteroskedastic Normal errors
U <- rnorm(n = 1000, mean = 0, sd = 0.0002 * X_dummy$income^2) 

# This is X * beta
X_beta <- X_dummy$income * true_beta[1] + X_dummy$vote_share * true_beta[2] + 
  X_dummy$age * true_beta[3] + X_dummy$race_b * true_beta[4] + 
  X_dummy$race_c * true_beta[5] + X_dummy$race_d * true_beta[6]

# The true equation has delta 0.
wage = exp(X_beta + true_delta * X_beta^2 + U)
df <- cbind(X_dummy, wage)
head(df)

#-----------------------------------------------------------------------

# The true equation has delta 0.
wage = exp(X_beta + true_delta * X_beta^2 + U)
df <- cbind(X_dummy, wage)
head(df)

nls_results <- nls(
  wage ~ exp(b1 * income + b2 * vote_share + 
               b3 * age + b4 * race_b + b5 * race_c + b6 * race_d),
  data=df,
  start=list(b1 = 0.01, b2 = 0.01, b3 = 0.01, 
             b4 = 0.01, b5 = 0.01, b6 = 0.01)
)

coef_nls1 <- coef(nls_results)

#------------------------------------------------------------------------

# Define a função

min.nls <- function(data, par) {
  with(df, sum((wage - exp(par[1] * income + par[2] * vote_share + 
                           par[3] * age + par[4] * race_b + par[5] * race_c + par[6] * race_d))^2))
}

# Minimiza a função

result <- optim(par = rep(0.001, 6), fn = min.nls, data = df, hessian = TRUE)


# Redefine matrizes

Y <- df$wage

X <- X_dummy %>% select(-race, - race_a) %>% as.matrix()

beta <- result$par

#------------------------------------------------------------------------

# Calcula o score e a hessiana

score <- Y - exp(X %*% beta) %*% t(exp(X %*% coef_nls1)) %*% X

A <- result$hessian

B <- t(score) %*% score
      
#------------------------------------------------------------------------

# Variância assintótica

Avar <- (inv(A) %*% B %*% inv(A))/1000
  


