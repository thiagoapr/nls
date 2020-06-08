library(mfx)
library(tidyverse)
library(lmtest)
library(ggplot2)
library(matlib)

#-----------------------------------------------------------------------------------------

rm(list=ls())

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
#------------------------------------------------------------------------

## Estima a função não linear

# Objective function to be minimized
modelo <- function(b) {
  m <- exp(b[1]*df$income + 
           b[2]*df$vote_share + 
           b[3]*df$age + 
           b[4]*df$race_b +
           b[5]*df$race_c + 
           b[6]*df$race_d)
  e <- wage - m
  return(sum(e^2))
}

# Gradient of function sse
grad<-function(b) {
  m <- exp(b[1]*df$income + 
             b[2]*df$vote_share + 
             b[3]*df$age + 
             b[4]*df$race_b +
             b[5]*df$race_c + 
             b[6]*df$race_d)
  e <- wage - m
  c(sum(-2*e*df$income*m),sum(-2*e*df$vote_share*m),sum(-2*e*df$age*m),sum(-2*e*df$race_b*m),
    sum(-2*e*df$race_c*m),sum(-2*e*df$race_d*m))
}


#Optimization process
M0 <- optim(c(0.01,0.01,0.01,0.01,0.01,0.01), modelo, grad, method = "BFGS")
#Checking convergence
M0$convergence
#Parameter values
par <- M0$par
M0$value

# Avar

#Hessian matrix
h<-optimHess(M0$par,modelo, grad)
h

sigma<-M0$value/(1000-6)
varcov<-sigma*inv(0.5*h)
varcov
d <- varcov
d <- diag(d)^(1/2)
d

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

## Calculo do gradiente e da hessiana (sem pré-programação)

## Atenção: para o Wooldridge é um vetor coluna 

# Ajusta base de dados

x <- df %>% select(-race, -race_a, -wage) 
y <- c(df$wage)
beta <- as.numeric(M0$par)


# Cálculo da hessiana

H <- matrix(0, nrow = 6, ncol = 6)

for(i in 1:1000) {
  
  h <- as.numeric(exp(2 %*% as.numeric(x[i,]) %*% beta)) * (as.numeric(x[i,]) %*% t(as.numeric(x[i,])))
  
  H <- H + h
  

}



