
# 1(a)

set.seed(243)
x = runif(10000)
mean(x^2)

# 1(b)

set.seed(243)
x = runif(10000)*4-2
y = runif(10000)
4*mean(x^2 * cos(x*y))

# 1(c)

set.seed(243)
x = rweibull(n = 10000, shape = 3, scale = 4^(1/3))
mean(x^2)

# 2

pnorm(2) - pnorm(1) # standard normal phi(2) - phi(1)

x = runif(100000)+1
mean(exp(-x^2/2)/sqrt(2*pi))

h = function(x, nu) {
  return((1/sqrt(2*pi*(nu^2)) * exp(-(x-1.5)^2/(2*nu^2))))
}

f = function(x) {
  return(((1)/(sqrt(2*pi))) * exp(-x^2/2))
}

g = function(x, nu) {
  return(((1)/(sqrt(2*pi*nu^2))) * exp(-(x - 1.5)^2/(2*nu^2)))
}

w = function(x, nu) {
  return(f(x)/g(x,nu = nu))  
}

set.seed(243)
# try to compare the 3 histograms using ggplot
nu = 0.1
x1 = rnorm(100000, 1.5, nu)
y1 = h(x1, nu)*w(x1, nu)
mean(y1)
hist(y1, main = "nu = 0.1")
plot(x1, y1)

nu = 1
x2 = rnorm(100000, 1.5, nu)
y2 = h(x2, nu)*w(x2, nu)
mean(y2)
hist(y2, main = "nu = 1")
plot(x2, y2)

nu = 10
x3 = rnorm(100000, 1.5, nu)
y3 = h(x3, nu)*w(x3, nu)
mean(y3)
hist(y3, main = "nu = 10")
plot(x3, y3)

# as the standard deviation grows larger, the distribution becomes more off-centered, so it's difficult
# to obtain the correct estimation using importance sampling
# looking for a g that has a higher probability of covering that point within the integral
# as the std. deviation grows larger, the prob. of covering that point becomes lower

# 3(a)

h2 = function(x) {
  return(1/(1 + x))
}

u = function(n) {
  return(runif(n))
}

IhatMC = function(n) {
  return(mean(h2(u(n))))
}

IhatMC(1500) # 0.6947216
log(2) # 0.6931472

# 3(b)

# Analytically calculate
# E[C(U)] = 1 + E(U) = 2

C = function(x) {
  return(1 + x)
}

bhat = function(u) {
  return((sum((h2(u) - mean(h2(u))) * (C(u) - mean(C(u))))) / sum((C(u) - mean(C(u)))^2))
}

IhatCV = function(n) {
  U = u(n)
  b = bhat(U)
  return(mean(h2(U)) - b*(mean(C(U)) - 1.5))
}

IhatCV(n = 1500)

u1 = u(1500)

IhatCV2 = function(u) {
  b = bhat(u)
  return(mean(h2(u)) - b*(mean(C(u)) - 1.5))
}

IhatMC2 = function(u) {
  return(mean(h2(u)))
}

IhatCV2(u1)
IhatMC2(u1)

set.seed(243)
cv = rep(0,10000)
mc = rep(0,10000)
for (i in 1:10000) {
  u = runif(1500)
  cv[i] = IhatCV2(u)
  mc[i] = IhatMC2(u)
}

var(cv)
var(mc)

# question 5

set.seed(243)

# 5(a)

Y = rpois(100, 2)
R = rbinom(100, 1, 0.3)
X = R*Y

# 5(b)

a = 1
b = 1

p0 = runif(1)
lambda0 = rgamma(1, shape = a, scale = 1/b)
r0 = rbinom(100, 1, p0)
xsum = sum(X)
n = length(X)

B = 1000

p = rep(0, B)
lambda = rep(0, B)
r = r0

for (i in 1:B){
  rsum = sum(r)
  lambda[i] = rgamma(1,shape = a+xsum, scale = 1/(b+rsum))
  p[i] = rbeta(1, shape1 = 1+rsum, shape2 = n+1-rsum)
  for (j in 1:length(X)){
    if (X[j]==0){
      r[j] = rbinom(1, 1, p[i]*exp(-lambda[i])/(p[i]*exp(-lambda[i])+1-p[i]))
    } else {
      r[j] = rbinom(1, 1, p[i]*exp(-lambda[i])/(p[i]*exp(-lambda[i])))
    }
  }
}

quantile(p, probs = c(0.025, 0.975)) # the true is p = 0.3
#      2.5%     97.5% 
# 0.1408957 0.3361286

quantile(lambda, probs = c(0.025, 0.975)) # the true is lambda = 2
#     2.5%    97.5% 
# 1.403276 2.870141

## Now, make it as function

gibbs = function(a, b, X, seed, B = 1000){
  
  p0 = runif(1)
  lambda0 = rgamma(1, shape = a, scale = 1/b)
  r0 = rbinom(100, 1, p0)
  xsum = sum(X)
  n = length(X)
  
  p = rep(0, B)
  lambda = rep(0, B)
  r = r0
  
  for (i in 1:B){
    rsum = sum(r)
    lambda[i] = rgamma(1,shape = a+xsum, scale = 1/(b+rsum))
    p[i] = rbeta(1, shape1 = 1+rsum, shape2 = n+1-rsum)
    for (j in 1:length(X)){
      if (X[j]==0){
        r[j] = rbinom(1, 1, p[i]*exp(-lambda[i])/(p[i]*exp(-lambda[i])+1-p[i]))
      } else {
        r[j] = rbinom(1, 1, p[i]*exp(-lambda[i])/(p[i]*exp(-lambda[i])))
      }
    }
  }
  
  cat("95% Bayesian confidence intervals for p is \n")
  cat(quantile(p, probs = c(0.025, 0.975))) 
  cat("\n\n")
  cat("95% Bayesian confidence intervals for lambda is \n")
  cat(quantile(lambda, probs = c(0.025, 0.975)))
  
  out = list(p=p, lambda=lambda)
  
  return(out)
  
}

# and let's try different a and b 

out1 = gibbs(a = 1, b = 10, X = X, seed=243, B=1000); p=out1$p; lambda=out1$lambda
# 95% Bayesian confidence intervals for p is 0.1699037 0.4199715
# 95% Bayesian confidence intervals for lambda is 0.8635715 1.823345

out2 = gibbs(a = 10, b = 1, X = X, seed=243, B=1000); p=out2$p; lambda=out2$lambda
# 95% Bayesian confidence intervals for p is 0.1325896 0.3056973
# 95% Bayesian confidence intervals for lambda is 1.925456 3.45623

out3 = gibbs(a = 1/2, b = 2, X = X, seed=243, B=1000); p=out3$p; lambda=out3$lambda
# 95% Bayesian confidence intervals for p is 0.1453662 0.3438351
# 95% Bayesian confidence intervals for lambda is 1.308008 2.703497

out4 = gibbs(a = 2, b = 1/2, X = X, seed=243, B=1000); p=out4$p; lambda=out4$lambda
# 95% Bayesian confidence intervals for p is 0.1411438 0.3306524
# 95% Bayesian confidence intervals for lambda is 1.372459 2.968282
