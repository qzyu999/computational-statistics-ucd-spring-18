# data points
x = c(-13.87, -2.53, -2.44, -2.40, -1.75, -1.34, -1.05, -0.23, -0.07, 
       0.27, 1.77, 2.76, 3.29, 3.47, 3.71, 3.80, 4.24, 4.53, 43.21, 56.75)

# test vector for theta
theta = seq(-100,100,1)

# function for log of Cauchy distribution
mylog = function(x, theta){
  n = length(x)
  return(-n*log(pi) - sum(log(1+(theta-x)^2)))
}

# create vector applying the function to the theta
w = sapply(theta, function(y) mylog(x, y))

# plot
plot(theta, w, main = "log likelihood function")


# 1d MLE using Newton-Raphson method
# create a function for the first log derivative
mylog1 = function(x, theta){
  n = length(x)
  return(-2*sum((theta-x)/(1+(theta-x)^2)))
}

# create a function for the second log derivative
mylog2 = function(x, theta){
  n = length(x)
  return(-2*sum((1-(theta-x)^2)/((1+(theta-x)^2)^2)))
}

# function h

h = function(x, theta){
  return(mylog1(x, theta)/mylog2(x, theta))
}

# create the Newton-Raphson function using the previous functions
nr = function(x, theta_int, cutoff){
  theta_new = theta_int
  theta_old = theta_int + 9999
  while(abs(theta_old-theta_new)>cutoff){
  theta_old = theta_new
  theta_new = theta_old + h(x, theta_old)
  }
  return(theta_new)
}

nr(x, 4, 0.0001)

# second set of data points
x2 = c(-11, -1, 0, 1.4, 4.1, 4.8, 7, 8, 38)

sapply(x2, function(y) nr(x,y,0.0001))

# 1e

# function h2, for fisher scoring method

h2 = function(x, theta){
  n = length(x)
  return(mylog1(x, theta)/(n/2))
}

nr2 = function(x, theta_int, cutoff){
  theta_new = theta_int
  theta_old = theta_int + 9999
  while(abs(theta_old-theta_new)>cutoff){
    theta_old = theta_new
    theta_new = theta_old + h2(x, theta_old)
  }
  return(theta_new)
}

theta_fsc = sapply(x2, function(y) nr2(x,y,0.0001))
theta_finalized = sapply(theta_fsc, function(y) nr(x,y,0.0001))
theta_finalized

#2

x3 = c(0.52, 1.96, 2.22, 2.28, 2.28, 2.46, 2.50, 2.53, 2.54, 2.99, 3.47, 3.53, 3.70, 3.88, 3.91, 4.04, 4.06, 4.82, 4.85, 5.46)

# log likelihood for the pdf
logfnc = function(x, theta){
    n = length(x)
    return(sum(log(1-cos(x-theta)))-n*log(2*pi))
}

theta = seq(-10, 10, 0.01)

lld = sapply(theta, function(y) logfnc(x3,y))
plot(theta, lld, main = "log likelihood")

# 2c
dlpdf = function(x, theta){
  return(-sum(sin(x-theta)/(1-cos(x-theta))))
}

d2lpdf = function(x, theta){
  return(sum(cos(x-theta)/(1-cos(x-theta))^2))
}

# function h2

h2 = function(x, theta){
  return(dlpdf(x, theta)/d2lpdf(x, theta))
}

# create the Newton-Raphson function using the previous functions
nrmle = function(x, theta_int, cutoff){
  theta_new = theta_int
  theta_old = theta_int + 9999
  while(abs(theta_old-theta_new)>cutoff){
    theta_old = theta_new
    theta_new = theta_old + h2(x, theta_old)
  }
  return(theta_new)
}

nrmle(x3, asin(mean(x3)-pi), 0.0001)

# 2d

nrmle(x3, -2.7, 0.0001)
nrmle(x3, 2.7, 0.0001)

# 2e

theta_test = c(-pi, ppoints(198)*2*pi-pi, pi)
theta_select = sapply(theta_test, function(y) nrmle(x3, y, 0.0001))

split(theta_test, round(theta_select, 1))

# 3
x4 = c(0.02, 0.06, 0.11, 0.22, 0.56, 1.10, 0.02, 0.06, 0.11, 0.22, 0.56, 1.10)
y = c(47, 97, 123, 152, 191, 200, 76, 107, 139, 159, 201, 207)
ystar = 1/y
u = 1/x4
fit1 = lm(ystar~u)
fit1
beta0 = 0.0051072
beta1 = 0.0002472
theta1_int = 1/beta0
theta2_int = beta1*theta1_int
theta1_int
theta2_int

# 3b

fnc0 = function(x, y, theta1, theta2){
  return(sum((y-theta1*x/(x+theta2))^2))
}

fnc1 = function(x,y,theta1,theta2){
  return(-2*sum(x*y)+2*sum((theta1*x^2)/(x+theta2)))
}

fnc2 = function(x,y,theta1,theta2){
  return(sum(2*(y-((theta1*x)/(x+theta2)))*((theta1*x)/(x+theta2)^2)))
}

fnc3 = function(x,y,theta1,theta2){
  return(c(fnc1(x,y,theta1,theta2),fnc2(x,y,theta1,theta2)))
}

fnc4 = function(x,theta2){
  return(2*sum(x^2/(x+theta2)))
}

fnc5 = function(x,y,theta1,theta2){
  return(-2*sum((2*y*theta1*x)/(x+theta2)^3)+3*sum((2*theta1^2*x^2)/(x+theta2)^4))
}

fnc6 = function(x,theta1,theta2){
  return(-2*sum(theta1*x^2/(x+theta2)^2))
}

fnc7 = function(x,y,theta1,theta2){
  return(matrix(c(fnc4(x,theta2), fnc6(x,theta1,theta2), fnc6(x,theta1,theta2), fnc5(x,y,theta1,theta2)), ncol = 2))
}

h = function(x, y, theta1, theta2){
  return(solve(fnc7(x, y, theta1, theta2)) %*% fnc3(x, y, theta1, theta2))
}

nr_mult = function(x, y, theta1_int, theta2_int, cutoff){
  theta1_new = theta1_int
  theta2_new = theta2_int
  theta1_old = theta1_int + 9999
  theta2_old = theta2_int + 9999
  while(abs(theta1_old-theta1_new)+abs(theta2_old-theta2_new)>cutoff){
    theta1_old = theta1_new
    theta2_old = theta2_new
    theta_new = c(theta1_old, theta2_old) - h(x, y, theta1_old, theta2_old)
    theta1_new = theta_new[1]
    theta2_new = theta_new[2]
  }
  return(theta_new)
}

nr_mult(x4, y, theta1_int, theta2_int, 0.0001)

# 3c

SA = function(x, y, theta1_int, theta2_int, alpha, cutoff){
  theta1_new = theta1_int
  theta2_new = theta2_int
  theta1_old = theta1_int + 9999
  theta2_old = theta2_int + 9999
  while(abs(theta1_old-theta1_new)>cutoff | abs(theta2_old-theta2_new)>cutoff){
    theta1_old = theta1_new
    theta2_old = theta2_new
    theta_new = c(theta1_old, theta2_old) - alpha*fnc3(x, y, theta1_old, theta2_old)
    theta1_new = theta_new[1]
    theta2_new = theta_new[2]
    Fval_old = fnc0(x, y, theta1_old, theta2_old)
    Fval_new = fnc0(x, y, theta1_new, theta2_new)
    #if (Fval_new<Fval_old) alpha = alpha/2
    alpha = alpha/2
    print(c(theta1_new, theta2_new))
  }
  return(theta_new)
}

SA(x4, y, theta1_int, theta2_int, 0.000001, 0.1)

# 3d

A = function(x, theta1, theta2){
  return(cbind(x/(x+theta2), -1*theta1*x/(x+theta2)^2))
}

z = function(x, y, theta1, theta2){
  return(y-theta1*x/(x+theta2))
}

GN = function(x, y, theta1_int, theta2_int, cutoff){
  theta1_new = theta1_int
  theta2_new = theta2_int
  theta1_old = theta1_int + 9999
  theta2_old = theta2_int + 9999
  while(abs(theta1_old-theta1_new)>cutoff | abs(theta2_old-theta2_new)>cutoff){
    theta1_old = theta1_new
    theta2_old = theta2_new
    theta_new = c(theta1_old, theta2_old) + solve(t(A(x, theta1_old, theta2_old))%*%A(x, theta1_old, theta2_old))%*%t(A(x, theta1_old, theta2_old))%*%z(x, y, theta1_old, theta2_old)
    theta1_new = theta_new[1]
    theta2_new = theta_new[2]
    print(c(theta1_new, theta2_new))
  }
  return(theta_new)
}

GN(x4, y, theta1_int, theta2_int, 0.0001)



