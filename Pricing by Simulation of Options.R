## calculating price of European Call options using Heston Model


library ('MASS')
params = c(1.18, 0.034, 3.52, -0.77, 0.052,265,0.015,1,0.0177)

euroCallSimulation = function (N, T, M, S0, K, rho)
{
  sigma = params[1]
  kappa = params[3]
  theta = params[5]
  rate  = params[7]
  q     = params[9]
  dt    = T/M
  covar = matrix(c(1, rho, rho,1), ncol = 2)
  v     = c()
  S     = c()
  S[1]  = S0
  v[1]  = params[2]
  
  price = c()
  
  
  for (i in 1:N)
  {
    rnd_v = mvrnorm(n = M, mu = c(0,0), Sigma = covar)
    W_1   = rnd_v[,1]
    W_2   = rnd_v[,2]
    dW_1  = sqrt(dt) * W_1 
    dW_2  = sqrt(dt) * W_2
    
    for (j in 2:M)
    {
      v[j] = v[j-1] +  kappa *(theta - v[j-1])*dt + sigma * sqrt(max(v[j-1],0)) * dW_2[j-1]
      S[j] = S[j-1] + (rate - q) * S[j-1] * dt + sigma * sqrt(max(v[j-1],0))  * S[j-1] *dW_1[j-1]
    }
    price[i] = max(S[M] - K, 0) * exp(-rate * T)
  }
  
  return(sum(price, na.rm = TRUE)/N)
  
}

euroCallSimulation(1000, 1, 252, 265, 285, params[4])
euroCallSimulation(10000, 1, 252, 265, 285, params[4])




## calculating price of Up and Out Call options using Heston Model


upAndOutCallSimulation = function (N, T, M, S0, K1, K2, rho)
{
  sigma = params[1]
  kappa = params[3]
  theta = params[5]
  rate  = params[7]
  q     = params[9]
  dt    = T/M
  covar = matrix(c(1, rho, rho,1), ncol = 2)
  v     = c()
  S     = c()
  S[1]  = S0
  v[1]  = params[2]
  
  price = c()
  
  
  for (i in 1:N)
  {
    rnd_v = mvrnorm(n = M, mu = c(0,0), Sigma = covar)
    W_1   = rnd_v[,1]
    W_2   = rnd_v[,2]
    dW_1  = sqrt(dt) * W_1 
    dW_2  = sqrt(dt) * W_2
    
    for (j in 2:M)
    {
      v[j] = v[j-1] +  kappa *(theta - v[j-1])*dt + sigma * sqrt(max(v[j-1],0)) * dW_2[j-1]
      S[j] = S[j-1] + (rate - q) * S[j-1] * dt + sigma * sqrt(max(v[j-1],0))  * S[j-1] *dW_1[j-1]
    }
    price[i] = ifelse(S[M] < K2, max(S[M] - K1, 0), 0) * exp(-rate * T)
  }
  
  return(sum(price, na.rm = TRUE)/N)
  
}

upAndOutCallSimulation(1000, 1, 252, 265, 285, 315, params[4]) # Answer: 3.8825
upAndOutCallSimulation(2000, 1, 252, 265, 285, 315, params[4]) # Answer: 3.998377
upAndOutCallSimulation(3000, 1, 252, 265, 285, 315, params[4]) # Answer: 3.713694
upAndOutCallSimulation(5000, 1, 252, 265, 285, 315, params[4]) # Answer: 3.733655



## calculating price of Up and Out Call options using Heston Model;
## control variate: European Call

upAndOutCallSimulation_Control = function (N, T, M, S0, K1, K2, rho)
{
  sigma = params[1]
  kappa = params[3]
  theta = params[5]
  rate  = params[7]
  q     = params[9]
  dt    = T/M
  covar = matrix(c(1, rho, rho,1), ncol = 2)
  v     = c()
  S     = c()
  S[1]  = S0
  v[1]  = params[2]
  
  price_euro = c()
  price_upAndOut = c()
  
  
  for (i in 1:N)
  {
    rnd_v = mvrnorm(n = M, mu = c(0,0), Sigma = covar)
    W_1   = rnd_v[,1]
    W_2   = rnd_v[,2]
    dW_1  = sqrt(dt) * W_1 
    dW_2  = sqrt(dt) * W_2
    
    for (j in 2:M)
    {
      v[j] = v[j-1] +  kappa *(theta - v[j-1])*dt + sigma * sqrt(max(v[j-1],0)) * dW_2[j-1]
      S[j] = S[j-1] + (rate - q) * S[j-1] * dt + sigma * sqrt(max(v[j-1],0))  * S[j-1] *dW_1[j-1]
    }
    price_upAndOut[i] = ifelse(S[M] < K2, max(S[M] - K1, 0), 0) * exp(-rate * T)
    price_euro[i] = max(S[M] - K1, 0) * exp(-rate * T)
  }
  
  exp_euro     = sum(price_euro, na.rm = TRUE)/N
  exp_upAndOut = sum(price_upAndOut, na.rm = TRUE)/N
  
  var_upAndOut = var(price_upAndOut)
  
  covar        = cov(price_upAndOut, price_euro)
  
  c            = -covar/var_upAndOut
  
  return (sum(exp_upAndOut + c * (price_euro-exp_euro)) / N)
}

upAndOutCallSimulation_Control(1000, 1, 252, 265, 285, 315, params[4]) # Answer: 3.74298
upAndOutCallSimulation_Control(2000, 1, 252, 265, 285, 315, params[4]) # Answer: 3.757449
upAndOutCallSimulation_Control(3000, 1, 252, 265, 285, 315, params[4]) # Answer: 3.710746
upAndOutCallSimulation_Control(5000, 1, 252, 265, 285, 315, params[4]) # Answer: 3.539985

