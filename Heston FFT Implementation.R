
## For furhter details on implementation of FFT see Chapter 2 of Computational Methods in Finance by Ali Hirsa.
## Characteristic function of Heston Stochastic Model can in found in same book.

lambda = function(u, parameters)
{
  i      = complex(real = 0, imaginary = 1)
  sigma  = parameters[1]
  kaappa = parameters[3]
  rho    = parameters[4]
  
  out_pu = sqrt(sigma^2 * ( u^2 + i*u) + (kaappa - i*rho*sigma*u)^2)
  return(out_pu)
}


omega = function( u, parameters )
{
  ii     = complex(real = 0, imaginary = 1)
  sigma  = parameters[1]
  kaappa = parameters[3]
  rho    = parameters[4]
  theta  = parameters[5]
  S0     = parameters[6]
  r      = parameters[7]
  t      = parameters[8]
  q      = parameters[9]
  
  num = exp(ii * u * log(S0) + ii*u *(r - q)*t + (kaappa * theta * t * (kaappa - ii*rho*sigma*u ))/(sigma^2))
  
  denom = (cosh(lambda(u, parameters) * (t/2)) + ((kaappa - ii*rho*sigma*u )/lambda(u, parameters)) * sinh(lambda(u, parameters) * (t/2))) ^ ((2 * kaappa * theta)/(sigma^2))
  
  return (num/denom)
}


psi = function( u, parameters )
{
  iii = complex(real = 0, imaginary = 1)
  sigma  = parameters[1]
  nhu_0  = parameters[2]
  kaappa = parameters[3]
  rho    = parameters[4]
  t      = parameters[8]
  
  num1   = -(u^2 + iii*u) * nhu_0
  
  denom1 = lambda(u, parameters) * (1/(tanh(lambda(u, parameters) *(t/2)))) + (kaappa - iii*rho*sigma*u )
  
  val1   = omega( u, parameters ) * exp(num1/denom1)
  
  return(val1)
}



phi = function(u, alpha, parameters)
{
  r  = parameters[7]
  t  = parameters[8]
  i1 = complex(real = 0, imaginary = 1)
  
  C  = exp(-r*t)
  
  phi_num = C * psi((u - (alpha + 1)*i1), parameters)
  phi_denom = (alpha + u*i1) * (alpha + u*i1 + 1)
  
  return (phi_num/phi_denom)
  
}



indicator = function(N) ifelse( N == 1, 1, 0)

# Problem a.i.:
# heston_call_price(1, 1000, 100, 100, parameters) Answer: 11.40794
heston_call_price = function(alpha, N, B, K, parameters)
{
  km    = log(K)
  d_nhu = B/N
  j     = 1:(N+1)
  vj    = (j - 1) * d_nhu
  i2    = complex(real = 0, imaginary = 1)
  c_p   = rep(0, N)
  for (x in 1:N)
  {
    part_1 = exp(-vj[x] * km * i2) * phi(vj[x], alpha, parameters)
    
    part_2 = exp(-vj[x+1] * km * i2) * phi(vj[x+1], alpha, parameters)
    
    c_p[x] = 0.5 * (part_1 + part_2) * d_nhu

  }
  
  call_p = ((exp(-alpha * km))/pi)*sum(Re(c_p), na.rm = TRUE)
  return(call_p)
}

test_alphas = c(0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2)
parameters = c(0.25, 0.06, 0.8, -0.25, 0.09, 100, 2.5/100, 1, 0)
x = heston_call_price(test_alphas, 1000, 100, 100, parameters)
plot(test_alphas, x, xlab = "Alphas", ylab = " Euro Call Price", type = "l", main = "Euro Call Price with Diff. Alphas")





heston_fft_call_price = function ( alpha, n, K, B, parameters)
{
  S0     = parameters[6]
  r      = parameters[7]
  t      = parameters[8]
  q      = parameters[9]
  
  N           = 2^n
  d_nhu       = B/N
  d_nhu_kappa = 2 * pi/N
  d_kappa     = d_nhu_kappa/d_nhu
  
  j    = 1:N
  vj   = (j - 1) * d_nhu
  m    = 1:N
  beta = log(S0) - (d_kappa * N)/2
  k_m  = beta + (m - 1) * d_kappa
  
  iiii = complex(real = 0, imaginary = 1)
  
  psi_vals = rep(0, N)
  for (xx in 1:N)
  {
    u        = vj[xx] - (alpha + 1)*iiii
    num2     = psi(u, parameters)
    denom2   = (alpha + iiii * vj[xx] ) * (1+ alpha + iiii*vj[xx])
    psi_vals[xx] = (num2/denom2)
  }
  
  wght = (d_nhu/3) * (3 + (-1)^j - indicator(1:N))
  
  
  x = psi_vals * wght * exp(-iiii * beta * vj) 
  
  y = fft(x)
  
  call_price = (exp(-alpha * k_m)/pi) * Re(y)
  
  mid_point = N/2
  
  
  return(call_price[mid_point + 1])
  
}

