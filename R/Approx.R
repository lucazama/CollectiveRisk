#'@import actuar

#' @title Aggregate distribution approximation based on moments
#' @description Function to approximate the Cumulative Distribution Function (CDF) of the aggregate claim size distribution.
#' @details The function takes as input the first raw moments of the aggregate claim size distribution S and returns
#' the desired approximation to P(S>m) based on moments. The available approximations are c("Normal","Gamma", "NormalPower", "NormalToGamma","GramCharlier", "Bowers").
#' For GramCharlier and Bowers, the approximations are available up to the 5th order (n=c(3,4,5))
#' @examples MomentsApprox(900,moments=AggregateRawMoments(Moments("Poisson",30,"raw"),Moments("Gamma",c(40,2),"raw")), "Bowers",n=3)
#' MomentsApprox(900,moments=AggregateRawMoments(Moments("Poisson",30,"raw"),Moments("Gamma",c(40,2),"raw")), "NormalToGamma")
#' @usage MomentsApprox(m,moments,type,n)
#'@param m probability to approximate:  P(S>m)
#'@param moments vector of raw moments of the aggregate distribution S. Besides the "Normal" approximation, all other
#'@param approximations require at least the first 3 raw moments of S. "GramCharlier" and "Bowers" require as many raw moments
#'as their chosen approximation order.
#'@param type the available approximation methods are c("Normal","Gamma", "NormalPower", "NormalToGamma","GramCharlier", "Bowers")

#'@param n approximation order, only needed for GramCharlier and Bowers c(3,4,5).
#'@export
MomentsApprox <- function(m,moments,type,n){
  ES <- moments[1]
  VarS <- rawtocentral(moments)[2]
  SkewS <- rawtostandard(moments)[3]
  KurtS <-  rawtostandard(moments)[4]
  z <-  (m-ES)/sqrt(VarS)
  alpha <- 4/SkewS^2
  beta <- 2/(SkewS*sqrt(VarS))
  x0 <-  ES- 2*sqrt(VarS)/SkewS

  if (type=="Normal"){
    return(1- pnorm(m,ES,sqrt(VarS)))
  }
  if (type=="Gamma"){
    return(1-pgamma(m-x0,alpha,beta))
  }
  if (type=="NormalToGamma"){
    return(1- pnorm(sqrt(16/SkewS^2 + 8*z/SkewS)- sqrt(16/SkewS^2 -1) ))
  }
  if (type=="NormalPower"){
    return(1- pnorm(sqrt(9/SkewS^2 + 6*z/SkewS+1)- 3/SkewS))
  }
  if (type=="GramCharlier"){
    if (n==3){
      val=integrate(function(x) dnorm(x)*(x^3-3*x)*SkewS/6, lower=-Inf,upper=(m-ES)/sqrt(VarS))
      return(1- (val$value + pnorm((m-ES)/sqrt(VarS))))
    }
    if (n==4){
      val=integrate(function(x) dnorm(x)*((x^3-3*x)*SkewS/6 + (x^4-6*x^2+3)* (KurtS -3)/24 ), lower=-Inf,upper=(m-ES)/sqrt(VarS))
      return(1- (val$value + pnorm((m-ES)/sqrt(VarS))))
    }
    if (n==5){
      val=integrate(function(x) dnorm(x)*((x^3-3*x)*SkewS/6 + (x^4-6*x^2+3)* (KurtS -3)/24  + (x^5-10*x^3+15*x) *(rawtostandard(moments)[5]-10*SkewS )/120 ), lower=-Inf,upper=(m-ES)/sqrt(VarS))
      return(1- (val$value + pnorm((m-ES)/sqrt(VarS))))
    }
    else return("Wrong input. The approximation orders supported are n= c(3,4,5)")
  }
  if (type=="Bowers"){
    a <- ES^2/VarS
    EZ <-  a
    EZ2 <- moments[2]* ES^2/VarS^2
    EZ3 <- ES^3*moments[3]/VarS^3
    EZ4 <- ES^4*moments[4]/VarS^4
    EZ5 <- ES^5*moments[5]/VarS^5
    mu3 <- rawtocentral(c(EZ,EZ2,EZ3,EZ4,EZ5))[3]
    mu4 <- rawtocentral(c(EZ,EZ2,EZ3,EZ4,EZ5))[4]
    mu5 <- rawtocentral(c(EZ,EZ2,EZ3,EZ4,EZ5))[5]
    A3 <- gamma(a)* (mu3-2*a)/(6*gamma(a+3))
    A4 <- gamma(a)* (mu4-12*mu3-3*a^2+18*a)/(24*gamma(a+4))
    A5 <- gamma(a)* (mu5 - 20*mu4 - (10*a -120)*mu3 +60*a^2-144*a)/(120*gamma(a+4))
    L3 <- function(x) x^3-3*(a+2)*x^2 +3*(a+2)*(a+1)*x - (a+2)*(a+1)*a
    L4 <- function(x) x^4-4*(a+3)*x^3 + 6*(a+3)*(a+2)*x^2 -4*(a+3)*(a+2)*(a+1)*x +(a+3)*(a+2)*(a+1)*a
    L5 <- function(x) x^5 - 5*(a+4)*x^4 +10*(a+4)*(a+3)*x^3 -10*(a+4)*(a+3)*(a+2)*x^2 + 5*(a+4)*(a+3)*(a+2)*(a+1)*x -(a+4)*(a+3)*(a+2)*(a+1)*a

    if (n==3){
      val=integrate(function(x) dgamma(x,a,1)* (1+A3*L3(x)), lower=0, upper=m*ES/VarS)
      return(1-val$value)
    }
    if (n==4){
      val=integrate(function(x) dgamma(x,a,1)* (1+A3*L3(x) +A4*L4(x)), lower=0, upper=m*ES/VarS)
      return(1-val$value)
    }
    if (n==5){
      val=integrate(function(x) dgamma(x,a,1)* (1+A3*L3(x) +A4*L4(x) +A5*L5(x)), lower=0, upper=m*ES/VarS)
      return(1-val$value)
    }
    else return("Wrong input. The approximation orders supported are n= c(3,4,5)")
  }

}

#'@title Panjer approximation
#'@description Computes the Panjer approximation to P(S>m) where S is the aggregate claim distribution resulting
#'from a claim frequency distribution in the (a,b,0) class and any continuous claim intensity distribution.
#'@usage panjer(fx, dist, par, step, m)

#'@param fx continuous Probability Density Function (PDF) for the claim intensities.
#'@param dist c("Poisson","NegBin", "Binomial")
#'@param par vector of parameters of the claim frequency distributions (as ordered as in the built in PDF function)
#'@param step width of the discretization (for a computing time lower than a couple of minutes choose step such that m/step <20000)
#'@param m probability to approximate P(S>m)
#'@examples panjer(fx=function(x) pgamma(x,40,2), dist="Poisson", par=30, step=0.1, 900)
#'@return the function returns the upper and lower bound probabilities. The thinner the grid (lower the step) the closer the two bounds will get.
#'@export
panjer= function(fx, dist, par, step, m) {
  M <- m/step
  Hu <-  discretize(fx, from=0,to=m+step ,step= step, method="upper")
  Hl <-  discretize(fx, from=0,to=m  ,step= step, method="lower")
  gu <- numeric(M)
  gl <- numeric(M)

  if (dist=="Poisson"){
    lambda <-  par
    a <- 0
    b <- lambda
    gu[1] <- exp(-lambda*(1-Hu[1])) #in the upper discretization we have a non zero probability mass in 0 (check notes)
    gl[1] <- exp(-lambda*(1-Hl[1])) #this is equivalent to dpois(0,lambda)

    #upper
    for (i in 1:M){
      steparray=numeric(i)
      for (k in 1:i){
        steparray[k]= k*Hu[k+1]* gu[i-k+1]
      }
      gu[i+1]=(b/i) * sum(steparray[1:i])
    }

    # lower
    for (i in 1:M){
      steparray=numeric(i)
      for (k in 1:i){
        steparray[k]= k*Hl[k+1]* gl[i-k+1]
      }
      gl[i+1]=(b/i) * sum(steparray[1:i])
    }
    return(c(1-sum(gu), 1-sum(gl)))
  }

  if (dist=="NegBin"){
    r <-  par[1]
    p <-  par[2]
    a <- 1-p
    b <- (1-p)*(r-1)
    gu[1] <- ((p)/(1-(1-p)*Hu[1]))^r #in the upper discretization we have a non zero probability mass in 0 (check notes)
    gl[1] <- ((p)/(1-(1-p)*Hl[1]))^r #this is equivalent to dnbinom(0,r,p)

    #upper

    for (i in 1:M){
      steparray=numeric(i)
      for (k in 1:i){
        steparray[k]= (a+b*k/i)*Hu[k+1]* gu[i-k+1]
      }
      gu[i+1]= 1/(1-a*Hu[1])*sum(steparray[1:i])
    }

    # lower
    for (i in 1:M){
      steparray=numeric(i)
      for (k in 1:i){
        steparray[k]= (a+ b*k/i)*Hl[k+1]* gl[i-k+1]
      }
      gl[i+1]= sum(steparray[1:i])
    }

    return(c( 1-sum(gu),1-sum(gl)))
  }

  if (dist=="Binomial"){
    n <- par[1]
    p <- par[2]
    a <- -p/(1-p)
    b <- p*(n+1)/(1-p)
    gu[1] <- (p*Hu[1]+ (1-p))^n #in the upper discretization we have a non zero probability mass in 0 (check notes)
    gl[1] <- (p*Hl[1]+ (1-p))^n#this is equivalent to dbinom(0,n,p)

    #upper

    for (i in 1:M){
      steparray=numeric(i)
      for (k in 1:i){
        steparray[k]= (a+b*k/i)*Hu[k+1]* gu[i-k+1]
      }
      gu[i+1]= 1/(1-a*Hu[1])*sum(steparray[1:i])
    }

    # lower
    for (i in 1:M){
      steparray=numeric(i)
      for (k in 1:i){
        steparray[k]= (a+ b*k/i)*Hl[k+1]* gl[i-k+1]
      }
      gl[i+1]= sum(steparray[1:i])
    }
    return(c( 1-sum(gu),1-sum(gl)))
  }
  else
    return("wrong frequency distribution input.")
}



#'@title Fast Fourier Transform
#'@description Computes the Fast Fourier Approximation to P(S>m), given the Characteristic function of S.

#'@usage FFT (phi, n,lower,upper, m )
#'@param phi Characteristic function of S.
#'@param n number of points used in the discretization (possibly a power of 2)
#'@param lower starting point of the discretization of S
#'@param upper end point of the discretization of S
#'@param m approximation to P(S>m)
#'@details The function is useful when the characteristic function of a random variable is known, while the CDF is not.
#'this is the case in the collective risk model with light tailed claim intensities. For the case in which the intensities are heavy tailed,
#'refer to FFT2.
#'@seealso FFT2

#'@examples FFT(phi= function(t) exp(30*((2/ (2- (0+1i)*t))^40 -1 )),n=2^19, lower=0,upper=3000,m=900)
#'@return the function returns P(S>m)
#'@export
FFT <- function(phi, n,lower,upper, m ){
  characteristic_function_to_density <- function(
    phi, # characteristic function;
    n,   # Number of points, ideally a power of 2
    a=lower, b=upper # Evaluate the density on [a,b[
  ) {
    i <- 0:(n-1)            # Indices
    dx <- (b-a)/n           # Step size, for the density
    x <- a + i * dx         # Grid, for the density
    dt <- 2*pi / ( n * dx ) # Step size, frequency space
    c <- -n/2 * dt          # Evaluate the characteristic function on [c,d]
    d <-  n/2 * dt          # (center the interval on zero)
    t <- c + i * dt         # Grid, frequency space
    phi_t <- phi(t)
    X <- exp( -(0+1i) * i * dt * a ) * phi_t
    Y <- fft(X)
    density <- dt / (2*pi) * exp( - (0+1i) * c * x ) * Y
    return(Re(density))
  }
  density <- characteristic_function_to_density( phi, n,lower,upper)
  k= m/upper * n
  return(1 - sum(density[1:floor(k)]) / sum(density))
}



#'@title Fast Fourier Transform
#'@description Computes the Fast Fourier Approximation to P(S>m), where S is the aggregate claim size.

#'@usage FFT2 (fx,fm, n,lower,upper, m )
#'
#'@param fx PDF of the claim intesities.
#'@param fm PDF of the claim frequency.
#'@param n number of points used in the discretization (possibly a power of 2)
#'@param lower starting point of the discretization of X
#'@param upper end point of the discretization of X
#'@param m approximation to P(S>m)
#'
#'@details The function works for any aggregate claim distribution. It is not that stable when the lower and upper bounds are not selected with accuracy.
#' (be sure to consider a great part of the domain of S)

#'@examples FFT2(fx=function(x) pgamma(x,40,2),fM= function(t) exp(30*(exp(t)-1)),n=2^19, lower=0, upper=3000,m=900)
#'@return the function returns P(S>m)
#'@seealso FFT
#'@note the discretization bounds (lower and upper) are generally different from those of FFT. In FFT the discretization
#'regards S, in FFT2 it concerns X.
#'@references "Panjer Recursion versus FFT For Compound Distributions", Paul Embrechets and Marco Frei (2009), Mathematical Methods of Operations Research.
#'@export


FFT2 <- function(fx, fM,n, lower,upper,m){
  step= upper/n
  f <- discretize(fx, from=lower, to=upper, by=step,
                  method="rounding")
  fhat <- fft(f, inverse=FALSE)
  P <- fM(log(fhat))
  g <- 1/n*fft(P, inverse=TRUE)

  return(1-sum(Re(g[1:floor(m/step)])))
}


#'@title Simulator of random numbers from some not so common distributions
#'@description The function simulates random numbers from some distributions
#'@param distN For now distN= c("NegBin","Binomial","Poisson")
#'@param distY For now distY= c("Gamma", "Exponential","LogNormal","Weibul")
#'@param parN the parameters of the chosen counting distribution
#'@param parY the parameters of the chosen claim intensities distribution.
#'@param n the number of simulations.
#'@return The function returns the simulations of S. If you want to compute the probability P(S<m)
#'you can simply run the following command mean(S<m)
#'@example S=simulation(distN= "Poisson", distY="Gama",parN=30,parY=c(3,2),n=10000)
#'mean(S>50)
#'@note Confidence intervals still to be added (along with more distributions, and truncated versions as well)
#'@export
simulator <- function(dist, par,n){
  if (dist=="Pareto"){
    alpha <- par[1]
    theta <- par[2]
    U <- runif(n)
    return(theta/(1-U)^(1/alpha))
  }
  if (dist=="TruncWeibull"){
    lambda <- par[1]
    k <- par[2]
    alpha <- par[3]
    U <- runif(n)
    return( lambda* (alpha/lambda -log(1-U))^(1/k))
  }
  if (dist=="TruncLogNormal"){
    mu <- par[1]
    sigma <- par[2]
    k <- par[3]
    t= (log(k)-mu)/sigma
    if (t>0){
      alfa=(t +sqrt(t^2+4))/2
    } else
      alfa=t

    m=runif(n*1.5)
    z=(-1/alfa)* (log(1-m)-t*alfa)
    p=exp(-(alfa-z)^2/2)
    u=runif(length(p),0,1)
    x=numeric(length(p))
    i=1
    count=0
    while (count<n){
      if (u[i]<p[i]){
        x[i]=z[i]
        count=count+1
      } else
        x[i]=0
      i=i+1
    }
    return(exp(x[x>0]*sigma+mu))
  }
}


# Moments("Pareto",c(3,2),"raw")





#'@title Simulation of Aggregate Claim Size
#'@description The function simulates the aggregate claim size distribution
#'@param distN For now distN= c("NegBin","Binomial","Poisson")
#'@param distY For now distY= c("Gamma", "Exponential","LogNormal","Weibull", "Pareto",
#'"TruncLogNormal","TruncWeibull").
#'@param parN the parameters of the chosen counting distribution
#'@param parY the parameters of the chosen claim intensities distribution.
#'@param n the number of simulations.
#'@return The function returns the simulations of S. If you want to compute the probability P(S<m)
#'you can simply run the following command mean(S<m)
#'@example S=simulation(distN= "Poisson", distY="Gama",parN=30,parY=c(3,2),n=10000)
#'mean(S>50)
#'@note Confidence intervals still to be added (along with more distributions, and truncated versions as well)
#'@export
simulation <- function(distN,distY, parN,parY,n){
  if (distN=="Poisson"){
    lambda=parN
    N=rpois(n,lambda)
  }
  else if (distN=="NegBin"){
   r <- parN[1]
   p <- parN[2]
    N=rnbinom(n,r,p)
  }
  else if (distN=="Binomial"){
    n2 <- parN[1]
    p <- parN[2]
    N=rbinom(n,n2,p)
  }
  else {
    return("wrong counting distribution")
  }
  S=numeric(n)
  if (distY=="Gamma"){
    for (i in 1:n){
      if (N[i]!=0){
        S[i]=sum(rgamma(N[i],parY[1],parY[2]))
      }
      else
        S[i]=0
    }
  }
  else if (distY=="Exponential"){
    for (i in 1:n){
      if (N[i]!=0){
        S[i]=sum(rexp(N[i],parY[1]))
      }
      else
        S[i]=0
    }
  }
  else if (distY=="LogNormal"){
    for (i in 1:n){
      if (N[i]!=0){
        S[i]=sum(rlnorm(N[i],parY[1],parY[2]))
      }
      else
        S[i]=0
    }
  }
  else if (distY=="Weibull"){
    for (i in 1:n){
      if (N[i]!=0){
        S[i]=sum(rweibull(N[i],parY[1],parY[2]))
      }
      else
        S[i]=0
    }
  }
  else if (distY=="Pareto"){
    for (i in 1:n){
      if (N[i]!=0){
        S[i]=sum(simulator("Pareto",par=parY,n=N[i]))
      }
      else
        S[i]=0
    }
  }
  else if (distY=="TruncLogNormal"){
    for (i in 1:n){
      if (N[i]!=0){
        S[i]=sum(simulator("TruncLogNormal",par=parY,n=N[i]))
      }
      else
        S[i]=0
    }
  }
  else if (distY=="TruncWeibull"){
    for (i in 1:n){
      if (N[i]!=0){
        S[i]=sum(simulator("TruncWeibull",par=parY,n=N[i]))
      }
      else
        S[i]=0
    }
  }
  else {
    return("wrong claim instensity distribution")
  }
  return(S)
}

