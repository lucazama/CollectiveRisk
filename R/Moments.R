#'@title Function to convert raw Moments to Central ones.
#'@description The function converts the raw moments of a distribution to the central ones.
#' Instead of the first central moment (which is always 0), the function displays the first raw moment instead.
#' @param moments is a vector of raw moments
#' @return A vector of central moments (the first one is the first raw moment.) the nth element of the vector
#' is the nth central moment.
#' @note Maximum length supported is n=5.
#' @export
rawtocentral= function(moments){
  EX <-  moments[1]
  EX2 <-  moments[2]
  EX3 <-  moments[3]
  EX4 <-  moments[4]
  EX5 <- moments[5]

  EC1 <- EX #the central moment would be 0, but using EX simplifies things when going back from central moments to raw ones.
  EC2 <- EX2-EX^2
  EC3 <- EX3 -3*EX2*EX + 2*EX^3
  EC4 <- EX4 - 4*EX3*EX + 6*EX2*EX^2 -3*EX^4
  EC5 <- EX5 -5*EX4*EX +10*EX3*EX^2 -10*EX2*EX^3 +4*EX^5

  if (length(moments)==2){
    return(c(EC1,EC2))
  }
  if (length(moments)==3){
    return(c(EC1,EC2,EC3))
  }
  if (length(moments)==4){
    return(c(EC1,EC2,EC3,EC4))
  }
  if (length(moments)==5){
    return(c(EC1,EC2,EC3,EC4,EC5))
  }
}

#'@title Function to convert Central Moments to raw ones.
#'@description The function converts the central moments of a distribution to the raw ones.
#' @param moments is a vector of central moments. The first element must be the first raw moment and the nth element is the nth central moment.
#' @note Maximum length supported is n=5.
#' @return A vector of raw moments. the nth element of the vector
#' is the nth raw moment.
#' @export

centraltoraw <- function(moments){
  EC1 <- moments[1]
  EC2 <- moments[2]
  EC3 <- moments[3]
  EC4 <- moments[4]
  EC5 <- moments[5]

  EX <- EC1
  EX2 <- EC2+EC1^2
  EX3 <- EC3 +3* EX2*EC1-2*EC1^3
  EX4 <- EC4 + 4*EX3*EC1 - 6*EX2*EC1^2 +3*EC1^4
  EX5 <- EC5 +5*EX4*EX -10*EX3*EX^2 +10*EX2*EX^3 -4*EX^5

  if (length(moments)==2){
    return(c(EX,EX2))
  }
  if (length(moments)==3){
    return(c(EX,EX2,EX3))
  }
  if (length(moments)==4){
    return(c(EX,EX2,EX3,EX4))
  }
  if (length(moments)==5){
    return(c(EX,EX2,EX3,EX4,EX5))
  }
}


#'@title Function to convert raw Moments to Standard ones.
#'@description The function converts the raw moments of a distribution to the standardized ones.
#' @param moments is a vector of raw moments.
#' @return A vector of standardized moments. the nth element of the vector
#' is the nth standardized moment.
#' @note Maximum length supported is n=5.
#' @export
rawtostandard <- function(moments){
  EX  <-  moments[1]
  EX2 <- moments[2]
  EX3 <-  moments[3]
  EX4 <- moments[4]
  EX5 <- moments[5]

  ES1 <- 0
  ES2 <- 1
  ES3 <- (EX3 -3*EX2*EX + 2*EX^3)/(EX2-EX^2)^(3/2)
  ES4 <- (EX4 - 4*EX3*EX + 6*EX2*EX^2 -3*EX^4)/(EX2-EX^2)^(2)
  ES5 <- (EX5 -5*EX4*EX +10*EX3*EX^2 -10*EX2*EX^3 +4*EX^5)/(EX2-EX^2)^(5/2)

  if (length(moments)==3){
    return(c(ES1,ES2,ES3))
  }
  if (length(moments)==4){
    return(c(ES1,ES2,ES3,ES4))
  }
  if (length(moments)==5){
    return(c(ES1,ES2,ES3,ES4,ES5))
  }
}


#'@title Aggregate Moments
#'@description the function returns the raw moments of the aggregate distribution of S in the collective risk model.
#'@param momentsN a vector of raw moments of the claim frequency N.
#'@param momentsY a vector (possibly of the same length of momentsN) of raw moments of the claim intensities Y.
#'@return The function returns the first raw moments of the aggregate distribution of S in the collective risk model.
#'@examples AggregateRawMoments(momentsN=Moments("Poisson",10,"raw"),momentsY=Moments("Gamma",c(3,4),"raw"))
#'@note The raw moments of S are obtained using the probability generating function of S.
#' @export

AggregateRawMoments <- function(momentsN, momentsY){
  EN  <-  momentsN[1]
  EN2 <- momentsN[2]
  EN3 <-  momentsN[3]
  EN4 <- momentsN[4]
  EN5 <- momentsN[5]

  EY  <-  momentsY[1]
  EY2 <- momentsY[2]
  EY3 <-  momentsY[3]
  EY4 <- momentsY[4]
  EY5 <- momentsY[5]

  ES  <- EN*EY
  ES2 <- (EN2-EN)*EY^2 +EN*EY2
  ES3 <-  (EN3-3*EN2+2*EN)*EY^3 + 3*(EN2-EN)*EY2*EY + EN*EY3
  ES4 <- (EN4-6*EN3+11*EN2-6*EN)* EY^4 + 6* ( (EN3-3*EN2 +2*EN)*EY2*EY^2)+3* ( (EN2-EN)* (EY3*EY+ (EY2)^2) ) + (EN2-EN)*(EY3*EY) + EN*EY4
  ES5 <- ((EN5-10*EN4+35*EN3-50*EN2+24*EN)*EY^5 +
            10*((EN4-6*EN3+11*EN2-6*EN)*EY^3*EY2) +
            10*((EN3-3*EN2+2*EN)*EY^2*EY3)+
            15*((EN3-3*EN2+2*EN)*EY*EY2^2)+
            10*((EN2-EN)*EY2*EY3)+
            5*((EN2-EN)*EY*EY4)+
            EN*EY5)
  if (length(momentsN)==3){
    return(c(ES,ES2,ES3))
  }
  if (length(momentsN)==4){
    return(c(ES,ES2,ES3,ES4))
  }
  if (length(momentsN)==5){
    return(c(ES,ES2,ES3,ES4,ES5))
  }
}
#'@title Raw Moments of common distribution functions
#'@description The function computes the first five moments of a distribution among the ones supported.
#'@param dist distribution for which you want to compute the moments: c("Binomial","Negbinomial"
#',"Geometric", "Poisson", "Pareto", "Gamma", "Exponential", "LogNormal")
#'@param par vector of parameters of the chosen distribution.
#'@param type the function can return c("central", "raw", "standardized") moments.
#'@return The  function returns the first 5 moments of the distribution chosen.
#'@note The function returns only the finite moments.
#' @export
Moments <- function(dist,par,type){
  if(dist=="Binomial"){
    n <- par[1]
    p <- par[2]
    E1 <- n*p
    E2 <- n*p*(1-p)
    E3 <-  n*p*(1-p)*(1-2*p)
    E4 <- n*p*(1-p)*(1+ (3*n-6)*p*(1-p))
    E5 <- n*p*(1-p)*(1-2*p)*(1+(10*n-12)*p*(1-p))
    rawmoments <- centraltoraw(c(E1,E2,E3,E4,E5))
    if (type=="central"){
      return(c(E1,E2,E3,E4,E5))
    }
    if (type== "raw"){
      return(rawmoments)
    }
    if (type== "standardized"){
      return(rawtostandard(rawmoments))
    }
  }
  if (dist=="Geometric"){
    p <- par
    E1 <- 1/p
    E2 <- (2-p)/p^2
    E3 <- ((p^2-2*p+1)+ 4-4*p+1)/p^3
    E4 <- -(p^3-3*p^2+3*p-1 -11*p^2+22*p-11+11*p-11-1)/p^4
    E5 <- (p^4-4*p^3+6*p^2-4*p+1-26*p^3+78*p^2-78*p+26 +66*p^2-132*p+66 +26-26*p+1)/p^5

    rawmoments <-  c(E1,E2,E3,E4,E5)
    if (type=="raw"){
      return(rawmoments)
    }
    if (type== "central"){
      return(rawtocentral(rawmoments))
    }
    if (type== "standardized"){
      return(rawtostandard(rawmoments))
    }
  }
  if (dist=="NegBin"){
    r <- par[1]
    p <- par[2]
    E1 <- (1-p)*r/p
    E2 <- (p-1)*r*((p-1)*r-1)/p^2
    E3 <- -(p-1)*r*((p-1)*(r*((p-1)*r-3)-1)+1)/p^3
    E4 <- (p-1)*r*((p-1)*((p-1)*(r*(r*((p-1)*r-6)-4)-1)+7*r+4)-1)/p^4
    E5 <- -(p-1)*r*((p-1)*((p-1)^2 *(r*(r*(r*((p-1)*r-10)-10)-5)-1)+(p-1)*(25*r^2+30*r+11)-15*r-11)+1)/p^5
    rawmoments <- c(E1,E2,E3,E4,E5)
    if (type=="raw"){
      return(rawmoments)
    }
    if (type== "central"){
      return(rawtocentral(rawmoments))
    }
    if (type== "standardized"){
      return(rawtostandard(rawmoments))
    }
  }
  if (dist=="Poisson"){
    lambda <- par
    E1 <- lambda
    E2 <- lambda*(lambda+1)
    E3 <- lambda*(lambda^2+ 3*lambda+1)
    E4 <- lambda*(lambda^3 +6*lambda^2 + 7*lambda +1)
    E5 <- lambda*(lambda^4 +10*lambda^3 +25* lambda^2 +15*lambda +1)
    rawmoments <- c(E1,E2,E3,E4,E5)
    if (type=="raw"){
      return(rawmoments)
    }
    if (type== "central"){
      return(rawtocentral(rawmoments))
    }
    if (type== "standardized"){
      return(rawtostandard(rawmoments))
    }
  }
  if (dist=="Pareto"){
    alpha <- par[1]
    x0 <- par[2]
    if (alpha>1){
      E1 <-  integrate(function(x) alpha*x0^alpha * x^(-alpha), lower=x0, upper=Inf)$val
      rawmoments <- E1
    }
    if (alpha>2){
      E2 <-  integrate(function(x) alpha*x0^alpha * x^(-alpha+1), lower=x0, upper=Inf)$val
      rawmoments <- c(E1,E2)
    }
    if (alpha>3){
      E3 <-  integrate(function(x)alpha*x0^alpha * x^(-alpha+2), lower=x0, upper=Inf)$val
      rawmoments <- c(E1,E2,E3)
    }
    if (alpha>4){
      E4 <-  integrate(function(x) alpha*x0^alpha * x^(-alpha+3), lower=x0, upper=Inf)$val
      rawmoments <-  c(E1,E2,E3,E4)
    }
    if (alpha>5){
      E5 <-  integrate(function(x) alpha*x0^alpha * x^(-alpha+4), lower=x0, upper=Inf)$val
      rawmoments <- c(E1,E2,E3,E4,E5)
    }
    if (type=="raw"){
      return(rawmoments)
    }
    if (type== "central"){
      return(rawtocentral(rawmoments))
    }
    if (type== "standardized"){
      return(rawtostandard(rawmoments))
    }
  }
  if(dist=="Exponential"){
    lambda <- par
    E1= integrate(function(x) x*dexp(x,lambda), lower=0,upper=Inf)$value
    E2= integrate(function(x) x^2*dexp(x,lambda), lower=0,upper=Inf)$value
    E3= integrate(function(x) x^3*dexp(x,lambda), lower=0,upper=Inf)$value
    E4= integrate(function(x) x^4*dexp(x,lambda), lower=0,upper=Inf)$value
    E5= integrate(function(x) x^5*dexp(x,lambda), lower=0,upper=Inf)$value
    rawmoments <- c(E1,E2,E3,E4,E5)
    if (type=="raw"){
      return(rawmoments)
    }
    if (type== "central"){
      return(rawtocentral(rawmoments))
    }
    if (type== "standardized"){
      return(rawtostandard(rawmoments))
    }
  }
  if (dist=="Gamma"){
    alpha <- par[1]
    lambda <- par[2]
    E1= integrate(function(x) x*dgamma(x,alpha,lambda), lower=0,upper=Inf)$value
    E2= integrate(function(x) x^2*dgamma(x,alpha,lambda), lower=0,upper=Inf)$value
    E3= integrate(function(x) x^3*dgamma(x,alpha,lambda), lower=0,upper=Inf)$value
    E4= integrate(function(x) x^4*dgamma(x,alpha,lambda), lower=0,upper=Inf)$value
    E5= integrate(function(x) x^5*dgamma(x,alpha,lambda), lower=0,upper=Inf)$value
    rawmoments <- c(E1,E2,E3,E4,E5)
    if (type=="raw"){
      return(rawmoments)
    }
    if (type== "central"){
      return(rawtocentral(rawmoments))
    }
    if (type== "standardized"){
      return(rawtostandard(rawmoments))
    }
  }
  if (dist=="LogNormal"){
    mu <- par[1]
    sigma <- par[2]
    MG <- function(n){
      exp(mu*n +0.5*sigma^2*n^2)
    }
    E1 <- MG(1)
    E2 <- MG(2)
    E3 <- MG(3)
    E4 <- MG(4)
    E5 <- MG(5)
    rawmoments <- c(E1,E2,E3,E4,E5)
    if (type=="raw"){
      return(rawmoments)
    }
    if (type== "central"){
      return(rawtocentral(rawmoments))
    }
    if (type== "standardized"){
      return(rawtostandard(rawmoments))
    }
  }
}

















