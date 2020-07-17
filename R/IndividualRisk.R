#'@title Aggregate distribution in the Individual Risk Model
#'@description The function returns the cumulative distribution function of the aggregate claim size.
#'@param x the independent variable P(S<x).
#'@param n The number of policies
#'@param dist the distribution of each single policy (only "Gamma" and "InvGauss" are supported)
#'@param par the parameters of the distribution of each single policy.
#'@return The function returns the cumulative distribution function of the aggregate claim size P(S<x)
#'@examples 1-IndividualAggregateDist(100,10,dist="InvGauss",par= c( 0.05,2.5))
#'@note Other distributions than the Gamma and the Inverse Gaussian might be considered. Nonetheless
#' It might be necessary to use convolutions or approximation methods to compute the CDF.
#' @export
IndividualAggregateDist <- function(x,n, dist,par){
  pInvGauss <- function(x,alpha, mu){
    return(pnorm(sqrt((alpha*x)/mu)- sqrt((alpha*mu)/x),0,1) + exp(2*alpha)* pnorm(-sqrt((alpha*x)/mu)- sqrt((alpha*mu)/x),0,1))
  }
  if (dist== "Gamma"){
    alpha=par[1]
    lambda=par[2]
    return(pgamma(x,alpha*n,lambda))
  }
  if (dist== "InvGauss"){
    alpha= par[1]
    mu= par[2]
    return(pInvGauss(x,alpha*n, mu*n))
  }
}
#'@title Aggregate distribution in the Individual Risk Model
#'@description The function returns the cumulative distribution function of the aggregate claim size.
#'in this model specification there are n policies X and each policy is in the form X=I*B. Where I is the
#'indicator function (bernoulli(q)) and B is either a Gamma or an Inverse Gaussian random variable
#'(describing claim intensity given occurrence).
#'@param x the independent variable P(S<x).
#'@param n The number of policies
#'@param dist the distribution of each single policy (only "Gamma" and "InvGauss" are supported)
#'@param par the parameters of the distribution of each single policy.
#'@return The function returns the cumulative distribution function of the aggregate claim size P(S<x)
#'@examples IndividualAggregateDist2(50,200,dist="Gamma",par=c(1,0.5),q=0.1)
#'@note Other distributions than the Gamma and the Inverse Gaussian might be considered. Nonetheless
#' It might be necessary to use convolutions or approximation methods to compute the CDF.
#' @export
IndividualAggregateDist2 <- function(x,n,dist,par, q){
  pInvGauss <- function(x,alpha, mu){
    return(pnorm(sqrt((alpha*x)/mu)- sqrt((alpha*mu)/x),0,1) + exp(2*alpha)* pnorm(-sqrt((alpha*x)/mu)- sqrt((alpha*mu)/x),0,1))
  }
  if (dist== "Gamma"){
    alpha=par[1]
    lambda=par[2]
    res=numeric(n)
    for (k in 1:n){
      res[k]= pgamma(x,k*alpha,lambda)*dbinom(k,n,q)
    }
    return(sum(res))
  }
  if (dist== "InvGauss"){
    alpha=par[1]
    mu=par[2]
    res=numeric(n)
    for (k in 1:n){
      res[k]= pInvGauss(x,k*alpha,mu*k)*dbinom(k,n,q)
    }
    return(sum(res))
  }
}

