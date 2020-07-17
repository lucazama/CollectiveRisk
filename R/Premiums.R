#you need the moment generating function of S. beta is the risk aversion coefficient.

ExponentialPremium <- function(beta, MGM){
  return(1/beta * MGM(beta))
}
