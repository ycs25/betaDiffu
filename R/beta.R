# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Simulate a diffusion process with beta invariant distribution
#'
#' @param dt Equal subintervals
#' @param alpha The beta distribution parameter alpha
#' @param beta The beta distribution parameter beta
#' @param theta The beta distribution parameter theta
#'
#' @return A simulation of diffusion process given beta invariant distribution
#' @export
#'
#' @examples
#' betaDiff()
#' betaDiff(alpha=2,beta=4,theta=0.5)
betaDiff <- function(T=10,dt=0.01,alpha=3,beta=5,theta=1) {
  t = seq(0,T,dt)
  #Wiener process
  w = rnorm(n = length(t) - 1, sd = sqrt(dt))
  #diffusion process
  N =T/dt
  x = alpha/(alpha+beta)
  X = rep(0,N)
  #Euler-Maruyama
  for (i in 1:N) {
    mu = theta*(alpha/(alpha+beta)-x)
    sigma = sqrt(2*theta/(alpha+beta)*x*(1-x))
    x = x + mu*dt + sigma*w[i]
    X[i]=x
  }
  X=c(alpha/(alpha+beta),X)
}

#' Title
#'
#' @param X A diffusion process with beta invariant distribution
#' @param dt Equal subintervals
#'
#' @return Esimates of parameters in diffusion model with beta invariant distribution
#' @export
#'
#' @examples
betaEstimating <- function(X,dt=0.01) {

  est<-function(Y) {

    g1=0;g2=0;g3=0;N=length(X)-1;Z=X-1;
    for (t in 1:N) {
      x1=Z[t]
      x2=Z[t+1]

      h1=sqrt((Y[1]+Y[2]+1)*Y[2]/Y[1])*(1+(Y[1]+Y[2])/Y[2]*x2)-exp(-Y[3]*dt)*sqrt((Y[1]+Y[2]+1)*Y[2]/Y[1])*(1+(Y[1]+Y[2])/Y[2]*x1)

      g1=g1+h1

      h2=sqrt((Y[1]+Y[2]+3)*Y[2]*(Y[2]+1)*(Y[1]+Y[2])/Y[1]/(Y[1]+1)/2)*
        (1+2*(Y[1]+Y[2]+1)/Y[2]*x2+(Y[1]+Y[2]+1)*(Y[1]+Y[2]+2)/Y[2]/(Y[2]+1)*x2^2)-
        exp(-Y[3]*2*dt*(Y[1]+Y[2]+1)/(Y[1]+Y[2]))*sqrt((Y[1]+Y[2]+3)*Y[2]*(Y[2]+1)*(Y[1]+Y[2])/Y[1]/(Y[1]+1)/2)*
        (1+2*(Y[1]+Y[2]+1)/Y[2]*x1+(Y[1]+Y[2]+1)*(Y[1]+Y[2]+2)/Y[2]/(Y[2]+1)*x1^2)

      g2=g2+h2

      h3=sqrt((Y[1]+Y[2]+5)*Y[2]*(Y[2]+1)*(Y[2]+2)*(Y[1]+Y[2])*(Y[1]+Y[2]+1)/Y[1]/(Y[1]+1)/(Y[1]+2)/6)*
        (1+3*(Y[1]+Y[2]+2)/Y[2]*x2+3*(Y[1]+Y[2]+2)*(Y[1]+Y[2]+3)/Y[2]/(Y[2]+1)*x2^2+
           (Y[1]+Y[2]+2)*(Y[1]+Y[2]+3)*(Y[1]+Y[2]+4)/Y[2]/(Y[2]+1)/(Y[2]+2)*x2^3)-
        exp(-Y[3]*3*dt*(Y[1]+Y[2]+2)/(Y[1]+Y[2]))*sqrt((Y[1]+Y[2]+5)*Y[2]*(Y[2]+1)*(Y[2]+2)*(Y[1]+Y[2])*(Y[1]+Y[2]+1)/Y[1]/(Y[1]+1)/(Y[1]+2)/6)*
        (1+3*(Y[1]+Y[2]+2)/Y[2]*x1+3*(Y[1]+Y[2]+2)*(Y[1]+Y[2]+3)/Y[2]/(Y[2]+1)*x1^2+
           (Y[1]+Y[2]+2)*(Y[1]+Y[2]+3)*(Y[1]+Y[2]+4)/Y[2]/(Y[2]+1)/(Y[2]+2)*x1^3)

      g3=g3+h3
    }


    return(c(F1=g1,F2=g2,F3=g3))
  }

  library(rootSolve)
  s=multiroot(est,c(3,5,1))
  return(s)
}
