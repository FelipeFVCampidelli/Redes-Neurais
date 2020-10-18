rm(list=ls())
library('plot3D')
library('roccv')
library('corpcor')

# Fun??es:

ELMmodel <- function(X, Y, nos) {
  Z <- replicate(nos, runif((dim(X)[2]+1),-0.5,0.5))
  Xaug <- cbind(replicate(dim(X)[1], 1), X)
  H <- tanh(Xaug %*% Z)
  W <- pseudoinverse(H) %*% Y
  parametros <- list(W, Z)
  return(parametros)
}

ELMpredict <- function(Xtest, Ytest, parametros){
  Xaug <- cbind(replicate(dim(Xtest)[1], 1), Xtest)
  H <- tanh(Xaug %*% parametros[[2]])
  Yhat <- sign(H %*% parametros[[1]])
  erro <- sum((Ytest-Yhat)^2)/4
  resultados <- list(Yhat, erro)
  return(resultados)
}

predict <- function(X, pesos, z){
  H <- tanh(X %*% t(z))
  Yhat <- sign(H %*% t(pesos))
  return(Yhat)
}
YELM<-function(xin, Z, W, par){
  n<-dim(xin)[2]
  
  if(par==1){
    xin<-cbind(xin)
  }
  
  H<-tanh(xin %*% Z)
  yhat<-sign(t(H) %*% W)
  
  return(yhat)
}