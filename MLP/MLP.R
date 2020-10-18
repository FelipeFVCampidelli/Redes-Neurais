rm(list=ls())
dev.off()

MLPerceptron <- function(xin, yd, eta, tol, maxepocas, neuronios, xtest, ytest, fold) {
  dimxin <- dim(xin)
  N <- dimxin[1]
  n <- dimxin[2]
  
  wo <- matrix( runif( (n+1)*neuronios, -0.5, 0.5), nrow =neuronios, ncol=n+1 )
  wt <- matrix(runif(neuronios+1)-0.5, nrow = 1)
  xin <- cbind(1, xin)
  xtest <- cbind(1, xtest)
  
  nepocas <- 0
  eepoca <- tol + 1
  
  evec <- matrix(0, nrow = 1, ncol = maxepocas)
  eTestvec <- matrix(0, nrow = 1, ncol = maxepocas)
  while((nepocas < maxepocas) && (eepoca > tol)) {
    erro <- errotest <- 0
    xseq <- sample(N)
    
    for(i in 1:N) {
      irand <- xseq[i]
      
      z1 <- wo %*% xin[irand, ]
      a1 <- rbind(1, tanh(z1))
      
      z2 <- wt %*% a1
      #yhati <- tanh(z2)
      yhati <- z2
      
      e <- yd[irand]-yhati
      deltaE2 <- -1*e
      dwt <- eta*deltaE2 %*% t(a1)
      
      dwo <- matrix(0,dim(wo)[1], dim(wo)[2])
      for(i in 1:dim(wo)[1]) {
        dwo[i,] <- ( eta*deltaE2*wt[,i+1]*( 1/cosh(z1[i,])^2 ) ) %*% t(xin[irand, ])
      }
      
      wt <- wt - dwt
      wo <- wo - dwo
      erro <- erro + e*e
    }
    
    xtestseq <- sample(dim(xtest)[1])
    for(i in 1:dim(xtest)[1]) {
      irandtest <- xtestseq[i]
      Z1test <- wo %*% xtest[irandtest, ]
      A1test <- tanh(Z1test)
      Yhattest <- wt %*% rbind(1,A1test)
      Predict <- Yhattest
      etest <- ytest[irandtest] - Predict
      errotest <- errotest + etest*etest
    }
    
    nepocas <- nepocas + 1
    
    evec[nepocas] <- erro/N
    eTestvec[nepocas] <- errotest/dim(xtest)[1]
    
    eepoca <- evec[nepocas]
    
    if(nepocas %% 100 == 0) cat("Erro[", fold, ",", nepocas,"]:", evec[nepocas], "\n")
  }
  retlist <- list(wo, wt, evec[1:nepocas], eTestvec[1:nepocas])
  return(retlist)
}

MLPredict <- function(xin, model) {
  W1 <- model[[1]]
  W2 <- model[[2]]
  X <- cbind(1,xin)
  Predict <- matrix(0, dim(xin)[1])
  
  for(i in 1:dim(X)[1]) {
    Z <- W1 %*% X[i,]
    A <- tanh(Z)
    Yhat <- tanh(W2 %*% rbind(1,A))
    Predict[i] <- Yhat
  }
  return(Predict)
}