rm(list=ls())
cat("\014")
library(matlib)

geraMatrizB <- function(sist, nc, nb) {
  matB <- matrix(0, nrow=nb, ncol=nb)
  for (ibarra in 1:length(indbarras)) {
    for (iconexao in 1:nc) {
      if ((ibarra == sist[iconexao,1]) || (ibarra == sist[iconexao,2])) {
        #soma todas as susceptancias ligadas a uma barra
        matB[ibarra,ibarra] <- matB[ibarra,ibarra] + sist[iconexao,"susceptancia"]
        #descobre a conexao da barra
        conexao <- sist[iconexao,which(sist[iconexao,1:2]!=ibarra)]
        #preenche as susceptancias em relacao as conexoes da barra
        matB[ibarra,conexao] <- -sist[iconexao,"susceptancia"]
        #print(paste("barra", ibarra, "aparece na conexao", iconexao, "conectada a", conexao, sep=" "))
      }
    }
  }
  
  return(matB)
}

calculaPerdasAtivas <- function(sist, nc, nb, theta) {
  Gkl <- sist[,"resistencia"]/(sist[,"resistencia"]^2 + sist[,"reatancia"]^2)
  
  perdasAtivas <- Gkl * ((theta[sist[,1]]-theta[sist[,2]])^2)
  perdasAtivas <- perdasAtivas/2
  perdasAtivas <- cbind(sist[,1:2],perdasAtivas)
  
  perdasAtivasBarra <- rep(0, nb)
  for (ibarra in 1:nb) {
    for (iconexao in 1:nc) {
      if ((ibarra == perdasAtivas[iconexao,1]) || (ibarra == perdasAtivas[iconexao,2])) {
        conexao <- perdasAtivas[iconexao,which(sist[iconexao,1:2]!=ibarra)]
        perdasAtivasBarra[ibarra] <- perdasAtivasBarra[ibarra] + perdasAtivas[iconexao,3]
      }
    }
  }
  
  return(perdasAtivasBarra)
}

calculaPoteAtivas <- function(sist, nc, nb, thetaComPerdas) {
  poteAtivasBarra <- matrix(0, nrow=nc, ncol=3)
  for (ibarra in 1:nb) {
    for (iconexao in 1:nc) {
      if ((ibarra == sist[iconexao,1]) || (ibarra == sist[iconexao,2])) {
        conexao <- sist[iconexao,which(sist[iconexao,1:2]!=ibarra)]
        
        poteAtivasBarra[iconexao,1] <- ibarra
        poteAtivasBarra[iconexao,2] <- conexao
        poteAtivasBarra[iconexao,3] <- sist[iconexao,"susceptancia"] * 
          (thetaComPerdas[ibarra] - thetaComPerdas[conexao])
        
      }
    }
  }
  
  return(poteAtivasBarra)
}

somaPotenciaNaBarra <- function(M, nb) {
  potTotal <- rep(NA, nb)
  for (i in 1:nb) {
    potTotal[i] <- sum(M[which(M[,1]==i | M[,2]==i),3])
  }
  return(potTotal)
}

geraMatPotenciaNaBarra <- function(MPot, nb) {
  intervalos <- length(MPot[1,])-2
  
  m <- matrix()[-1]
  for (i in 3:(intervalos+2)) {
    m <- rbind(m, somaPotenciaNaBarra(MPot[,c(1:2,i)], nb))
  }
  
  return(m)
}

#potencia base
pbase <- 100

#numero de barras conectadas pelo menos uma vez
nconexoes <- 7

#intervalo de variacao da potencia
linf <- 1
lsup <- 1.5

#MATRIZ de impedancias (JÁ ESTÁ EM PU)
sist <- matrix(c(1,2,0.02,0.06,
                 1,3,0.08,0.24,
                 2,3,0.06,0.18,
                 2,4,0.06,0.18,
                 2,5,0.04,0.12,
                 3,4,0.01,0.03,
                 4,5,0.08,0.24), nrow=nconexoes, ncol=4, byrow=TRUE)

sist <- cbind(sist, 1/sist[,4])
colnames(sist) <- c("p","q","resistencia","reatancia","susceptancia")

indbarras <- data.frame(table(sist[,1:2]))
indbarras <- as.numeric(as.vector(indbarras[,1]))
nbarras <- length(indbarras)

#MATRIX de potencias
poteOrig <- matrix(c(1,0.06,0,0,0,0,0,
                 2,1,0,40,30,20,10,
                 3,1,0,0,0,45,15,
                 4,1,0,0,0,40,5,
                 5,1,0,0,0,60,10), nrow=nbarras, ncol=7, byrow=TRUE)

colnames(poteOrig) <- c("p","v_real","v_imag","pa_gen","pr_gen","pa_carga","pr_carga")

#MATRIX de susceptancias B nas barras
matB <- geraMatrizB(sist, nconexoes, nbarras)

#MATRIX que armazena as potencias a cada experimento
matPot <- matrix()[-1]
#MATRIX que armazena os angulos de tensao sem perdas
matThetaSemPerdas <- matrix()[-1]
#MATRIX que armazena os angulos de tensao com perdas
matThetaComPerdas <- matrix()[-1]

cont <- 0
for(iMultPot in seq(linf, lsup, by = 0.1)) {
  #calculo dos angulos SEM PERDAS
  matBAjus <- matB[-1,-1]
  
  poteOrigAjus <- cbind(poteOrig[,1:3],poteOrig[,4:7]*iMultPot)
  
  poteAtivas <- (poteOrigAjus[,"pa_gen"] - poteOrigAjus[,"pa_carga"])/pbase
  poteAtivasAjus <- poteAtivas[-1]
  thetaSemPerdas <- poteAtivasAjus %*% Ginv(matBAjus)
  thetaSemPerdas <- cbind(0, thetaSemPerdas)
  
  #calculo das PERDAS
  perdasAtivasBarra <- calculaPerdasAtivas(sist, nconexoes, nbarras, thetaSemPerdas)
  
  #calculo dos angulos COM PERDAS
  perdasAtivasBarraAjus <- perdasAtivasBarra[-1]
  thetaComPerdas <- (poteAtivasAjus - perdasAtivasBarraAjus) %*% Ginv(matBAjus)
  
  #calculo das POTENCIAS
  thetaComPerdas <- cbind(0, thetaComPerdas)
  poteAtivasBarra <- calculaPoteAtivas(sist, nconexoes, nbarras, thetaComPerdas)
  colnames(poteAtivasBarra) <- c("p","q","potencia")
  
  matPot <- cbind(matPot, poteAtivasBarra[,"potencia"])
  matThetaSemPerdas <- rbind(matThetaSemPerdas, thetaSemPerdas)
  matThetaComPerdas <- rbind(matThetaComPerdas, thetaComPerdas)
}

matPot <- cbind(poteAtivasBarra[,c("p","q")], matPot)
#matPotBarras <- geraMatPotenciaNaBarra(matPot, nbarras)
