rm(list=ls())
cat("\014")
library(matlib)

#-----------------------------------------------------------------------------#
#                                                                             #
#                FUNCOES PARA CALCULO DO FLUXO DE POTENCIA AC                 #
#                                                                             #
#                     CINIRO APARECIDO LEITE NAMETALA                         #
#-----------------------------------------------------------------------------#

#conversao de graus-rad e rad-graus
rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}

#funcao que retorna uma matriz com os valores escalares de g e b separados por linha
calcAdmitancias <- function(mLinhas) {
  admitancias <- matrix(0, nrow=nlinhas, ncol=2)
  colnames(admitancias) <- c("g","b")
  
  for (ilinha in 1:nlinhas) {
    #calcula a admitancia
    admitancias[ilinha,"g"] <- mLinhas[ilinha,"r"]/(mLinhas[ilinha,"r"]^2 + mLinhas[ilinha,"x"]^2)
    admitancias[ilinha,"b"] <- - mLinhas[ilinha,"x"]/(mLinhas[ilinha,"r"]^2 + mLinhas[ilinha,"x"]^2)
  }
  
  return(admitancias)
}

#funcao que retorna as matrizes G e B
calcMatrizesGB <- function(mLinhas, mBarras) {
  
  matG <- matrix(0, nrow=nbarras, ncol=nbarras)
  matB <- matrix(0, nrow=nbarras, ncol=nbarras)
  
  for (ibarra in 1:nbarras) {
    
    #soma todas as susceptancias shunt ligadas a barra
    if (!is.na(mBarras[ibarra,"bshk"]))
      matB[ibarra,ibarra] <- matB[ibarra,ibarra] + mBarras[ibarra,"bshk"]
    
    for (ilinha in 1:nlinhas) {
      if ((ibarra == mLinhas[ilinha,"k"]) || (ibarra == mLinhas[ilinha,"m"])) {
        #soma todas as condutancias e susceptancias das linhas ligadas a uma barra
        matG[ibarra,ibarra] <- matG[ibarra,ibarra] + mLinhas[ilinha,"g"]
        matB[ibarra,ibarra] <- matB[ibarra,ibarra] + mLinhas[ilinha,"b"]
        
        #soma todas as susceptancias shunt ligadas as linhas onde a barra esta conectada
        if (!is.na(mLinhas[ilinha,"bshkm"]))
          matB[ibarra,ibarra] <- matB[ibarra,ibarra] + mLinhas[ilinha,"bshkm"]
        
        #descobre a conexao da barra
        conexao <- mLinhas[ilinha, which(mLinhas[ilinha,1:2]!=ibarra)]
        
        #preenche as condutancias e susceptancias em relacao as conexoes da barra
        #com sinal invertido
        matG[ibarra,conexao] <- -mLinhas[ilinha,"g"]
        matB[ibarra,conexao] <- -mLinhas[ilinha,"b"]
      }
    }
  }
  
  return(list(matG, matB))
}

#funcao que retorna vetores de identificadores de equacoes e incognitas do problema,
#alem de potencias a serem calculadas (faltantes)
identificaEqIncPot <- function(mBarras) {
  vEquacoes <- vector()
  vIncognitas <- vector()
  vPotencias <- vector()
  
  #identificando equacoes de potencia ativa
  for (ibarra in 1:nbarras) {
    if (mBarras[ibarra,"tipo"] == "PQ" || mBarras[ibarra,"tipo"] == "PV") {
      vEquacoes <- c(vEquacoes, paste("P", ibarra, sep=""))
      vIncognitas <- c(vIncognitas, paste("T", ibarra, sep=""))
    }
  }
  
  #identificando equacoes de potencia reativa
  for (ibarra in 1:nbarras) {
    if (mBarras[ibarra,"tipo"] == "PQ") {
      vEquacoes <- c(vEquacoes, paste("Q", ibarra, sep=""))
      vIncognitas <- c(vIncognitas, paste("V", ibarra, sep=""))
      
    }
  }
  
  #identificando equacoes de potencia ativa faltante
  for (ibarra in 1:nbarras) {
    if (mBarras[ibarra,"tipo"] == "VT") {
      vPotencias<- c(vPotencias, paste("P", ibarra, sep=""))
    }
  }
  
  #identificando equacoes de potencia reativa faltante
  for (ibarra in 1:nbarras) {
    if (mBarras[ibarra,"tipo"] == "VT" || mBarras[ibarra,"tipo"] == "PV") {
      vPotencias<- c(vPotencias, paste("Q", ibarra, sep=""))
    }
  }
  
  return(list(vEquacoes, vIncognitas, vPotencias))
}

#funcao que retorna um vetor com valores 0 para theta e 1 para tensao (chutes iniciais)
geraVetorXInicial <- function(vetIncognitas) {
  x <- vector()
  
  for (iinc in vetIncognitas) {
    tipo <- substr(iinc, 1, 1)
    k <- substr(iinc, 2, 3)
    
    if (tipo == "T")
      x <- c(x, 0)
    else
      x <- c(x, 1)
  }
  
  return(x)
}

#funcao que retorna um vetor com as potencias calculadas
calculaPotencias <- function(vetEquacoes, vetIncognitas, x, barras, linhas, matG, matB) {
  
  x <- cbind(vetIncognitas, x)
  colnames(x) <- c("incognita","valor")
  
  vetF <- vector()
  
  for (ieq in vetEquacoes) {
    tipo <- substr(ieq, 1, 1)
    k <- substr(ieq, 2, 3)
    
    vetF <- c(vetF, calcDeltaPotencia(tipo, k, x, barras, linhas, matG, matB))
  }
  
  return(vetF)
}

#funcao que retorna o delta das potencias ativa e reativa, por tipo "P" ou "Q"
calcDeltaPotencia <- function(tipo, k, x, mBarras, mLinhas, matG, matB) {
  
  k <- as.numeric(k)
  
  #pega os dados da barra K
  thetav <- pegaDadosBarra(k, mBarras, x)
  tk <- as.numeric(thetav[1])
  vk <- as.numeric(thetav[2])
  rm(thetav)
  
  #captura o valor de potencia especificada para P ou Q
  if (tipo == "P") {
    potEspk <- as.numeric(mBarras[k, "p"])
  } else {
    potEspk <- as.numeric(mBarras[k, "q"]) 
  }
  
  #cria o vetor de barras m conectas a uma dada barra k
  ms <- pegaBarrasConectas(k, mLinhas)
  
  #inclui a propria barra k
  ms <- c(k, ms)
  
  #inicia o calculo da potencia
  somatorio <- 0
  
  for (m in ms) {
    
    m <- as.numeric(m)
    
    #captura dados da barra m
    thetav <- pegaDadosBarra(m, mBarras, x)
    tm <- as.numeric(thetav[1])
    vm <- as.numeric(thetav[2])
    rm(thetav)
    
    #calcula o termo do somatorio da equacao de potencia (ativa ou reativa)
    if (tipo == "P") {
      somatorio <- somatorio + ( vm * ( matG[k,m]*cos(tk - tm) + matB[k,m]*sin(tk - tm) ) )
    } else {
      somatorio <- somatorio + ( vm * ( matG[k,m]*sin(tk - tm) - matB[k,m]*cos(tk - tm) ) )
    }
    
  }
  
  #finaliza o calculo da potencia
  potK <- potEspk - vk*somatorio
  
  return (potK)
}

#funcao que retorna os valores de Theta e V para barras PQ ou PV ou VT
pegaDadosBarra <- function(indBarra, mBarras, x) {
  if (mBarras[indBarra, "tipo"] == "PQ") {
    tk <- x[which(x[,1] == paste("T",indBarra,sep="")),"valor"]
    vk <- x[which(x[,1] == paste("V",indBarra,sep="")),"valor"]
  }
  else if (mBarras[indBarra, "tipo"] == "PV") {
    tk <- x[which(x[,1] == paste("T",indBarra,sep="")),"valor"]
    vk <- mBarras[indBarra, "v"]
  }
  else {
    tk <- mBarras[indBarra, "theta"]
    vk <- mBarras[indBarra, "v"]
  }
  
  return (c(tk, vk))
}

#funcao que retornas as barras conectadas a barra m
pegaBarrasConectas <- function(k, mLinhas) {
  ms <- vector()
  
  for (ilinha in 1:nlinhas) {
    if (mLinhas[ilinha, "k"] == k)
      ms <- c(ms, as.numeric(mLinhas[ilinha, "m"]))
    
    if (mLinhas[ilinha, "m"] == k)
      ms <- c(ms, as.numeric(mLinhas[ilinha, "k"]))
  }
  
  return(ms)
}

#funcao que retorna TRUE se uma dada barra m estiver conectada a uma barra k
verificaConexaoKM <- function(k, m, mLinhas) {
  ms <- vector()
  
  for (ilinha in 1:nlinhas) {
    if (mLinhas[ilinha, "k"] == k)
      ms <- c(ms, as.numeric(mLinhas[ilinha, "m"]))
    
    if (mLinhas[ilinha, "m"] == k)
      ms <- c(ms, as.numeric(mLinhas[ilinha, "k"]))
  }
  
  conexao <- FALSE
  
  if (length(ms[which(ms == m)])>0 || k==m)
    conexao <- TRUE
  
  return(conexao)
}

#funcao que retorna os elementos textuais para se calcular na matriz J
montaJacobianaNominal <- function(vetEquacoes, vetIncognitas, mLinhas) {
  
  jacN <- matrix(0, nrow=length(vetEquacoes), ncol=length(vetIncognitas))
  i <- 0
  j <- 0
  for (eq in vetEquacoes) {
    i <- i + 1
    for (inc in vetIncognitas) {
      j <- j + 1
      
      tipo <- paste(substr(eq, 1, 1), substr(inc, 1, 1), sep="")
      k <- substr(eq, 2, 3)
      m <- substr(inc, 2, 3)
      indicesBarras <- paste(k, m, sep=",")
      
      if (verificaConexaoKM(k, m, mLinhas) == TRUE) {
        switch(tipo,
               #quadrante H
               PT={
                 jacN[i,j] <- paste("H",indicesBarras,sep="")
               },
               #quadrante N
               PV={
                 jacN[i,j] <- paste("N",indicesBarras,sep="")
               },
               #quadrante M
               QT={
                 jacN[i,j] <- paste("M",indicesBarras,sep="")
               },
               #quadrante L
               QV={
                 jacN[i,j] <- paste("L",indicesBarras,sep="")
               },
               stop("Valor invalido!")
        )
      }
    }
    j <- 0
  }
  
  return(jacN)
}

#funcao que calcula a matriz Jacobiana
calculaJacobiana <- function(jacNominal, mBarras, mLinhas, matG, matB, x) {
  x <- cbind(vetIncognitas, x)
  colnames(x) <- c("incognita","valor")
  
  matJ <- jacNominal
  
  for (i in 1:dim(matJ)[1]) {
    for (j in 1:dim(matJ)[2]) {
      if (matJ[i,j] != 0) {
        
        #captura os dados da equacao que ser calculada
        quadrante <- substr(matJ[i,j], 1, 1)
        k <- as.numeric(substr(matJ[i,j],2,(which(strsplit(matJ[i,j], "")[[1]]==",")-1)))
        m <- as.numeric(substr(matJ[i,j],(which(strsplit(matJ[i,j], "")[[1]]==",")+1),nchar(matJ[i,j])))
        
        #pega os dados da barra K
        thetav <- pegaDadosBarra(k, mBarras, x)
        tk <- as.numeric(thetav[1])
        vk <- as.numeric(thetav[2])
        rm(thetav)
        
        switch(quadrante,
               #quadrante H
               H={
                 if (k==m) {
                   #cria o vetor de barras m conectas a uma dada barra k
                   ms <- pegaBarrasConectas(k, mLinhas)
                   
                   hkk <- 0
                   
                   for (m in ms) {
                     
                     m <- as.numeric(m)
                     
                     #captura dados da barra m
                     thetav <- pegaDadosBarra(m, mBarras, x)
                     tm <- as.numeric(thetav[1])
                     vm <- as.numeric(thetav[2])
                     rm(thetav)
                     
                     #calcula o termo do somatorio da equacao
                     hkk <- hkk + ( vm * ( - matG[k,m]*sin(tk - tm) + matB[k,m]*cos(tk - tm) ) )
                   }
                   
                   #finaliza o calculo da potencia
                   matJ[i,j] <- vk*hkk
                 } else {
                   #captura dados da barra m
                   thetav <- pegaDadosBarra(m, mBarras, x)
                   tm <- as.numeric(thetav[1])
                   vm <- as.numeric(thetav[2])
                   rm(thetav)
                   
                   matJ[i,j] <- vk*vm*(matG[k,m]*sin(tk - tm) - matB[k,m]*cos(tk - tm))
                 }
               },
               #quadrante N
               N={
                 if (k==m) {
                   #cria o vetor de barras m conectas a uma dada barra k
                   ms <- pegaBarrasConectas(k, mLinhas)
                   
                   nkk <- 0
                   
                   for (m in ms) {
                     
                     m <- as.numeric(m)
                     
                     #captura dados da barra m
                     thetav <- pegaDadosBarra(m, mBarras, x)
                     tm <- as.numeric(thetav[1])
                     vm <- as.numeric(thetav[2])
                     rm(thetav)
                     
                     #calcula o termo do somatorio da equacao
                     nkk <- nkk + ( vm * ( matG[k,m]*cos(tk - tm) + matB[k,m]*sin(tk - tm) ) )
                   }
                   
                   #finaliza o calculo da potencia
                   matJ[i,j] <- (2*vk*matG[k,k]) + nkk
                 } else {
                   #captura dados da barra m
                   thetav <- pegaDadosBarra(m, mBarras, x)
                   tm <- as.numeric(thetav[1])
                   vm <- as.numeric(thetav[2])
                   rm(thetav)
                   
                   matJ[i,j] <- vk*(matG[k,m]*cos(tk - tm) + matB[k,m]*sin(tk - tm))
                 }
               },
               #quadrante M
               M={
                 if (k==m) {
                   #cria o vetor de barras m conectas a uma dada barra k
                   ms <- pegaBarrasConectas(k, mLinhas)
                   
                   mkk <- 0
                   
                   for (m in ms) {
                     
                     m <- as.numeric(m)
                     
                     #captura dados da barra m
                     thetav <- pegaDadosBarra(m, mBarras, x)
                     tm <- as.numeric(thetav[1])
                     vm <- as.numeric(thetav[2])
                     rm(thetav)
                     
                     #calcula o termo do somatorio da equacao
                     mkk <- mkk + ( vm * ( matG[k,m]*cos(tk - tm) + matB[k,m]*sin(tk - tm) ) )
                   }
                   
                   #finaliza o calculo da potencia
                   matJ[i,j] <- vk*mkk
                 } else {
                   #captura dados da barra m
                   thetav <- pegaDadosBarra(m, mBarras, x)
                   tm <- as.numeric(thetav[1])
                   vm <- as.numeric(thetav[2])
                   rm(thetav)
                   
                   matJ[i,j] <- -vk*vm*(matG[k,m]*cos(tk - tm) + matB[k,m]*sin(tk - tm))
                 }
               },
               #quadrante L
               L={
                 if (k==m) {
                   #cria o vetor de barras m conectas a uma dada barra k
                   ms <- pegaBarrasConectas(k, mLinhas)
                   
                   lkk <- 0
                   
                   for (m in ms) {
                     
                     m <- as.numeric(m)
                     
                     #captura dados da barra m
                     thetav <- pegaDadosBarra(m, mBarras, x)
                     tm <- as.numeric(thetav[1])
                     vm <- as.numeric(thetav[2])
                     rm(thetav)
                     
                     #calcula o termo do somatorio da equacao
                     lkk <- lkk + ( vm * ( matG[k,m]*sin(tk - tm) - matB[k,m]*cos(tk - tm) ) )
                   }
                   
                   #finaliza o calculo da potencia
                   matJ[i,j] <- (-2*vk*matB[k,k]) + lkk
                 } else {
                   #captura dados da barra m
                   thetav <- pegaDadosBarra(m, mBarras, x)
                   tm <- as.numeric(thetav[1])
                   vm <- as.numeric(thetav[2])
                   rm(thetav)
                   
                   matJ[i,j] <- vk*(matG[k,m]*sin(tk - tm) - matB[k,m]*cos(tk - tm))
                 }
               },
               stop("Valor invalido!")
        )
        
      }
    }
  }
  
  matJ <- mapply(matJ, FUN=as.numeric)
  return( matrix(matJ, nrow=dim(jacNominal)[1], ncol=dim(jacNominal)[2]))
}

#funcao que retorna o valor final com base em um angulo e tensao calculados finais
#para potencia P ou Q
calcPotenciaFinal <- function(tipo, k, x, mBarras, mLinhas, matG, matB) {
  
  k <- as.numeric(k)
  
  #pega os dados da barra K
  thetav <- pegaDadosBarra(k, mBarras, x)
  tk <- as.numeric(thetav[1])
  vk <- as.numeric(thetav[2])
  rm(thetav)
  
  #cria o vetor de barras m conectas a uma dada barra k
  ms <- pegaBarrasConectas(k, mLinhas)
  
  #inclui a propria barra k
  ms <- c(ms, k)
  
  #inicia o calculo da potencia
  somatorio <- 0
  
  for (m in ms) {
    
    m <- as.numeric(m)
    
    #captura dados da barra m
    thetav <- pegaDadosBarra(m, mBarras, x)
    tm <- as.numeric(thetav[1])
    vm <- as.numeric(thetav[2])
    rm(thetav)
    
    #calcula o termo do somatorio da equacao de potencia (ativa ou reativa)
    if (tipo == "P") {
      somatorio <- somatorio + ( vm * ( matG[k,m]*cos(tk - tm) + matB[k,m]*sin(tk - tm) ) )
    } else {
      somatorio <- somatorio + ( vm * ( matG[k,m]*sin(tk - tm) - matB[k,m]*cos(tk - tm) ) )
    }
    
  }
  
  #finaliza o calculo da potencia
  potK <- vk*somatorio
  
  return (potK)
}

#funcao que calcula as potencias faltantes no sistema com vetX encontrado final
calculaPotenciaBarrasFinal <- function(vetPotencias, x, barras, linhas, matG, matB) {
  
  x <- cbind(vetIncognitas, x)
  colnames(x) <- c("incognita","valor")
  
  vetP <- vector()
  
  for (ieq in vetPotencias) {
    tipo <- substr(ieq, 1, 1)
    k <- substr(ieq, 2, 3)
    
    vetP <- c(vetP, calcPotenciaFinal(tipo, k, x, barras, linhas, matG, matB))
  }
  
  return(vetP)
}

#funcao que atualiza a matriz de barras com os dados faltantes calculados finais
atualizaMatrizBarras <- function(barras, vetPotencias, vetP, vetIncognitas, vetX) {
  
  novoBarras <- barras
  
  i <- 1
  for (ipot in vetPotencias) {
    tipo <- substr(ipot, 1, 1)
    
    switch(tipo,
           "P"={tipo <- "p"},
           "Q"={tipo <- "q"}
    )
    
    k <- as.numeric(substr(ipot, 2, 3))
    
    novoBarras[k, tipo] <- vetP[i]
    i <- i + 1
  }
  
  i <- 1
  for (iinc in vetIncognitas) {
    tipo <- substr(iinc, 1, 1)
    
    switch(tipo,
           "T"={tipo <- "theta"},
           "V"={tipo <- "v"}
    )
    
    k <- as.numeric(substr(iinc, 2, 3))
    
    novoBarras[k, tipo] <- vetX[i]
    i <- i + 1
  }
  
  # for (i in 1:dim(barras)[1]) {
  #   if (barras[i, "limites"] == 1) {
  #     novoBarras[i, "q"] <- barras[i, "q"]
  #     novoBarras[i, "v"] <- barras[i, "v"]
  #   }
  # }
  
  return(novoBarras)
}

#funcao que atualiza a matriz de linhas com os dados das correntes entre barras
calculaCorrentesEntreBarras <- function(matG, matB, mBarras, mLinhas) {
  novaLinhas <- mLinhas
  
  i <- vector()
  
  for (ilinha in 1:dim(mLinhas)[1]) {
    k <- as.numeric(mLinhas[ilinha, "k"])
    m <- as.numeric(mLinhas[ilinha, "m"])
    
    ykm <- complex(real = mLinhas[ilinha, "g"], imaginary = mLinhas[ilinha, "b"])
    
    vk <- as.numeric(mBarras[k, "v"])
    tk <- as.numeric(mBarras[k, "theta"])
    vk <- complex(real = (vk*cos(tk)), imaginary = (vk*sin(tk)))
    
    vm <- as.numeric(mBarras[m, "v"])
    tm <- as.numeric(mBarras[m, "theta"])
    vm <- complex(real = (vm*cos(tm)), imaginary = (vm*sin(tm)))
    
    bshkm <- complex(real = 0, imaginary = 0)
    if (is.na(mLinhas[ilinha, "bshkm"]) == FALSE)
      bshkm <- complex(real = 0, imaginary = as.numeric(mLinhas[ilinha, "bshkm"]))
    
    ikm <- (ykm*(vk - vm)) + (bshkm*vk)
    
    i <- c(i, ikm)
  }
  
  return(cbind(novaLinhas, i))
}

#funcao que atualiza a matriz de linhas com os dados das potencias entre barras
calculaPotenciasEntreBarras <- function(mLinhas, mBarras) {
  novaLinhas <- mLinhas
  
  s <- vector()
  
  for (ilinha in 1:dim(mLinhas)[1]) {
    k <- as.numeric(mLinhas[ilinha, "k"])
    
    ikm <- as.complex(mLinhas[ilinha, "i"])
    
    vk <- as.numeric(mBarras[k, "v"])
    tk <- as.numeric(mBarras[k, "theta"])
    vk <- complex(real = (vk*cos(tk)), imaginary = (vk*sin(tk)))
    
    skm <- vk*Conj(ikm)
    
    s <- c(s, skm)
  }
  
  return(cbind(novaLinhas, s))
}

#funcao que averigua a violacao dos limites de reativos
violaReativos <- function(barras) {
  #captura somente as barras PV
  matPV <- barras[barras[, "tipo"] == "PV",]
  violouTotal <- FALSE
  
  for(i in 1:dim(matPV)[1]) {
    violou <- "NA"
    
    #confere maximo e minimo
    if (matPV[i, "q"] >= matPV[i, "qmax"]) {
      violou <- "max"
    }
    else if (matPV[i, "q"] <= matPV[i, "qmin"]) {
      violou <- "min"
    }
    
    #se houver violacao entao transforma a PV para PQ
    if (violou == "max") {
      indbarra <- matPV[i, "barra"]
      barras[indbarra, "tipo"] <- "PQ"
      barras[indbarra, "q"] <- matPV[i, "qmax"]
      barras[indbarra, "limites"] <- 1
    } else if (violou == "min") {
      indbarra <- matPV[i, "barra"]
      barras[indbarra, "tipo"] <- "PQ"
      barras[indbarra, "q"] <- matPV[i, "qmin"]
      barras[indbarra, "limites"] <- 1
    }
    
    if (violou != "NA")
      violouTotal <- TRUE
  }
  
  return(list(violouTotal,barras))
}

#funcao que averigua a violacao dos limites de tensao para barras convertidas
violaTensao <- function(barras, vmax, vmin) {
  #captura somente as barras convertidas
  matPQ <- barras[barras[, "limites"] == 1,]
  violouTotal <- FALSE
  
  if (dim(matPQ)[1] == 0)
    return(violouTotal)
  
  for(i in 1:dim(matPQ)[1]) {
    violou <- "NA"
    
    #confere maximo e minimo
    if (matPQ[i, "v"] > vmax) {
      violou <- "max"
    }
    else if (matPQ[i, "v"] < vmin)  {
      violou <- "min"
    }
    
    #se houver violacao entao transforma a PQ de volta para PV
    if (violou == "max") {
      indbarra <- matPQ[i, "barra"]
      barras[indbarra, "tipo"] <- "PV"
      barras[indbarra, "v"] <- vmax
      barras[indbarra, "limites"] <- 0
    } else if (violou == "min") {
      indbarra <- matPQ[i, "barra"]
      barras[indbarra, "tipo"] <- "PV"
      barras[indbarra, "v"] <- vmin
      barras[indbarra, "limites"] <- 0
    }
    
    if (violou != "NA")
      violouTotal <- TRUE
  }
  
  return(list(violouTotal,barras))
}

#funcao que atualiza o problema baseado na violacao de reativos
atualizaProblema <- function(barras, vetX, vetIncognitas, vetEquacoes) {
  #atualiza equacoes e incognitas do problema
  eqincpot <- identificaEqIncPot(barras)
  vetEquacoesNovo <- eqincpot[[1]]
  vetIncognitasNovo <- eqincpot[[2]]
  vetPotenciasNovo <- eqincpot[[3]]
  rm(eqincpot)
  
  #atualiza o vetX (incognitas)
  vetXNovo <- geraVetorXInicial(vetIncognitasNovo)
  
  kNovo <- 1
  for (incNovo in vetIncognitasNovo) {
    kAntigo <- 1
    for (incAntigo in vetIncognitas) {
      if (incNovo == incAntigo) {
        vetXNovo[kNovo] <- vetX[kAntigo]
        break()
      }
      kAntigo <- kAntigo + 1
    }
    kNovo <- kNovo + 1
  }
  
  return(list(vetXNovo, vetEquacoesNovo, vetIncognitasNovo, vetPotenciasNovo))
}

# RESOLUÇÃO -------------------------------------------------------------------------

#Escolhe o sistema para leitura
#sistema <- "2barras"
#sistema <- "3barrasteste"
#sistema <- "3barrastrabalho"
sistema <- "14barras"

#erro aceitavel
tolerancia <- 0.0001

#executa controle de reativos
controleReativo <- FALSE

#executa controle de tensao
controleTensao <- FALSE
vmax <- 1.2
vmin <- 0.8

#Executa a leitura dos arquivos de barras e linhas do sistema
barras <- read.table(paste("dados/",sistema,"_barras.csv",sep=""), header=T, sep=";")

linhas <- read.table(paste("dados/",sistema,"_linhas.csv",sep=""), header=T, sep=";")

nbarras <- dim(barras)[1]
nlinhas <- dim(linhas)[1]

#Adiciona coluna de flag para controle de reativos
if (controleReativo == TRUE) {
  barras <- cbind(barras, rep(0, nbarras))
  colnames(barras) <- c("barra","tipo","p","q","v","theta","bshk","pg","qg","pc","qc","qmin","qmax","limites")
}

#PASSO 1: CALCULAR AS ADMITANCIAS DO SISTEMA
linhas <- cbind(linhas, calcAdmitancias(linhas))

#PASSO 2: CALCULAR AS MATRIZES DE CONDUTANCIAS(G) E SUSCEPTANCIAS(B)
matY <- calcMatrizesGB(linhas, barras)
matG <- matY[[1]]
matB <- matY[[2]]
rm(matY)

#PASSO 3: IDENTIFICAR EQUACOES E INCOGNITAS DO PROBLEMA
eqincpot <- identificaEqIncPot(barras)
vetEquacoes <- eqincpot[[1]]
vetIncognitas <- eqincpot[[2]]
vetPotencias <- eqincpot[[3]]
rm(eqincpot)

#PASSO 4: CALCULAR AS EQUACOES DE POTENCIAS ATIVA E REATIVA
#preenchendo o vetor x com os chutes iniciais
vetX <- geraVetorXInicial(vetIncognitas)

#calculando o valor das potencias P e Q (vetor F)
vetF <- calculaPotencias(vetEquacoes, vetIncognitas, vetX, barras, linhas, matG, matB)

#AVALIA A CONVERGENCIA
cont <- 1
while (max(abs(vetF), na.rm = TRUE)>tolerancia) {
  
  #PASSO 5A: MONTANDO AS EQUACOES DA MATRIZ JACOBIANA (J) COM CONTROLE DE REATIVOS
  if (controleReativo == FALSE) {
    jacN <- montaJacobianaNominal(vetEquacoes, vetIncognitas, linhas)
  } else {
    #averigua violacao de reativos
    violacaoReativos <- violaReativos(barras)
    
    #no caso de haver violacao entao atualiza a matriz de barras
    if (violacaoReativos[[1]] == TRUE) {
      barras <- violacaoReativos[[2]]
      
      #atualiza o problema com base na violacao de reativos
      atualizacao <- atualizaProblema(barras, vetX, vetIncognitas, vetEquacoes)
      vetX <- atualizacao[[1]]
      vetEquacoes <- atualizacao[[2]]
      vetIncognitas <- atualizacao[[3]]
      vetPotencias <- atualizacao[[4]]

      #atualiza o vetF (potencias)
      vetF <- calculaPotencias(vetEquacoes, vetIncognitas, vetX, barras, linhas, matG, matB)
    }
    
    #atualiza a jacobiana
    jacN <- montaJacobianaNominal(vetEquacoes, vetIncognitas, linhas)
  }
  
  #PASSO 5B: CALCULAR A MATRIZ JACOBIANA (J)
  matJ <- calculaJacobiana(jacN, barras, linhas, matG, matB, vetX)

  #PASSO 6: APLICAR O MÉTODO DE NEWTON-RAPHSON
  deltaX <- Ginv(matJ)%*%vetF
  vetX <- vetX + deltaX
  
  #REcalculando o valor das potencias P e Q com novos valores de X (vetor F)
  vetF <- calculaPotencias(vetEquacoes, vetIncognitas, vetX, barras, linhas, matG, matB)
  
  #controle de tensao
  if (controleTensao == TRUE) {
    #averigua violacao de tensao
    violacaoTensao <- violaTensao(barras, vmax, vmin)
    
    #no caso de haver violacao entao atualiza a matriz de barras
    if (violacaoTensao[[1]] == TRUE) {
      barras <- violacaoTensao[[2]]
      
      #atualiza o problema com base na violacao de reativos
      atualizacao <- atualizaProblema(barras, vetX, vetIncognitas, vetEquacoes)
      vetX <- atualizacao[[1]]
      vetEquacoes <- atualizacao[[2]]
      vetIncognitas <- atualizacao[[3]]
      vetPotencias <- atualizacao[[4]]
      
      #atualiza o vetF (potencias)
      vetF <- calculaPotencias(vetEquacoes, vetIncognitas, vetX, barras, linhas, matG, matB)
    }
    
    #atualiza a jacobiana
    jacN <- montaJacobianaNominal(vetEquacoes, vetIncognitas, linhas)
  }
  
  #atualiza o contador de iteracoes
  cont <- cont + 1
}

#PASSO 7: CALCULAR AS POTENCIAS P e Q NAS BARRAS VT e PV (faltantes)
vetP <- calculaPotenciaBarrasFinal(vetPotencias, vetX, barras, linhas, matG, matB)

#atualiza as matrizes de linhas e barras
barrasFinal <- atualizaMatrizBarras(barras, vetPotencias, vetP, vetIncognitas, vetX)

#PASSO 8: CALCULAR AS CORRENTES ENTRE BARRAS
linhasFinal <- calculaCorrentesEntreBarras(matG, matB, barrasFinal, linhas)

#PASSO 9: CALCULAR AS POTENCIAS ENTRE BARRAS
linhasFinal <- calculaPotenciasEntreBarras(linhasFinal, barrasFinal)

#PASSSO 10: MOSTRA RESULTADOS
print("MATRIZ G:")
print(matG)
print("MATRIZ B:")
print(matB)
print("EQUACOES DO PROBLEMA:")
print(vetEquacoes)
print("INCOGNITAS DO PROBLEMA:")
print(vetIncognitas)
print("POTENCIAS FALTANTES:")
print(vetPotencias)
print("RAIZES ENCONTRADAS NA ULTIMA ITERACAO:")
print(vetF)
print("ELEMENTOS DA JACOBIANA:")
print(jacN)
print("ULTIMA JACOBIANA CALCULADA:")
print(matJ)
print("ULTIMO DELTA:")
print(deltaX)
print("INCOGNITAS RESOLVIDAS:")
print(vetX)
print("POTENCIAS NAS BARRAS FALTANTES:")
print(vetP)
print("MATRIZ DE BARRAS FINAL:")
print(barrasFinal)
print("MATRIZ DE LINHAS FINAL:")
print(linhasFinal)
print("TOTAL DE ITERACOES:")
print(cont)

