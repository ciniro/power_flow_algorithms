#-----------------------------------------------------------------------------#
#                                                                             #
#                FUNCOES PARA CALCULO DO FLUXO DE POTENCIA AC                 #
#          DISCIPLINA: ANÁLISE ESTÁTICA DE SISTEMAS DE ENERGIA ELÉTRICA       #
#                                                                             #
#                     CINIRO APARECIDO LEITE NAMETALA                         #
#-----------------------------------------------------------------------------#

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
    k <- substr(iinc, 2, 2)
    
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
    k <- substr(ieq, 2, 2)
    
    vetF <- c(vetF, calcDeltaPotencia(tipo, k, x, barras, linhas, matG, matB))
  }
  
  return(vetF)
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
      k <- substr(eq, 2, 2)
      m <- substr(inc, 2, 2)
      indicesBarras <- paste(k, m, sep="")
      
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
        k <- as.numeric(substr(matJ[i,j], 2, 2))
        m <- as.numeric(substr(matJ[i,j], 3, 3))
        
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
    k <- substr(ieq, 2, 2)
    
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
             #quadrante H
             "P"={tipo <- "p"},
             "Q"={tipo <- "q"}
             )
               
      k <- as.numeric(substr(ipot, 2, 2))
      
      novoBarras[k, tipo] <- vetP[i]
      i <- i + 1
  }
  
  i <- 1
  for (iinc in vetIncognitas) {
    tipo <- substr(iinc, 1, 1)
    
    switch(tipo,
           #quadrante H
           "T"={tipo <- "theta"},
           "V"={tipo <- "v"}
    )
    
    k <- as.numeric(substr(iinc, 2, 2))
    
    novoBarras[k, tipo] <- vetX[i]
    i <- i + 1
  }
  
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
