rm(list=ls())
cat("\014")
library(matlib)
source(file="funcoes_fluxoAC.R")

#Escolhe o sistema para leitura
#sistema <- "2barras"
#sistema <- "3barrasteste"
sistema <- "3barrastrabalho"

#erro aceitavel
tolerancia <- 0.003
  
#Executa a leitura dos arquivos de barras e linhas do sistema
barras <- read.table(paste("dados/",sistema,"_barras.csv",sep=""), header=T, sep=";")
linhas <- read.table(paste("dados/",sistema,"_linhas.csv",sep=""), header=T, sep=";")

nbarras <- dim(barras)[1]
nlinhas <- dim(linhas)[1]

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

#PASSO 5A: MONTANDO AS EQUACOES DA MATRIZ JACOBIANA (J)
jacN <- montaJacobianaNominal(vetEquacoes, vetIncognitas, linhas)

#AVALIA A CONVERGENCIA
cont <- 1
while (max(abs(vetF))>tolerancia) {
  
  #PASSO 5B: CALCULAR A MATRIZ JACOBIANA (J)
  matJ <- calculaJacobiana(jacN, barras, linhas, matG, matB, vetX)

  #PASSO 6: APLICAR O MÉTODO DE NEWTON-RAPHSON
  deltaX <- Ginv(matJ)%*%vetF
  vetX <- vetX + deltaX
  
  #REcalculando o valor das potencias P e Q com novos valores de X (vetor F)
  vetF <- calculaPotencias(vetEquacoes, vetIncognitas, vetX, barras, linhas, matG, matB)
  
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