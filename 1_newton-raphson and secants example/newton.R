rm(list=ls())
cat("\014")
library(compiler)
library(numDeriv)

#-----------------------------------------------------------------------------#
#                                                                             #
#                           MÉTODO DE NEWTON-RAPHSON                          #
#                                                                             #
#                     CINIRO APARECIDO LEITE NAMETALA                         #
#-----------------------------------------------------------------------------#

newtonx <- function(f, a, b, tol = 1e-10, n = 1000)
{
  #avalia se a raiz esta no limite inferior
  fa <- f(a)
  if (fa == 0.0) {
    return(a)
  }
  
  #avalia se a raiz esta no limite superior
  fb <- f(b)
  if (fb == 0.0) {
    return(b)
  }
  
  xcurrent <- a
  k <- n
  
  for (i in 1:n) {
    #Derivada de primeira ordem f'(xcurrent)
    dx <- genD(func = f, x = xcurrent)$D[1]
    #Calcula o proximo valor x (xnext)
    xnext <- xcurrent - (f(xcurrent) / dx)
    #Guarda xnext
    k[i] <- xnext
    
    #avalia a convergencia por erro
    #(abs(xnext - xcurrent) < tol) 
    if (abs(f(xnext)) < tol) {
      res <- list('Raiz' = tail(k, n=1), 'Iteracoes' = k)
      return(res)
    }
    
    #Se nao convergir atualiza o xcurrent
    xcurrent <- xnext
  }
  
  print('Nao convergiu!')
}

secantex <- function(f, a, b, tol = 1e-5, n = 1000)
{
  #avalia se a raiz esta no limite inferior
  fa <- f(a)
  if (fa == 0.0) {
    return(a)
  }
  
  #avalia se a raiz esta no limite superior
  fb <- f(b)
  if (fb == 0.0) {
    return(b)
  }
  
  xprevious <- a
  xcurrent <- b
  
  k <- n
  
  for (i in 1:n) {
    #Calcula o proximo valor x (xnext)
    xnext <- xcurrent - ((f(xcurrent)*(xcurrent-xprevious))/(f(xcurrent)-f(xprevious)))
    #Guarda xnext
    k[i] <- xnext
    
    #avalia a convergencia por erro
    if (abs(f(xnext)) < tol) {
      res <- list('Raiz' = tail(k, n=1), 'Iteracoes' = k)
      return(res)
    }
    
    #Se nao convergir atualiza o xcurrent
    xprevious <- xcurrent
    xcurrent <- xnext
  }
  
  print('Nao convergiu!')
}

funcaox <- function(theta)
{
  return(0.1923 * (1-cos(theta)) + 0.9615 * sin(theta))
}

funcao2x <- function(x) {
  x^3 - 2* x - 5
}

funcao3x <- function(x) {
  exp(2 * x) - x - 6
}

funcao <- compiler::cmpfun(funcaox)
funcao2 <- compiler::cmpfun(funcao2x)
funcao3 <- compiler::cmpfun(funcao3x)
newton <- compiler::cmpfun(newtonx)
secante <- compiler::cmpfun(secantex)

#-------------
f <- funcao
linf <- -1
lsup <- 1

#---------------------------
#CALCULO COM NEWTON-RAPHSON
#---------------------------

raiz <- newton(f, linf, lsup)

print("RESULTADO COM NEWTON-RAPHSON")
print(raiz)

x <- seq(-4,4,0.1)
plot(x, f(x), pch=19, type="l", col="blue", xlim=c(-4,4), main="Raiz da Função por Newton-Raphson", cex.lab=0.7, cex.main=0.9, cex.sub=0.7, xaxt="n", yaxt="n")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
axis(side=1, at=c(-4:4), cex.axis=0.7)
axis(side=2, at=c(-1:1), cex.axis=0.7)
abline(h=0, lty=5)
abline(v=raiz$Iteracoes)
points(raiz$Iteracoes, rep(0,length(raiz$Iteracoes)), pch=19, col="red")
points(raiz$Raiz, 0, pch=19, col="forestgreen")

#---------------------------
#CALCULO COM MÉTODO DAS SECANTES
#---------------------------
raiz <- secante(f, linf, lsup)

print("RESULTADO COM SECANTES")
print(raiz)

x <- seq(-4,4,0.1)
plot(x, f(x), pch=19, type="l", col="blue", xlim=c(-4,4), main="Raiz da Função por Secantes", cex.lab=0.7, cex.main=0.9, cex.sub=0.7, xaxt="n", yaxt="n")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
axis(side=1, at=c(-4:4), cex.axis=0.7)
axis(side=2, at=c(-1:1), cex.axis=0.7)
abline(h=0, lty=5)
abline(v=raiz$Iteracoes)
points(raiz$Iteracoes, rep(0,length(raiz$Iteracoes)), pch=19, col="red")
points(raiz$Raiz, 0, pch=19, col="forestgreen")

