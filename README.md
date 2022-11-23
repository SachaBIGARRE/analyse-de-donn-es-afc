# analyse-de-donn-es-afc
---
title: "AFC_dm"
author: "Sacha BIgarré"
date: "16/11/2022"
output: html_document
---

```{r setup, include=FALSE}
getwd()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
bordeaux <- read.csv("bordeaux.csv", sep = ";")
AFC <- function(x){
  Nind = sum(X)
  #calcul des métriques
  Di <- diag(rowSums(prop.table(X)))
  Dj <- diag(colSums(prop.table(X)))
  #tableau des écarts à l'indépendance
  vec_P <- unlist(prop.table(X))
  P <- matrix(vec_P, nrow(X), ncol(X))
  Ecart_indp <- (solve(Di) %*% P %*% solve(Dj)) - 1
  #recherche de X
  X <- Nind*(Di %*% (Ecart_indp+1) %*% Dj)
  #calcul de la matrice V
  V <- t(X) %*% Di %*% X
  #vecteurs propres et mise en place de la matrice normée A
  VDj <- V %*% Dj
  eigen(VDj) # decomposition en vecteurs et valeurs propres
  # on norme avec (Dj'Dj)-1/2
  normDj=(t(Dj)%*%Dj)^(-1/2)
  A <- eigen(VDj)$vectors%*% normDj
  #vecteurs propres DI-normés à l'unité de WDI
  W <- X %*% Dj %*% t(X)
  WDi <- W %*% Di
  eigen(WDi) # decomposition en vecteurs et valeurs propres
  C <- eigen(WDi)$vectors/norm(Di,type="O")
  #vecteur du pourcentage d'inertie de chaque axe
  inertie <- Ecart_indp %*% Dj %*% t(Ecart_indp) %*% Di
  vec_prct_inertie <- as.vector(eigen(inertie)$values*100)
  return(vec_prct_inertie)
}
barplot(AFC(bordeaux))
```
#################################################






# Analyse de données, devoir maison - SUJET 1
```{r selection du repertoire de travail, include=F}
setwd("D:/master 1 - MAS/M1_MAS_analyse_de_donnees/projet")
```

on crée une fonction qui prend un taleau de contingeance en entrée et

(a) calcul les metriques Di et Dj
```{r centrage réduction}
cr=function(X){
  n=nrow(X)
  
  # centrage, on soustrait la moyenne de la colonne de l'élément
  Xmean=apply(X,MARGIN=2,FUN=mean)
  Xmatmean=matrix(rep(Xmean,n),nrow = n,byrow = T)
  Xc=X- Xmatmean
  
  # reduction, on divise par l'écart type de la colonne de l'élément
  Xcet=sqrt(apply(Xc, MARGIN=2, var))
  Xc_matmean=matrix(rep(Xcet,n),nrow = n,byrow = T)
  Xcr=Xc/Xc_matmean  
  return(Xcr)
}
```


```{r}

```

```{r}

f = function(N){
  # N=as.matrix(N)
  #(a) calcul des métriques
  n=sum(N) # nombre d'element dans le tableau de contingeance
  
  Pi.=rowSums(N)/n # proba théorique de la ligne/modalité i
  P.j=colSums(N)/n # proba théorique de la colonne/modalité j
  P = N/n
  Di.=diag(nrow(N))*Pi. # calcul de la métrique Di.
  D.j=diag(ncol(N))*P.j # calcul de la métrique D.j

  #(b) tableau des écarts à l'indépendance


  Xtild = solve(Di.) %*% P %*% solve(D.j) # tableau dérivé
  Z=Xtild-1

  # centrage réduction de Z : ( on a créé la fonction centrage-réduction avant : cr() )
  Zcr=scale(Z)
  
#   # calcul de VDj et des valeurs propres
  V = t(Zcr)%*%Di.%*%Zcr
  VD.j=D.j%*%V
  
#   
  eigV=eigen(VD.j)
  
  # on norme par D.j 
  # pourtant on a deja normé V par D.J ??
  A=eigV$vectors/D.j# non !!
  
  print(eigV$vectors)
  print(D.j)
  C=(1/eigV$values)%*%Zcr%*%A
  sigma= eigV$values/tr(VD.j)
  ind_coord=Zcr%*%VD.j
  var_coord=C%*%D.j%*%Zcr
#   
#   return( sigma, A, ind_coord, C, var_coord)

}
f(mbor)

```


```{r TEST import, include=F}
readLines("bordeaux.csv",10)
bord=read.csv("bordeaux.csv",sep=";",header = T,row.names = 1)
```
