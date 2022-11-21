# analyse-de-donn-es-afc
---
title: "AFC_dm"
author: "Sacha BIgarré"
date: "16/11/2022"
output: html_document
---
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
