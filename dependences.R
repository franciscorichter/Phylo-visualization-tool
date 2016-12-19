list.of.packages <- c('phytools', "shiny", 'ape', 'subplex')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(shiny)
library(phytools)
library(ape)
library(subplex)

#' Run a model
#' @param model the type of model, must
#'   be 'Diversity-dependence' or 'Protracted'
#' @param tt crown age time, in million years
#' @param lambda0 ?birth rate, per million years per lineage
#' @param K carrying capacity, in number of species
#' @param draw ?no indea
#' @param printEv ?no idea
#' @param seed random number generator seed, can be any integer number
phyl2 <- function(
  model,
  tt = 15,
  lambda0=0.8,
  mu0=0.1,
  K=40,
  draw=TRUE,
  printEv=FALSE,
  seed=1
){
  if (model != "Diversity-dependence") {
    stop("Incorrect model")
  }
  # Phylogenetic tree simulation
  set.seed(seed)
  reboot = 0
  N = 1 # Number of species
  i = 1
  Tm = NULL
  sumt = 0
  sigma = 0
  E = NULL # vector with 0 if extinction and 1 if speciation
  n = NULL # vector with number of species at time t_i
  newick = paste(sl[1],";",sep="")  # Newick tree
  identf = data.frame(Spec="aa",Time=0)
  while (sumt<tt){
    if (model == "Diversity-dependence"){
      lambda = max(0,lambda0 - (lambda0-mu0)*N/K)
      mu = mu0
      lambda = rep(lambda,N)
      mu = rep(mu,N)
    }
    s = sum(lambda)+sum(mu)
    if (s == 0){break}
    tm = rexp(1,s)  # waiting time of iteration i
    if(tm+sumt>tt){break}
    sumt = tm + sumt
    prob = c(lambda,mu)/s  # Probability of extinctions and speciations
    BD = sample(2*N,1,prob=prob)  # speciation/extinction & identification of the species.
    n[i] = N
    if(BD > N){   # Extinction
      E[i] = 0
      ## for newick output
      species = identf[BD-N,1]
      ind = regexpr(species,newick)[1] + 2
      atm=sumt-identf[which(identf[,1]==species),2]
      identf = identf[-(BD-N),]
      newick = paste(substr(newick,1,ind),as.character(atm),substring(newick,ind+2),sep="")
      #
      N = N-1
      if(printEv){print(paste("extinction in time",sumt, sep=" "))}
    }else{  # Speciation
      E[i] = 1
      ## for newick output
      species = as.character(identf[BD,1])
      ind = regexpr(species,newick)[1]-1
      atm=sumt-identf[which(identf[,1]==species),2]
      newick = paste(substr(newick,1,ind),"(",substr(newick,ind+1,ind+4),",",sl[i+1],"):",as.character(atm),substring(newick,ind+5),sep="")
      identf = rbind(identf,data.frame(Spec=substr(sl[i+1],1,2),Time=sumt))
      identf[identf$Spec == species,2] = sumt
      #
      N = N+1
      if(printEv){print(paste("speciation in time",sumt,sep=" "))}
    }
    if (N==0){ # In case all species got extinct: restart
      reboot = reboot + 1
      N = 1 # Number of species
      i = 1
      Tm = NULL
      sumt = 0
      sigma = 0
      E = NULL # vector with 0 if extinction and 1 if speciation
      n = NULL # vector with number of species at time t_i
      newick = paste(sl[1],";",sep="")  # Newick tree
      identf = data.frame(Spec="aa",Time=0)
    }else { # Otherwise, update values and go to next iteration
      Tm[i] = tm
      sigma[i] = s
      i<-i+1
    }
  }
  vals = data.frame(time=cumsum(Tm),n=n)
  newick = compphyl(newi=newick,identf=identf,sumt=sumt)
  newick = read.tree(text=newick)
  treeD = list(t=Tm, E=E, r=reboot, i=i, n=n, vals=vals, newick=newick)
}

llik = function(b,n,E,t){
  # lo-likelihood function
  sigma = n*(b[1]-b[2]*n + b[3]) #n-dimentional
  rho = pmax(b[1]*E-b[2]*n*E+b[3]*(1-E),0)
  l = -sum(-sigma*t+log(rho))
  if(min(b)<0){l = Inf}
  return(l)
}

# this is just a verctor with labels
sl = paste(letters[1],letters,":0",sep="")
for (i in 2:26){
  ll = paste(letters[i],letters,":0",sep="")
  sl = c(sl,ll)
}

compphyl <- function(newi,identf,sumt){
  #test if the newick phylo is consistent with the identif matrix
  identf[,1]=as.character(identf[,1])
  identf[,2]=sumt-identf[,2]
  for(i in 1:length(identf[,1])){
    ind = regexpr(identf[i,1],newi)[1] + 2
    newi = paste(substr(newi,1,ind),as.character(identf[i,2]),substring(newi,ind+2),sep="")
  }
  return(newi)
}
