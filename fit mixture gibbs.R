rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(11)

setwd('U:\\GIT_models\\mixture_binom')
source('mixture gibbs functions.R')
sourceCpp('aux1.cpp')
source('mixture gibbs main function.R')

#get data
nl=read.csv('fake data mixture nl.csv',as.is=T)$x
y=read.csv('fake data mixture y.csv',as.is=T)
dat=data.matrix(y)

ngibbs=1000
burnin=ngibbs/2
ngroup=50
a.prior=b.prior=1

res=mixture.gibbs.main.func(dat=dat,nl=nl,ngroup=ngroup,ngibbs=ngibbs,burnin=burnin,
                            a.prior=a.prior,b.prior=b.prior)

    
