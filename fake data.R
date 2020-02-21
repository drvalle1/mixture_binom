rm(list=ls(all=TRUE))
set.seed(3)

#general settings
nloc=1000
nspp=50
ngroup=8

#set parameters
z.true=z=sample(1:ngroup,size=nloc,replace=T)
phi.true=phi=matrix(rbeta(ngroup*nspp,1,1),ngroup,nspp)
nl=rpois(nloc,lambda=1)+1; table(nl)

#generate data
y=matrix(NA,nloc,nspp)
for (i in 1:nloc){
  phi1=phi[z[i],]
  y[i,]=rbinom(nspp,size=nl[i],prob=phi1)
}
colnames(y)=paste0('spp',1:nspp)
rownames(y)=paste0('loc',1:nloc)

#export results
setwd('U:\\GIT_models\\mixture_binom')
write.csv(y,'fake data mixture y.csv',row.names=F)
write.csv(nl,'fake data mixture nl.csv',row.names=F)

