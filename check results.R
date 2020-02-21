#look at theta
theta=res$theta[nrow(res$theta),]
plot(theta,type='h')
sum(theta>0.01)

#look at z
z=res$z[nrow(res$z),]
k=table(z.true,z);k
ind=numeric()
for (i in 1:nrow(k)){
  tmp=which(k[i,]==max(k[i,]))
  ind=c(ind,tmp)
}
ind
k[,ind]

#look at phi
nspp=ncol(dat)
phi=matrix(res$phi[nrow(res$phi),],50,nspp)
plot(phi.true,phi[ind,])