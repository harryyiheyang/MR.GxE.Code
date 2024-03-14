library(foreach)
library(doParallel)
library(parallel)
source("basicfunction.R")
(numCores <- detectCores()/2)
cl <- makeCluster(numCores)
registerDoParallel(cl)
(Miter=floor(100/numCores))

h=.1;
N1=c(1:160000);nn1=length(N1);N=nn1;
N2=c(1:80000)
M=200
mecc=0.1
ME=seq(0.05,0.5,by=0.05)
P=array(0,c(Miter*numCores,3,10))

for(k in 10:1){

for(iters in 1:Miter){
results <- foreach(i=c(1:(numCores)),.errorhandling="pass",.combine="c") %dopar% {
MAF=runif(M,0.05,0.5)
MAF[1]=0.3
G=matrix(NA,N,M);
for(m in 1:M){
g=rbinom(N,2,MAF[m])
G[,m]=g/sqrt(2*MAF[m]*(1-MAF[m]))
} 
beta=rnorm(M,0,sqrt(h/M))
e=rbinom(length(N2),1,ME[k])
e=e/sqrt(mean(e)*(1-mean(e)))
mug=G%*%beta;mug=as.vector(mug)
rhog=h/var(mug)
beta=beta*sqrt(rhog);
mug=G%*%beta;mug=as.vector(mug)
mue=e*0.1;mui1=G[,1]*e;
beta3=0.005
mui1=beta3*mui1;
mu=mug+mue+mui1
sig=var(as.matrix(G%*%beta))/h-var(mu);sig=sig[1,1]
y=mu+rnorm(N,0,1)*sqrt(sig)
B=matrix(NA,M,6)
y1=y[N1];G1=G[N1,];e1=e[N1]
y2=y[N2];G2=G[N2,];e2=e[N2]
remove(G)
fit0=biggwas(y1,G1)
B[,1]=as.vector(fit0$est)
B[,4]=as.vector(fit0$std)
for(i in 1:M){
a2=lm(y2~G2[,M-i+1]*e2)
B[M-i+1,2]=a2$coefficients[2]
B[M-i+1,3]=a2$coefficients[4]
B[M-i+1,5]=summary(a2)$coefficients[2,2]
B[M-i+1,6]=summary(a2)$coefficients[4,2]
}
fit1=lm(B[3:M,1]~B[3:M,2])
hattheta=fit1$coefficients[-1]
hatbeta3=B[2,3]
hatpleio=(B[2,1]-B[2,3]*hattheta)/mean(e)
return(c(hattheta,hatbeta3,hatpleio))
}
B=matrix(results,3,numCores)
B=t(B)
ind=c(((iters-1)*numCores+1):(iters*numCores))
P[ind,,k]=B
print(c(iters,k,colMeans(P[1:max(ind),,k])))
}
#save.image(filename)
}

