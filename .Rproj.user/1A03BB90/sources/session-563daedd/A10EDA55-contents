library(foreach)
library(doParallel)
library(parallel)
source("basicfunction.R")
(numCores <- detectCores()/2)
cl <- makeCluster(numCores)
registerDoParallel(cl)
(Miter=floor(100/numCores))

N1=c(1:160000);N2=c(1:80000);n1=length(N1);n2=length(N2);no=80000;N=max(N1,N2)
h=.1;ME=seq(0.05,0.5,by=0.05)
POW1=POW2=matrix(NA,10,2)
M=200
mecc=0.1

ssign=1 ### It can also be -1 and 0
for(k in 10:1){
P=matrix(NA,Miter*numCores,4)
for(iters in 1:Miter){
results <- foreach(i=c(1:(numCores)),.errorhandling="pass",.combine = "c") %dopar% {
MAF=runif(M,0.05,0.5)
MAF[1]=0.3
G=matrix(NA,N,M);
for(m in 1:M){
g=rbinom(N,2,MAF[m])
G[,m]=g/sd(g)
} 
beta=rnorm(M,0,sqrt(h/M))
beta[1]=2*sqrt(h/M)*ssign
#######################################
e=rnorm(N,0,1)*0
e[N2]=rbinom(length(N2),1,mecc)
e[setdiff(N1,N2)]=rbinom(length(setdiff(N1,N2)),1,(ME[k]*2-mecc))
e=e/sd(e)
#######################################
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
y=y/sqrt(sig)
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
pv1=1-pchisq(abs(B[1,3])^2/B[1,6]^2,1)
pv2=1-pchisq(abs(B[2,3])^2/B[2,6]^2,1)

rho0=no/sqrt(n1)/sqrt(n2)/sqrt(1+mean(e2)^2)
fit1=lm(B[3:M,1]~B[3:M,2])
hattheta=fit1$coefficients[2]
hatthetase=summary(fit1)$coefficients[2,2]
betahat=var(B[,2])-mean(B[,5]^2)
r1=B[1,1]-hattheta*B[1,2]
r2=B[2,1]-hattheta*B[2,2]
pvv=t_pleio(bx=B[,2],by=B[,1],bx_se=B[,5],by_se=B[,4],sig=betahat,betahat=hattheta,betahat_se=hatthetase,rho=rho0)
pv3=pvv[1];pv4=pvv[2]
a=c(as.numeric(pv1<0.05),as.numeric(pv3<0.05),as.numeric(pv2<0.05),as.numeric(pv4<0.05))
return(a)
}
B=t(matrix(results,4,numCores))
ind=c(((iters-1)*numCores+1):(iters*numCores))
P[ind,]=B
print(c(iters,k,colMeans(P[1:max(ind),])))
}
POW1[k,]=colMeans(P[,1:2])
POW2[k,]=colMeans(P[,3:4])
#save.image("~/BD/fix1/var1-positive-mecc2-001-2020-sub.RData")
}
