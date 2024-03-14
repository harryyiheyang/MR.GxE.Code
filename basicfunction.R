
biggwas=function(x,G){
  x=as.vector(x)
  ux=mean(x)
  vx=var(x);vx=as.numeric(vx)
  ug=colMeans(G)
  G=t(t(G)-ug)
  vg=colSums(G^2)
  b=(t(G)%*%(x-ux))/vg
  sdb=(vx-b^2*vg/length(x))/length(x)
  A=list()
  A$est=b
  A$std=sqrt(sdb)
  return(A)
}

demean=function(A){
  p=dim(A)[2]
  for(i in 1:p){
    A[,i]=A[,i]-mean(A[,i])  
  }
  return(A)
}

simulationp=function(mu,S,alpha=0.05,k=10000){
  a=mvrnorm(k,mu,S);a=a^2/diag(S)
  pv=pchisq(a,1)
  pv1=apply(pv,MARGIN=1,FUN=min)
  thres=quantile(pv1,alpha)
}

optimalshrinkage=function(a,b,C,k){
  w=c(0,c(1:(k-1))/(k-1))
  pv=c(1:k)
  for(i in 1:k){
    s=a*w[i]+b*(1-w)[i]
    v=t(c(w[i],(1-w[i])))%*%C%*%c(w[i],(1-w[i]))
    pv[i]=pchisq(s^2/v,1,lower.tail=F)
  }
  return(min(pv))
}

rlm=function(X,y,yse,q,k=5){
  w=1/yse^2
  for(i in 1:k){
    fit=lm(y~X,weight=w)
    sig=1.48*median(abs(fit$residuals))
    w=pnorm(fit$residuals/sqrt(sig*yse));w=w^(1-q)/yse^2
  }
  fit$weight=w
  return(fit)
}

t_pleio=function(bx, by, bx_se, by_se, betahat, betahat_se, rho,sig) {
  T_var=by_se^2 + betahat^2*bx_se^2 - 2*betahat*rho*by_se*bx_se+(1-betahat)^2*sig
  T_var=abs(T_var)
  .T=(by-betahat*bx)^2/T_var
  T_p=1-pchisq(.T,1)
  T_p}

estimateh=function(n1,n2,no,mue,theta){
  theta=1-theta
  x=mue^2+1-no/n1
  y=mue^2+1
  h=x/t-y
  h=h/n2
  return(h)
}
