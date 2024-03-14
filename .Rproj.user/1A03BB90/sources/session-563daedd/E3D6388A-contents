library(ggplot2)
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
library(ggpubr)
## This program is a simulation for generating the Fig 2.E and Fig 2.F, as well as Supplementary Fig. S6
####### no mediation, environment contribution 1%
n2=20000
B1=matrix(0,4,3)
Theta1=matrix(0,1000,8)
### Fig 1, mue=1, FigS7, mue=0.5
mue=1    
p=0.3
mug=2*p/0.648
rho1=0

for(i in 1:4){
  if (i==1) n1=20000
  else n1= 100000+(i-2)*100000
  N1=c(1:n1)
  N2=c(1:n2)
  PV=matrix(1000,1000,3)
  PM=matrix(20,20,2)
  for(j in 1:1000){
    for (k in 1:20){
      g=(rbinom(n1,2,p))/0.648
      e=rho1*g+rnorm(n1,2*mue,sqrt(1-rho1^2))
      e[N2]=rho1*g[N2]+rnorm(n2,mue,sqrt(1-rho1^2))
      y=0.1*g+e+rnorm(n1,0,10)
      fit1=lm(y~g)
      y2=y[N2];g2=g[N2];e2=e[N2]
      fit2=lm(y2~g2*e2)
      pv1=summary(fit2)$coefficient[4,4]
      r=summary(fit1)$coefficient[2,1]-summary(fit2)$coefficient[2,1]
      rho=sqrt(n2/n1)/sqrt(1./(1-rho1^2)+(rho1*(1-mug)+mue)^2)
      varr=summary(fit1)$coefficient[2,2]^2+summary(fit2)$coefficient[2,2]^2-2*rho*summary(fit1)$coefficient[2,2]*summary(fit2)$coefficient[2,2]
      pv2=1-pchisq(r^2/varr,1)
      PM[k,]=c(pv1,pv2)
    }
    a=apply(PM,2,min)
    PV[j,1:2]=c(a[1]<0.0025,a[2]<0.0025)
    
    k0=sum(PM[,2]<0.0025)
    
    PV[j,3]=0
    if (k0>0) {
      if (k0==1){
        a1=PM[PM[,2]<0.0025,]
        PV[j,3]=c(a[1]< 0.05/k0)
      } else {
        a1=PM[PM[,2]<0.0025,]
        a=apply(a1,2,min)
        PV[j,3]=c(a[1]<0.05/k0)
      }
    }
    
    if(j %% 100==0){print(c(i,j))}
    Theta1[j,i]=r
  }
  B1[i,]=colMeans(PV)
}


####### environment contribution 1%
n2=20000
B1=matrix(0,4,3)
Theta1=matrix(0,1000,8)


for(i in 1:4){
  if (i==1) n1=20000
  else n1= 100000+(i-2)*100000
  N1=c(1:n1)
  N2=c(1:n2)
  PV=matrix(1000,1000,3)
  PM=matrix(20,20,2)
  for(j in 1:1000){
    
    rho1=0.01
    g=(rbinom(n1,2,0.3))/0.648
    
    
    e=rho1*g+rnorm(n1,2*mue,sqrt(1-rho1^2))
    e[N2]=rho1*g[N2]+rnorm(n2,mue,sqrt(1-rho1^2))
    
    y=0.2*g+e+rnorm(n1,0,10)
    fit1=lm(y~g)
    y2=y[N2];g2=g[N2];e2=e[N2]
    fit2=lm(y2~g2*e2)
    pv1=summary(fit2)$coefficient[4,4]
    r=summary(fit1)$coefficient[2,1]-summary(fit2)$coefficient[2,1]
    rho=sqrt(n2/n1)/sqrt(1./(1-rho1^2)+(rho1*(1-mug)+mue)^2)
    #rho=sqrt(n2/n1)/sqrt(1.0+0.5^2)
    varr=summary(fit1)$coefficient[2,2]^2+summary(fit2)$coefficient[2,2]^2-2*rho*summary(fit1)$coefficient[2,2]*summary(fit2)$coefficient[2,2]
    pv2=1-pchisq(r^2/varr,1)
    PM[1,]=c(pv1,pv2)
    
    rho1=0
    for (k in 2:20){
      g=(rbinom(n1,2,0.3))/0.648
      
      
      e=rho1*g+rnorm(n1,2*mue,sqrt(1-rho1^2))
      e[N2]=rho1*g[N2]+rnorm(n2,mue,sqrt(1-rho1^2))
      
      y=0.2*g+e+rnorm(n1,0,10)
      fit1=lm(y~g)
      y2=y[N2];g2=g[N2];e2=e[N2]
      fit2=lm(y2~g2*e2)
      pv1=summary(fit2)$coefficient[4,4]
      r=summary(fit1)$coefficient[2,1]-summary(fit2)$coefficient[2,1]
      rho=sqrt(n2/n1)/sqrt(1./(1-rho1^2)+(rho1*(1-mug)+mue)^2)
      #rho=sqrt(n2/n1)/sqrt(1.+0.5^2)
      
      varr=summary(fit1)$coefficient[2,2]^2+summary(fit2)$coefficient[2,2]^2-2*rho*summary(fit1)$coefficient[2,2]*summary(fit2)$coefficient[2,2]
      pv2=1-pchisq(r^2/varr,1)
      PM[k,]=c(pv1,pv2)
    }
    
    
    a=apply(PM,2,min)
    PV[j,1:2]=c(a[1]<0.0025,a[2]<0.0025)
    
    k0=sum(PM[,2]<0.0025)
    
    PV[j,3]=0
    if (k0>0) {
      if (k0==1){
        a1=PM[PM[,2]<0.0025,]
        PV[j,3]=c(a[1]< 0.05/k0)
      } else {
        a1=PM[PM[,2]<0.0025,]
        a=apply(a1,2,min)
        PV[j,3]=c(a[1]<0.05/k0)
      }
    }
    
    if(j %% 100==0){print(c(i,j))}
    Theta1[j,i]=r
  }
  B1[i,]=colMeans(PV)
}

####### environment contribution 5%
n2=20000
B1=matrix(0,4,3)
Theta1=matrix(0,1000,8)

for(i in 1:4){
  if (i==1) n1=20000
  else n1= 100000+(i-2)*100000
  N1=c(1:n1)
  N2=c(1:n2)
  PV=matrix(1000,1000,3)
  PM=matrix(20,20,2)
  for(j in 1:1000){
    
    rho1=0.05
    g=(rbinom(n1,2,0.3))/0.648
    
    e=rho1*g+rnorm(n1,2*mue,sqrt(1-rho1^2))
    e[N2]=rho1*g[N2]+rnorm(n2,mue,sqrt(1-rho1^2))
    
    y=0.2*g+sqrt(5)*e+rnorm(n1,0,10)
    fit1=lm(y~g)
    y2=y[N2];g2=g[N2];e2=e[N2]
    fit2=lm(y2~g2*e2)
    pv1=summary(fit2)$coefficient[4,4]
    r=summary(fit1)$coefficient[2,1]-summary(fit2)$coefficient[2,1]
    rho=sqrt(n2/n1)/sqrt(1./(1-rho1^2)+(rho1*(1-mug)+mue)^2)
    
    
    varr=summary(fit1)$coefficient[2,2]^2+summary(fit2)$coefficient[2,2]^2-2*rho*summary(fit1)$coefficient[2,2]*summary(fit2)$coefficient[2,2]
    pv2=1-pchisq(r^2/varr,1)
    PM[1,]=c(pv1,pv2)
    
    rho1=0
    for (k in 2:20){
      g=(rbinom(n1,2,0.3))/0.648
      
      e=rho1*g+rnorm(n1,2*mue,sqrt(1-rho1^2))
      e[N2]=rho1*g[N2]+rnorm(n2,mue,sqrt(1-rho1^2))
      
      y=0.2*g+sqrt(5)*e+rnorm(n1,0,10)
      fit1=lm(y~g)
      y2=y[N2];g2=g[N2];e2=e[N2]
      fit2=lm(y2~g2*e2)
      pv1=summary(fit2)$coefficient[4,4]
      r=summary(fit1)$coefficient[2,1]-summary(fit2)$coefficient[2,1]
      rho=sqrt(n2/n1)/sqrt(1./(1-rho1^2)+(rho1*(1-mug)+mue)^2)
      #rho=sqrt(n2/n1)/sqrt(1.+0.5^2)
      varr=summary(fit1)$coefficient[2,2]^2+summary(fit2)$coefficient[2,2]^2-2*rho*summary(fit1)$coefficient[2,2]*summary(fit2)$coefficient[2,2]
      pv2=1-pchisq(r^2/varr,1)
      PM[k,]=c(pv1,pv2)
    }
    
    
    a=apply(PM,2,min)
    PV[j,1:2]=c(a[1]<0.0025,a[2]<0.0025)
    
    k0=sum(PM[,2]<0.0025)
    
    PV[j,3]=0
    if (k0>0) {
      if (k0==1){
        a1=PM[PM[,2]<0.0025,]
        PV[j,3]=c(a[1]< 0.05/k0)
      } else {
        a1=PM[PM[,2]<0.0025,]
        a=apply(a1,2,min)
        PV[j,3]=c(a[1]<0.05/k0)
      }
    }
    
    if(j %% 100==0){print(c(i,j))}
    Theta1[j,i]=r
  }
  B1[i,]=colMeans(PV)
}

################## Power analysis
####### no mediation, environment contribution 1%
n2=20000
B1=matrix(0,4,3)
Theta1=matrix(0,1000,8)
rho1=0

for(i in 1:4){
  if (i==1) n1=20000
  else n1= 100000+(i-2)*100000
  N1=c(1:n1)
  N2=c(1:n2)
  PV=matrix(1000,1000,3)
  PM=matrix(20,20,2)
  for(j in 1:1000){
    
    g=(rbinom(n1,2,0.3))/0.648
    
    e=rho1*g+rnorm(n1,2*mue,sqrt(1-rho1^2))
    e[N2]=rho1*g[N2]+rnorm(n2,mue,sqrt(1-rho1^2))
    
    y=0.1*g+e+g*e*0.1+rnorm(n1,0,10)
    fit1=lm(y~g)
    y2=y[N2];g2=g[N2];e2=e[N2]
    fit2=lm(y2~g2*e2)
    pv1=summary(fit2)$coefficient[4,4]
    r=summary(fit1)$coefficient[2,1]-summary(fit2)$coefficient[2,1]
    rho=sqrt(n2/n1)/sqrt(1./(1-rho1^2)+(rho1*(1-mug)+mue)^2)
    varr=summary(fit1)$coefficient[2,2]^2+summary(fit2)$coefficient[2,2]^2-2*rho*summary(fit1)$coefficient[2,2]*summary(fit2)$coefficient[2,2]
    pv2=1-pchisq(r^2/varr,1)
    PM[1,]=c(pv1,pv2)
    
    rho1=0
    for (k in 2:20){
      g=(rbinom(n1,2,0.3))/0.648
      
      e=rho1*g+rnorm(n1,2*mue,sqrt(1-rho1^2))
      e[N2]=rho1*g[N2]+rnorm(n2,mue,sqrt(1-rho1^2))
      
      
      y=0.1*g+e+rnorm(n1,0,10)
      fit1=lm(y~g)
      y2=y[N2];g2=g[N2];e2=e[N2]
      fit2=lm(y2~g2*e2)
      pv1=summary(fit2)$coefficient[4,4]
      r=summary(fit1)$coefficient[2,1]-summary(fit2)$coefficient[2,1]
      rho=sqrt(n2/n1)/sqrt(1./(1-rho1^2)+(rho1*(1-mug)+mue)^2)
      varr=summary(fit1)$coefficient[2,2]^2+summary(fit2)$coefficient[2,2]^2-2*rho*summary(fit1)$coefficient[2,2]*summary(fit2)$coefficient[2,2]
      pv2=1-pchisq(r^2/varr,1)
      PM[k,]=c(pv1,pv2)
    }
    
    
    a=apply(PM,2,min)
    PV[j,1:2]=c(a[1]<0.0025,a[2]<0.0025)
    
    k0=sum(PM[,2]<0.0025)
    
    PV[j,3]=0
    if (k0>0) {
      if (k0==1){
        a1=PM[PM[,2]<0.0025,]
        PV[j,3]=c(a[1]< 0.05/k0)
      } else {
        a1=PM[PM[,2]<0.0025,]
        a=apply(a1,2,min)
        PV[j,3]=c(a[1]<0.05/k0)
      }
    }
    
    if(j %% 100==0){print(c(i,j))}
    Theta1[j,i]=r
  }
  B1[i,]=colMeans(PV)
}

####### environment contribution 1%
n2=20000
B1=matrix(0,4,3)
Theta1=matrix(0,1000,8)


for(i in 1:4){
  if (i==1) n1=20000
  else n1= 100000+(i-2)*100000
  N1=c(1:n1)
  N2=c(1:n2)
  PV=matrix(1000,1000,3)
  PM=matrix(20,20,2)
  for(j in 1:1000){
    
    rho1=0.01
    g=(rbinom(n1,2,0.3))/0.648
    
    e=rho1*g+rnorm(n1,2*mue, sqrt(1-rho1^2))
    e[N2]=rho1*g[N2]+rnorm(n2,mue, sqrt(1-rho1^2))
    
    y=0.2*g+e+g*e*0.2+rnorm(n1,0,10)
    fit1=lm(y~g)
    y2=y[N2];g2=g[N2];e2=e[N2]
    fit2=lm(y2~g2*e2)
    pv1=summary(fit2)$coefficient[4,4]
    r=summary(fit1)$coefficient[2,1]-summary(fit2)$coefficient[2,1]
    
    rho=sqrt(n2/n1)/sqrt(1./(1-rho1^2)+(rho1*(1-mug)+mue)^2)
    varr=summary(fit1)$coefficient[2,2]^2+summary(fit2)$coefficient[2,2]^2-2*rho*summary(fit1)$coefficient[2,2]*summary(fit2)$coefficient[2,2]
    pv2=1-pchisq(r^2/varr,1)
    PM[1,]=c(pv1,pv2)
    
    rho1=0
    for (k in 2:20){
      g=(rbinom(n1,2,0.3))/0.648
      
      e=rho1*g+rnorm(n1,2*mue, sqrt(1-rho1^2))
      e[N2]=rho1*g[N2]+rnorm(n2,mue, sqrt(1-rho1^2))
      
      
      y=0.2*g+e+rnorm(n1,0,10)
      fit1=lm(y~g)
      y2=y[N2];g2=g[N2];e2=e[N2]
      fit2=lm(y2~g2*e2)
      pv1=summary(fit2)$coefficient[4,4]
      r=summary(fit1)$coefficient[2,1]-summary(fit2)$coefficient[2,1]
      
      rho=sqrt(n2/n1)/sqrt(1./(1-rho1^2)+(rho1*(1-mug)+mue)^2)
      
      varr=summary(fit1)$coefficient[2,2]^2+summary(fit2)$coefficient[2,2]^2-2*rho*summary(fit1)$coefficient[2,2]*summary(fit2)$coefficient[2,2]
      pv2=1-pchisq(r^2/varr,1)
      PM[k,]=c(pv1,pv2)
    }
    
    a=apply(PM,2,min)
    PV[j,1:2]=c(a[1]<0.0025,a[2]<0.0025)
    
    k0=sum(PM[,2]<0.0025)
    
    PV[j,3]=0
    if (k0>0) {
      if (k0==1){
        a1=PM[PM[,2]<0.0025,]
        PV[j,3]=c(a[1]< 0.05/k0)
      } else {
        a1=PM[PM[,2]<0.0025,]
        a=apply(a1,2,min)
        PV[j,3]=c(a[1]<0.05/k0)
      }
    }
    
    if(j %% 100==0){print(c(i,j))}
    Theta1[j,i]=r
  }
  B1[i,]=colMeans(PV)
}

####### environment contribution 5%
n2=20000
B1=matrix(0,4,3)
Theta1=matrix(0,1000,8)


for(i in 1:4){
  if (i==1) n1=20000
  else n1= 100000+(i-2)*100000
  N1=c(1:n1)
  N2=c(1:n2)
  PV=matrix(1000,1000,3)
  PM=matrix(20,20,2)
  for(j in 1:1000){
    
    rho1=0.05
    g=(rbinom(n1,2,0.3))/0.648
    
    e=rho1*g+rnorm(n1,2*mue, sqrt(1-rho1^2))
    
    e[N2]=rho1*g[N2]+rnorm(n2,mue, sqrt(1-rho1^2))
    mue1=mean(e[N2])
    
    y=0.1*g+sqrt(5)*e+g*e*0.2+rnorm(n1,0,10)
    fit1=lm(y~g)
    y2=y[N2];g2=g[N2];e2=e[N2]
    fit2=lm(y2~g2*e2)
    pv1=summary(fit2)$coefficient[4,4]
    r=summary(fit1)$coefficient[2,1]-summary(fit2)$coefficient[2,1]
    
    rho=sqrt(n2/n1)/sqrt(1./(1-rho1^2)+(rho1*(1-mug)+mue1)^2)
    
    varr=summary(fit1)$coefficient[2,2]^2+summary(fit2)$coefficient[2,2]^2-2*rho*summary(fit1)$coefficient[2,2]*summary(fit2)$coefficient[2,2]
    pv2=1-pchisq(r^2/varr,1)
    PM[1,]=c(pv1,pv2)
    rho1=0
    for (k in 2:20){
      g=(rbinom(n1,2,0.3))/0.648
      
      e=rho1*g+rnorm(n1,2*mue, sqrt(1-rho1^2))
      e[N2]=rho1*g[N2]+rnorm(n2,mue, sqrt(1-rho1^2))
      mue1=mean(e[N2])
      
      y=0.2*g+sqrt(5)*e+rnorm(n1,0,10)
      fit1=lm(y~g)
      y2=y[N2];g2=g[N2];e2=e[N2]
      fit2=lm(y2~g2*e2)
      pv1=summary(fit2)$coefficient[4,4]
      r=summary(fit1)$coefficient[2,1]-summary(fit2)$coefficient[2,1]
      
      rho=sqrt(n2/n1)/sqrt(1./(1-rho1^2)+(rho1*(1-mug)+mue1)^2)
      
      varr=summary(fit1)$coefficient[2,2]^2+summary(fit2)$coefficient[2,2]^2-2*rho*summary(fit1)$coefficient[2,2]*summary(fit2)$coefficient[2,2]
      pv2=1-pchisq(r^2/varr,1)
      PM[k,]=c(pv1,pv2)
    }
    
    
    a=apply(PM,2,min)
    PV[j,1:2]=c(a[1]<0.0025,a[2]<0.0025)
    
    k0=sum(PM[,2]<0.0025)
    
    PV[j,3]=0
    if (k0>0) {
      if (k0==1){
        a1=PM[PM[,2]<0.0025,]
        PV[j,3]=c(a[1]< 0.05/k0)
      } else {
        a1=PM[PM[,2]<0.0025,]
        a=apply(a1,2,min)
        PV[j,3]=c(a[1]<0.05/k0)
      }
    }
    
    if(j %% 100==0){print(c(i,j))}
    Theta1[j,i]=r
  }
  B1[i,]=colMeans(PV)
}
