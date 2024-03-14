library(egg)
library(ggplot2)
library(RColorBrewer)
getPalette=colorRampPalette(brewer.pal(9, "RdBu"))
se2=2*sqrt(0.95*0.05/3000)
setwd("~/BD")

load("Binom-MEAN-heter-positive-005.RData")
A=as.data.frame(as.vector(Theta))
names(A)="theta"
B=Theta
for(i in 1:10) B[,i]=glue::glue("{ME[i]}")
A$sample=as.vector(B)
A$sample=ordered(A$sample,levels=B[1,])
plot1=ggplot(data=A,aes(x=sample,y=theta,fill=sample))+geom_boxplot(fill="#45d9fd",alpha=0.5)+theme_bw()+theme(panel.grid = element_blank())+
  scale_y_continuous(limits=c(0.9,1.03),breaks=seq(0.9,1,by=0.025))+
  theme(legend.position = "none")+ylab(bquote(hat(theta)))+xlab(bquote("GWAS environmental factor mean"~mu[E]^mar))+ggtitle(bquote("(A). Estimation of"~theta))+
  geom_hline(yintercept=1,linetype="dashed")+theme(text = element_text(size = 22),panel.border = element_rect(size = 1.5))
################################################################################################################################################################

load("Binom-MEAN-heter-positive-005.RData")
BETA=BETA[1:10,,]

A=NULL
A=as.data.frame(as.vector(BETA[,,c(3:4)]))
names(A)="est"
B=BETA[,,c(3:4)]
for(i in 1:10) B[i,,]=glue::glue("{ME[i]}")
A$sample=as.vector(B)
A$sample=ordered(A$sample,levels=B[,1,1])
B[,,1]="Direct estimate";B[,,2]="MR-GxE estimate"
A$method=as.vector(B)
A$s="s = +1"
A1=A

load("Binom-MEAN-heter-1608080-negative-005.RData")
POW1=POW1[1:10,];POW2=POW2[1:10,]
BETA=BETA[1:10,,]
Theta=Theta[,1:10]

A=NULL
A=as.data.frame(as.vector(BETA[,,c(3:4)]))
names(A)="est"
B=BETA[,,c(3:4)]
for(i in 1:10) B[i,,]=glue::glue("{ME[i]}")
A$sample=as.vector(B)
A$sample=ordered(A$sample,levels=B[,1,1])
B[,,1]="Direct estimate";B[,,2]="MR-GxE estimate"
A$method=as.vector(B)
A$s="s = -1"
A1=rbind(A,A1)

load("Binom-MEAN-heter-1608080-zero-005.RData")
POW1=POW1[1:10,];POW2=POW2[1:10,]
BETA=BETA[1:10,,]
Theta=Theta[,1:10]

A=NULL
A=as.data.frame(as.vector(BETA[,,c(3:4)]))
names(A)="est"
B=BETA[,,c(3:4)]
for(i in 1:10) B[i,,]=glue::glue("{ME[i]}")
A$sample=as.vector(B)
A$sample=ordered(A$sample,levels=B[,1,1])
B[,,1]="Direct estimate";B[,,2]="MR-GxE estimate"
A$method=as.vector(B)
A$s="s = 0"
A1=rbind(A,A1)

plot2=ggplot(data=A1,aes(x=sample,y=est,fill=sample))+geom_boxplot(fill="#45d9fd",alpha=0.5)+theme_bw()+theme(panel.grid = element_blank(),panel.border = element_rect(size = 1.5))+
  facet_grid(method~s)+
  scale_y_continuous(limits=c(-0.1,0.1),breaks=seq(-0.1,0.1,by=0.05))+  scale_x_discrete(label=c("","0.1","","0.2","","0.3","","0.4","","0.5"))+
  theme(legend.position = "none")+ylab(bquote(hat(beta)[3]))+xlab(bquote("GWAS environmental factor mean"~mu[E]^mar))+ggtitle(bquote("B. Estimation of interaction effect"~beta[3]))+geom_hline(yintercept=0,linetype="dashed")+
  theme(text = element_text(size = 22))
################################################################################################################################################################################################################################

load("Binom-MEAN-heter-positive-005.RData")

A=NULL
A=as.data.frame(as.vector(POW1))
names(A)="power"
B=POW1
for(i in 1:10) B[i,]=B[i,]=ME[i]
A$sample=as.vector(B)
B[,1]="Direct test";B[,2]="MR-GxE test"
A$method=as.vector(B)
A$s="s = +1"
A$typeI=as.vector(POW2)
A1=A

load("Binom-MEAN-heter-negative-005.RData")
POW1=POW1[1:10,];POW2=POW2[1:10,]
BETA=BETA[1:10,,]
Theta=Theta[,1:10]

A=NULL
A=as.data.frame(as.vector(POW1))
names(A)="power"
B=POW1
for(i in 1:10) B[i,]=B[i,]=ME[i]
A$sample=as.vector(B)
B[,1]="Direct test";B[,2]="MR-GxE test"
A$method=as.vector(B)
A$s="s = -1"
A$typeI=as.vector(POW2)
A1=rbind(A,A1)

load("Binom-MEAN-heter-zero-005.RData")

A=NULL
A=as.data.frame(as.vector(POW1))
names(A)="power"
B=POW1
for(i in 1:10) B[i,]=ME[i]
A$sample=as.vector(B)
B[,1]="Direct test";B[,2]="MR-GxE test"
A$method=as.vector(B)
A$s="s = 0"
A$typeI=as.vector(POW2)
A1=rbind(A,A1)

se2=2*sqrt(0.05*0.95/3000)
plot4=ggplot(data=A1,aes(x=sample,y=power,color=method))+geom_line(size=3)+theme_bw()+theme(panel.grid = element_blank())+
  facet_grid(~s)+scale_x_continuous(limits=c(0,0.55),breaks=seq(0,0.5,by=0.1))+ theme(legend.title=element_blank())+
  theme(legend.position=c(0.02,0.98), legend.justification=c(0.02,0.98))+ylab("power")+xlab(bquote("GWAS environmental factor mean"~mu[E]^mar))+ggtitle("(D).  Power of the direct and MR-GxE tests, no mediation")+
  scale_color_manual(values=c('#45d9fd', '#fe4365'))+  theme(text = element_text(size = 22),panel.border = element_rect(size = 1.5))

plot3=ggplot(data=A1,aes(x=sample,y=typeI,color=method))+geom_line(size=3)+theme_bw()+theme(panel.grid = element_blank())+
  facet_grid(~s)+scale_x_continuous(limits=c(0,0.55),breaks=seq(0,0.5,by=0.1))+ theme(legend.title=element_blank())+
  theme(legend.position=c(0.02,0.98), legend.justification=c(0.02,0.98))+ylab("type-I error")+xlab(bquote("GWAS environmental factor mean"~mu[E]^mar))+ggtitle("C.  Type I error of the direct and MR_GxE tests, no mediation")+
  scale_y_continuous(limits=c(0.02,0.08),breaks=seq(0.02,0.08,by=0.01))+geom_hline(yintercept=0.05)+geom_hline(yintercept=0.05+se2,linetype="dashed")+geom_hline(yintercept=0.05-se2,linetype="dashed")+
  scale_color_manual(values=c('#45d9fd', '#fe4365'))+theme(text = element_text(size = 22),panel.border = element_rect(size = 1.5))
################################################################################################################################################################################################################################

load("~/BD/GESimu1.RData.RData")
se3=2*sqrt(0.05*0.95/1000)

SAMPLE=C1;SAMPLE[1,]="20K";SAMPLE[2,]="100K";SAMPLE[3,]="200K";SAMPLE[4,]="300K"
METHOD=C1;METHOD[,1]="Direct test";METHOD[,2]="MR-GxE test";METHOD[,3]="Two-step test"
K=as.data.frame(as.vector(C1));names(K)="power"
K$method=as.vector(METHOD);K$sample=as.vector(SAMPLE)
K$Condition="No mediation, E contributes 1% phenotypic variance"
K1=K

K=as.data.frame(as.vector(C2));names(K)="power"
K$method=as.vector(METHOD);K$sample=as.vector(SAMPLE)
K$Condition="Mediation, E contributes 1% phenotypic variance"
K1=rbind(K,K1)

K=as.data.frame(as.vector(C3));names(K)="power"
K$method=as.vector(METHOD);K$sample=as.vector(SAMPLE)
K$Condition="Mediation, E contributes 5% phenotypic variance"
K1=rbind(K,K1)

K1$Condition=ordered(K1$Condition,levels=c("No mediation, E contributes 1% phenotypic variance","Mediation, E contributes 1% phenotypic variance","Mediation, E contributes 5% phenotypic variance"))

plot5=ggplot(K1, aes(fill=method, y=power, x=sample))+theme_bw()+theme(panel.grid = element_blank())+
  geom_bar(position="dodge", stat="identity",color="black")+facet_grid(~Condition)+theme(legend.title=element_blank(),panel.border = element_rect(size = 1.5))+
  xlab(bquote(n[1]))+ylab("type-I error")+theme(legend.position=c(0.02,0.98), legend.justification=c(0.02,0.98))+scale_y_continuous(limits=c(0,0.5),breaks=seq(0,0.5,by=0.1))+
  geom_hline(yintercept=0.05)+geom_hline(yintercept=0.05+se3,linetype="dashed")+geom_hline(yintercept=0.05-se3,linetype="dashed")+
  ggtitle("E. Type I error of the direct, MR_GxE and two-step tests, 20 variants, no mediation or mediation")+theme(text = element_text(size = 22))+
  scale_fill_manual(values=c('#45d9fd', '#fe4365','#fff1b9'))
###################################################################################################################################################################################

K=as.data.frame(as.vector(C4));names(K)="power"
K$method=as.vector(METHOD);K$sample=as.vector(SAMPLE)
K$Condition="No mediation, E contributes 1% phenotypic variance"
K1=K

K=as.data.frame(as.vector(C5));names(K)="power"
K$method=as.vector(METHOD);K$sample=as.vector(SAMPLE)
K$Condition="Mediation, E contributes 1% phenotypic variance"
K1=rbind(K,K1)

K=as.data.frame(as.vector(C6));names(K)="power"
K$method=as.vector(METHOD);K$sample=as.vector(SAMPLE)
K$Condition="Mediation, E contributes 5% phenotypic variance"
K1=rbind(K,K1)

K1$Condition=ordered(K1$Condition,levels=c("No mediation, E contributes 1% phenotypic variance","Mediation, E contributes 1% phenotypic variance","Mediation, E contributes 5% phenotypic variance"))

plot6=ggplot(K1, aes(fill=method, y=power, x=sample))+theme_bw()+theme(panel.grid = element_blank())+
  geom_bar(position="dodge", stat="identity",color="black")+facet_grid(~Condition)+theme(legend.title=element_blank(),panel.border = element_rect(size = 1.5))+
  xlab(bquote(n[1]))+ylab("power")+theme(legend.position=c(0.02,0.98), legend.justification=c(0.02,0.98))+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.2))+
  ggtitle("F. Power of the direct, MR_GxE and two-step tests, 20 variants, no mediation or mediation")+theme(text = element_text(size = 22))+
  scale_fill_manual(values=c('#45d9fd', '#fe4365','#fff1b9'))
#####################################################################################################################################################################################################

pdf("mainplotbinomE.pdf",width=25,height=30)
pushViewport(viewport(layout = grid.layout(5, 3)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(plot1, vp = vplayout(1, 1))
print(plot2, vp = vplayout(1, c(2,3)))
print(plot3, vp = vplayout(2, c(1,3)))
print(plot4, vp = vplayout(3, c(1,3)))
print(plot5, vp = vplayout(4, c(1,3)))
print(plot6, vp = vplayout(5, c(1,3)))
dev.off()
