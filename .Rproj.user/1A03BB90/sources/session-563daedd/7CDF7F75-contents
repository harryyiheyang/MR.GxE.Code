library(CMplot)
library(tidyverse)
library(dplyr)
setwd("~/plott")
bottom1 <- readRDS("~/plott/bottom1.rds")
top1 <- readRDS("~/plott/top1.rds")
joint1=merge(top1,bottom1,by="SNP")
joint1=dplyr::inner_join(
  top1%>%tidyr::drop_na(SNP)%>%distinct(SNP,.keep_all=TRUE)%>%dplyr::select(-CHR,-BP),
  bottom1%>%tidyr::drop_na(SNP)%>%distinct(SNP,.keep_all=TRUE),by="SNP")

names(joint1)=c("SNP","TG-CD","Chromosome","Position","TG-ED")
joint1= joint1 %>% select("SNP","Chromosome","Position","TG-CD","TG-ED")

CMplot(joint1,type="p",plot.type="c",chr.labels=paste("Chr",c(1:22),sep=""),r=0.4,cir.legend=T,
       outward=F,cir.legend.col="black",cir.chr.h=1.3,chr.den.col="black",signal.line=0.5,threshold.col="black",
       amplify=T,signal.col="#f1404b",signal.cex=0.5,col=matrix(c("#80d4f6","#03a6ff","#e3dede","#7C7877"),2,2),threshold=5e-8,
       file="jpg",dpi=300)

scatter1 <- readRDS("~/plott/scatter1.rds")
scatter2 <- readRDS("~/plott/scatter2.rds")
#add1=scatter1 %>% filter((CHR=="2"|CHR=="19")&IsPleiotropic==1)
#add1=add1[c(3,4),]
#add2=scatter2 %>% filter((CHR=="2"|CHR=="19")&IsPleiotropic==1)
#add2=add2[7,]
MR1 <- readRDS("~/plott/MR1.rds")
MR2 <- readRDS("~/plott/MR2.rds")
plink1= scatter1 %>% select(SNP,CHR,BP,PleioP_MR)
names(plink1)[4]="P"
plink2= scatter2 %>% select(SNP,CHR,BP,PleioP_MR)
names(plink2)[4]="P"
write.table(plink1,file="~/plott/plinkfile/TGDrink1.txt", quote=F, sep="\t", row.name=F)
write.table(plink2,file="~/plott/plinkfile/TGDrink2.txt", quote=F, sep="\t", row.name=F)
# plink --bfile /mnt/rstor/SOM_EPBI_XXZ10/njl96/data/1000G/1kg_phase3_EUR_only --clump ~/plott/plinkfile/TGDrink1.txt   --clump-field P  --clump-kb 500 --clump-p1 8e-1 --clump-p2 8e-1 --clump-r2 0.01 --out ~/plott/plinkfile/TGDrink1
# plink --bfile /mnt/rstor/SOM_EPBI_XXZ10/njl96/data/1000G/1kg_phase3_EUR_only --clump ~/plott/plinkfile/TGDrink2.txt   --clump-field P  --clump-kb 500 --clump-p1 8e-1 --clump-p2 8e-1 --clump-r2 0.01 --out ~/plott/plinkfile/TGDrink2

ind1=fread("plinkfile/TGDrink1.clumped")
ind2=fread("plinkfile/TGDrink2.clumped")
indclump1=which(scatter1$SNP %in% ind1$SNP)
clump1=scatter1[indclump1,]
#clump1=rbind(clump1,add1)
clump1$CHR=as.numeric(clump1$CHR)
clump1$BP=as.numeric(clump1$BP)
clump1=dplyr::arrange(clump1,CHR,BP)
(TG_Current_Smoke=clump1[which(clump1$IsPleiotropic==1),])
clump2=clump2[-c(which(clump2$SNP=="rs483082"),which(clump2$SNP=="rs438811"))]
gene1=c("LPL","ZPR1","APOE")
TG_Current_Smoke$gene=gene1
clump1$Interaction[clump1$IsPleiotropic==1]=gene1
clump1$Interaction[clump1$IsPleiotropic==0]="No Interaction"

indclump2=which(scatter2$SNP %in% ind2$SNP)
clump2=scatter2[indclump2,]
#clump2=rbind(clump2,add2)
clump2$CHR=as.numeric(clump2$CHR)
clump2$BP=as.numeric(clump2$BP)
clump2=dplyr::arrange(clump2,CHR,BP)
#clump2=clump2[-which(clump2$SNP=="rs4420638"),]
(TG_Ever_Smoke=clump2[which(clump2$IsPleiotropic==1),])
gene2=c("LPL","LINC00861","ZPR1","APOE")
TG_Ever_Smoke$gene=gene2
clump2$Interaction[clump2$IsPleiotropic==1]=gene2
clump2$Interaction[clump2$IsPleiotropic==0]="No Interaction"

clump2$Interaction=ordered(clump2$Interaction,levels=c(gene2,"No Interaction"))
clump1$Interaction=ordered(clump1$Interaction,levels=c(gene1,"No Interaction"))
color2=c("#fe4365","#FADAD8","#a3daff","#113285")
color1=color2[-4]
color1[4]="#BDC0BA";color2[5]="#BDC0BA";

q0=qnorm(1-0.05/2)

plot1=ggplot(clump1, aes(x=x1,y=x2, fill=Interaction,color=Interaction)) +
  geom_point(size=8,shape=21) +theme_bw()+
  scale_fill_manual(values=color1)+theme(panel.grid=element_blank())+
  scale_color_manual(values=c(rep("black",3),rep("#BDC0BA",615-3)))+
  guides(color=NULL)+theme(legend.position=c(0.125,0.75),legend.title=element_text(size=0),legend.background=element_rect(linetype="solid",color="black"))+
  geom_vline(xintercept=0,linetype="dashed") +
  geom_hline(yintercept=0,linetype="dashed") +
  scale_x_continuous(limits=c(-0.032,0.061),breaks=seq(-0.03,0.06,by=0.03))+
  scale_y_continuous(limits=c(-0.065,0.12),breaks=seq(-0.06,0.12,by=0.06))+
  geom_abline(intercept=0, slope=MR1$CausalEstimate,color="black",size=1.5) +
  theme(legend.text=element_text(size=19),axis.title.x=element_text(size=22),axis.title.y=element_text(size=22))+
  theme(plot.title=element_text(size=25,hjust=0.5))+theme(axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))+
  labs(x=bquote("Main Effect"~beta[1]),y=bquote("Marginal Effect"~alpha),title="TG x Current Drinking")

plot2=ggplot(clump2, aes(x=x1,y=x2, fill=Interaction,color=Interaction)) +
  geom_point(size=8,shape=21) +theme_bw()+
  scale_fill_manual(values=color2)+theme(panel.grid=element_blank())+
  scale_color_manual(values=c(rep("black",4),rep("#BDC0BA",616-4)))+
  guides(color=NULL)+theme(legend.position=c(0.125,0.75),legend.title=element_text(size=0),legend.background=element_rect(linetype="solid",color="black"))+
  geom_vline(xintercept=0,linetype="dashed") +
  geom_hline(yintercept=0,linetype="dashed") +
  scale_x_continuous(limits=c(-0.032,0.061),breaks=seq(-0.03,0.06,by=0.03))+
  scale_y_continuous(limits=c(-0.065,0.12),breaks=seq(-0.06,0.12,by=0.06))+
  geom_abline(intercept=0, slope=MR2$CausalEstimate,color="black",size=1.5) +
  theme(legend.text=element_text(size=19),axis.title.x=element_text(size=22),axis.title.y=element_text(size=22))+
  theme(plot.title=element_text(size=25,hjust=0.5))+theme(axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))+
  labs(x=bquote("Main Effect"~beta[1]), y=bquote("Marginal Effect"~alpha),title="TG x Regular Drinking")
