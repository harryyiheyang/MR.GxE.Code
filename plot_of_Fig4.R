library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(glue)

HDL_EUR<-fread("/mnt/vstor/SOM_EPBI_XXZ10/njl96/MR/data/gene_exposure/lipids_new/ancestry_specific/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz",header=T)%>%mutate(pvalue=as.numeric(pvalue),pvalue_GC=as.numeric(pvalue_GC),lipid="HDL")
HDL_AFR<-fread("/mnt/vstor/SOM_EPBI_XXZ10/njl96/MR/data/gene_exposure/lipids_new/ancestry_specific/HDL_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz",header=T)%>%mutate(pvalue=as.numeric(pvalue),pvalue_GC=as.numeric(pvalue_GC),lipid="HDL")
HDL_HIS<-fread("/mnt/vstor/SOM_EPBI_XXZ10/njl96/MR/data/gene_exposure/lipids_new/ancestry_specific/HDL_INV_HIS_1KGP3_ALL.meta.singlevar.results.gz",header=T)%>%mutate(pvalue=as.numeric(pvalue),pvalue_GC=as.numeric(pvalue_GC),lipid="HDL")
HDL_EAS<-fread("/mnt/vstor/SOM_EPBI_XXZ10/njl96/MR/data/gene_exposure/lipids_new/ancestry_specific/HDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz",header=T)%>%mutate(pvalue=as.numeric(pvalue),pvalue_GC=as.numeric(pvalue_GC),lipid="HDL")

snplist2=c("rs4970836","rs10779836","rs61400250","rs34497903","rs1820163",
           "rs1569209","rs28597716",
           "rs112875651","rs7846466",
           "rs11591480",
           "rs3741298","rs61905078","rs5128",
           "rs1601934","rs12708454","rs2070895","rs573922",
           "rs8045855","rs247616",
           "rs56097211","rs28685665",
           "rs10402112",
           "rs739846",
           "rs6857","rs3208856","rs4803759","rs4420638","rs1601934"
)

snplist3=readxl::read_xlsx("/mnt/vstor/SOM_EPBI_XXZ10/gxl352/home/snplist3.xlsx",sheet="Sheet1")%>%as.data.frame()
snplistunion=union(snplist3$rs_dbSNP150,snplist2)
HDL_AFR=HDL_AFR[which(HDL_AFR$rsID%in%snplistunion),]
HDL_EUR=HDL_EUR[which(HDL_EUR$rsID%in%snplistunion),]
HDL_EAS=HDL_EAS[which(HDL_EAS$rsID%in%snplistunion),]
HDL_HIS=HDL_HIS[which(HDL_HIS$rsID%in%snplistunion),]
data=list(HDL_EUR,HDL_AFR,HDL_HIS,HDL_EAS)

#############################################################################
LDL_EUR<-fread("/mnt/vstor/SOM_EPBI_XXZ10/njl96/MR/data/gene_exposure/lipids_new/ancestry_specific/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz",header=T)%>%mutate(pvalue=as.numeric(pvalue),pvalue_GC=as.numeric(pvalue_GC),lipid="LDL")
LDL_AFR<-fread("/mnt/vstor/SOM_EPBI_XXZ10/njl96/MR/data/gene_exposure/lipids_new/ancestry_specific/LDL_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz",header=T)%>%mutate(pvalue=as.numeric(pvalue),pvalue_GC=as.numeric(pvalue_GC),lipid="LDL")
LDL_HIS<-fread("/mnt/vstor/SOM_EPBI_XXZ10/njl96/MR/data/gene_exposure/lipids_new/ancestry_specific/LDL_INV_HIS_1KGP3_ALL.meta.singlevar.results.gz",header=T)%>%mutate(pvalue=as.numeric(pvalue),pvalue_GC=as.numeric(pvalue_GC),lipid="LDL")
LDL_EAS<-fread("/mnt/vstor/SOM_EPBI_XXZ10/njl96/MR/data/gene_exposure/lipids_new/ancestry_specific/LDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz",header=T)%>%mutate(pvalue=as.numeric(pvalue),pvalue_GC=as.numeric(pvalue_GC),lipid="LDL")

snplist2=c("rs4970836","rs10779836","rs61400250","rs34497903","rs1820163",
           "rs1569209","rs28597716",
           "rs112875651","rs7846466",
           "rs11591480",
           "rs3741298","rs61905078","rs5128",
           "rs1601934","rs12708454","rs2070895","rs573922",
           "rs8045855","rs247616",
           "rs56097211","rs28685665",
           "rs10402112",
           "rs739846",
           "rs6857","rs3208856","rs4803759","rs4420638"
)
snplist3=readxl::read_xlsx("/mnt/vstor/SOM_EPBI_XXZ10/gxl352/home/snplist3.xlsx",sheet="Sheet1")%>%as.data.frame()
snplistunion=union(snplist3$rs_dbSNP150,snplist2)
LDL_AFR=LDL_AFR[which(LDL_AFR$rsID%in%snplistunion),]
LDL_EUR=LDL_EUR[which(LDL_EUR$rsID%in%snplistunion),]
LDL_EAS=LDL_EAS[which(LDL_EAS$rsID%in%snplistunion),]
LDL_HIS=LDL_HIS[which(LDL_HIS$rsID%in%snplistunion),]
data=list(LDL_EUR,LDL_AFR,LDL_HIS,LDL_EAS)

###########################################################################################
TG_EUR<-fread("/mnt/vstor/SOM_EPBI_XXZ10/njl96/MR/data/gene_exposure/lipids_new/ancestry_specific/logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz",header=T)%>%mutate(pvalue=as.numeric(pvalue),pvalue_GC=as.numeric(pvalue_GC),lipid="TG")
TG_AFR<-fread("/mnt/vstor/SOM_EPBI_XXZ10/njl96/MR/data/gene_exposure/lipids_new/ancestry_specific/logTG_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz",header=T)%>%mutate(pvalue=as.numeric(pvalue),pvalue_GC=as.numeric(pvalue_GC),lipid="TG")
TG_HIS<-fread("/mnt/vstor/SOM_EPBI_XXZ10/njl96/MR/data/gene_exposure/lipids_new/ancestry_specific/logTG_INV_HIS_1KGP3_ALL.meta.singlevar.results.gz",header=T)%>%mutate(pvalue=as.numeric(pvalue),pvalue_GC=as.numeric(pvalue_GC),lipid="TG")
TG_EAS<-fread("/mnt/vstor/SOM_EPBI_XXZ10/njl96/MR/data/gene_exposure/lipids_new/ancestry_specific/logTG_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz",header=T)%>%mutate(pvalue=as.numeric(pvalue),pvalue_GC=as.numeric(pvalue_GC),lipid="TG")

snplist2=c("rs4970836","rs10779836","rs61400250","rs34497903","rs1820163",
           "rs1569209","rs28597716",
           "rs112875651","rs7846466",
           "rs11591480",
           "rs3741298","rs61905078","rs5128",
           "rs1601934","rs12708454","rs2070895","rs573922",
           "rs8045855","rs247616",
           "rs56097211","rs28685665",
           "rs10402112",
           "rs739846",
           "rs6857","rs3208856","rs4803759","rs4420638"
)
snplistunion=union(snplist3$rs_dbSNP150,snplist2)
TG_AFR=TG_AFR[which(TG_AFR$rsID%in%snplistunion),]
TG_EUR=TG_EUR[which(TG_EUR$rsID%in%snplistunion),]
TG_EAS=TG_EAS[which(TG_EAS$rsID%in%snplistunion),]
TG_HIS=TG_HIS[which(TG_HIS$rsID%in%snplistunion),]

HDL_AFR1=HDL_AFR%>%dplyr::select(rsID,CHR=CHROM,BP=POS_b37,BETA=EFFECT_SIZE,SE,lipid,P=pvalue_GC,AF=POOLED_ALT_AF)
HDL_EUR1=HDL_EUR%>%dplyr::select(rsID,CHR=CHROM,BP=POS_b37,BETA=EFFECT_SIZE,SE,lipid,P=pvalue_GC,AF=POOLED_ALT_AF)
HDL_EAS1=HDL_EAS%>%dplyr::select(rsID,CHR=CHROM,BP=POS_b37,BETA=EFFECT_SIZE,SE,lipid,P=pvalue_GC,AF=POOLED_ALT_AF)
HDL_HIS1=HDL_HIS%>%dplyr::select(rsID,CHR=CHROM,BP=POS_b37,BETA=EFFECT_SIZE,SE,lipid,P=pvalue_GC,AF=POOLED_ALT_AF)
LDL_AFR1=LDL_AFR%>%dplyr::select(rsID,CHR=CHROM,BP=POS_b37,BETA=EFFECT_SIZE,SE,lipid,P=pvalue_GC,AF=POOLED_ALT_AF)
LDL_EUR1=LDL_EUR%>%dplyr::select(rsID,CHR=CHROM,BP=POS_b37,BETA=EFFECT_SIZE,SE,lipid,P=pvalue_GC,AF=POOLED_ALT_AF)
LDL_EAS1=LDL_EAS%>%dplyr::select(rsID,CHR=CHROM,BP=POS_b37,BETA=EFFECT_SIZE,SE,lipid,P=pvalue_GC,AF=POOLED_ALT_AF)
LDL_HIS1=LDL_HIS%>%dplyr::select(rsID,CHR=CHROM,BP=POS_b37,BETA=EFFECT_SIZE,SE,lipid,P=pvalue_GC,AF=POOLED_ALT_AF)
TG_AFR1=TG_AFR%>%dplyr::select(rsID,CHR=CHROM,BP=POS_b37,BETA=EFFECT_SIZE,SE,lipid,P=pvalue_GC,AF=POOLED_ALT_AF)
TG_EUR1=TG_EUR%>%dplyr::select(rsID,CHR=CHROM,BP=POS_b37,BETA=EFFECT_SIZE,SE,lipid,P=pvalue_GC,AF=POOLED_ALT_AF)
TG_EAS1=TG_EAS%>%dplyr::select(rsID,CHR=CHROM,BP=POS_b37,BETA=EFFECT_SIZE,SE,lipid,P=pvalue_GC,AF=POOLED_ALT_AF)
TG_HIS1=TG_HIS%>%dplyr::select(rsID,CHR=CHROM,BP=POS_b37,BETA=EFFECT_SIZE,SE,lipid,P=pvalue_GC,AF=POOLED_ALT_AF)

EAS=rbind(LDL_EAS1,HDL_EAS1,TG_EAS1)
EUR=rbind(LDL_EUR1,HDL_EUR1,TG_EUR1)
AFR=rbind(LDL_AFR1,HDL_AFR1,TG_AFR1)
HIS=rbind(LDL_HIS1,HDL_HIS1,TG_HIS1)

data=list(EUR,AFR,HIS,EAS)
k=0
coll=c('#ffc952', '#47b8e0', '#ff7473')
for(i in 1:3){
  for(j in (i+1):4){
    
    data1=data[[i]]%>%dplyr::filter(rsID%in%snplist3$rs_dbSNP150&P<5e-8)%>%dplyr::select(rsID,X=BETA,Lipid=lipid)
    data2=data[[j]]%>%dplyr::filter(rsID%in%snplist3$rs_dbSNP150&P<5e-8)%>%dplyr::select(rsID,X=BETA,Lipid=lipid)
    plotdata=merge(data1,data2,by="rsID")
    ss=which(plotdata$Lipid.x==plotdata$Lipid.y)
    plotdata=plotdata[ss,]
    rho1=cor(plotdata$X.x[which(plotdata$Lipid.x=="HDL")],plotdata$X.y[which(plotdata$Lipid.x=="HDL")])
    rho2=cor(plotdata$X.x[which(plotdata$Lipid.x=="LDL")],plotdata$X.y[which(plotdata$Lipid.x=="LDL")])
    rho3=cor(plotdata$X.x[which(plotdata$Lipid.x=="TG")],plotdata$X.y[which(plotdata$Lipid.x=="TG")])
    rho=c(rho1,rho2,rho3)
    plotdata$Lipid.x[which(plotdata$Lipid.x=="HDL")]=paste0("HDL-C: rho=",round(rho1,3))
    plotdata$Lipid.x[which(plotdata$Lipid.x=="LDL")]=paste0("LDL-C: rho=",round(rho2,3))
    plotdata$Lipid.x[which(plotdata$Lipid.x=="TG")]=paste0("TG: rho=",round(rho3,3))
    names(plotdata)[3]="Lipid"
    plot1=ggplot(data=plotdata,aes(x=X.x,y=X.y,color=Lipid,fill=Lipid))+theme_bw()+theme(panel.grid=element_blank())+
      geom_point(shape=21,size=5)+
      geom_smooth(aes(color=Lipid),method="lm",formula = 'y ~ x',size=2)+
      geom_vline(xintercept = 0,size=1.5,linetype="dashed")+
      geom_hline(yintercept = 0,size=1.5,linetype="dashed")+
      xlab((Ancestry[i]))+
      ylab(Ancestry[j])+
      scale_fill_manual(values=coll)+scale_color_manual(values=coll)+
      theme(text=element_text(size=22),panel.border = element_rect(size = 2),legend.text = element_text(face = "bold"))+
      ggtitle(glue("({LETTERS[k+1]})"))+theme(legend.position =c(0.01,.99),legend.justification = c(0.01,.99),legend.key.size = unit(1, 'cm'),legend.title=element_blank())
    
    k=k+1
    plot_list[[k]]=plot1
    
    
    data1=data[[i]]
    data1=data1%>%dplyr::filter(rsID%in%snplist2|(rsID=="rs1601934"&AF>0.05))%>%dplyr::select(rsID,X=BETA,Lipid=lipid)
    data2=data[[j]]
    data2=data2%>%dplyr::filter(rsID%in%snplist2|(rsID=="rs1601934"&AF>0.05))%>%dplyr::select(rsID,X=BETA,Lipid=lipid)
    plotdata=merge(as.data.frame(data1),as.data.frame(data2),by="rsID")
    ss=which(plotdata$Lipid.x==plotdata$Lipid.y)
    plotdata=plotdata[ss,]
    rho1=cor(plotdata$X.x[which(plotdata$Lipid.x=="HDL")],plotdata$X.y[which(plotdata$Lipid.x=="HDL")])
    rho2=cor(plotdata$X.x[which(plotdata$Lipid.x=="LDL")],plotdata$X.y[which(plotdata$Lipid.x=="LDL")])
    rho3=cor(plotdata$X.x[which(plotdata$Lipid.x=="TG")],plotdata$X.y[which(plotdata$Lipid.x=="TG")])
    rho=c(rho1,rho2,rho3)
    plotdata$Lipid.x[which(plotdata$Lipid.x=="HDL")]=paste0("HDL-C: rho=",round(rho1,3))
    plotdata$Lipid.x[which(plotdata$Lipid.x=="LDL")]=paste0("LDL-C: rho=",round(rho2,3))
    plotdata$Lipid.x[which(plotdata$Lipid.x=="TG")]=paste0("TG: rho=",round(rho3,3))
    names(plotdata)[3]="Lipid"
    plot2=ggplot(data=plotdata,aes(x=X.x,y=X.y,color=Lipid,fill=Lipid))+theme_bw()+theme(panel.grid=element_blank())+
      geom_point(shape=21,size=5)+
      geom_smooth(aes(color=Lipid),method="lm",formula = 'y ~ x',size=2)+
      geom_vline(xintercept = 0,size=1.5,linetype="dashed")+
      geom_hline(yintercept = 0,size=1.5,linetype="dashed")+
      xlab((Ancestry[i]))+
      ylab(Ancestry[j])+
      scale_fill_manual(values=coll)+scale_color_manual(values=coll)+
      theme(text=element_text(size=22),panel.border = element_rect(size = 2),legend.text = element_text(face = "bold"))+
      ggtitle(glue("({LETTERS[k+1]})"))+theme(legend.position =c(0.01,.99),legend.justification = c(0.01,.99),legend.key.size = unit(1, 'cm'),legend.title=element_blank())
    
    k=k+1
    plot_list[[k]]=plot2
  }
}
layout <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow = 3)
png("figure_0818.png",width=36*0.9,height=24*0.9,unit="in",res=300)
egg::ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],nrow=3)
dev.off()
pdf("figure_0818.pdf",width=36*0.9,height=24*0.9)
egg::ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],nrow=3)
dev.off()