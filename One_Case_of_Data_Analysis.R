library(data.table)
library(IMRP)
library(dplyr)
SearchPleio=function (BetaOutcome, BetaExposure, SdOutcome, SdExposure, data, 
                      rho, CausalEstimate, SdCausalEstimate) 
{
  X <- data[, c(BetaOutcome)] - CausalEstimate * data[, c(BetaExposure)]
  Y <- sqrt(data[, c(SdOutcome)]^2 + CausalEstimate^2 * data[, 
                                                             c(SdExposure)]^2 + data[, c(BetaExposure)]^2 * SdCausalEstimate^2 - 
              2 * CausalEstimate * rho * data[, c(SdExposure)] * data[, 
                                                                      c(SdOutcome)])
  pval <- apply(cbind(2 * pnorm(-abs(X/Y)), 1), 1, min)
  return(list(pleio_p = pval,beta=X,SE=Y))
}

# LDL_Drinking
HDL_EUR_1<-fread("/mnt/rstor/SOM_EPBI_XXZ10/xxz10/gene-life/AlcoholLipids/META.LDL.CURDRINK.M1_2df1.txt.gz")
HDL_EUR_2<-fread("/mnt/rstor/SOM_EPBI_XXZ10/njl96/MR/data/gene_exposure/lipids_new/trans_ancestry/meta-analysis_AFR_EAS_EUR_HIS_SAS_LDL_INV_ALL_with_N_1.gz")
# Interaction data QCs
HDL_EUR_1$Allele1=toupper(HDL_EUR_1$Allele1)
HDL_EUR_1$Allele2=toupper(HDL_EUR_1$Allele2)
A=strsplit( HDL_EUR_1$MarkerName , ":" ) 
B=unlist(lapply(A,length))
HDL_EUR_1$VarType=B    # VarType: 1=rs number, 2= Chr+Pos, 3= I/D 
HDL_EUR_1_rs=HDL_EUR_1[HDL_EUR_1$VarType==1,]
HDL_EUR_1_ChrPos=HDL_EUR_1[HDL_EUR_1$VarType==2,]
HDL_EUR_1_ID=HDL_EUR_1[HDL_EUR_1$VarType==3,]
A=unlist(strsplit(HDL_EUR_1_ID$MarkerName, ":" )) 
B=matrix(A,ncol=3, byrow=TRUE)
B=data.frame(B)
names(B) <- c("CHROM","POS_b37","Type")
B=B[,c("CHROM","POS_b37")]
HDL_EUR_1_ID=cbind(B,HDL_EUR_1_ID)
A=unlist(strsplit(HDL_EUR_1_ChrPos$MarkerName, ":" )) 
B=matrix(A,ncol=2, byrow=TRUE)
B=data.frame(B)
names(B) <- c("CHROM","POS_b37")
HDL_EUR_1_ChrPos=cbind(B,HDL_EUR_1_ChrPos)
HDL_EUR_1_ID=HDL_EUR_1_ID[,c("CHROM", "POS_b37", "MarkerName", "Allele1", "Allele2", "Freq1", "Effect","StdErr","IntEffect", "IntStdErr","IntCov","ChiSq2df", "P-value","N","VarType")]
HDL_EUR_1_ChrPos=HDL_EUR_1_ChrPos[,c("CHROM", "POS_b37", "MarkerName", "Allele1", "Allele2", "Freq1", "Effect","StdErr","IntEffect", "IntStdErr","IntCov","ChiSq2df", "P-value","N","VarType")]
HDL_EUR_1_rs=HDL_EUR_1_rs[,c("MarkerName", "Allele1", "Allele2", "Freq1", "Effect","StdErr","IntEffect", "IntStdErr","IntCov","ChiSq2df", "P-value","N", "VarType")]
HDL_EUR_ID=merge(HDL_EUR_1_ID,HDL_EUR_2,by=c("CHROM","POS_b37"))
HDL_EUR_ChrPos=merge(HDL_EUR_1_ChrPos,HDL_EUR_2,by=c("CHROM","POS_b37"))
HDL_EUR_rs=merge(HDL_EUR_1_rs,HDL_EUR_2,by.x="MarkerName", by.y="rsID")
HDL_EUR_rs$rsID=HDL_EUR_rs$MarkerName
HDL_EUR=rbind(HDL_EUR_ChrPos,HDL_EUR_ID,HDL_EUR_rs) # curently there are 8344809 rows.
names(HDL_EUR)[14] <- c("TotalSampleSize")
names(HDL_EUR)[19] <- c("N")
names(HDL_EUR)[13]<-c("Pvalue")
HDL_EUR$FreqMin=ifelse(HDL_EUR$Freq1<0.5,HDL_EUR$Freq1, 1-HDL_EUR$Freq1)              
HDL_EUR$POOLED_Min_AF=ifelse(HDL_EUR$POOLED_ALT_AF<0.5,HDL_EUR$POOLED_ALT_AF, 1-HDL_EUR$POOLED_ALT_AF)
HDL_EUR=HDL_EUR[abs(HDL_EUR$FreqMin-HDL_EUR$POOLED_Min_AF)<0.15,]  # 8272473 rows
HDL_EUR=HDL_EUR[HDL_EUR$VarType==3 | (HDL_EUR$Allele1==HDL_EUR$REF & HDL_EUR$Allele2==HDL_EUR$ALT) | (HDL_EUR$Allele1==HDL_EUR$ALT & HDL_EUR$Allele2==HDL_EUR$REF),] #8233711 rows
A=HDL_EUR[HDL_EUR$VarType==3 | (HDL_EUR$Allele1==HDL_EUR$REF & HDL_EUR$Allele2==HDL_EUR$ALT) | (HDL_EUR$Allele1==HDL_EUR$ALT & HDL_EUR$Allele2==HDL_EUR$REF),] 
HDL_EUR=HDL_EUR[HDL_EUR$VarType<3 | (HDL_EUR$VarType==3 & !((HDL_EUR$REF %in% c("A","C","T","G")) & (HDL_EUR$ALT %in% c("A","C","T","G")))),]    # there are 8215940 rows
A=HDL_EUR[duplicated(HDL_EUR$MarkerName),]  # there 6169 duplicate rows
HDL_EUR=HDL_EUR[!(HDL_EUR$MarkerName %in% A$MarkerName) | ((HDL_EUR$MarkerName %in% A$MarkerName) & abs(HDL_EUR$FreqMin-HDL_EUR$POOLED_Min_AF)<0.01),] # there are 8209791 rows
HDL_EUR=HDL_EUR[!duplicated(HDL_EUR$MarkerName),]  # drop further 1651 duplicate rows. Now there are 8208140 rows finally. These variants may still be problematic and need to be further checked if there are any findings.
A=HDL_EUR[HDL_EUR$VarType < 3,]
B=HDL_EUR[HDL_EUR$VarType == 3,]
A$Effect_new=ifelse(A$Allele1 == A$ALT, A$Effect, -A$Effect)
B$Effect_new=ifelse(nchar(B$ALT) < nchar(B$REF), B$Effect, -B$Effect)
HDL_EUR=rbind(A,B)
"testA"
HDL_EUR$x1 <- HDL_EUR$Effect_new/HDL_EUR$StdErr/sqrt(HDL_EUR$TotalSampleSize)
HDL_EUR$x2 <- HDL_EUR$METAL_Effect/HDL_EUR$METAL_StdErr/sqrt(HDL_EUR$N)
HDL_EUR$x1_se <- 1/sqrt(HDL_EUR$TotalSampleSize)
HDL_EUR$x2_se <- 1/sqrt(HDL_EUR$N)
names(HDL_EUR)[1]<-c("CHR")
names(HDL_EUR)[2]<-c("BP")
names(HDL_EUR)[16]<-c("SNP")
HDL_EUR=HDL_EUR[as.numeric(HDL_EUR$N) >100000,]
HDL_EUR=HDL_EUR[as.numeric(HDL_EUR$TotalSampleSize) >50000,]
Alcloci=fread("/mnt/rstor/SOM_EPBI_XXZ10/xxz10/gene-life/SmkAlcGWAS/AlcSNPs.txt")
HDL_EUR=HDL_EUR[!HDL_EUR$SNP %in% Alcloci$RSID,]
A=HDL_EUR[as.numeric(HDL_EUR$pvalue_GC) < 5e-8,]
HDL_EUR_MR<-fread("/mnt/rstor/SOM_EPBI_XXZ10/xxz10/gene-life/AlcoholLipids/dataanalysis/LDL_TrEtn_CurDring_1000GPlink.clumped")
HDL_EUR_MR<-merge(A, HDL_EUR_MR, by="SNP")
A=HDL_EUR[as.numeric(HDL_EUR$pvalue_GC) > 0.05 & as.numeric(HDL_EUR$Pvalue) > 0.05,]
rho=cor(A$Effect_new/A$StdErr,A$METAL_Effect/A$METAL_StdErr)
MR1<-MR_pleio("x2","x1","x2_se","x1_se",as.data.frame(HDL_EUR_MR),SignifThreshold=0.05,rho=rho)  
A <-SearchPleio("x2","x1","x2_se","x1_se",as.data.frame(HDL_EUR),rho,MR1$CausalEstimate,MR1$SdCausalEstimate)
HDL_EUR$PleioP_MR=A$pleio_p
HDL_EUR$Pleiobeta=A$beta
HDL_EUR$PleioSE=A$SE
A=HDL_EUR[HDL_EUR$PleioP_MR<5e-8 | HDL_EUR$PleioP_MR_rho1<5e-8,] 
HDL_EUR$CHR=as.numeric(as.character(HDL_EUR$CHR))
HDL_EUR$BP=as.numeric(as.character(HDL_EUR$BP))
top=HDL_EUR%>%dplyr::select(SNP,CHR,BP,PleioP_MR,Pleiobeta,PleioSE,Allele1,Allele2,Freq1,IntEffect,IntStdErr,TotalSampleSize,METAL_Effect,METAL_StdErr)
#Sig_1=HDL_EUR%>%filter(PleioP_MR<5e-8)
"testB"
HDL_EUR$IsPleiotropic=ifelse(HDL_EUR$PleioP_MR<5e-8,"1","0")
HDL_EUR$IsPleiotropic=ordered(HDL_EUR$IsPleiotropic,levels=c("1","0"))
upper_highlight=HDL_EUR%>%filter(PleioP_MR<5e-8)%>%filter(SNP!="NA")
HDL_EUR$pvalue_GC=as.numeric(HDL_EUR$pvalue_GC)
HDL_EUR$Pvalue=as.numeric(HDL_EUR$Pvalue)
HDL_EUR_scatter=HDL_EUR%>%filter(PleioP_MR<5e-8|pvalue_GC < 5e-8|Pvalue<5e-8)
HDL_EUR_scatter=HDL_EUR_scatter %>% dplyr::select(x1,x2,x1_se,x2_se,IsPleiotropic,PleioP_MR,SNP,CHR,BP)
saveRDS(HDL_EUR_scatter,"~/plott/newrds/scatter1.rds")
saveRDS(MR1,"~/plott/newrds/MR1.rds")
"testC"
HDL_EUR_1<-fread("/mnt/rstor/SOM_EPBI_XXZ10/xxz10/gene-life/AlcoholLipids/META.LDL.REGDRINK.M1_2df1.txt.gz")
# Interaction data QCs
HDL_EUR_1$Allele1=toupper(HDL_EUR_1$Allele1)
HDL_EUR_1$Allele2=toupper(HDL_EUR_1$Allele2)
A=strsplit( HDL_EUR_1$MarkerName , ":" ) 
B=unlist(lapply(A,length))
HDL_EUR_1$VarType=B    # VarType: 1=rs number, 2= Chr+Pos, 3= I/D 
HDL_EUR_1_rs=HDL_EUR_1[HDL_EUR_1$VarType==1,]
HDL_EUR_1_ChrPos=HDL_EUR_1[HDL_EUR_1$VarType==2,]
HDL_EUR_1_ID=HDL_EUR_1[HDL_EUR_1$VarType==3,]
A=unlist(strsplit(HDL_EUR_1_ID$MarkerName, ":" )) 
B=matrix(A,ncol=3, byrow=TRUE)
B=data.frame(B)
names(B) <- c("CHROM","POS_b37","Type")
B=B[,c("CHROM","POS_b37")]
HDL_EUR_1_ID=cbind(B,HDL_EUR_1_ID)
A=unlist(strsplit(HDL_EUR_1_ChrPos$MarkerName, ":" )) 
B=matrix(A,ncol=2, byrow=TRUE)
B=data.frame(B)
names(B) <- c("CHROM","POS_b37")
HDL_EUR_1_ChrPos=cbind(B,HDL_EUR_1_ChrPos)

HDL_EUR_1_ID=HDL_EUR_1_ID[,c("CHROM", "POS_b37", "MarkerName", "Allele1", "Allele2", "Freq1", "Effect","StdErr","IntEffect", "IntStdErr","IntCov","ChiSq2df", "P-value","N","VarType")]
HDL_EUR_1_ChrPos=HDL_EUR_1_ChrPos[,c("CHROM", "POS_b37", "MarkerName", "Allele1", "Allele2", "Freq1", "Effect","StdErr","IntEffect", "IntStdErr","IntCov","ChiSq2df", "P-value","N","VarType")]
HDL_EUR_1_rs=HDL_EUR_1_rs[,c("MarkerName", "Allele1", "Allele2", "Freq1", "Effect","StdErr","IntEffect", "IntStdErr","IntCov","ChiSq2df", "P-value","N", "VarType")]
HDL_EUR_ID=merge(HDL_EUR_1_ID,HDL_EUR_2,by=c("CHROM","POS_b37"))
HDL_EUR_ChrPos=merge(HDL_EUR_1_ChrPos,HDL_EUR_2,by=c("CHROM","POS_b37"))
HDL_EUR_rs=merge(HDL_EUR_1_rs,HDL_EUR_2,by.x="MarkerName", by.y="rsID")
HDL_EUR_rs$rsID=HDL_EUR_rs$MarkerName
HDL_EUR=rbind(HDL_EUR_ChrPos,HDL_EUR_ID,HDL_EUR_rs) # curently there are 8344809 rows.
names(HDL_EUR)[14] <- c("TotalSampleSize")
names(HDL_EUR)[19] <- c("N")
names(HDL_EUR)[13]<-c("Pvalue")
HDL_EUR$FreqMin=ifelse(HDL_EUR$Freq1<0.5,HDL_EUR$Freq1, 1-HDL_EUR$Freq1)              
HDL_EUR$POOLED_Min_AF=ifelse(HDL_EUR$POOLED_ALT_AF<0.5,HDL_EUR$POOLED_ALT_AF, 1-HDL_EUR$POOLED_ALT_AF)
HDL_EUR=HDL_EUR[abs(HDL_EUR$FreqMin-HDL_EUR$POOLED_Min_AF)<0.15,]  # 8272473 rows
HDL_EUR=HDL_EUR[HDL_EUR$VarType==3 | (HDL_EUR$Allele1==HDL_EUR$REF & HDL_EUR$Allele2==HDL_EUR$ALT) | (HDL_EUR$Allele1==HDL_EUR$ALT & HDL_EUR$Allele2==HDL_EUR$REF),] #8233711 rows
A=HDL_EUR[HDL_EUR$VarType==3 | (HDL_EUR$Allele1==HDL_EUR$REF & HDL_EUR$Allele2==HDL_EUR$ALT) | (HDL_EUR$Allele1==HDL_EUR$ALT & HDL_EUR$Allele2==HDL_EUR$REF),] 
HDL_EUR=HDL_EUR[HDL_EUR$VarType<3 | (HDL_EUR$VarType==3 & !((HDL_EUR$REF %in% c("A","C","T","G")) & (HDL_EUR$ALT %in% c("A","C","T","G")))),]    # there are 8215940 rows
A=HDL_EUR[duplicated(HDL_EUR$MarkerName),]  # there 6169 duplicate rows
HDL_EUR=HDL_EUR[!(HDL_EUR$MarkerName %in% A$MarkerName) | ((HDL_EUR$MarkerName %in% A$MarkerName) & abs(HDL_EUR$FreqMin-HDL_EUR$POOLED_Min_AF)<0.01),] # there are 8209791 rows
HDL_EUR=HDL_EUR[!duplicated(HDL_EUR$MarkerName),]  # drop further 1651 duplicate rows. Now there are 8208140 rows finally. These variants may still be problematic and need to be further checked if there are any findings.
A=HDL_EUR[HDL_EUR$VarType < 3,]
B=HDL_EUR[HDL_EUR$VarType == 3,]
A$Effect_new=ifelse(A$Allele1 == A$ALT, A$Effect, -A$Effect)
B$Effect_new=ifelse(nchar(B$ALT) < nchar(B$REF), B$Effect, -B$Effect)
HDL_EUR=rbind(A,B)

"testA"
HDL_EUR$x1 <- HDL_EUR$Effect_new/HDL_EUR$StdErr/sqrt(HDL_EUR$TotalSampleSize)
HDL_EUR$x2 <- HDL_EUR$METAL_Effect/HDL_EUR$METAL_StdErr/sqrt(HDL_EUR$N)
HDL_EUR$x1_se <- 1/sqrt(HDL_EUR$TotalSampleSize)
HDL_EUR$x2_se <- 1/sqrt(HDL_EUR$N)
names(HDL_EUR)[1]<-c("CHR")
names(HDL_EUR)[2]<-c("BP")
names(HDL_EUR)[16]<-c("SNP")
#names(HDL_EUR)[13]<-c("Pvalue")
HDL_EUR=HDL_EUR[as.numeric(HDL_EUR$N) >100000,]
HDL_EUR=HDL_EUR[as.numeric(HDL_EUR$TotalSampleSize) >50000,]
Alcloci=fread("/mnt/rstor/SOM_EPBI_XXZ10/xxz10/gene-life/SmkAlcGWAS/AlcSNPs.txt")
HDL_EUR=HDL_EUR[!HDL_EUR$SNP %in% Alcloci$RSID,]
A=HDL_EUR[as.numeric(HDL_EUR$pvalue_GC) < 5e-8,]

HDL_EUR_MR<-fread("/mnt/rstor/SOM_EPBI_XXZ10/xxz10/gene-life/AlcoholLipids/dataanalysis/LDL_TrEtn_CurDring_1000GPlink.clumped")
HDL_EUR_MR<-merge(A, HDL_EUR_MR, by="SNP")
A=HDL_EUR[as.numeric(HDL_EUR$pvalue_GC) > 0.05 & as.numeric(HDL_EUR$Pvalue) > 0.05,]
rho=cor(A$Effect_new/A$StdErr,A$METAL_Effect/A$METAL_StdErr)
MR1<-MR_pleio("x2","x1","x2_se","x1_se",as.data.frame(HDL_EUR_MR),SignifThreshold=0.05,rho=rho)  
A <-SearchPleio("x2","x1","x2_se","x1_se",as.data.frame(HDL_EUR),rho,MR1$CausalEstimate,MR1$SdCausalEstimate)
HDL_EUR$PleioP_MR=A$pleio_p
HDL_EUR$Pleiobeta=A$beta
HDL_EUR$PleioSE=A$SE
A=HDL_EUR[HDL_EUR$PleioP_MR<5e-8 | HDL_EUR$PleioP_MR_rho1<5e-8,] 
HDL_EUR$CHR=as.numeric(as.character(HDL_EUR$CHR))
HDL_EUR$BP=as.numeric(as.character(HDL_EUR$BP))
bottom=HDL_EUR%>%dplyr::select(SNP,CHR,BP,PleioP_MR,Pleiobeta,PleioSE,Allele1,Allele2,Freq1,IntEffect,IntStdErr,TotalSampleSize,METAL_Effect,METAL_StdErr)
HDL_EUR$IsPleiotropic=ifelse(HDL_EUR$PleioP_MR<5e-8,"1","0")
HDL_EUR$IsPleiotropic=ordered(HDL_EUR$IsPleiotropic,levels=c("1","0"))
lower_highlight=HDL_EUR%>%filter(PleioP_MR<5e-8)%>%filter(SNP!="NA")
HDL_EUR$pvalue_GC=as.numeric(HDL_EUR$pvalue_GC)
HDL_EUR$Pvalue=as.numeric(HDL_EUR$Pvalue)
HDL_EUR_scatter=HDL_EUR%>%filter(PleioP_MR<5e-8|pvalue_GC < 5e-8 | Pvalue<5e-8)

HDL_EUR_scatter=HDL_EUR_scatter %>% dplyr::select(x1,x2,x1_se,x2_se,IsPleiotropic,PleioP_MR,SNP,CHR,BP)
names(top)[1]=c("rsid")
names(bottom)[1]=c("rsid")
names(top)[2]=c("chr")
names(bottom)[2]=c("chr")
names(top)[3]=c("pos")
names(bottom)[3]=c("pos")
names(top)[4]=c("pvalue")
names(bottom)[4]=c("pvalue")

saveRDS(top,"top1.rds")
saveRDS(bottom,"bottom1.rds")
saveRDS(HDL_EUR_scatter,"scatter1.rds")
saveRDS(MR1,"MR1.rds")
