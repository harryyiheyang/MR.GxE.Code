library(bigsnpr)
library(ldscR)
library(glue)
library(data.table)
data("EURLDSC")
data("AFRLDSC")
data("EASLDSC")
data("AMRLDSC")
HISLDSC=AMRLDSC
CrossLDSC=EURLDSC
Env=c("CD","RD","CS","ES")
Trait=c("HDL","LDL","TG")

iii=1
POP="EUR"
A=list()
for(i in 1:3){
for(j in 1:4){
intercept=NULL
w1=fread(glue("Interaction/{POP}_INT/{POP}_{Env[j]}_{Trait[i]}.txt"))
LDSC=get(glue("{POP}LDSC"))
w1=merge(w1,LDSC,by="SNP")
fit1=snp_ldsc(ld_score=w1$LDSC,ld_size=dim(w1)[1],chi2=w1$BETA^2/w1$SE^2,sample_size=w1$N,intercept=intercept,blocks=200)
fit1=data.frame(int=fit1[1],intse=fit1[2],h2=fit1[3],h2se=fit1[4],environment=Env[j],lipid=Trait[i])
A[[iii]]=fit1
iii=iii+1
print(fit1)
}
}

saveRDS(A,"Interaction/EURH2.rds")