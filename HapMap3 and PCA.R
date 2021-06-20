library(R.utils)
library(tidyr)
library(data.table)
library(tidyverse)
library(snpReady)
library(impute)
library(dplyr)
library(kernlab)
library(ggfortify)
library(stringr)
library(ggplot2)
library(lubridate)
library(generics)


#Load HapMap data for each population 

CEU=fread("CEU.hmap.gz")     
MKK=fread("MKK.hmap.gz")     
ASW=fread("ASW.hmap.gz")     
CHB=fread("CHB.hmap.gz")     
CHD=fread("CHD.hmap.gz")     
JPT=fread("JPT.hmap.gz")     
LWK=fread("LWK.hmap.gz")     
TSI=fread("TSI.hmap.gz")     
YRI=fread("YRI.hmap.gz")     
MEX=fread("MEX.hmap.gz")     
GIH=fread("GIH.hmap.gz")

#remove unnecessary columns, and flip matrix (columns: SNPs, rows: Individuals)

CEU=data.matrix(t(CEU[,-(2:11)]))
MKK=data.matrix(t(MKK[,-(2:11)]))
JPT=data.matrix(t(JPT[,-(2:11)]))
CHB=data.matrix(t(CHB[,-(2:11)]))
MEX=data.matrix(t(MEX[,-(2:11)]))
GIH=data.matrix(t(GIH[,-(2:11)]))
ASW=data.matrix(t(ASW[,-(2:11)]))
CHD=data.matrix(t(CHD[,-(2:11)]))
TSI=data.matrix(t(TSI[,-(2:11)]))
YRI=data.matrix(t(YRI[,-(2:11)]))
LWK=data.matrix(t(LWK[,-(2:11)]))

#name columns by SNP name(rs#)

colnames(CEU) <- CEU[1,]
colnames(MKK) <- MKK[1,]
colnames(JPT) <- JPT[1,]
colnames(CHB) <- CHB[1,]
colnames(MEX) <- MEX[1,]
colnames(GIH) <- GIH[1,]
colnames(ASW) <- ASW[1,]
colnames(CHD) <- CHD[1,]
colnames(TSI) <- TSI[1,]
colnames(YRI) <- YRI[1,]
colnames(LWK) <- LWK[1,]

#reduce all populations to all common snps 

x = Reduce(intersect, list(ASW,MEX,LWK,MKK,YRI,CEU,TSI,GIH,CHB,JPT,CHD))

MKK = MKK[-1,colnames(MKK)%in%x]
ASW = ASW[-1,colnames(ASW)%in%x]
LWK = LWK[-1,colnames(LWK)%in%x]
YRI = YRI[-1,colnames(YRI)%in%x]
MEX = MEX[-1,colnames(MEX)%in%x]
CEU = CEU[-1,colnames(CEU)%in%x]
TSI = TSI[-1,colnames(TSI)%in%x]
GIH = GIH[-1,colnames(GIH)%in%x]
CHB = CHB[-1,colnames(CHB)%in%x]
JPT = JPT[-1,colnames(JPT)%in%x]
CHD = CHD[-1,colnames(CHD)%in%x]

#transform (raw) SNPs to 0,1,2 and the mean for NAs and save transformed data

ASW[ASW=="NN"]<-NA
ASW=raw.data(data = ASW,frame="wide",base=TRUE,imput = TRUE,
             imput.type = "mean", call.rate = 0, maf=0, sweep.sample = 1)
ASW=ASW$M.clean
saveRDS(ASW, file = "ASW_r")

MKK[MKK=="NN"]<-NA
MKK=raw.data(data = MKK,frame="wide",base=TRUE,imput = TRUE,
             imput.type = "mean", call.rate = 0, maf=0, sweep.sample = 1)
MKK=MKK$M.clean
saveRDS(MKK, file = "MKk_r")

LWK[LWK=="NN"]<-NA
LWK=raw.data(data = LWK,frame="wide",base=TRUE,imput = TRUE,
             imput.type = "mean", call.rate = 0, maf=0, sweep.sample = 1)

LWK=LWK$M.clean
saveRDS(LWK, file = "LWK_r")

YRI[YRI=="NN"]<-NA
YRI=raw.data(data = YRI,frame="wide",base=TRUE,imput = TRUE,
             imput.type = "mean", call.rate = 0, maf=0, sweep.sample = 1)
YRI=YRI$M.clean
saveRDS(YRI, file = "YRI_r")

MEX[MEX=="NN"]<-NA
MEX=raw.data(data = MEX,frame="wide",base=TRUE,imput = TRUE,
             imput.type = "mean", call.rate = 0, maf=0, sweep.sample = 1)
MEX=MEX$M.clean
saveRDS(MEX, file = "MEX_r")

CEU[CEU=="NN"]<-NA
CEU=raw.data(data = CEU,frame="wide",base=TRUE,imput = TRUE,
             imput.type = "mean", call.rate = 0, maf=0, sweep.sample = 1)
CEU=CEU$M.clean
saveRDS(CEU, file = "CEU_r")

TSI[TSI=="NN"]<-NA
TSI=raw.data(data = TSI,frame="wide",base=TRUE,imput = TRUE,
             imput.type = "mean", call.rate = 0, maf=0, sweep.sample = 1)
TSI=TSI$M.clean
saveRDS(TSI, file = "TSI_r")

GIH[GIH=="NN"]<-NA
GIH=raw.data(data = GIH,frame="wide",base=TRUE,imput = TRUE,
             imput.type = "mean", call.rate = 0, maf=0, sweep.sample = 1)
GIH=GIH$M.clean
saveRDS(GIH, file = "GIH_r")

CHB[CHB=="NN"]<-NA
CHB=raw.data(data = CHB,frame="wide",base=TRUE,imput = TRUE,
             imput.type = "mean", call.rate = 0, maf=0, sweep.sample = 1)
CHB=CHB$M.clean
saveRDS(CHB, file = "CHB_r")

JPT[JPT=="NN"]<-NA
JPT=raw.data(data = JPT,frame="wide",base=TRUE,imput = TRUE,
             imput.type = "mean", call.rate = 0, maf=0, sweep.sample = 1)
JPT=JPT$M.clean
saveRDS(JPT, file = "JPT_r")

CHD[CHD=="NN"]<-NA
CHD=raw.data(data = CHD,frame="wide",base=TRUE,imput = TRUE,
             imput.type = "mean", call.rate = 0, maf=0, sweep.sample = 1)
CHD=CHD$M.clean
saveRDS(CHD, file = "CHD_r")

#---------------------------------------------------------------------------

CEU = readRDS('CEU_r')
JPT = readRDS('JPT_r')
CHB = readRDS('CHB_r')
MEX = readRDS('MEX_r')
GIH = readRDS('GIH_r')
ASW = readRDS('ASW_r')
CHD = readRDS('CHD_r')
TSI = readRDS('TSI_r')
YRI = readRDS('YRI_r')
LWK = readRDS('LWK_r')
MKK = readRDS('MKk_r')

#combine all populations in one dataset

All=rbind(ASW,MEX,LWK,MKK,YRI,CEU,TSI,GIH,CHB,JPT,CHD)
saveRDS(All, file = "All")
All = readRDS("All")

#create population names column

pops=c(rownames(CEU),rownames(CHB),rownames(ASW),rownames(LWK),rownames(CHD),
       rownames(GIH),rownames(JPT),rownames(MEX),rownames(MKK),rownames(TSI),rownames(YRI))

#add column to full dataset and assign each individul to its population
All1=cbind(All, pops)
All1[rownames(All1)%in%rownames(ASW),1047056] <- 'ASW'
All1[rownames(All1)%in%rownames(CEU),1047056] <- 'CEU'
All1[rownames(All1)%in%rownames(CHB),1047056] <- 'CHB'
All1[rownames(All1)%in%rownames(LWK),1047056] <- 'LWK'
All1[rownames(All1)%in%rownames(CHD),1047056] <- 'CHD'
All1[rownames(All1)%in%rownames(GIH),1047056] <- 'GIH'
All1[rownames(All1)%in%rownames(JPT),1047056] <- 'JPT'
All1[rownames(All1)%in%rownames(MEX),1047056] <- 'MEX'
All1[rownames(All1)%in%rownames(MKK),1047056] <- 'MKK'
All1[rownames(All1)%in%rownames(TSI),1047056] <- 'TSI'
All1[rownames(All1)%in%rownames(YRI),1047056] <- 'YRI'

All1 = as.data.frame(All1)
All1$pops=as.factor(as.character(All1$pops))

#save full data set

saveRDS(All1, file = "All1")
All1 = readRDS("All1")

#Apply PCA on full dataset (without population variable)
pca.all = prcomp(All,scale. = FALSE, rank. = 10)
saveRDS(pca.all, file = "pca.all")

All1 = readRDS("All1")
pca.all = readRDS("pca.all")

All1 %>% 
  rename(populations = pops)

#Plot PC1,2
library(ggfortify)
autoplot(pca.all, data= All1,colour= "population",label.size=3, main = "PCA on the HAPMAP3 data")
legend("topright", legend=levels(All1$populations),col = 1:11, pch = 19, bty = "n", cex = 0.6)

#Plot PC2,3
autoplot(pca.all, data= All1, x=2, y=3, colour="population", label.size=3)
legend("topright", legend=levels(All1$populations),col = 1:11, pch = 19, bty = "n", cex = 0.6)