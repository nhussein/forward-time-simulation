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

CEU= fread("CEU.hmap.gz")   
CEU= data.matrix(t(CEU))
colnames(CEU) <- CEU[1,]

CEU = CEU[-(4:10),]
CEU = t(CEU)
CEU = CEU[CEU[,3]=="chr22",]
CEU = CEU[,-(2:4)]

MKK= fread("MKK.hmap.gz")   
MKK= data.matrix(t(MKK))
colnames(MKK) <- MKK[1,]

MKK = MKK[-(4:10),]
MKK = t(MKK)
MKK = MKK[MKK[,3]=="chr22",]
MKK = MKK[,-(2:4)]

ASW= fread("ASW.hmap.gz")   
ASW= data.matrix(t(ASW))
colnames(ASW) <- ASW[1,]

ASW = ASW[-(4:10),]
ASW = t(ASW)
ASW = ASW[ASW[,3]=="chr22",]
ASW = ASW[,-(2:4)]


CHB= fread("CHB.hmap.gz")   
CHB= data.matrix(t(CHB))
colnames(CHB) <- CHB[1,]

CHB = CHB[-(4:10),]
CHB = t(CHB)
CHB = CHB[CHB[,3]=="chr22",]
CHB = CHB[,-(2:4)]


CHD= fread("CHD.hmap.gz")   
CHD= data.matrix(t(CHD))
colnames(CHD) <- CHD[1,]

CHD = CHD[-(4:10),]
CHD = t(CHD)
CHD = CHD[CHD[,3]=="chr22",]
CHD = CHD[,-(2:4)]


JPT= fread("JPT.hmap.gz")   
JPT= data.matrix(t(JPT))
colnames(JPT) <- JPT[1,]

JPT = JPT[-(4:10),]
JPT = t(JPT)
JPT = JPT[JPT[,3]=="chr22",]
JPT = JPT[,-(2:4)]

LWK= fread("LWK.hmap.gz")   
LWK= data.matrix(t(LWK))
colnames(LWK) <- LWK[1,]

LWK = LWK[-(4:10),]
LWK = t(LWK)
LWK = LWK[LWK[,3]=="chr22",]
LWK = LWK[,-(2:4)]

TSI= fread("TSI.hmap.gz")   
TSI= data.matrix(t(TSI))
colnames(TSI) <- TSI[1,]

TSI = TSI[-(4:10),]
TSI = t(TSI)
TSI = TSI[TSI[,3]=="chr22",]
TSI = TSI[,-(2:4)]

YRI= fread("YRI.hmap.gz")   
YRI= data.matrix(t(YRI))
colnames(YRI) <- YRI[1,]

YRI = YRI[-(4:10),]
YRI = t(YRI)
YRI = YRI[YRI[,3]=="chr22",]
YRI = YRI[,-(2:4)]

MEX= fread("MEX.hmap.gz")   
MEX= data.matrix(t(MEX))
colnames(MEX) <- MEX[1,]

MEX = MEX[-(4:10),]
MEX = t(MEX)
MEX = MEX[MEX[,3]=="chr22",]
MEX = MEX[,-(2:4)]

GIH= fread("GIH.hmap.gz")   
GIH= data.matrix(t(GIH))
colnames(GIH) <- GIH[1,]

GIH = GIH[-(4:10),]
GIH = t(GIH)
GIH = GIH[GIH[,3]=="chr22",]
GIH = GIH[,-(2:4)]

CEU = t(CEU)
MKK = t(MKK)
ASW = t(ASW)
CHB = t(CHB)
CHD = t(CHD)
JPT = t(JPT)
LWK = t(LWK)
TSI = t(TSI)
YRI = t(YRI)
MEX = t(MEX)
GIH = t(GIH)

#common snps
x = Reduce(intersect, list(ASW,MEX,LWK,MKK,YRI,
                           CEU,TSI,GIH,CHB,JPT,CHD))

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

#combine All observations and add subpopulation vector
All=rbind(ASW,MEX,LWK,MKK,YRI,CEU,TSI,GIH,CHB,JPT,CHD)

pops=c(rownames(CEU),rownames(CHB),rownames(ASW),rownames(LWK),rownames(CHD),
       rownames(GIH),rownames(JPT),rownames(MEX),rownames(MKK),rownames(TSI),rownames(YRI))

All=cbind(All, pops)
All[rownames(All)%in%rownames(ASW),14200] <- 'ASW'
All[rownames(All)%in%rownames(CEU),14200] <- 'CEU'
All[rownames(All)%in%rownames(CHB),14200] <- 'CHB'
All[rownames(All)%in%rownames(LWK),14200] <- 'LWK'
All[rownames(All)%in%rownames(CHD),14200] <- 'CHD'
All[rownames(All)%in%rownames(GIH),14200] <- 'GIH'
All[rownames(All)%in%rownames(JPT),14200] <- 'JPT'
All[rownames(All)%in%rownames(MEX),14200] <- 'MEX'
All[rownames(All)%in%rownames(MKK),14200] <- 'MKK'
All[rownames(All)%in%rownames(TSI),14200] <- 'TSI'
All[rownames(All)%in%rownames(YRI),14200] <- 'YRI'

write.csv(All[,14200], "subpopulations.csv")

subpopNames = c("ASW","MEX","LWK","MKK","YRI","CEU","TSI","GIH","CHB","JPT","CHD")
indPerSubpops = vector()

#how many individuals in each subpop.
for(i in subpopNames) {
  indPerSubpops[i] = sum(All[,14200] == i)             
}

write.csv(All, "All_chr22.csv")

All = read.csv('All_chr22.csv')

All[All == "NN"] <- NA
All = All[-(1:7),]
rownames(All) <- All[,1]
All = All[,-1]

chr22= raw.data(data = t(All[,-14200]),frame="wide",base=TRUE,imput = TRUE,
                imput.type = "mean", call.rate = 0.95, maf=0.05, sweep.sample = 1)
chr22 = chr22$M.clean
chr22 = t(chr22)

pca.all = prcomp(chr22,scale. = FALSE, rank. = 10)

colnames(All)[14200] <- "Populations" 

library(ggfortify)
autoplot(pca.all, data= All,colour = "Populations",label.size=3, 
         main = "PCA on the HapMap3 Populations chromosome22")

#european american 
TSI = subset(All, Populations == "TSI")
CEU = subset(All, Populations == "CEU")
MEX = subset(All, Populations == "MEX")
EA = rbind(TSI, CEU, MEX)
chr22.EA= raw.data(data = t(EA[,-14200]),frame="wide",base=TRUE,imput = TRUE,
                imput.type = "mean", call.rate = 0.95, maf=0.05, sweep.sample = 1)
chr22.EA = chr22.EA$M.clean
chr22.EA = t(chr22.EA)

pca.EA = prcomp(chr22.EA,scale. = FALSE, rank. = 10)
autoplot(pca.EA, data= EA,colour = "Populations",label.size=3, 
         main = "PCA on European American group")
autoplot(pca.EA,x=3,y=4,data= EA,colour = "Populations",label.size=3, 
         main = "PCA on European American group")

#African american 
YRI = subset(All, Populations == "YRI")
LWK = subset(All, Populations == "LWK")
ASW = subset(All, Populations == "ASW")
MKK = subset(All, Populations == "MKK")
GIH = subset(All, Populations == "GIH")
AA = rbind(YRI, LWK, ASW, MKK, GIH)
chr22.AA= raw.data(data = t(AA[,-14200]),frame="wide",base=TRUE,imput = TRUE,
                   imput.type = "mean", call.rate = 0.95, maf=0.05, sweep.sample = 1)
chr22.AA = chr22.AA$M.clean
chr22.AA = t(chr22.AA)

pca.AA = prcomp(chr22.AA,scale. = FALSE, rank. = 10)
autoplot(pca.AA, data= AA,colour = "Populations",label.size=3, 
         main = "PCA on African American group")
autoplot(pca.AA,x=3,y=4,data= AA,colour = "Populations",label.size=3, 
         main = "PCA on African American group")



All1 = raw.data(data = t(All[,14200]),frame="wide",base=TRUE,imput = TRUE, outfile = "structure",
                imput.type = "mean", call.rate = 0.95, maf=0.05, sweep.sample = 1)
All1 = All1$M.clean
All1 = t(All1)

write.csv(All1, "raw_chr22.csv")

#-------------------------------------------- Python
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "miniconda3/bin/python")

use_condaenv("r_reticulate")

repl_python()

import simuOpt
simuOpt.setOptions(alleleType='short', quiet=True)
import simuPOP as sim
import pandas as pd
import numpy as np
import collections
np.set_printoptions(suppress=True, precision=3)
from simuPOP.utils import saveCSV
import random

#make hapmap (initial population) a simupop object 

initPop = sim.Population(size=[71, 71, 83, 171, 163, 162, 77, 83, 82, 82, 70], 
                        subPopNames=["ASW","MEX","LWK","MKK","YRI","CEU",
                                     "TSI","GIH","CHB","JPT","CHD"],
                        ploidy=2, loci=14199)

df = pd.read_csv("raw_chr22.csv", index_col=0)
df.to_csv("raw_chr22.csv", index=False)
genotypes = np.array(pd.read_csv('raw_chr22.csv'))
genotypes = genotypes.tolist()

for i, ind in enumerate(initPop.individuals()):
  ind.setGenotype(genotypes[i])

initPop.save('initPop.pop')

saveCSV(initPop, "initPop.csv")
