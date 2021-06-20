
#------------------python in R------------------------
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "miniconda3/bin/python")
reticulate::py_config()
use_condaenv("r_reticulate")

#run python code line by line in R console
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

#create an initial population
pop = sim.Population(size = [70,64,50,56,80],
                    ploidy = 2, loci = [6939,5171,3570],
                    chromNames = ['chr1','chr2','chr3'],
                    alleleNames = ['A','C','G','T'], subPopNames = ['A','B','C','D','E'],
                    infoFields = ['ind_id','father_id','mother_id'])


#set genotypes by allele frequency 
for i in range(5):
  fr = [random.random() for x in range(4)] #alleles 0,1,2 and3.
  sim.initGenotype(pop, freq=[0, fr[1], 1-fr[1], 0], subPops=[i],
                   loci=range(2613)) #alleles C, G
  sim.initGenotype(pop, freq=[fr[0], 1-fr[0], 0, 0], subPops=[i],
                   loci=range(2613,5226)) #alleles A, C
  sim.initGenotype(pop, freq=[fr[0], 0, 1-fr[0], 0], subPops=[i],
                   loci=range(5226,7839)) #alleles A,G
  sim.initGenotype(pop, freq=[0,0, fr[2], 1-fr[2]], subPops=[i],
                   loci=range(7839, 10452)) #alleles G, T
  sim.initGenotype(pop, freq=[0,fr[1],0,1-fr[1]], subPops=[i],
                   loci=range(10452, 13065)) #alleles C, T      
  sim.initGenotype(pop, freq=[fr[0],0, 0, 1-fr[0]], subPops=[i],
                   loci=range(13065, 15680)) #alleles A, T

saveCSV(pop, filename='initialPopulation.csv',
infoFields=['ind_id', 'father_id', 'mother_id'], sep= ",")

#save as a simupop object
pop.save('initialPopulation.pop')


#evolve simulated data using random mating (mendelian laws), mating happens in each subpopualtion by default 
#5th generation
sim.IdTagger().reset(startID=1)
pop.evolve(
  initOps=[
    sim.InitSex(),
    sim.IdTagger(),
    sim.PedigreeTagger(output='>>pedigree.txt'),
    ],
  matingScheme=sim.RandomMating(ops=[
    # give new born an ID
    sim.IdTagger(), 
    # track parents of each individual
    sim.PedigreeTagger(output='>>pedigree.txt'),
    sim.MendelianGenoTransmitter()]
  ),
  gen=5
)

#10th generation
saveCSV(pop, filename='5thGen.csv',
        infoFields=['ind_id', 'father_id', 'mother_id'], sep= ",")
pop.evolve(
  initOps=[
    sim.InitSex(),
    sim.IdTagger(),
    sim.PedigreeTagger(output='>>pedigree.txt'),
    ],
  matingScheme=sim.RandomMating(ops=[
    # give new born an ID
    sim.IdTagger(), 
    # track parents of each individual
    sim.PedigreeTagger(output='>>pedigree.txt'),
    sim.MendelianGenoTransmitter()]
  ),
  gen=5
)

saveCSV(pop, filename='10thGen.csv',
        infoFields=['ind_id', 'father_id', 'mother_id'], sep= ",")

#20th generation
pop.evolve(
  initOps=[
    sim.InitSex(),
    sim.IdTagger(),
    sim.PedigreeTagger(output='>>pedigree.txt'),
    ],
  matingScheme=sim.RandomMating(ops=[
    # give new born an ID
    sim.IdTagger(), 
    # track parents of each individual
    sim.PedigreeTagger(output='>>pedigree.txt'),
    sim.MendelianGenoTransmitter()]
  ),
  gen=10
)

saveCSV(pop, filename='20thGen.csv',
        infoFields=['ind_id', 'father_id', 'mother_id'], sep= ",")

#---------------------------back to R-----------------------------
library(snpReady)
library(impute)
library(dplyr)
library(tidyverse)

initialPopulation= read.csv("initialPopulation.csv", sep =",")
#name each subpopulation (by default subpopulation mate within each subpop and stay in the same order)
subPopulations = rep(c("A","B","C","D","E"),times=c(70,64,50,56,80))
initialPopulation = cbind(subPopulations, initialPopulation)

initialPopulation.clean = initialPopulation[,-(1:6)]
#combine adjacent columns(diplid), to create SNPs (e.g AC,CA,AA and CC) 
initialPopulation.clean = data.frame( initialPopulation.clean[0],
                                    mapply( paste0, initialPopulation.clean[][c(T,F)], 
                                            initialPopulation.clean[][c(F,T)], 
                                            MoreArgs= list(sep = "") ))
initialPopulation.clean = as.matrix(initialPopulation.clean)

#transform raw SNPs to 0,1,2
initialPopulation.clean=raw.data(data = initialPopulation.clean,frame="wide",base=TRUE,imput = TRUE,
                               imput.type = "mean", call.rate = 0, maf=0, sweep.sample = 1)
initialPopulation.clean = initialPopulation.clean$M.clean

#PCA on initial population
pca.initialPop= prcomp(initialPopulation.clean,scale. = FALSE, rank. = 10)
library(ggfortify)

#plot PC1,2 of PCA on initial population
autoplot(pca.initialPop, data= initialPopulation, colour= "subPopulations",
         label.size=3, main = "PCA on the simulated initial Pop")


#do same process for 5th,10th and 20th generation 
fifthGen= read.csv("5thGen.csv", sep =",")
fifthGen = cbind(subPopulations, fifthGen)

fifthGen.clean = fifthGen[,-(1:6)]

fifthGen.clean = data.frame( fifthGen.clean[0],
                                      mapply( paste0, fifthGen.clean[][c(T,F)], 
                                              fifthGen.clean[][c(F,T)], 
                                              MoreArgs= list(sep = "") ))
fifthGen.clean = as.matrix(fifthGen.clean)

fifthGen.clean=raw.data(data = fifthGen.clean,frame="wide",base=TRUE,imput = TRUE,
                                 imput.type = "mean", call.rate = 0, maf=0, sweep.sample = 1)
fifthGen.clean = fifthGen.clean$M.clean

#PCA on 5th generation
pca.5thGen= prcomp(fifthGen.clean,scale. = FALSE, rank. = 10)
autoplot(pca.5thGen, data= fifthGen, colour= "subPopulations",
         label.size=3, main = "PCA on the 5th Generation (simulated)")


tenthGen= read.csv("10thGen.csv", sep =",")
tenthGen = cbind(subPopulations, tenthGen)

tenthGen.clean = tenthGen[,-(1:6)]

tenthGen.clean = data.frame( tenthGen.clean[0],
                             mapply( paste0, tenthGen.clean[][c(T,F)], 
                                     tenthGen.clean[][c(F,T)], 
                                     MoreArgs= list(sep = "") ))
tenthGen.clean = as.matrix(tenthGen.clean)

tenthGen.clean=raw.data(data = tenthGen.clean,frame="wide",base=TRUE,imput = TRUE,
                        imput.type = "mean", call.rate = 0, maf=0, sweep.sample = 1)
tenthGen.clean = tenthGen.clean$M.clean

#PCA on 10th generation
pca.10thGen= prcomp(tenthGen.clean,scale. = FALSE, rank. = 10)
autoplot(pca.10thGen, data= tenthGen, colour= "subPopulations",
         label.size=3, main = "PCA on the 10th Generation (simulated)")


twentiethGen= read.csv("20thGen.csv", sep =",")
twentiethGen = cbind(subPopulations, twentiethGen)

twentiethGen.clean = twentiethGen[,-(1:6)]

twentiethGen.clean = data.frame( twentiethGen.clean[0],
                             mapply( paste0, twentiethGen.clean[][c(T,F)], 
                                     twentiethGen.clean[][c(F,T)], 
                                     MoreArgs= list(sep = "") ))
twentiethGen.clean = as.matrix(twentiethGen.clean)

twentiethGen.clean=raw.data(data = twentiethGen.clean,frame="wide",base=TRUE,imput = TRUE,
                        imput.type = "mean", call.rate = 0, maf=0, sweep.sample = 1)
twentiethGen.clean = twentiethGen.clean$M.clean

#PCA on 20th generation
pca.20thGen= prcomp(twentiethGen.clean,scale. = FALSE, rank. = 10)
autoplot(pca.20thGen, data= twentiethGen, colour= "subPopulations",
         label.size=3, main = "PCA on the 20th Generation (simulated)")
