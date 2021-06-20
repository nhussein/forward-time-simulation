#--------------------------Python in R---------------------------------------------
library(reticulate)
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

#initial population (HapMap chromosome 22)
initPop = sim.loadPopulation("initPop.pop")

#we have no special info assigned to the members of our initial population. so next we will add infofields.
initPop.setInfoFields(['ind_id','parent_idx'])

sim.IdTagger().reset(startID=1)
initPop.evolve(
  initOps=[
    sim.InitSex(),
    sim.IdTagger(),
    ],
  matingScheme=sim.RandomSelection(ops=[
    sim.ParentsTagger(infoFields='parent_idx'),
    sim.CloneGenoTransmitter(),
    ]),
  gen=2
)

saveCSV(initPop, filename='2ndGenAsexual.csv',
        infoFields=['ind_id','parent_idx'])

#use esc key to go back to R
secondGenAsexual= read.csv('2ndGenAsexual.csv')
subPopulations = read.csv('subpopulations.csv')
subPopulations = subPopulations[-(1:7),]

secondGenAsexual = cbind(subPopulations, secondGenAsexual)
names(secondGenAsexual)[names(secondGenAsexual) == "x"] <- "SubPopulation"
secondGenAsexual = secondGenAsexual[,-1]
secondGenAsexual[secondGenAsexual == 1] <- "A"
secondGenAsexual[secondGenAsexual == 2] <- "C"
secondGenAsexual[secondGenAsexual == 3] <- "G"
secondGenAsexual[secondGenAsexual == 4] <- "T"
secondGenAsexual[secondGenAsexual == 247] <- NA

secondGenAsexual.clean = secondGenAsexual[,-(1:5)]

secondGenAsexual.clean = data.frame( secondGenAsexual.clean[0],
                                    mapply( paste0, secondGenAsexual.clean[][c(T,F)], 
                                            secondGenAsexual.clean[][c(F,T)] ))
secondGenAsexual.clean = as.matrix(secondGenAsexual.clean)

'%ni%' <- Negate('%in%')
X = c("AA","AC","AT","AG","CC","CA","CT","CG","GG","GA","GT","GC","TT","TC","TA","TG")
secondGenAsexual.clean[secondGenAsexual.clean %ni% X  ] <- NA

secondGenAsexual.clean=raw.data(data = secondGenAsexual.clean,frame="wide",base=TRUE,imput = TRUE,
                               imput.type = "mean", call.rate = 0, maf=0, sweep.sample = 1)
secondGenAsexual.clean = secondGenAsexual.clean$M.clean

#PCA
pca.2ndGenAsexual = prcomp(secondGenAsexual.clean,scale. = FALSE, rank. = 10)
library(ggfortify)
autoplot(pca.2ndGenAsexual, data= secondGenAsexual,colour= "SubPopulation",
         label.size=3, main = "2nd Generation Asexual Mating")
