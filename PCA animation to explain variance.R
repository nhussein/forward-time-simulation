library(ggplot2)
library(stats)

set.seed(10)             
x <- 1:100              
erx <- rnorm(100, 0, 10) 
ery <- rnorm(100, 0, 10) 
y <- 30 + 2*x     
X <- x + erx         
Y <- y + ery

toyData <- data.frame(X = X, Y = Y)

PC = prcomp(toyData)$rotation
PofVE = prcomp(toyData)$sdev^2

m1 = PC[2,1]/PC[1,1]
b1 = mean(Y)-(PC[2,1]/PC[1,1])*mean(X)

m2 = PC[2,2]/PC[1,2]
b2 = mean(Y)-(PC[2,2]/PC[1,2])*mean(X)

g <- ggplot(toyData, mapping = aes(X,Y))+
  geom_point()+
  geom_point(aes(x=mean(X),y=mean(Y), col = "red" ))+
  geom_abline(slope = m1, intercept = b1, col = "blue")+
  geom_abline(slope = m2, intercept = b2, col = "blue")+
  coord_fixed() +
  theme(legend.position = "none")+
  ggtitle("PCA Maximizing Variance")

  
tab <- as.data.frame(c(
  "Variance by PC1" = PofVE[1],
  "Variance by PC2" = PofVE[2])
  )

# using gridExtra
library(gridExtra)
p_tab <- tableGrob(unname(tab))
final <- grid.arrange(g, p_tab, ncol = 2)

ggsave(final, file=paste0("final_PCA",".png"))
  
  
M=c(mean(X), mean(Y))
distance = sqrt((M[2]-toyData[,2])^2+(M[1]-toyData[,1])^2) 

iteration = c(1/10*pi, 1/5*pi, 3/10*pi, 2/5*pi, 1/2*pi, 
              3/5*pi, 7/10*pi, 4/5*pi, 9/10*pi, pi, 
              11/10*pi, 6/5*pi, 13/10*pi, 7/5*pi, 3/2*pi, 
              8/5*pi, 17/10*pi, 9/5*pi, 19/10*pi, 2*pi)

for(i in iteration ){
  
  g <- ggplot(toyData, mapping = aes(X,Y))+
    geom_point()+
    geom_point(aes(x=mean(X),y=mean(Y), col = "red" ))+
    geom_abline(slope = tan(i), intercept = mean(Y)-mean(X)*tan(i), col = "purple")+
    geom_abline(slope = -1/tan(i), intercept = mean(Y)-mean(X)*(-1/tan(i)), col = "purple")+
    coord_fixed() +
    theme(legend.position = "none")+
    ggtitle("PCA Maximizing Variance")
  
  projection1 = distance*cos(atan(toyData[,2]/toyData[,1]) - i)
  var1 = mean(projection1^2)
  
  projection2 = distance*sin(atan(toyData[,2]/toyData[,1]) - i)
  var2 = mean(projection2^2)
  
  tab <- as.data.frame(c(
    "Variance by axis1" = var1,
    "Variance by axis2" = var2)
  )
  p_tab <- tableGrob(unname(tab))
  final <- grid.arrange(g, p_tab, ncol = 2)
  
  ggsave(final, file=paste0("final",i ,".png"))

  dev.off()
}


ggplot(toyData, mapping = aes(X,Y))+
  geom_point()+
  geom_point(aes(x=mean(X),y=mean(Y), col = "red" ))+
  geom_abline(slope = tan(0.5), intercept = mean(Y)-mean(X)*tan(0.5), col = "purple")+
  geom_abline(slope = -1/tan(0.5), intercept = mean(Y)-mean(X)*(-1/tan(0.5)), col = "purple")+
  coord_fixed() +
  theme(legend.position = "none")+
  ggtitle("PCA Maximizing Variance")
  

                
M=c(mean(X), mean(Y))
distance = sqrt((M[2]-toyData[,2])^2+(M[1]-toyData[,1])^2) 
projection1 = distance*cos(atan(toyData[,2]/toyData[,1]) - 0.5)
var1 = mean(projection1^2)

distance = sqrt((M[2]-toyData[,2])^2+(M[1]-toyData[,1])^2) 
projection2 = distance*sin(atan(toyData[,2]/toyData[,1]) - 0.5)
var2 = mean(projection2^2)

#---------------------python---------------------------

library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "miniconda3/bin/python")
reticulate::py_config()
use_condaenv("r_reticulate")

repl_python()

import imageio
import os
import matplotlib.pyplot as plt

images = []
dir = "/home/hnoor/PCA"

files = os.listdir(dir)
sorted_files =  sorted(files)

for image in sorted_files:
  images.append(imageio.imread(dir + "/" + image))

imageio.mimsave("/home/hnoor/PCA.gif", images, duration = 0.5)
