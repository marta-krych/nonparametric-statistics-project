## Useful libraries:
library(mvtnorm)
library(MVN)
library(rgl)
library(car)
library(dbscan)
library(cluster)
library(fields)

## Set the parameters:
B = 1e5
seed = 26111992

## Import the dataset:
dataset <- read.csv("parkinson.csv")
View(dataset)

colnames(dataset) <- dataset[1,]
dataset <- dataset[-1,]
rownames(dataset) <- 1:dim(dataset)[1]

sani <- which(startsWith(dataset$`Participant code`, "HC"))
malati <- which(startsWith(dataset$`Participant code`, "PD"))
rischio <- which(startsWith(dataset$`Participant code`, "RBD"))

dataset <- dataset[,-1]
dataset <- type.convert(dataset, as.is=T)
dataset <- dataset[,-c(7,8,10)]


## Let's start to use the depth measure:
parkinson_ill <- dataset[1:30,]
RBD <-  dataset[31:80,]
health <- dataset[81:130,]


## Extract the main differences between PD and HC:
#   • but first i build the test:
perm_t_test=function(x,y,iter=1e3){
  
  T0=abs(mean(x)-mean(y))  # define the test statistic
  T_stat=numeric(iter) # a vector to store the values of each iteration
  x_pooled=c(x,y) # pooled sample
  n=length(x_pooled)
  n1=length(x)
  
  for(perm in 1:iter){ # loop for conditional MC
    # permutation:
    permutation <- sample(1:n)
    x_perm <- x_pooled[permutation]
    x1_perm <- x_perm[1:n1]
    x2_perm <- x_perm[(n1+1):n]
    # test statistic:
    T_stat[perm] <- abs(mean(x1_perm) - mean(x2_perm))
    
  }
  
  # p-value
  p_val <- sum(T_stat>=T0)/iter
  return(p_val)
}
#   •let's go
# B2 = 1e4
p.value <- numeric(24)
#pb=progress_bar$new(total=B)
#pb$tick(0)
#set.seed(seed)
cont <- 1
for (i in 38:61){
  x.1 <- parkinson_ill[,i]
  x.2 = health[,i]
  #       x.2 = RBD[which(cluster.ew==2), i]
  p.value[cont] <- perm_t_test(x.1,x.2,iter=1e3)
  cont <- cont + 1
}
p.value 

# P-values: 0.798 0.000 0.867 0.000 0.031 0.680 0.223 0.412 0.677 0.274 0.802 0.685 0.521 0.003 0.538 0.000
# 0.055 0.521 0.010 0.245 0.321 0.108 0.088 0.115 . The main differences are in Rate of speech timing, Duration of pause ntervals, Duration of voiced intervals
# and "Rate of speech timing (-/min).1" , "Duration of pause intervals (ms).1","Duration of pause intervals (ms).1" , "Duration of voiced intervals (ms).1" ,
# "Duration of unvoiced stops (ms).1"  and "Rate of speech respiration(- /min)". 
