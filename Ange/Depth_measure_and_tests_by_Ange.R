## Useful libraries:
library(readr)
library(MASS)
library(rgl)
library(DepthProc)
library(hexbin)
library(aplpack)
library(robustbase)
library(MDBED)

## Set B and seed:
B = 1e5
seed = 26111992

## Import the dataset:
dataset <- read.csv("dataset_preprocessed.csv")
View(dataset)

dataset <- dataset[,-1]

# to fix and order if is possible.
# what about the speech categorical variable?

## Let's start to use the depth measure:
parkinson_ill <- dataset[1:30,]
RBD <-  dataset[31:80,]
health <- dataset[81:130,]

#  • i start look on the numerical variables:
#    - Parkinson vs RBD:
     ddPlot(x = parkinson_ill[,14:25],y = RBD[,14:25],depth_params = list(method='Tukey'))
     
#    - Parkinson vs health:
     ddPlot(x = parkinson_ill[,14:25],y = health[,14:25],depth_params = list(method='Tukey'))     
     
#    - RBD vs health:
     ddPlot(x = RBD[,14:25],y = health[,14:25],depth_params = list(method='Tukey'))
     
#  • what about the categorical ones? :
     #    - Parkinson vs RBD:
     ddPlot(x = sapply(parkinson_ill[,7:9], as.numeric),y = sapply(RBD[,7:9], as.numeric) ,depth_params = list(method='Tukey'))
     
# NB: for the categorical variables focused only to charachter, speech and facial expression because 
#     Total.Stability ,Total.Spasm, Total.Rigidity, Total.Tremor  are constructed after permutational Anova tests
#      that highlight some differences between the group of the patients!
     

## I can try to do some tests:
# • Two multivariate populations test:
#  - Parinson vs RBD:
     t1 =  parkinson_ill[,14:25]
     t2 =  RBD[,14:25]
     
     t1.mean = colMeans(t1)
     t2.mean = colMeans(t2)
     
     n1 = dim(t1)[1]
     n2 = dim(t2)[1]
     n  = n1 + n2
     
     T20 = as.numeric(t(t1.mean-t2.mean) %*% (t1.mean-t2.mean))  # matrix (vector) product
     T20     
     
     T2 = numeric(B)
     set.seed(seed)
     for(perm in 1:B){
       # Random permutation of indexes
       # When we apply permutations in a multivariate case, we keep the units together
       # i.e., we only permute the rows of the data matrix
       t_pooled = rbind(t1,t2)
       permutation = sample(n)
       t_perm = t_pooled[permutation,]
       t1_perm = t_perm[1:n1,]
       t2_perm = t_perm[(n1+1):n,]
       
       # Evaluation of the test statistic on permuted data
       t1.mean_perm = colMeans(t1_perm)
       t2.mean_perm = colMeans(t2_perm)
       T2[perm]  = t(t1.mean_perm-t2.mean_perm) %*% (t1.mean_perm-t2.mean_perm) 
     }
     
     hist(T2,xlim=range(c(T2,T20)))
     abline(v=T20,col=3,lwd=4)
     
     plot(ecdf(T2))
     abline(v=T20,col=3,lwd=4)
     
     p_val = sum(T2>=T20)/B
     p_val # 0.33684. Absolutely not satisfing
     
# - parkinson vs health:
     t1 =  parkinson_ill[,14:25]
     t2 =  health[,14:25]
     
     t1.mean = colMeans(t1)
     t2.mean = colMeans(t2)
     
     n1 = dim(t1)[1]
     n2 = dim(t2)[1]
     n  = n1 + n2
     
     T20 = as.numeric(t(t1.mean-t2.mean) %*% (t1.mean-t2.mean))  # matrix (vector) product
     T20     
     
     T2 = numeric(B)
     set.seed(seed)
     for(perm in 1:B){
       # Random permutation of indexes
       # When we apply permutations in a multivariate case, we keep the units together
       # i.e., we only permute the rows of the data matrix
       t_pooled = rbind(t1,t2)
       permutation = sample(n)
       t_perm = t_pooled[permutation,]
       t1_perm = t_perm[1:n1,]
       t2_perm = t_perm[(n1+1):n,]
       
       # Evaluation of the test statistic on permuted data
       t1.mean_perm = colMeans(t1_perm)
       t2.mean_perm = colMeans(t2_perm)
       T2[perm]  = t(t1.mean_perm-t2.mean_perm) %*% (t1.mean_perm-t2.mean_perm) 
     }
     
     hist(T2,xlim=range(c(T2,T20)))
     abline(v=T20,col=3,lwd=4)
     
     plot(ecdf(T2))
     abline(v=T20,col=3,lwd=4)
     
     p_val = sum(T2>=T20)/B
     p_val #  0.00023, here is ok.
     
  # - RBD vs health:
     t1 =  RBD[,14:25]
     t2 =  health[,14:25]
     
     t1.mean = colMeans(t1)
     t2.mean = colMeans(t2)
     
     n1 = dim(t1)[1]
     n2 = dim(t2)[1]
     n  = n1 + n2
     
     T20 = as.numeric(t(t1.mean-t2.mean) %*% (t1.mean-t2.mean))  # matrix (vector) product
     T20     
     
     T2 = numeric(B)
     set.seed(seed)
     for(perm in 1:B){
       # Random permutation of indexes
       # When we apply permutations in a multivariate case, we keep the units together
       # i.e., we only permute the rows of the data matrix
       t_pooled = rbind(t1,t2)
       permutation = sample(n)
       t_perm = t_pooled[permutation,]
       t1_perm = t_perm[1:n1,]
       t2_perm = t_perm[(n1+1):n,]
       
       # Evaluation of the test statistic on permuted data
       t1.mean_perm = colMeans(t1_perm)
       t2.mean_perm = colMeans(t2_perm)
       T2[perm]  = t(t1.mean_perm-t2.mean_perm) %*% (t1.mean_perm-t2.mean_perm) 
     }
     
     hist(T2,xlim=range(c(T2,T20)))
     abline(v=T20,col=3,lwd=4)
     
     plot(ecdf(T2))
     abline(v=T20,col=3,lwd=4)
     
     p_val = sum(T2>=T20)/B
     p_val # 0.00047
     
# This tests show that the numerical variables pattern show important differences between 
# the health patients with RBD and Parkinson. But it seem the this test is not enough powerful to show important differences 
# between RBD and Parkison patients.
# Honestly i want to go deeper to see some differences between RBD and Parkinson patiens.
     
# • Permutational one way Manova:
     reduced_dataset <- dataset[1:80,]
     reduced_dataset$Category <- as.factor(c( rep("Ill", 30) , rep("At Risk", 50) ))
     n1 <- 30
     n2 <- 50
     n_test  <- n1+n2
     
     g  <- 2
     p  <- 11
     
     fit <- manova(as.matrix(reduced_dataset[,14:25]) ~ reduced_dataset$Category)
     summary.manova(fit,test="Wilks") 
     
     T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
     T0
     
     set.seed(seed)
     T_stat <- numeric(B)
     
     for(perm in 1:B){
       # choose random permutation
       permutation <- sample(1:n_test)
       category.perm <- reduced_dataset$Category[permutation]
       fit.perm <- manova(as.matrix(reduced_dataset[,14:25]) ~ category.perm)
       T_stat[perm] <- -summary.manova(fit.perm,test="Wilks")$stats[1,2]
     }
     
     hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
     abline(v=T0,col=3,lwd=2)
     
     plot(ecdf(T_stat),xlim=c(-2,1))
     abline(v=T0,col=3,lwd=4)
     
     p_val <- sum(T_stat>=T0)/B
     p_val # 0.3994.
     
# In conclusion Test of multivariate data and one-way Manova bring the same result: non difference between RBD an Parkinson.
     