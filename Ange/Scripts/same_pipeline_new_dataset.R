## Useful libraries:
library(readr)
library(MASS)
library(rgl)
library(DepthProc)
library(hexbin)
library(aplpack)
library(robustbase)
library(MDBED)
library(mvtnorm)
library(MVN)
library(rgl)
library(car)
library(dbscan)
library(cluster)
library(fields)

## Set B and seed:
B = 1e5
seed = 26111992

## Import the dataset:
dataset <- read.csv("speech_no_outliers.csv")
View(dataset)

# dataset <- dataset[,-1]
# dataset <- dataset[-1,]


## Let's start to use the depth measure:
parkinson_ill <- dataset[1:23,]
health <-  dataset[24:66,]
RBD<- dataset[67:116,]

#  • i start look on the numerical variables:
#    - Parkinson vs RBD:

ddPlot(x = parkinson_ill[,c(13,13+ 12)],y = RBD[,c(13,13+12)],depth_params = list(method='Tukey'))


#    - Parkinson vs health:
ddPlot(x = parkinson_ill[,14:25],y = health[,14:25],depth_params = list(method='Tukey'))     

#    - RBD vs health:
ddPlot(x = RBD[,14:25],y = health[,14:25],depth_params = list(method='Tukey'))



# NB: for the categorical variables focused only to charachter, speech and facial expression because 
#     Total.Stability ,Total.Spasm, Total.Rigidity, Total.Tremor  are constructed after permutational Anova tests
#      that highlight some differences between the group of the patients!


## I can try to do some tests:
# • Two multivariate populations test:
#  - def of the test (to run in parallel):
multi_perm_t_test=function(x,y,iter=B){
t1 <- x
t2 <- y

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
 return(p_val)
}

#  - Parinson vs RBD:
 multi_perm_t_test(parkinson_ill[,c(2:25)],RBD[,c(2:25)]) # p-value: 0.6614

# - parkinson vs health:
 multi_perm_t_test(parkinson_ill[,c(2:25)], health[,c(2:25)]) # p-value: 0.00817

# - RBD vs health:
 multi_perm_t_test(RBD[,c(2:25)], health[,c(2:25)]) # p-value: 0.00114


# This tests show that the numerical variables pattern show important differences between 
# the health patients with RBD and Parkinson. But it seem the this test is not enough powerful to show important differences 
# between RBD and Parkison patients.
# Honestly i want to go deeper to see some differences between RBD and Parkinson patiens.

# • Permutational one way Manova:
reduced_dataset <- dataset[c(1:23,67:116),]
category <- factor(dataset[c(1:23,67:116),26])
n1 <- 23
n2 <- 50
n_test  <- n1+n2

g  <- 2
p  <- 24

fit <- manova(as.matrix(reduced_dataset[,2:25]) ~ category)
summary.manova(fit,test="Wilks") 

T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
T0

set.seed(seed)
T_stat <- numeric(B)

for(perm in 1:B){
  # choose random permutation
  permutation <- sample(1:n_test)
  category.perm <- category[permutation]
  fit.perm <- manova(as.matrix(reduced_dataset[,2:25]) ~ category.perm)
  T_stat[perm] <- -summary.manova(fit.perm,test="Wilks")$stats[1,2]
}

hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

plot(ecdf(T_stat),xlim=c(-2,1))
abline(v=T0,col=3,lwd=4)

p_val <- sum(T_stat>=T0)/B
p_val # 0.74978

# In conclusion Test of multivariate data and one-way Manova bring the same result: non difference between RBD an Parkinson.

## Hierarchical clustering (Only for RBD patients):
dataset.e <- dist(RBD[,2:25], method='euclidean')
dataset.m <- dist(RBD[,2:25], method='manhattan')
dataset.c <- dist(RBD[,2:25], method='canberra')

x11()
image(1:50,1:50,as.matrix(dataset.e), main='metrics: Euclidean', asp=1, xlab='i', ylab='j' )

x11()
image(1:50,1:50,as.matrix(dataset.c), main='metrics: Canberra', asp=1, xlab='i', ylab='j' )

x11()
image(1:50,1:50,as.matrix(dataset.m), main='metrics: Manhattan', asp=1, xlab='i', ylab='j' )

# proceed with the Euclidean distance.

dataset.es <- hclust(dataset.e, method='single')
dataset.ea <- hclust(dataset.e, method='average')
dataset.ec <- hclust(dataset.e, method='complete')
dataset.ew <- hclust(dataset.e, method='ward.D2')

x11()
plot(dataset.es, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(dataset.es, k=2)

x11()
plot(dataset.ec, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(dataset.ec, k=2)

x11()
plot(dataset.ea, main='euclidean-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(dataset.ea, k=2)

x11()
plot(dataset.ew, main='euclidean-ward', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(dataset.ew, k=2)


cluster.ew <- cutree(dataset.ew, k=2)

# number on each cluster:
length(which(cluster.ew==1))

# confirm the the clustering is good throught a permutational MANOVA:
n1 <- 31
n2 <- 19
n_test  <- n1+n2

g  <- 2
p  <- 24

fit <- manova(as.matrix(RBD[,2:25]) ~ cluster.ew)
summary.manova(fit,test="Wilks") 

T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
T0

set.seed(seed)
T_stat <- numeric(B)

for(perm in 1:B){
  # choose random permutation
  permutation <- sample(1:n_test)
  category.perm <- cluster.ew[permutation]
  fit.perm <- manova(as.matrix(RBD[,14:25]) ~ category.perm)
  T_stat[perm] <- -summary.manova(fit.perm,test="Wilks")$stats[1,2]
}

hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

plot(ecdf(T_stat),xlim=c(-2,1))
abline(v=T0,col=3,lwd=4)

p_val <- sum(T_stat>=T0)/B
p_val # 0, the clustering is good.

# Now which of the two cluster is more similar to PD?
# • Start with cluster 1:
#   - try to build the intersted dataset:
my_dataset <- rbind(parkinson_ill[,2:25],RBD[which(cluster.ew==1), 2:25])
PD_RBD_category <- as.factor(c( rep("0", 23) , rep("1", 31) ))

#   - Manova test:    
n1 <- 23
n2 <- 31
n_test  <- n1+n2

g  <- 2
p  <- 11

fit <- manova(as.matrix(my_dataset) ~ PD_RBD_category)
summary.manova(fit,test="Wilks") 

T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
T0

set.seed(seed)
T_stat <- numeric(B)

for(perm in 1:B){
  # choose random permutation
  permutation <- sample(1:n_test)
  category.perm <- PD_RBD_category[permutation]
  fit.perm <- manova(as.matrix(my_dataset) ~ category.perm)
  T_stat[perm] <- -summary.manova(fit.perm,test="Wilks")$stats[1,2]
}

hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

plot(ecdf(T_stat),xlim=c(-2,1))
abline(v=T0,col=3,lwd=4)

p_val <- sum(T_stat>=T0)/B
p_val # 0.69374

# • Continue with cluster 2:
#   - try to build the intersted dataset:
my_dataset <- rbind(parkinson_ill[,2:25],RBD[which(cluster.ew==2), 2:25])
PD_RBD_category <- as.factor(c( rep("0", 23) , rep("2", 19) ))

#   - Manova test:    
n1 <- 23
n2 <- 19
n_test  <- n1+n2

g  <- 2
p  <- 24

fit <- manova(as.matrix(my_dataset) ~ PD_RBD_category)
summary.manova(fit,test="Wilks") 

T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
T0

set.seed(seed)
T_stat <- numeric(B)

for(perm in 1:B){
  # choose random permutation
  permutation <- sample(1:n_test)
  category.perm <- PD_RBD_category[permutation]
  fit.perm <- manova(as.matrix(my_dataset) ~ category.perm)
  T_stat[perm] <- -summary.manova(fit.perm,test="Wilks")$stats[1,2]
}

hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

plot(ecdf(T_stat),xlim=c(-2,1))
abline(v=T0,col=3,lwd=4)

p_val <- sum(T_stat>=T0)/B
p_val # 0.09972 at lthe limit

# - try to see the means to detect the possible differences:
colMeans(parkinson_ill[,2:25]) - colMeans(RBD[which(cluster.ew==2), 2:25]) 

# - patiently i try to see which variable highlights the difference:
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
for (i in 2:25){
  x.1 <- parkinson_ill[,i]
  x.2 = RBD[which(cluster.ew==2),i]
  p.value[cont] <- perm_t_test(x.1,x.2,iter=1e3)
  cont <- cont + 1
}
p.value

# RBD(cluster 2) vs PD:
# p-value: 0.401 0.511 0.349 0.337 0.300 0.293 0.593 0.908 0.179 0.156 0.235 0.891 0.112 0.000 0.046 0.027 0.043 0.018 0.011 0.582 0.357
# 0.006 0.214 0.028
# Differences in RST_m, AST_m, DPI_m, DVI_m, GVI_m, DUS_m, PIR_m and LRE_m


# Health vs PD:
# p-values: 0.115 0.043 0.827 0.009 0.314 0.630 0.464 0.826 0.444 0.647 0.669 0.960 0.865 0.031 0.323
# 0.001 0.716 0.956 0.020 0.785 0.684 0.056 0.092 0.111.
# Differences in: RST_r, DPI_r, RST_m, DPI_m, DUS_m, PIR_m and RSR_m