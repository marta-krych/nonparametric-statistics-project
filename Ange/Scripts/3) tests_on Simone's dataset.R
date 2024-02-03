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

colnames(dataset) <- pdataset[1,]
dataset <- dataset[-1,]
rownames(dataset) <- 1:dim(dataset)[1]

dataset <- dataset[,-1]

dataset <- type.convert(dataset, as.is=T)

dataset <- dataset[,-c(7,8,10)]
## Let's start to use the depth measure:
parkinson_ill <- dataset[1:30,]
RBD <-  dataset[31:80,]
health <- dataset[81:130,]

## I can try to do some tests:
# • Two multivariate populations test:
#  - Parinson vs RBD:
t1 =  parkinson_ill[,15:28]
t2 =  RBD[,15:28]

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
p_val # 0.25614. Absolutely not satisfing

# - parkinson vs health:
t1 =  parkinson_ill[,15:28]
t2 =  health[,15:28]

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
p_val #  0.00789 here is ok.

# - RBD vs health:
t1 =  RBD[,15:28]
t2 =  health[,15:28]

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
p_val # 0.00232


# • Permutational one way Manova:
reduced_dataset <- dataset[1:80,]
reduced_dataset$Category <- as.factor(c( rep("Ill", 30) , rep("At Risk", 50) ))
n1 <- 30
n2 <- 50
n_test  <- n1+n2

g  <- 2
p  <- 14

fit <- manova(as.matrix(reduced_dataset[,15:28]) ~ reduced_dataset$Category)
summary.manova(fit,test="Wilks") 

T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
T0

set.seed(seed)
T_stat <- numeric(B)

for(perm in 1:B){
  # choose random permutation
  permutation <- sample(1:n_test)
  category.perm <- reduced_dataset$Category[permutation]
  fit.perm <- manova(as.matrix(reduced_dataset[,15:28]) ~ category.perm)
  T_stat[perm] <- -summary.manova(fit.perm,test="Wilks")$stats[1,2]
}

hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

plot(ecdf(T_stat),xlim=c(-2,1))
abline(v=T0,col=3,lwd=4)

p_val <- sum(T_stat>=T0)/B
p_val # 0.54311. Try another strategy.

## Hierarchical clustering:
dataset.e <- dist(RBD[,15:28], method='euclidean')
dataset.m <- dist(RBD[,15:28], method='manhattan')
dataset.c <- dist(RBD[,15:28], method='canberra')

x11()
image(1:50,1:50,as.matrix(dataset.e), main='metrics: Euclidean', asp=1, xlab='i', ylab='j' )

x11()
image(1:50,1:50,as.matrix(dataset.c), main='metrics: Canberra', asp=1, xlab='i', ylab='j' )

x11()
image(1:50,1:50,as.matrix(dataset.m), main='metrics: Manhattan', asp=1, xlab='i', ylab='j' )

dataset.es <- hclust(dataset.e, method='single')
dataset.ea <- hclust(dataset.e, method='average')
dataset.ec <- hclust(dataset.e, method='complete')
dataset.ew<- hclust(dataset.e, method='ward.D2')

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

cluster.ec <- cutree(dataset.ec, k=2)


# confirm the the clustering is good throught a permutational MANOVA:
n1 <- 44
n2 <- 6
n_test  <- n1+n2

g  <- 2
p  <- 14

fit <- manova(as.matrix(RBD[,15:28]) ~ cluster.ec)
summary.manova(fit,test="Wilks") 

T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
T0

set.seed(seed)
T_stat <- numeric(B)

for(perm in 1:B){
  # choose random permutation
  permutation <- sample(1:n_test)
  category.perm <- cluster.ec[permutation]
  fit.perm <- manova(as.matrix(RBD[,15:28]) ~ category.perm)
  T_stat[perm] <- -summary.manova(fit.perm,test="Wilks")$stats[1,2]
}

hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

plot(ecdf(T_stat),xlim=c(-2,1))
abline(v=T0,col=3,lwd=4)

p_val <- sum(T_stat>=T0)/B
p_val # 0


# Now which of the two cluster is more similar to PD?
# • Start with cluster 1:
#   - try to build the intersted dataset:
my_dataset <- rbind(parkinson_ill[,15:28],RBD[which(cluster.ec==1), 15:28])
PD_RBD_category <- as.factor(c( rep("0", 30) , rep("1", 44) ))

#   - Manova test:    
n1 <- 30
n2 <- 44
n_test  <- n1+n2

g  <- 2
p  <- 14

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
p_val # 0.46648


# • Start with cluster 1:
#   - try to build the intersted dataset:
my_dataset <- rbind(parkinson_ill[,15:28],RBD[which(cluster.ec==2), 15:28])
PD_RBD_category <- as.factor(c( rep("0", 30) , rep("2", 6) ))

#   - Manova test:    
n1 <- 30
n2 <- 6
n_test  <- n1+n2

g  <- 2
p  <- 14

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
p_val # 0.16331


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
p.value <- numeric(14)
#pb=progress_bar$new(total=B)
#pb$tick(0)
#set.seed(seed)
cont <- 1
for (i in 14:27){
  x.1 <- parkinson_ill[,i]
  x.2 = health[,i]
  #       x.2 = RBD[which(cluster.ew==2), i]
  p.value[cont] <- perm_t_test(x.1,x.2,iter=1e3)
  cont <- cont + 1
}
p.value

# p-values:  0.800 0.003 0.889 0.691 0.217 0.419 0.685 0.296 0.674 0.523 0.004 0.503 0.275 0.112.
# The variables in which there are the differences are Rate of speech timing and Rate of speech timing.1

