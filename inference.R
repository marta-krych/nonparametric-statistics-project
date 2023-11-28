library(readr)
library(ggplot2)
library(tidyr)

parkinson <- read_csv("data/parkinson_preprocessed.csv")
parkinson <- parkinson[,-1]

pd <- 1:30
rbd <- 31:80
hc <- 81:130

B <- 1000
alpha <- 0.05
set.seed(2023)

############ PD vs RBD

#### tremor, rigidity, spasm, stability
#we want to test if the 2 populations PD and RBD have the same distribution in each of these variables

pooled <- sapply(parkinson[1:80,11:14], as.numeric)

group <- as.factor(c(rep("pd",30), rep("rbd",50)))


for (i in 1:ncol(pooled))
  boxplot(pooled[,i] ~ group, main = colnames(pooled)[i])

n1 <- length(pd)
n2 <- length(rbd)
n <- n1 + n2

#Mann-Whitney U test: H0: X =^d Y vs H1: X !=^d Y
p.values <- numeric(4)
for (i in 1:4) {
  p.values[i] <- wilcox.test(x=pooled[pd,i], y=pooled[rbd,i], 
                             paired=F, conf.level = 1-alpha, conf.int = TRUE)$p.value
}
p.values
#low p-values in every variable. We reject the null hypothesis that the 2 populations are equally distributed

#Mann-Whitney U test: H0: X =^d Y+mu vs H1: X !=^d Y+mu
pd.median <- rbd.median <-numeric(4)
for (i in 1:4) {
  pd.median[i] <- median(pooled[pd,i])
  rbd.median[i] <- median(pooled[rbd,i])
}
mu <- pd.median - rbd.median

p.values <- numeric(4)
for (i in 1:4) {
  p.values[i] <- wilcox.test(x=pooled[pd,i], y=pooled[rbd,i], 
                             paired=F, conf.level = 1-alpha, mu=mu[i], conf.int=TRUE)$p.value
}
p.values
#we cannot reject the null hypothesis

# Permutation test

pd.mean <- colMeans(pooled[pd,])
rbd.mean <- colMeans(pooled[rbd,])
pd.cov  <-  cov(pooled[pd,])
rbd.cov  <-  cov(pooled[rbd,])

Sp      <- ((n1-1)*pd.cov + (n2-1)*rbd.cov)/(n1+n2-2)  # pooled cov matrix
Spinv   <- solve(Sp)

T2 <- as.numeric( t(pd.mean-rbd.mean) %*% Spinv %*% (pd.mean-rbd.mean))
T.stat <- numeric(B)

for(perm in 1:B)
{
  permutation <- sample(n) 
  pooled.perm <- pooled[permutation,]
  pd.mean.perm <- colMeans(pooled.perm[1:n1, ])
  rbd.mean.perm <- colMeans(pooled.perm[(n1+1):n,])
  T.stat[perm] <- as.numeric( t(pd.mean.perm-rbd.mean.perm) %*% Spinv %*% (pd.mean.perm-rbd.mean.perm))
}

hist(T.stat,xlim=range(c(T2,T.stat)))
abline(v=T2,col=3,lwd=4)

#reject the null hypothesis, there is statistical evidence of the difference between PD and RBD in 
#tremor, rigidity, spasm and stability looked together.

# Permutation test for X, Y+mu

pd.mean <- colMeans(pooled[pd,])
rbd.mean <- colMeans(pooled[rbd,])+mu
pd.cov  <-  cov(pooled[pd,])
rbd.cov  <-  cov(pooled[rbd,])

Sp      <- ((n1-1)*pd.cov + (n2-1)*rbd.cov)/(n1+n2-2)  # pooled cov matrix
Spinv   <- solve(Sp)

T2 <- as.numeric( t(pd.mean-rbd.mean) %*% Spinv %*% (pd.mean-rbd.mean))
T.stat <- numeric(B)

for(perm in 1:B)
{
  permutation <- sample(n) 
  pooled.perm <- pooled[permutation,]
  pd.mean.perm <- colMeans(pooled.perm[1:n1, ])
  rbd.mean.perm <- colMeans(pooled.perm[(n1+1):n,])+mu
  T.stat[perm] <- as.numeric( t(pd.mean.perm-rbd.mean.perm) %*% Spinv %*% (pd.mean.perm-rbd.mean.perm))
}

hist(T.stat,xlim=range(c(T2,T.stat)))
abline(v=T2,col=3,lwd=4)

plot(ecdf(T.stat), xlim=range(c(T2,T.stat)))
abline(v=T2,col=3,lwd=4)

(p.value <- sum(T.stat>=T2)/B)

#we can confirm that we cannot reject the hypothesis that pd =^d rbd + mu 

####################### Speech data

speech.data <- parkinson[,15:28]
speech.data$Group <- c(rep("pd",30), rep("rbd",50), rep("hc",50))

# Reshape the data to long format for ggplot
data_long <- gather(speech.data, key = "Variable", value = "Value", -Group)

# Plot histograms using ggplot2 with facets
ggplot(data_long, aes(x = Value, fill = Group)) +
  geom_histogram(position = "identity", bins = 60) +
  facet_wrap(~ Variable, scales = "free") +
  labs(title = "Histograms for speech variables by Group")


speech.data.pd <- data.frame(speech.data[1:30,-15])
speech.data.rbd <- data.frame(speech.data[31:80,-15])
speech.data.hc <- data.frame(speech.data[81:130,-15])


############## PD vs NOPD
speech.data.nopd <- data.frame(speech.data[-(1:30),-15])

pd.median <- nopd.median <- numeric(14)
for (i in 1:14) {
  pd.median[i] <- median(speech.data.pd[,i])
  nopd.median[i] <- median(speech.data.nopd[,i])
}
mu <- pd.median - nopd.median

pd.mean <- colMeans(speech.data.pd)
nopd.mean.shifted <- colMeans(speech.data.nopd)+mu

T2 <- as.numeric( t(pd.mean - nopd.mean.shifted) %*% (pd.mean - nopd.mean.shifted))
T.stat <- numeric(B)
n1 <- dim(speech.data.pd)[1]
n2 <- dim(speech.data.nopd)[1]
n <- n1+n2

for(perm in 1:B) {
  permutation <- sample(n) 
  speech.data.perm <- speech.data[permutation,-15]
  pd.mean.perm <- colMeans(speech.data.perm[1:n1, ])
  nopd.mean.perm <- colMeans(speech.data.perm[(n1+1):n,])+mu
  T.stat[perm] <- as.numeric( t(pd.mean.perm-nopd.mean.perm) %*% (pd.mean.perm-nopd.mean.perm))
}

hist(T.stat,xlim=range(c(T2,T.stat)))
abline(v=T2,col=3,lwd=4)

(p.value <- sum(T.stat>=T2)/B)

#we cannot reject the null hypothesis. We can then assume that mu is the shifting difference
#between pd and non-pd patients.

################# RBD VS HC

rbd.median <- hc.median <- numeric(14)
for (i in 1:14) {
  rbd.median[i] <- median(speech.data.rbd[,i])
  hc.median[i] <- median(speech.data.hc[,i])
}
mu <- rbd.median - hc.median

rbd.mean <- colMeans(speech.data.rbd)
hc.mean.shifted <- colMeans(speech.data.hc)+mu

T2 <- as.numeric( t(rbd.mean - hc.mean.shifted) %*% (rbd.mean - hc.mean.shifted))
T.stat <- numeric(B)
n1 <- dim(speech.data.rbd)[1]
n2 <- dim(speech.data.hc)[1]
n <- n1+n2

for(perm in 1:B) {
  permutation <- sample(n) 
  speech.data.perm <- speech.data[permutation,-15]
  rbd.mean.perm <- colMeans(speech.data.perm[1:n1, ])
  hc.mean.perm <- colMeans(speech.data.perm[(n1+1):n,])+mu
  T.stat[perm] <- as.numeric( t(rbd.mean.perm-hc.mean.perm) %*% (rbd.mean.perm-hc.mean.perm))
}

hist(T.stat,xlim=range(c(T2,T.stat)))
abline(v=T2,col=3,lwd=4)

(p.value <- sum(T.stat>=T2)/B)

#we cannot reject the null hypothesis. We can then assume that mu is the shifting difference
#between rbd and pc patients




######################## paired speech data
original <- read_csv("data/parkinson.csv")
colnames(original) <- original[1,]
original<-original[-1,-1]

task1 <- sapply(original[,41:52],as.numeric)
task2 <- sapply(original[,53:64],as.numeric)

diff <- data.frame(task1-task2)
colnames(diff) <- colnames(task1)
diff <- as.data.frame(sapply(diff, as.numeric))

diff$Group <- c(rep("pd",30), rep("rbd",50), rep("hc",50))

# Reshape the data to long format for ggplot
data_long <- gather(diff, key = "Variable", value = "Value", -Group)

# Plot histograms using ggplot2 with facets
ggplot(data_long, aes(x = Value, fill = Group)) +
  geom_histogram(position = "identity", bins = 60) +
  facet_wrap(~ Variable, scales = "free") +
  labs(title = "Histograms for speech task differences by Group")
diff$Group <-NULL

t1.mean <- colMeans(task1)
t2.mean <- colMeans(task2)
t1.cov  <-  cov(task1)
t2.cov  <-  cov(task2)
Sp      <- (t1.cov + t2.cov)/2  # pooled cov matrix
Spinv   <- solve(Sp)

delta.0 <- rep(0,12)

diff.mean <- colMeans(diff)
diff.cov <- cov(diff)
diff.invcov <- solve(diff.cov)

T20 <- as.numeric(t(diff.mean-delta.0)  %*% diff.invcov %*% (diff.mean-delta.0))

T2 <- numeric(B)
n <- dim(diff)[1]
p <- dim(diff)[2]

for(perm in 1:B) {
 
  signs.perm <- rbinom(n, 1, 0.5)*2 - 1
  
  diff_perm <- diff * matrix(signs.perm,nrow=n,ncol=p,byrow=FALSE)
  diff.mean_perm <- colMeans(diff_perm)
  diff.cov_perm <- cov(diff_perm)
  diff.invcov_perm <- solve(diff.cov_perm)
  
  T2[perm] <- as.numeric(t(diff.mean_perm-delta.0) %*%
                           diff.invcov_perm %*% (diff.mean_perm-delta.0))
}

hist(T2,xlim=range(c(T2,T20)),breaks=100)
abline(v=T20,col=3,lwd=4)

(p_val <- sum(T2>=T20)/B)








