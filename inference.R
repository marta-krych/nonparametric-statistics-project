library(readr)
library(ggplot2)
library(tidyr)
library(DepthProc)

parkinson <- read_csv("data/parkinson_preprocessed.csv")

B <- 1000
alpha <- 0.05
seed <- 2024

############ Exploratory analysis on categorical variables

#we want to test if the 2 populations PD and RBD have the same distribution in each 
#of these variables

categoric <- sapply(parkinson[1:80,9:35], as.numeric)

category <- as.factor(parkinson$Category[1:80])

for (i in 1:ncol(categoric))
  boxplot(categoric[,i] ~ category, main = colnames(categoric)[i])

pd <- which(category=="PD")
rbd <- which(category=="RBD")

n1 <- length(pd)
n2 <- length(rbd)
n <- n1 + n2

#Mann-Whitney U test: H0: X =^d Y vs H1: X !=^d Y
p.values <- numeric(dim(categoric)[2])
names(p.values) <- colnames(categoric)

for (i in 1:27) {
  p.values[i] <- wilcox.test(x=categoric[pd,i], y=categoric[rbd,i], 
                             paired=F, conf.level = 1-alpha, conf.int = TRUE)$p.value
}
names(which(p.values >= alpha))
#low p-values in every variable except "Tremor at Rest - LLE" and "Arising from Chair"
#We can state that PD and RBD are equally distributed only in this variables


##### tremor 
#We study if the overall effect (their sum) of the tremor variables underlines
#a difference between PD and RBD, via permutational ANOVA

tremor <- categoric[,3:9]
summary(tremor)

sum.tremor <- rowSums(tremor)

boxplot(sum.tremor ~ category)

fit <- aov(sum.tremor ~ category)
T0 <- summary(fit)[[1]][1,4]

set.seed(seed)

T_stat <- numeric(B) 
n <- length(sum.tremor)

for(perm in 1:B){
  # Permutation:
  permutation <- sample(1:n)
  sum.tremor_perm <- sum.tremor[permutation]
  fit_perm <- aov(sum.tremor_perm ~ category)
  
  # Test statistic:
  T_stat[perm] <- summary(fit_perm)[[1]][1,4]
}

hist(T_stat,xlim=range(c(T_stat,T0)))
abline(v=T0,col=3,lwd=2)

p.val <- sum(T_stat>=T0)/B
p.val  #the sum of the tremor variables is significant, as expected


##### rigidity 
#We study if the overall effect (their sum) of the rigidity variables underlines
#a difference between PD and RBD, via permutational ANOVA

rigidity <- categoric[,10:14]

sum.rigidity <- rowSums(rigidity)

boxplot(sum.rigidity ~ category)

fit <- aov(sum.rigidity ~ category)
T0 <- summary(fit)[[1]][1,4]

set.seed(seed)

T_stat <- numeric(B) 
n <- length(sum.rigidity)

for(perm in 1:B){
  # Permutation:
  permutation <- sample(1:n)
  sum.rigidity_perm <- sum.rigidity[permutation]
  fit_perm <- aov(sum.rigidity_perm ~ category)
  
  # Test statistic:
  T_stat[perm] <- summary(fit_perm)[[1]][1,4]
}

hist(T_stat,xlim=range(c(T_stat,T0)))
abline(v=T0,col=3,lwd=2)

p.val <- sum(T_stat>=T0)/B
p.val  #the sum of the rigidity variables is significant, as expected 


##### spasm 
#We study if the overall effect (their sum) of the spasm variables underlines
#a difference between PD and RBD, via permutational ANOVA

spasm <- categoric[,15:20]

sum.spasm <- rowSums(spasm)

boxplot(sum.spasm ~ category)

fit <- aov(sum.spasm ~ category)
T0 <- summary(fit)[[1]][1,4]

set.seed(seed)

T_stat <- numeric(B) 
n <- length(sum.spasm)

for(perm in 1:B){
  # Permutation:
  permutation <- sample(1:n)
  sum.spasm_perm <- sum.spasm[permutation]
  fit_perm <- aov(sum.spasm_perm ~ category)
  
  # Test statistic:
  T_stat[perm] <- summary(fit_perm)[[1]][1,4]
}

hist(T_stat,xlim=range(c(T_stat,T0)))
abline(v=T0,col=3,lwd=2)

p.val <- sum(T_stat>=T0)/B
p.val  #the sum of the spasm variables is significant, as expected 


##### stability 
#We study if the overall effect (their sum) of the stability variables underlines
#a difference between PD and RBD, via permutational ANOVA

stability <- categoric[,21:27]

sum.stability <- rowSums(stability)

boxplot(sum.stability ~ category)

fit <- aov(sum.stability ~ category)
T0 <- summary(fit)[[1]][1,4]

set.seed(seed)

T_stat <- numeric(B) 
n <- length(sum.stability)

for(perm in 1:B){
  # Permutation:
  permutation <- sample(1:n)
  sum.stability_perm <- sum.stability[permutation]
  fit_perm <- aov(sum.stability_perm ~ category)
  
  # Test statistic:
  T_stat[perm] <- summary(fit_perm)[[1]][1,4]
}

hist(T_stat,xlim=range(c(T_stat,T0)))
abline(v=T0,col=3,lwd=2)

p.val <- sum(T_stat>=T0)/B
p.val  #the sum of the stability variables is significant, as expected 

## PD and RBD patients are clearly distinguishable from these variables




####################### Speech data

data <- read.csv("data/speech.csv")
speech.data <- sapply(data[,-c(25:27)], as.numeric)
speech.data <- as.data.frame(speech.data)

age <- data$Age
gender <- data$Gender
Category <- factor(data$Category)

# Reshape the data to long format for ggplot
data_long <- gather(cbind(speech.data,Category), 
                    key = "Variable", value = "Value", -Category)

# Plot histograms using ggplot2 with facets
ggplot(data_long, aes(x = as.numeric(Value), fill = Category)) +
  geom_histogram(position = "identity") +
  facet_wrap(~ Variable, scales = "free") +
  labs(title = "Histograms for speech variables by Group")


speech.data.pd <- speech.data[Category == "PD",]
speech.data.rbd <- speech.data[Category == "RBD",]
speech.data.hc <- speech.data[Category == "Control",]


#####  PD vs (HC + RBD)
#I apply a global 2 populations permutation test to assess if PD patients 
#are equally distributed to the other categories

speech.data.nopd <- speech.data[Category != "PD",]

med.pd <- depthMedian(speech.data.pd, depth_params = list(method='Tukey'))
med.nopd <- depthMedian(speech.data.nopd, depth_params = list(method='Tukey'))

n1 <- dim(speech.data.pd)[1]
n2 <- dim(speech.data.nopd)[1]
n  <- n1 + n2

T0 <- max(abs(med.pd-med.nopd))  

T2 <- numeric(B)

set.seed(seed)

for(perm in 1:B){
  
  permutation <- sample(n)
  df.perm <- speech.data[permutation,]
  df1.perm <- df.perm[1:n1,]
  df2.perm <- df.perm[(n1+1):n,]
  
  x.med.perm <- depthMedian(df1.perm, depth_params = list(method='Tukey'))
  y.med.perm <- depthMedian(df2.perm, depth_params = list(method='Tukey'))
  
  T2[perm] <- max(abs(x.med.perm-y.med.perm))  
  
}

plot(ecdf(T2), xlim=range(c(T2,T0)))
abline(v=T0,col=3,lwd=4)

p.val <- sum(T2>=T0)/B
p.val

#0.473, so we cannot reject the null hypothesis 
#(maybe outliers or the presence of HC patients influence the result)


################# HC vs (PD + RBD)
#I apply a global 2 populations permutation test to assess if HC patients 
#are equally distributed to the other categories
speech.data.nohc <- speech.data[Category != "Control",]

med.hc <- depthMedian(speech.data.hc, depth_params = list(method='Tukey'))
med.nohc <- depthMedian(speech.data.nohc, depth_params = list(method='Tukey'))

n1 <- dim(speech.data.hc)[1]
n2 <- dim(speech.data.nohc)[1]
n  <- n1 + n2

T0 <- max(abs(med.hc-med.nohc))  

T2 <- numeric(B)

set.seed(seed)

for(perm in 1:B){
  
  permutation <- sample(n)
  df.perm <- speech.data[permutation,]
  df1.perm <- df.perm[1:n1,]
  df2.perm <- df.perm[(n1+1):n,]
  
  x.med.perm <- depthMedian(df1.perm, depth_params = list(method='Tukey'))
  y.med.perm <- depthMedian(df2.perm, depth_params = list(method='Tukey'))
  
  T2[perm] <- max(abs(x.med.perm-y.med.perm))  
  
}

plot(ecdf(T2), xlim=range(c(T2,T0)))
abline(v=T0,col=3,lwd=4)

p.val <- sum(T2>=T0)/B
p.val

#0.035: we can reject the null hypothesis. We can then state that Control patients
#are not distributed as the others


#### PD vs RBD
#I apply a global 2 populations permutation test to assess if PD patients 
#are equally distributed to RBD patients

med.pd <- depthMedian(speech.data.pd, depth_params = list(method='Tukey'))
med.rbd <- depthMedian(speech.data.rbd, depth_params = list(method='Tukey'))

n1 <- dim(speech.data.pd)[1]
n2 <- dim(speech.data.rbd)[1]
n  <- n1 + n2

T0 <- max(abs(med.pd-med.rbd))  

T2 <- numeric(B)

set.seed(seed)

for(perm in 1:B){
  
  permutation <- sample(n)
  df.perm <- speech.data[permutation,]
  df1.perm <- df.perm[1:n1,]
  df2.perm <- df.perm[(n1+1):n,]
  
  x.med.perm <- depthMedian(df1.perm, depth_params = list(method='Tukey'))
  y.med.perm <- depthMedian(df2.perm, depth_params = list(method='Tukey'))
  
  T2[perm] <- max(abs(x.med.perm-y.med.perm))  
  
}

plot(ecdf(T2), xlim=range(c(T2,T0)))
abline(v=T0,col=3,lwd=4)

p.val <- sum(T2>=T0)/B
p.val

##0.005, we reject the null hypothesis. As expected RBD and PD are not equally distributed



####### One-way MANOVA
#Via a MANOVA test we want to assess whether speech variables differ in distribution
#among the 3 groups

fit <- manova(as.matrix(speech.data) ~ Category)

T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
T_stat <- numeric(B)

n <- dim(speech.data)[1]

set.seed(seed)

for(perm in 1:B){
  permutation <- sample(1:n)
  group.perm <- Category[permutation]
  fit.perm <- manova(as.matrix(speech.data) ~ group.perm)
  T_stat[perm] <- -summary.manova(fit.perm, test="Wilks")$stats[1,2]
}

p_val <- sum(T_stat>=T0)/B
p_val  # 0

# Indeed, there are some variables in which there is significant difference among
# the groups (PD, RBD and Control)
# To check specifically which are those variables, we will use a one-way permutational ANOVA


### One-way ANOVA

g <- nlevels(Category)

## Parametric approach:
# (we first check whether the assumption of parametric ANOVA hold, if that's
#  not the case, we will use a nonparametric approach instead)

# 1) univariate normality in each group
Ps <- NULL
p <- dim(speech.data)[2]

for (idx in 1:p) {
  shap.res <- c(shapiro.test(speech.data[Category == 'Control', idx])$p,
          shapiro.test(speech.data[Category == 'RBD', idx])$p,
          shapiro.test(speech.data[Category == 'PD', idx])$p)
  Ps <- rbind(Ps, shap.res)
}
rownames(Ps) <- 1:p

#For which columns can't we reject the hypothesis that all the groups 
#are normally distributed?

normal.cols <- which(apply(Ps, 1, function(x) (all(x > alpha))))
# 6  9 14 18 20 21 23 

# 2) Bartlett's test for homogeneity of variances
bartlett.pvs <- numeric(length(normal.cols))
for (idx in 1:length(normal.cols)) {
  bartlett.pvs[idx] <- bartlett.test(speech.data[,idx], Category)$p.value
}

#homogeneity of variances is verified in the following columns
param.cols <- which(bartlett.pvs > alpha) #2 3 6 7
param.cols

param.aov.pvs <- NULL
for (idx in param.cols) {
  fit <- aov(speech.data[,idx] ~ Category)
  param.aov.pvs <- c(param.aov.pvs, summary(fit)[[1]][1,5])
}
param.aov.pvs < alpha
#only in column 2 the parametric ANOVA underlines a difference among
#the categories


## Nonparametric approach:

p.vals <- NULL

for (idx in 1:p) {
  
  if (! (idx %in% param.cols) ) {
    
    fit <- aov(speech.data[,idx] ~ Category)
    T0 <- summary(fit)[[1]][1,4]
    T_stat <- numeric(B) 
    
    set.seed(seed)
    
    for(perm in 1:B){
      permutation <- sample(1:n)
      variable_perm <- speech.data[permutation,idx]
      fit_perm <- aov(variable_perm ~ Category)
      
      T_stat[perm] <- summary(fit_perm)[[1]][1,4]
    }
    
    p.vals <- c(p.vals, sum(T_stat>=T0)/B)
    
  }
  
}

names(p.vals) <- colnames(speech.data)[-param.cols]
p.vals

# significant idx: 2 (0.005), 4 (0.001), 5 (0.039), 12 (0.06), 14 (0), 16 (0.002), 
#                  17 (0.07), 19 (0.001), 22 (0.016), 23 (0.04)

colnames(speech.data)[c(2,4,5,12,14,16,17,19,22,23)]
# "RST_r" "DPI_r" "DVI_r" "LRE_r" "RST_m" "DPI_m" "DVI_m" "DUS_m" "PIR_m" "RSR_m"
