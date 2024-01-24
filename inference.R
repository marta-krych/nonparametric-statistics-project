library(readr)
library(ggplot2)
library(tidyr)
library(DepthProc)

parkinson <- read_csv("data/parkinson_preprocessed.csv")

B <- 1000
alpha <- 0.05
seed <- 2024

############ Exploratory analysis on categorical variables

#we want to test if the 2 populations PD and RBD have the same distribution in each of these variables

pooled <- sapply(parkinson[1:80,9:35], as.numeric)

category <- as.factor(parkinson$Category[1:80])

for (i in 1:ncol(pooled))
  boxplot(pooled[,i] ~ category, main = colnames(pooled)[i])

pd <- which(category=="PD")
rbd <- which(category=="RBD")

n1 <- length(pd)
n2 <- length(rbd)
n <- n1 + n2

#Mann-Whitney U test: H0: X =^d Y vs H1: X !=^d Y
p.values <- numeric(27)
names(p.values) <- colnames(pooled)

for (i in 1:27) {
  p.values[i] <- wilcox.test(x=pooled[pd,i], y=pooled[rbd,i], 
                             paired=F, conf.level = 1-alpha, conf.int = TRUE)$p.value
}
names(which(p.values >= alpha))
#low p-values in every variable except "Tremor LLE" and "Arising from Chair"
#We can state that PD and RBD are equally distributed only in this variables




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
#(maybe outliers or the presence of HC patitens influence the result)


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



