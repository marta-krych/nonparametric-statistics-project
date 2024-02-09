library(readr)
library(ggplot2)
library(tidyr)
library(mgcv)
library(progress)

B <- 1000
alpha <- 0.05
seed <- 2024
set.seed(seed)

###### VARIABLE SELECTION ON OUTLIERS-FREE DATASET

speech.df <- read.csv("speech_no_outliers.csv")

rbd <- speech.df[which(speech.df$Category == "RBD"),-27]
no.rbd <- speech.df[-which(speech.df$Category == "RBD"),-c(25:27)]

grouping <- speech.df$Category[-which(speech.df$Category == "RBD")]
grouping <- as.factor(grouping)


### ANOVA on the no.rbd dataset to highlight the difference among PD and Control

anova.test <- function(column, group, seed, nperm) {

  variable <- no.rbd[,column]
  fit <- aov(variable ~ grouping)
  T0 <- summary(fit)[[1]][1,4]

  T_stat <- numeric(nperm)
  n <- length(variable)
  set.seed(seed)

  for(perm in 1:nperm){
    # Permutation:
    permutation <- sample(1:n)
    variable_perm <- variable[permutation]
    fit_perm <- aov(variable_perm ~ grouping)

    # Test statistic:
    T_stat[perm] <- summary(fit_perm)[[1]][1,4]
  }

  hist(T_stat,xlim=range(c(T_stat,T0)), main=colnames(no.rbd)[column])
  abline(v=T0,col=3,lwd=2)

  p.val <- sum(T_stat>=T0)/nperm
  return(p.val)
}

p <- dim(no.rbd)[2]
p.values <- numeric(p)
names(p.values) <- colnames(no.rbd)

for (i in 1:p )
  p.values[i] <- anova.test(column=i, group=grouping, seed=seed, nperm=B)


p.values

##These are the columns in which the difference between the groups is noticeable.
names(p.values)[which(p.values < alpha)]
#"RST_r" "DPI_r" "RST_m" "DPI_m" "DUS_m"

c("RST_r", "DPI_r" ,"RST_m", "DPI_m", "DUS_m")

#dataset for PD+HC patients with only the significant variables, plus age and gender
reduced.no.rbd <- speech.df[speech.df$Category != "RBD",
            c(names(which(p.values < alpha)),"Age","Gender")] # here i do a runn with age and a run without it
# 
# reduced.no.rbd <- speech.df[speech.df$Category != "RBD",
#                             c("RST_r", "DPI_r" ,"RST_m", "DPI_m", "DUS_m","Age","Gender")]


###### LOGISTIC REGRESSION MODELS ON SPEECH VARIABLES

#original dataset to compute UPDRSIII* for RBD patients
original.df <- read.csv("parkinson.csv")
colnames(original.df) <- original.df[1,]
original.df <- original.df[-1,]
rownames(original.df) <- 1:dim(original.df)[1]


updrsIII_new <- as.numeric(original.df[31:80,14]) -  as.numeric(original.df[31:80,22]) - 
  as.numeric(original.df[31:80,23])


#BRUTE FORCE APPROACH: to build the best generalized linear model we build all the possible
#models using the variables selected and keep the best one in terms of LOOCV 
#performance on the training set (no.rbd patients)

cv.performance <- function(model, treshold) {
  y <- model$y
  perf <- numeric(length(y))
  
  for (i in 1:length(y)) {
    
    y.i <- y[-i]
    x <- data.frame(model$data, row.names = 1:length(y))
    test.x <- data.frame(x[i,])
    dataset.i <- data.frame(x[-i,])
    colnames(test.x) <- colnames(x) <- colnames(dataset.i) <- colnames(model$data)
    dataset.i$y <- y.i
    
    model.i <- glm(model$formula, family=model$family, data = dataset.i)
    
    pred.i <- as.numeric(predict(model.i, test.x, type="response") > treshold)
    
    perf[i] <- as.numeric( pred.i == y[i])
  }
  
  return(mean(perf))
}

### GENERALIZED LINEAR MODELS 
y <- as.numeric(grouping == "PD")
treshold <- 0.5
p <- dim(reduced.no.rbd)[2]

models.nvars <- function(nvars, treshold) {
  
  results <- as.data.frame(matrix(0, ncol = p+1, nrow = choose(p,nvars)))
  colnames(results) <- c(colnames(reduced.no.rbd), "cv.performance")
  
  combinations <- combn(p,nvars)
  
  pb <- progress_bar$new(total=ncol(combinations))
  pb$tick(0)
  
  max <- 0.5 
  for (j in 1:ncol(combinations)) {
    
    columns <- combinations[,j]
    newdata <- as.data.frame(reduced.no.rbd[,columns])
    colnames(newdata) <- colnames(reduced.no.rbd)[columns]
    model <- glm(y ~ ., family = 'binomial', data=newdata)
    
    perf <- cv.performance(model, treshold)
    if (perf >= max)  {
      results[j,dim(results)[2]] <- perf
      results[j, columns] <- 1
      max <- perf
    }
    
    
    pb$tick()
  }
  
  results.ordered <- results[order(results$cv.performance,decreasing=TRUE),]
  return(results.ordered)
  
}

get.max <- function(dataframe) {
  max.row <- which(dataframe$cv.performance == max(dataframe$cv.performance))
  return(dataframe[max.row,])
}

results.1.var <- models.nvars(1, treshold)
best.1.var <- get.max(results.1.var)

results.2.var <- models.nvars(2, treshold)
best.2.var <- get.max(results.2.var)

results.3.var <- models.nvars(3, treshold)
best.3.var <- get.max(results.3.var)

results.4.var <- models.nvars(4, treshold)
best.4.var <- get.max(results.4.var)

results.5.var <- models.nvars(5, treshold)
best.5.var <- get.max(results.5.var)

results.6.var <- models.nvars(6, treshold)
best.6.var <- get.max(results.6.var)

results.7.var <- models.nvars(7, treshold)
best.7.var <- get.max(results.7.var)


#the best models (2,3 variables) reach all 71.21% of LOOCV accuracy. 
#So we perform a ChiSq-test to understand if the least complex model is sufficient.

#models definition (without repeating all the loops)
model2 <- glm(y ~ DPI_m + Gender, family = 'binomial', data=reduced.no.rbd)

summary(model2)$coefficients

#               Estimate  Std. Error   z value    Pr(>|z|)
# (Intercept) -3.63391837 1.257730950 -2.889265 0.003861432
# DPI_m        0.01835659 0.006039795  3.039273 0.002371500
# GenderM     -1.33746409 0.709909111 -1.883993 0.059565866

model3 <- glm(y ~ DUS_m + DPI_m + Gender, family = 'binomial', data=reduced.no.rbd)

summary(model3)$coefficients

#               Estimate  Std. Error    z value    Pr(>|z|)
# (Intercept) -4.46773707 1.598333964 -2.7952463 0.005186019
# DUS_m        0.04912476 0.055641980  0.8828722 0.377305325
# DPI_m        0.01602978 0.006451571  2.4846323 0.012968533
# GenderM     -1.28458074 0.717531236 -1.7902785 0.073409152

anova(model2,model3, test = "Chisq") #0.3762, reduced model is sufficient
#Model 2 has all significant terms at level 0.06.


cv.performance(model2, treshold) #0.7121212


#Let's see how the model classify the RBD patients
model2.rbd.predictions <- predict(model2, rbd, type="response") #probabilities
model2.rbd.classifications <- as.numeric(model2.rbd.predictions>treshold)

table(true_label = as.numeric(updrsIII_new>3), prediction = model2.rbd.classifications)
#           prediction
# true_label  0  1
#           0 17 10
#           1 14  9

errors = (as.numeric(updrsIII_new>3) != model2.rbd.classifications)
ER =  sum(errors)/length(updrsIII_new)
ER




############GENERALIZED ADDITIVE MODELS

y <- as.factor(as.numeric(grouping == "PD"))
treshold <- 0.5
p <- dim(reduced.no.rbd)[2]



###Let's start from the variables of the best parametric model

parametric.best <- gam(y ~ DPI_m + Gender, family = 'binomial', data=reduced.no.rbd)
summary(parametric.best) #significant terms at level 6%
cv.performance(parametric.best, treshold) #71.12%
#           prediction
# true_label   0  1
#           0 17 10
#           1 14  9

gam0 <- gam(y ~ DPI_m + Gender + Age, family = 'binomial', data=reduced.no.rbd)
summary(gam0) #only DPI_m is significant 
cv.performance(gam0, treshold) #68.18%

gam1 <- gam(y ~ s(DPI_m,bs='cr') + Gender, family = 'binomial', data=reduced.no.rbd)
summary(gam1) # not significant terms
cv.performance(gam1, treshold) #71.12%

### Go back to the full models

gam2 <- gam(y ~ s(RST_r,bs='cr') + s(DPI_r,bs='cr') + s(RST_m,bs='cr') + 
              s(DPI_m,bs='cr') + s(DUS_m,bs='cr') + Age + 
              Gender, family = 'binomial', data=reduced.no.rbd)
summary(gam2) # not significant terms
cv.performance(gam2, treshold) #62.12%


gam3 <- gam(y ~ RST_r + s(DPI_r,bs='cr') + s(RST_m,bs='cr') + 
              s(DPI_m,bs='cr') + DUS_m + Age + 
              Gender, family = 'binomial', data=reduced.no.rbd)
summary(gam3) # not significant terms
cv.performance(gam3, treshold) #60.61%



gam4 <- gam(y ~ RST_r + s(DPI_r,bs='cr') + RST_m + 
              s(DPI_m,bs='cr') + DUS_m + 
              Gender , family = 'binomial', data=reduced.no.rbd)
summary(gam4) # not significant terms
cv.performance(gam4, treshold) #74.24%
table(true_label = as.numeric(updrsIII_new>3), 
      prediction = as.numeric(predict(gam4, rbd, "response")>treshold))
#           prediction
# true_label   0  1
#           0 21  6
#           1 16  7


gam5 <- gam(y ~ RST_r + s(DPI_m,bs='cr') + RST_m + DUS_m + 
              Gender , family = 'binomial', data=reduced.no.rbd)
summary(gam5) # some terms are significant
cv.performance(gam5, treshold) #66.67%


gam6 <- gam(y ~ RST_r + s(DPI_m,bs='cr') + RST_m + DUS_m + 
              Gender -1 , family = 'binomial', data=reduced.no.rbd)
summary(gam6) # not significant terms
cv.performance(gam6, treshold) #66.67%


gam7 <- gam(y ~ RST_r + s(DPI_m,bs='cr') + RST_m +
              Gender , family = 'binomial', data=reduced.no.rbd)
summary(gam7) # some terms are significant
cv.performance(gam7, treshold) #65.15%


gam8 <- gam(y ~ RST_r + s(DPI_m,bs='cr') +
              Gender , family = 'binomial', data=reduced.no.rbd)
summary(gam8) # parametric terms are significant at level 10%
cv.performance(gam8, treshold) #68.18%
table(true_label = as.numeric(updrsIII_new>3), 
      prediction = as.numeric(predict(gam8, rbd, "response")>treshold))
#           prediction
# true_label   0  1
#           0 22  5
#           1 13  10



plot(gam8)


gam9 <- gam(y ~ RST_r + DPI_m +
              Gender , family = 'binomial', data=reduced.no.rbd)
summary(gam9) #DPI_m and Gender are significant at level 5%
cv.performance(gam9, treshold) #65.15%





##the best model in terms of cv performance is model gam4 but its terms are not 
#significant and the confusion matrix is poor. 
#The best confusion matrix is given by model gam8, which terms are significant at 
#level 10% even if CV performance is not so high. 

concurvity(gam8)
#concurvity (colinaerity in smooth terms) is not a problem 


#by changing the threshold both the cv performance and the confusion matrix get 
#better
cv.performance(gam8, 0.6)  #72.73%
table(true_label = as.numeric(updrsIII_new>3), 
      prediction = as.numeric(predict(gam8, rbd, "response")>0.6))
#           prediction
# true_label   0  1
#           0 24  3
#           1 14  9

#sensitivity = 0.3913
#specificity = 0.8889

errors = (as.numeric(updrsIII_new>3) != as.numeric(predict(gam8, rbd, "response")>0.6))
ER =  sum(errors)/length(updrsIII_new)
ER # 0.34


### Let's see another approach: cluster analysis:
rbd_of_interest <- speech.df[speech.df$Category == "RBD",
                            c(names(which(p.values < alpha)),"Age","Gender")]

rbd_of_interest$Gender <- ifelse(rbd_of_interest$Gender == 'F',1,0)

reduced.no.rbd$Gender <- ifelse(reduced.no.rbd$Gender == 'F',1,0)

### Hierarchical clustering:
dataset.e <- dist(rbd_of_interest, method='euclidean')
dataset.ew<- hclust(dataset.e, method='ward.D2')


x11()
plot(dataset.ew, main='euclidean-ward', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(dataset.ew, k=2)

cluster.ew <- cutree(dataset.ew, k=2)
length(which(cluster.ew==1))

# Asses the goodness of clustering with permutational MANOVA:
n1 <- 26
n2 <- 24
n_test  <- n1+n2

g  <- 2
p  <- 6

fit <- manova(as.matrix(rbd_of_interest) ~ cluster.ew)
summary.manova(fit,test="Wilks") 

T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
T0

set.seed(seed)
T_stat <- numeric(B)

for(perm in 1:B){
  # choose random permutation
  permutation <- sample(1:n_test)
  category.perm <- cluster.ew[permutation]
  fit.perm <- manova(as.matrix(rbd_of_interest) ~ category.perm)
  T_stat[perm] <- -summary.manova(fit.perm,test="Wilks")$stats[1,2]
}

x11()
hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

x11()
plot(ecdf(T_stat),xlim=c(-2,1))
abline(v=T0,col=3,lwd=4)

p_val <- sum(T_stat>=T0)/B
p_val # 0

# Assess similarity with PD with permutational Manova:
# Start with cluster 1:
my_dataset <- rbind(reduced.no.rbd[1:23,], rbd_of_interest[which(cluster.ew==1), ])
category <- as.factor(c( rep("0", 23) , rep("1", 26) ))

n1 <- 23
n2 <- 26
n_test  <- n1+n2

g  <- 2
p  <- 8

fit <- manova(as.matrix(my_dataset) ~ category)
summary.manova(fit,test="Wilks") 

T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
T0

set.seed(seed)
T_stat <- numeric(B)

for(perm in 1:B){
  # choose random permutation
  permutation <- sample(1:n_test)
  category.perm <- category[permutation]
  fit.perm <- manova(as.matrix(my_dataset) ~ category.perm)
  T_stat[perm] <- -summary.manova(fit.perm,test="Wilks")$stats[1,2]
}

x11()
hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

x11()
plot(ecdf(T_stat),xlim=c(-2,1))
abline(v=T0,col=3,lwd=4)

p_val <- sum(T_stat>=T0)/B
p_val # 0.035 with Age variable (0.021 without)

# cluster 1: we are in limit case with Age variable while is clearly control-like  without it


#  Continue with cluster 2:
my_dataset <- rbind(reduced.no.rbd[1:23,],  rbd_of_interest[which(cluster.ew==2), ])
category <- as.factor(c( rep("0", 23) , rep("1", 24) ))

n1 <- 23
n2 <- 24
n_test  <- n1+n2

g  <- 2
p  <- 8

fit <- manova(as.matrix(my_dataset) ~ category)
summary.manova(fit,test="Wilks") 

T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
T0

set.seed(seed)
T_stat <- numeric(B)

for(perm in 1:B){
  # choose random permutation
  permutation <- sample(1:n_test)
  category.perm <- category[permutation]
  fit.perm <- manova(as.matrix(my_dataset) ~ category.perm)
  T_stat[perm] <- -summary.manova(fit.perm,test="Wilks")$stats[1,2]
}

x11()
hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

x11()
plot(ecdf(T_stat),xlim=c(-2,1))
abline(v=T0,col=3,lwd=4)

p_val <- sum(T_stat>=T0)/B
p_val # 0.032 with variable Age (0.018 without)


# adjust the labels according with the two previous p-values
cluster.ew[which(cluster.ew == 1)] <- 0
cluster.ew[which(cluster.ew == 2)] <- 1

# table and the error:
table(true_label=as.numeric(updrsIII_new>3), prediction=cluster.ew)
errors = (cluster.ew!= as.numeric(updrsIII_new>3))
ER =  sum(errors)/length(as.numeric(updrsIII_new>3))
ER # 0.38 with and without Age


### K-means:
result.k <- kmeans(rbd_of_interest, centers=2)
length(which(result.k$cluster==1))

# Assess the quality of clustering:
n1 <- 9 
n2 <- 41 
n_test  <- n1+n2

g  <- 2
p  <- 7

fit <- manova(as.matrix(rbd_of_interest) ~ result.k$cluster)
summary.manova(fit,test="Wilks") 

T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
T0

set.seed(seed)
T_stat <- numeric(B)

for(perm in 1:B){
  # choose random permutation
  permutation <- sample(1:n_test)
  category.perm <- result.k$cluster[permutation]
  fit.perm <- manova(as.matrix(rbd_of_interest) ~ category.perm)
  T_stat[perm] <- -summary.manova(fit.perm,test="Wilks")$stats[1,2]
}

x11()
hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

x11()
plot(ecdf(T_stat),xlim=c(-2,1))
abline(v=T0,col=3,lwd=4)

p_val <- sum(T_stat>=T0)/B
p_val # 0

# Verify cluster 1
my_dataset <- rbind(reduced.no.rbd[1:23,], rbd_of_interest[which(result.k$cluster==1), ])
category <- as.factor(c( rep("0", 23) , rep("1", 20) )) # rep("1", 9), with age variable

n1 <- 23
n2 <- 20 # with Age variable
n_test  <- n1+n2

g  <- 2
p  <- 6

fit <- manova(as.matrix(my_dataset) ~ category)
summary.manova(fit,test="Wilks") 

T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
T0

set.seed(seed)
T_stat <- numeric(B)

for(perm in 1:B){
  # choose random permutation
  permutation <- sample(1:n_test)
  category.perm <- category[permutation]
  fit.perm <- manova(as.matrix(my_dataset) ~ category.perm)
  T_stat[perm] <- -summary.manova(fit.perm,test="Wilks")$stats[1,2]
}

x11()
hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

x11()
plot(ecdf(T_stat),xlim=c(-2,1))
abline(v=T0,col=3,lwd=4)

p_val <- sum(T_stat>=T0)/B
p_val 



# cluster 2:
my_dataset <- rbind(reduced.no.rbd[1:23,], rbd_of_interest[which(result.k$cluster==2), ])
category <- as.factor(c( rep("0", 23) , rep("1", 30) ))

n1 <- 23
n2 <- 30
n_test  <- n1+n2

g  <- 2
p  <- 6

fit <- manova(as.matrix(my_dataset) ~ category)
summary.manova(fit,test="Wilks") 

T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
T0

set.seed(seed)
T_stat <- numeric(B)

for(perm in 1:B){
  # choose random permutation
  permutation <- sample(1:n_test)
  category.perm <- category[permutation]
  fit.perm <- manova(as.matrix(my_dataset) ~ category.perm)
  T_stat[perm] <- -summary.manova(fit.perm,test="Wilks")$stats[1,2]
}

x11()
hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

x11()
plot(ecdf(T_stat),xlim=c(-2,1))
abline(v=T0,col=3,lwd=4)

p_val <- sum(T_stat>=T0)/B
p_val 

# adjust the labels according to the two previous p-values
result.k$cluster[which(result.k$cluster == 1)] <- 0
result.k$cluster[which(result.k$cluster == 2)] <- 1

# table and the error:
table(true_label=as.numeric(updrsIII_new>3), prediction=result.k$cluster)
errors = (result.k$cluster!= as.numeric(updrsIII_new>3))
ER =  sum(errors)/length(as.numeric(updrsIII_new>3))
ER 

# > 0.6 with and without variable Age;
# It has happened that ER = 0.38, i think that the removal of the outlier is much sgnificant.


# Let's use a robust version of K-means: K-medoids:
if (!requireNamespace("cluster", quietly = TRUE)) {
  install.packages("cluster")
}

# Load the 'cluster' package
library(cluster)

result <- pam(rbd_of_interest, 2)
result$clustering
length(which(result$clustering==1))


# Assess the quality of clustering:
n1 <- 31# with Age variable
n2 <- 19 # " "
n_test  <- n1+n2

g  <- 2
p  <- 7

fit <- manova(as.matrix(rbd_of_interest) ~ result$clustering)
summary.manova(fit,test="Wilks") 

T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
T0

set.seed(seed)
T_stat <- numeric(B)

for(perm in 1:B){
  # choose random permutation
  permutation <- sample(1:n_test)
  category.perm <- result$clustering[permutation]
  fit.perm <- manova(as.matrix(rbd_of_interest) ~ category.perm)
  T_stat[perm] <- -summary.manova(fit.perm,test="Wilks")$stats[1,2]
}

x11()
hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

x11()
plot(ecdf(T_stat),xlim=c(-2,1))
abline(v=T0,col=3,lwd=4)

p_val <- sum(T_stat>=T0)/B
p_val # 0

# Verify cluster 1
my_dataset <- rbind(reduced.no.rbd[1:23,], rbd_of_interest[which(result$clustering==1), ])
category <- as.factor(c( rep("0", 23) , rep("1", 31) )) # rep("1", 9), with age variable

n1 <- 23
n2 <- 31 # with Age variable
n_test  <- n1+n2

g  <- 2
p  <- 6

fit <- manova(as.matrix(my_dataset) ~ category)
summary.manova(fit,test="Wilks") 

T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
T0

set.seed(seed)
T_stat <- numeric(B)

for(perm in 1:B){
  # choose random permutation
  permutation <- sample(1:n_test)
  category.perm <- category[permutation]
  fit.perm <- manova(as.matrix(my_dataset) ~ category.perm)
  T_stat[perm] <- -summary.manova(fit.perm,test="Wilks")$stats[1,2]
}

x11()
hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

x11()
plot(ecdf(T_stat),xlim=c(-2,1))
abline(v=T0,col=3,lwd=4)

p_val <- sum(T_stat>=T0)/B
p_val 



# cluster 2:
my_dataset <- rbind(reduced.no.rbd[1:23,], rbd_of_interest[which(result$clustering==2), ])
category <- as.factor(c( rep("0", 23) , rep("1", 19) ))

n1 <- 23
n2 <- 19
n_test  <- n1+n2

g  <- 2
p  <- 6

fit <- manova(as.matrix(my_dataset) ~ category)
summary.manova(fit,test="Wilks") 

T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
T0

set.seed(seed)
T_stat <- numeric(B)

for(perm in 1:B){
  # choose random permutation
  permutation <- sample(1:n_test)
  category.perm <- category[permutation]
  fit.perm <- manova(as.matrix(my_dataset) ~ category.perm)
  T_stat[perm] <- -summary.manova(fit.perm,test="Wilks")$stats[1,2]
}

x11()
hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

x11()
plot(ecdf(T_stat),xlim=c(-2,1))
abline(v=T0,col=3,lwd=4)

p_val <- sum(T_stat>=T0)/B
p_val 


# adjust the labels
result.k$cluster[which(result$clustering == 1)] <- 0
result.k$cluster[which(result$clustering == 2)] <- 0

# table and the error:
table(true_label=as.numeric(updrsIII_new>3), prediction=result$clustering)
errors = (result.k$cluster!= as.numeric(updrsIII_new>3))
ER =  sum(errors)/length(as.numeric(updrsIII_new>3))
ER # >= 0.6








 
