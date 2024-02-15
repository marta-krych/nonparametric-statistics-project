
seed <- 20012024
B <- 1000

### Pre-processing part

parkinson <- read.csv("data/parkinson.csv")

colnames(parkinson) <- parkinson[1,]
parkinson <- parkinson[-1,]
rownames(parkinson) <- 1:dim(parkinson)[1]

control <- which(startsWith(parkinson$`Participant code`, "HC"))
pd <- which(startsWith(parkinson$`Participant code`, "PD"))
rbd <- which(startsWith(parkinson$`Participant code`, "RBD"))

parkinson <- parkinson[,-1]

parkinson <- type.convert(parkinson, as.is=T)

parkinson$Category <- as.factor(c(rep("PD",30), rep("RBD",50), rep("Control",50)))


# Extracting 'UPDRS III Total' column and subtracting 'Action tremor' from it
updrsIII <- as.numeric(parkinson[rbd, 13])
updrsIII.new <- updrsIII - (as.numeric(parkinson[rbd, 21]) + as.numeric(parkinson[rbd, 22]))


# Dataset with only speech variables
speech_dataset <- read.csv("data/speech.csv")

### Boxplots

library(ggplot2)
data <- speech_dataset
features <- names(data)[ !(names(data) %in% c("Category", "Gender"))]

for (feature in features) {
  plot <- ggplot(data, aes_string(x = "Category", y = feature, color = "Category")) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = paste("Box Plot of", feature), x = "Category", y = feature)
  
  print(plot)
}



### Checking for speech variables with Gaussian distribution

filtered_data <- subset(data, Category %in% c("PD", "Control"))
filtered_data_RBD <- subset(data, Category %in% c("RBD"))
features <- names(filtered_data)[!(names(filtered_data) %in% c("Age","Category","Gender"))]

normality_test_results <- list()
normal_variables <- c()

set.seed(seed)

for (feature in features) {
  pd_data <- filtered_data[filtered_data$Category == "PD", feature]
  hc_data <- filtered_data[filtered_data$Category == "Control", feature]
  
  shapiro_pd <- shapiro.test(pd_data)
  shapiro_hc <- shapiro.test(hc_data)
  
  normality_test_results[[feature]] <- list(PD_Shapiro_p_value = shapiro_pd$p.value, 
                                            HC_Shapiro_p_value = shapiro_hc$p.value)
  
  if (shapiro_pd$p.value > 0.05 && shapiro_hc$p.value > 0.05) {
    normal_variables <- c(normal_variables, feature)
  }
}

print(normal_variables)
# The following variables are normally distributed:
# "AST_r" "GVI_r" "RLR_r" "RST_m" "AST_m" "GVI_m" "DUF_m" "RLR_m" "RSR_m"



### Checking for variance homogeneity of the Gaussian variables

library(car)

homogeneity_test_results <- list()
homogeneous_normal_variables <- c()

set.seed(seed)

for (feature in normal_variables) {
  
  homogeneity_test_results[[feature]] <- list(var.test(filtered_data[filtered_data$Category == 'PD', feature],
                                                       filtered_data[filtered_data$Category == 'Control', feature])$p.value)
  
  if (homogeneity_test_results[[feature]] > 0.05) {
    homogeneous_normal_variables <- c(homogeneous_normal_variables, feature)
  }
}

print(homogeneous_normal_variables)

# The following variables are normally distributed with homogeneous variance:
# "AST_r" "GVI_r" "RLR_r" "RST_m" "GVI_m" "DUF_m" "RLR_m" "RSR_m"



### Testing the normally distributed variables for the difference in mean between HC and PD

t_test_results <- list()
significant_normal_variables <- c()

set.seed(seed)

for (feature in homogeneous_normal_variables) {
  pd_data <- filtered_data[filtered_data$Category == "PD", feature]
  hc_data <- filtered_data[filtered_data$Category == "Control", feature]
  
  test <- t.test(pd_data, hc_data, alternative="two.sided", mu=0)
  
  t_test_results[[feature]] <- list(t_test_p_value = test$p.value)
  
  if (test$p.value < 0.05) {
    significant_normal_variables <- c(significant_normal_variables, feature)
  }
}

print(significant_normal_variables)

# Among the normally distributed speech variables, only "RST_m" is significant



### Testing the non-normally distributed variables for the difference in median between HC and PD
### (using a permutation test)

not_normal_variables <- setdiff(features, homogeneous_normal_variables)

n_permutations <- 1000
significant_p_values <- c()

set.seed(seed)

for (feature in not_normal_variables) {
  
  x1 <- filtered_data[filtered_data$Category == 'PD', feature]
  x2 <- filtered_data[filtered_data$Category == 'Control', feature]
  
  x_pooled <- c(x1,x2)
  n  <- length(x_pooled)
  n1 <- length(x1)
  
  T0 <- abs(median(x1) - median(x2)) # test statistic
  T_stat <- numeric(n_permutations)
  
  set.seed(seed)
  
  for(perm in 1:n_permutations) {
    permutation <- sample(1:n)
    x_perm <- x_pooled[permutation]
    x1_perm <- x_perm[1:n1]
    x2_perm <- x_perm[(n1+1):n]
    
    T_stat[perm] <- abs(median(x1_perm) - median(x2_perm))
  }
  
  p_value <- sum(T_stat>=T0)/n_permutations
  
  if (p_value < 0.05) {
    significant_p_values[feature] <- p_value
  }
}

print(significant_p_values)

# Among the non-normally distributed speech variables, "RST_r" "DPI_r" "DUS_r" "DPI_m" 
# have significant p-values



### Boxplots of significant features (HC vs PD)

library(ggplot2)
data <- speech_dataset[speech_dataset$Category != "RBD",] # speech dataset
features <- c("RST_r", "DPI_r", "DUS_r", "DPI_m", "RST_m")

for (feature in features) {
  plot <- ggplot(data, aes(x=Category, y=data[,feature], fill=Category)) +
    scale_fill_manual(values=c("lightskyblue","lightcyan4")) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = paste("Box Plot of", feature), x = "Category", y = feature)
  
  print(plot)
}



### Modeling step: semiparametric model (i dind't run this section for the cluster analysis, with some exceptions)

library(mgcv)
library(splines)


## Building a model with: Gender RST_r DPI_r DUS_r RST_m DPI_m 

df_model_g <- filtered_data[,c(26,2,4,7,14,16)] # this yes
df_model_g$Category <- ifelse(filtered_data$Category == 'PD',1,0)

df_model_g <- data.frame(Gender=parkinson[c(pd, control),2], df_model_g)
df_model_g$Gender <- ifelse(df_model_g$Gender == 'F',1,0) # this yes

with(df_model_g, scatterplotMatrix(data.frame(Category, Gender, RST_r, DPI_r, DUS_r, DPI_m, RST_m)))

# linear relationship: Gender, RST_r, RST_m
# non-linear relationship: DPI_r, DPI_m, DUS_r

cor(as.matrix(df_model_g[,c(3,6)]))
# cor(RST_r, RST_m) = 0.6322871 (highly correlated) 
# --> keeping only RST_m to avoid collinearity in the linear part of the model


### 1)

model_semipar <- gam(Category ~ Gender + s(DPI_r) + s(DUS_r) + RST_m + s(DPI_m), 
                     family=binomial, data=df_model_g)
summary(model_semipar)

concurvity(model_semipar)  # some issues with DPI_r and DPI_m
# (concurvity is the analogous of collinearity but for smoothed terms)
# --> removing DPI_m, as DPI_r is significant in the model


### 2)

model_semipar_red <- gam(Category ~ Gender + s(DPI_r) + s(DUS_r) + RST_m, 
                         family=binomial, data=df_model_g)
summary(model_semipar_red)

concurvity(model_semipar_red)
# --> removing DUS_r


### 3)

model_semipar_f <- gam(Category ~ Gender + s(DPI_r) + RST_m, 
                       family=binomial, data=df_model_g)
summary(model_semipar_f)

concurvity(model_semipar_f)



### 'Traditional' prediction using model_semipar_f

test <- speech_dataset[rbd, c(4,14)]
test <- data.frame(Gender=ifelse(parkinson[rbd, 2] == 'F',1,0),
                   test)

predictions <- predict(model_semipar_f, test, se=TRUE, type="response")
# predictions <- predict(model_semipar_g, test, se=TRUE)
# pfit <- exp(predictions$fit)/(1 + exp(predictions$fit))

model_RBD <- ifelse(predictions$fit > .5, 1, 0)
updrs_RBD <- ifelse(updrsIII.new > 3, 1, 0) # this yes

# table and the error:
table(true_label=updrs_RBD, prediction=model_RBD)

errors = (model_RBD != updrs_RBD)
ER =  sum(errors)/length(updrs_RBD)
ER


#           prediction
# true_label   0  1
#          0  18  9
#          1  11 12
# ER = 0.4


DPI_r.grid <- seq(range(df_model_g$DPI_r)[1] - 0.05 * diff(range(df_model_g$DPI_r)),
                  range(df_model_g$DPI_r)[2] + 0.05 * diff(range(df_model_g$DPI_r)),
                  length.out=20)
RST_m.grid <- seq(range(df_model_g$RST_m)[1] - 0.05 * diff(range(df_model_g$RST_m)),
                  range(df_model_g$RST_m)[2] + 0.05 * diff(range(df_model_g$RST_m)),
                  length.out=20)

gridM <- expand.grid(DPI_r = DPI_r.grid, RST_m = RST_m.grid, Gender = 0)
predsM <- predict(model_semipar_f, gridM, type="response")
zM <- matrix(predsM, length(DPI_r.grid))

gridF <- expand.grid(DPI_r = DPI_r.grid, RST_m = RST_m.grid, Gender = 1)
predsF <- predict(model_semipar_f, gridF, type="response")
zF <- matrix(predsF, length(DPI_r.grid))


par(mfrow=c(1,2))
persp(DPI_r.grid, RST_m.grid, zM, xlab="DPI_r", ylab = "RST_m", zlab="risk",  theta = 230, phi = 20, col="blue")
persp(DPI_r.grid, RST_m.grid, zF, xlab="DPI_r", ylab = "RST_m", zlab="risk",  theta = 230, phi = 20, col="pink")


### 3D Logistic surfaces for the best model

library(plotly)

plotM <- plot_ly(x=DPI_r.grid, y=RST_m.grid, z=zM) %>% layout(scene = list(xaxis = list(title = "DPI_r"), 
                                                                           yaxis = list(title = "RST_m"),
                                                                           zaxis = list(title = "risk"))) %>%
  add_surface(colorscale = list(list(0,"e6e6fa"), list(1,"00008b")))
plotM

plotF <- plot_ly(x=DPI_r.grid, y=RST_m.grid, z=zF) %>% layout(scene = list(xaxis = list(title = "DPI_r"), 
                                                                           yaxis = list(title = "RST_m"),
                                                                           zaxis = list(title = "risk"))) %>%
  add_surface(colorscale = list(list(0,"ffe4e1"), list(1,"8b0000")))
plotF


### Full conformal prediction

library(conformalInference)

train_gam <- function(x, y, out=NULL){
  colnames(x) <- c('Gender','DPI_r','RST_m')
  train_data <- data.frame(y, x)
  model_gam <- gam(as.factor(y) ~ Gender + s(DPI_r) + RST_m, family=binomial, data=train_data)
}

predict_gam <- function(obj, new_x){
  new_x <- data.frame(new_x)
  colnames(new_x) <- c('Gender','DPI_r','RST_m')
  predict.gam(obj, new_x, type="response")
}

model_gam <- gam(Category ~ Gender + s(DPI_r) + RST_m, family=binomial, data=df_model_g)

target <- df_model_g$Category
covariates <- df_model_g[,c(1,4,6)]

c_preds <- conformal.pred(covariates, target, as.matrix(test), alpha=0.05, verbose=TRUE, 
                          train.fun=train_gam, predict.fun=predict_gam, num.grid.pts=200)

data.frame(lower=c_preds$lo, 
           prediction=c_preds$pred, 
           upper=c_preds$up)

model_RBD <- ifelse(c_preds$pred > .5, 1, 0)
updrs_RBD <- ifelse(updrsIII.new > 3, 1, 0)

# table and the error:
table(true_label=updrs_RBD, prediction=model_RBD)
errors = (model_RBD != updrs_RBD)
ER =  sum(errors)/length(updrs_RBD)
ER


### Let's see another approach: cluster analysis:
filtered_data_RBD <- filtered_data_RBD[,c(26,2,4,7,14,16)]
filtered_data_RBD$Gender <- ifelse(filtered_data_RBD$Gender == 'F',1,0)

### Hierarchical clustering:
dataset.e <- dist(filtered_data_RBD, method='euclidean')
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

fit <- manova(as.matrix(filtered_data_RBD) ~ cluster.ew)
summary.manova(fit,test="Wilks") 

T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
T0

set.seed(seed)
T_stat <- numeric(B)

for(perm in 1:B){
  # choose random permutation
  permutation <- sample(1:n_test)
  category.perm <- cluster.ew[permutation]
  fit.perm <- manova(as.matrix(filtered_data_RBD) ~ category.perm)
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
my_dataset <- rbind(df_model_g[1:30,], filtered_data_RBD[which(cluster.ew==1), ])
category <- as.factor(c( rep("0", 30) , rep("1", 26) ))

n1 <- 30
n2 <- 26
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
p_val # 0.009

# cluster 1 is Control-like


#  Continue with cluster 2:
my_dataset <- rbind(df_model_g[1:30,], filtered_data_RBD[which(cluster.ew==2), ])
category <- as.factor(c( rep("0", 30) , rep("1", 24) ))

n1 <- 30
n2 <- 24
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
p_val # 0.025

# Cluster 1: health-like (label 0)
# Cluster 2: PD-like (label 1)

cluster.ew[which(cluster.ew == 1)] <- 0
cluster.ew[which(cluster.ew == 2)] <- 1

# table and the error:
table(true_label=updrs_RBD, prediction=cluster.ew)
errors = (cluster.ew!= updrs_RBD)
ER =  sum(errors)/length(updrs_RBD)
ER # 0.38


### K-means:
result.k <- kmeans(filtered_data_RBD, centers=2)
length(which(result.k$cluster==1))

# Assess the quality of clustering:
n1 <- 30
n2 <- 20
n_test  <- n1+n2

g  <- 2
p  <- 6

fit <- manova(as.matrix(filtered_data_RBD) ~ result.k$cluster)
summary.manova(fit,test="Wilks") 

T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
T0

set.seed(seed)
T_stat <- numeric(B)

for(perm in 1:B){
  # choose random permutation
  permutation <- sample(1:n_test)
  category.perm <- result.k$cluster[permutation]
  fit.perm <- manova(as.matrix(filtered_data_RBD) ~ category.perm)
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

# Verify cluster 1:
my_dataset <- rbind(df_model_g[1:30,], filtered_data_RBD[which(result.k$cluster==1), ])
category <- as.factor(c( rep("0", 30) , rep("1", 41) ))

n1 <- 30
n2 <- 41
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
my_dataset <- rbind(df_model_g[1:30,], filtered_data_RBD[which(result.k$cluster==2), ])
category <- as.factor(c( rep("0", 30) , rep("1", 9) ))

n1 <- 30
n2 <- 9
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

result.k$cluster[which(result.k$cluster == 1)] <- 0
result.k$cluster[which(result.k$cluster == 2)] <- 1

# table and the error:
table(true_label=updrs_RBD, prediction=result.k$cluster)

errors = (result.k$cluster != updrs_RBD)
ER =  sum(errors)/length(updrs_RBD)
ER # 0.38 (in one run i had 0.36)



 

