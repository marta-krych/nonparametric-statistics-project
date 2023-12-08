library(readr)
library(ggplot2)
library(tidyr)
library(progress)

B <- 1000
alpha <- 0.05
seed <- 5122023
set.seed(seed)

###### VARIABLE SELECTION ON OUTLIERS-FREE DATASET

speech.df <- read.csv("data/speech_no_outliers.csv")

ids <- speech.df$X 
speech.df$X <- NULL
rownames(speech.df) <- ids

rbd <- speech.df[which(speech.df$Category == "RBD"),-25]
no.rbd <- speech.df[-which(speech.df$Category == "RBD"),-25]

grouping <- speech.df$Category[-which(speech.df$Category == "RBD")]
grouping <- as.factor(grouping)


### ANOVA on the no.rbd dataset to highlight the difference among PD and Control

# anova.test <- function(column, group, seed, nperm) {
#   
#   variable <- no.rbd[,column]
#   fit <- aov(variable ~ grouping)
#   T0 <- summary(fit)[[1]][1,4]
#   
#   T_stat <- numeric(nperm) 
#   n <- length(variable)
#   set.seed(seed)
#   
#   for(perm in 1:nperm){
#     # Permutation:
#     permutation <- sample(1:n)
#     variable_perm <- variable[permutation]
#     fit_perm <- aov(variable_perm ~ grouping)
#     
#     # Test statistic:
#     T_stat[perm] <- summary(fit_perm)[[1]][1,4]
#   }
#   
#   hist(T_stat,xlim=range(c(T_stat,T0)), main=colnames(no.rbd)[column])
#   abline(v=T0,col=3,lwd=2)
#   
#   p.val <- sum(T_stat>=T0)/nperm
#   return(p.val)
# }
# 
# p <- dim(no.rbd)[2]
# p.values <- matrix(NA, ncol=p, nrow = 1)
# colnames(p.values) <- colnames(no.rbd)
# 
# for (i in 1:p ) 
#   p.values[i] <- anova.test(column=i, group=grouping, seed=seed, nperm=B)
# 
# 
# p.values
# 
# #These are the columns in which the difference between the groups is noticeable.
# #the value 0.15 is high in order to keep a sufficient number of variables
# colnames(p.values)[which(p.values < 0.15)]

# reduced.no.rbd <- no.rbd[, which(p.values < 0.15)]
# reduced.rbd <- rbd[,which(p.values < 0.15)]

#without repeating everything
columns <- c(1,2,4,14,16,19,22,23,24)
reduced.no.rbd <- data.frame(scale(no.rbd[,columns]))
reduced.rbd <- data.frame(scale(rbd[,columns]))

###### LOGISTIC REGRESSION MODELS ON SPEECH VARIABLES

#original dataset to retrieve age and gender and to compute UPDRSIII*
original.df <- read.csv("data/parkinson.csv")
colnames(original.df) <- original.df[1,]
original.df <- original.df[-1,]
rownames(original.df) <- 1:dim(original.df)[1]

rows <- as.numeric(rownames(reduced.no.rbd))
reduced.no.rbd$Age <- as.numeric(original.df$`Age (years)`[rows])
reduced.no.rbd$Gender <- original.df$Gender[rows]

rows <- as.numeric(rownames(reduced.rbd))
reduced.rbd$Age <- as.numeric(original.df$`Age (years)`[rows])
reduced.rbd$Gender <- original.df$Gender[rows]

updrsIII_new <- as.numeric(original.df[31:80,14]) -  as.numeric(original.df[31:80,22]) - 
  as.numeric(original.df[31:80,23])

#function for computing probability estimates
pred.model <- function(model, dataset) {
  preds <- predict(model, dataset, se=T)
  pfit <- exp(preds$fit )/(1+ exp( preds$fit ))
  return(pfit)
}

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
    
    pred.i <- as.numeric(pred.model(model.i, test.x) > treshold)
    
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

results.8.var <- models.nvars(8, treshold)
best.8.var <- get.max(results.8.var)

results.9.var <- models.nvars(9, treshold)
best.9.var <- get.max(results.9.var)

results.10.var <- models.nvars(10, treshold)
best.10.var <- get.max(results.10.var)

results.11.var <- models.nvars(11, treshold)
best.11.var <- get.max(results.11.var)

#the best models (3,4,5 variables) reach all 74.24% of LOOCV accuracy. 
#So we perform a ChiSq-test to understand if the least complex model is sufficient.

#models definition (without repeating all the loops)
model3 <- glm(y ~ EST_r + DPI_r + Gender, family = 'binomial', data=reduced.no.rbd)
model4 <- glm(y ~ EST_r + DPI_r + LRE_m + Gender, family = 'binomial', data=reduced.no.rbd)
model5 <- glm(y ~ EST_r + DPI_r + RST_m + DPI_m + Gender, family = 'binomial', data=reduced.no.rbd)
model5.1 <- glm(y ~ EST_r + DPI_r + LRE_m + Age + Gender, family = 'binomial', data=reduced.no.rbd)

anova(model3,model4, test = "Chisq") #0.8973, reduced model is sufficient
#anova(model3,model5.1, test = "Chisq") #0.989, reduced model is sufficient

#the simpler linear model in terms of number of covariates among all the best LOOCV-score models
#is model3
summary(model3)

#Let's see if we can consider a reduction of model3 not losing a lot of LOOCV-score ,
#a possible solution is removing the intercept
model3.reduced <- glm(y ~  EST_r + DPI_r + Gender -1 , family = 'binomial', data=reduced.no.rbd)
summary(model3.reduced)


summary(model5)
summary(model5.1)
#not very good p.values for the coefficients

cv.performance(model3, treshold) #74.2%
cv.performance(model3.reduced, treshold) #74.2%

#Let's see how both models classify the RBD patients
model3.rbd.predictions <- as.numeric(pred.model(model3, reduced.rbd)>treshold)
model3.rbd.predictions

model3.reduced.rbd.predictions <- as.numeric(pred.model(model3.reduced, reduced.rbd)>treshold)
model3.reduced.rbd.predictions
#same predictions as model3, let's stick with the reduced version then.

table(true_label = as.numeric(updrsIII_new>3), prediction = model3.reduced.rbd.predictions)




















 
