library(readr)
library(ggplot2)
library(tidyr)
library(progress)

B <- 1000
alpha <- 0.05
seed <- 2024

###### VARIABLE SELECTION ON OUTLIERS-FREE DATASET

speech.df <- read.csv("data/speech_no_outliers.csv")

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
            c(names(which(p.values < alpha)),"Age","Gender")]

# reduced.no.rbd <- speech.df[speech.df$Category != "RBD",
#                             c("RST_r", "DPI_r" ,"RST_m", "DPI_m", "DUS_m","Age","Gender")]


###### LOGISTIC REGRESSION MODELS ON SPEECH VARIABLES

#original dataset to compute UPDRSIII* for RBD patients
original.df <- read.csv("data/parkinson.csv")
colnames(original.df) <- original.df[1,]
original.df <- original.df[-1,]
rownames(original.df) <- 1:dim(original.df)[1]


updrsIII_new <- as.numeric(original.df[31:80,14]) -  as.numeric(original.df[31:80,22]) - 
  as.numeric(original.df[31:80,23])


#BRUTE FORCE APPROACH: to build the best generalized linear model I build all the possible
#models using the variables I selected and keep the best one in terms of Cross Validation 
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



















 
