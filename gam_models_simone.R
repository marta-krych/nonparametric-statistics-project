library(readr)
library(ggplot2)
library(tidyr)
library(mgcv)
library(progress)

B <- 1000
alpha <- 0.05
seed <- 2024
set.seed(seed)

###### SPEECH DATASET
speech.df <- read.csv("data/speech_no_outliers.csv")


rbd <- speech.df[which(speech.df$Category == "RBD"),-27]
no.rbd <- speech.df[-which(speech.df$Category == "RBD"),]

grouping <- speech.df$Category[-which(speech.df$Category == "RBD")]
grouping <- as.factor(grouping)

#reduced datasets with variables coming from the anova test
columns <- c("RST_r", "DPI_r" ,"RST_m", "DPI_m", "DUS_m","Age","Gender")
reduced.no.rbd <- no.rbd[,columns]
reduced.no.rbd$Gender <- factor(reduced.no.rbd$Gender )

#original dataset to retrieve age and gender and to compute UPDRSIII*
original.df <- read.csv("data/parkinson.csv")
colnames(original.df) <- original.df[1,]
original.df <- original.df[-1,]
rownames(original.df) <- 1:dim(original.df)[1]


updrsIII_new <- as.numeric(original.df[31:80,14]) -  as.numeric(original.df[31:80,22]) - 
  as.numeric(original.df[31:80,23])


#function for computing LOOCV perfromance with a given treshold 
cv.performance <- function(gam_model, treshold) {
  y <- gam_model$y
  perf <- numeric(length(y))
  
  for (i in 1:length(y)) {
    
    y.i <- y[-i]
    x <- data.frame(gam_model$model, row.names = 1:length(y))
    test.x <- data.frame(x[i,])
    dataset.i <- data.frame(x[-i,])
    colnames(test.x) <- colnames(x) <- colnames(dataset.i) <- colnames(gam_model$model)
    dataset.i$y <- y.i
    
    model.i <- gam(gam_model$formula, family=gam_model$family, data = dataset.i)
    
    pred.i <- as.numeric(predict(model.i, test.x, type="response") > treshold)
    
    perf[i] <- as.numeric( pred.i == y[i])
  }
  
  return(mean(perf))
}


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


