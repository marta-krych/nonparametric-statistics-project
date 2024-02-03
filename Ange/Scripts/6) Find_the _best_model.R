## Useful libraries:
library(ISLR2)
library(car)
library(carData)
library(kernlab)
library(locfit)
library(np)
library(splines)
library(fda)
library(magrittr)
library(KernSmooth)
library(glmnet)
library(splines)
library(caret)
library(mgcv)
library(rgl)
library(pbapply)

## set the seed and B:
B = 1e5
seed = 26111992

## Import the dataset:
dataset <- read.csv("speech_no_outliers.csv")
dataset <- dataset[,-1]
View(dataset)


## Settings to search the best model:
reduced_dataset <- dataset[1:66,]
parkinson_ill <- dataset[1:23,]
control <- dataset[24:66,]
RBD <- dataset[67:116,]
# reduced_dataset$Category <- as.numeric(factor(reduced_dataset[,25]))
# reduced_dataset$Category[which(reduced_dataset$Category==2)] <- 0
# 0 are the control patients
# 1 are the parkinson patients

## Throught the hierarchical clustering i estimatre the labels of the RBD patients:
dataset.e <- dist(RBD[,-25], method='euclidean')

# proceed with the Euclidean distance.
dataset.ew <- hclust(dataset.e, method='ward.D2')

x11()
plot(dataset.ew, main='euclidean-ward', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(dataset.ew, k=2)

cluster.ew <- cutree(dataset.ew, k=2)
cluster.ew

# number on cluster 1:
length(which(cluster.ew==1))

# Which cluster is more similar to PD? Let's discover with permutational Manova:
my_dataset <- rbind(parkinson_ill[,-25], RBD[which(cluster.ew==1), -25])
PD_RBD_category <- as.factor(c( rep("0", 23) , rep("2", 31) )) 

    
n1 <- 23
n2 <- 31
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
p_val # 0.69374 for cluster 1 and 0.09972 for cluster 2

#cluster.ew[which(cluster.ew == 1)] <- 0
cluster.ew[which(cluster.ew == 2)] <- 0
#cluster.ew

## attach the dataset:
attach(reduced_dataset)

## Try different models:

# Logistic regression:
model_1 <- glm(category ~ RST_m+AST_m+DPI_m+DVI_m+GVI_m+DUS_m+PIR_m+LRE_m , data = reduced_dataset, family = binomial)
summary(model_1)

model_2 <- glm(category ~ AST_m+DPI_m+GVI_m+DUS_m, data = reduced_dataset, family = binomial)
summary(model_2)

model_3 <- glm(category ~ DPI_m+GVI_m, data = reduced_dataset, family = binomial)
summary(model_3)

model_4 <- glm(category ~ DPI_m, data = reduced_dataset, family = binomial)
summary(model_4)


# Polinomial logistic regression:
DPI_m_2 <- (DPI_m)^2
model_5 <- glm(category ~ DPI_m + DPI_m_2 , data = reduced_dataset, family = binomial)
summary(model_5) # nononono

DPI_m_sqrt <- sqrt(DPI_m)
model_6 <- glm(category ~ DPI_m + DPI_m_sqrt , data = reduced_dataset, family = binomial)
summary(model_6)

# Local logistic regression (training+validation):      
# •  Lasso regression:

          #   - set the variables: 
          x <- model.matrix( ~ RST_m+AST_m+DPI_m+DVI_m+GVI_m+DUS_m+PIR_m+LRE_m)[,-1]
          y <- as.factor(reduced_dataset$Category)
        
        #   - Let's set a grid of candidate lambda's for the estimate:
        lambda.grid <- 10^seq(5,-3,length=100)
        fit.lasso <- glmnet(x,y, family = "binomial", lambda = lambda.grid) # default: alpha=1 -> lasso 
        #    NB: if alpha=0 -> ridge regression. What follows can be done also with ridge!
        
        #    plot the coefficients:
        plot(fit.lasso,xvar='lambda',label=TRUE, col = rainbow(dim(x)[2]))
        legend('topright', dimnames(x)[[2]], col =  rainbow(dim(x)[2]), lty=1, cex=1)
        
        #   - Let's set lambda via cross validation:
        cv.lasso <- cv.glmnet(x,y,family = "binomial", lambda=lambda.grid) # default: 10-fold CV # remeber to set alpha = 0 if want to perform Ridge
        
        bestlam.lasso <- cv.lasso$lambda.min
        bestlam.lasso
        
        plot(cv.lasso)
        abline(v=log(bestlam.lasso), lty=1)
        
        #   - Get the coefficients for the optimal lambda:
        coef.lasso <- predict(fit.lasso, s=bestlam.lasso, type = 'coefficients')
        coef.lasso 
        
#   • Number of rows in the dataset
     n <- nrow(reduced_dataset)

#   • Create placeholders for predictions and actual values
     predictions <- rep(NA, n)
     actual_values <- rep(NA, n)

#   • Perform leave-one-out cross-validation
     for (i in 1:n) {
       # Exclude the i-th observation
       train_data <- reduced_dataset[-i, ]
       test_data <- reduced_dataset[i, ]
      
      # Fit a local logistic regression model on the training data
      model <- locfit((as.numeric(as.factor(Category))- 1) ~ DPI_m+DUS_m, alpha = 0.2, family = "binomial", data = train_data)
      
      # Predict on the left-out data point
      prediction <- predict(model, newdata = test_data[-25], type = "response")
      predictions[i] <- prediction
      
      predicted_classes <- ifelse(predictions > 0.7, 1, 0)
      predicted_classes 
      
      # Store the actual value
      actual_values[i] <- ifelse(test_data[25] == "Control", 0, 1)
    }

#  • Evaluate the performance (e.g., accuracy, ROC curve, etc.) using predictions and actual values
#    For example, calculate accuracy
      correct_predictions <- sum(predicted_classes == actual_values)
      accuracy <- correct_predictions / n

#  • Output the accuracy
      print(paste("Accuracy:", accuracy)) 
      # 0.53 
      
#  • Retrain on the reduced_dataset:
      model <- locfit((as.numeric(as.factor(Category))- 1) ~ DPI_m+DUS_m, alpha = 0.2, family = "binomial", data = reduced_dataset)
      
#  • Predicting probabilities for new data (test set):
      predictions <- predict(model, newdata = RBD, type = "response")
      predictions

#  • Predicting classes (0 or 1) based on a threshold 
      predicted_classes <- ifelse(predictions > 0.7, 1, 0)
      predicted_classes 
      
#  • check the error (APER):
      prior <- c(0.6515152, 0.3484848)
      G <- 2 # number of the groups
      misc <- table(class.true=cluster.ew, class.assigned=predicted_classes)
       APER <- 0
       for(g in 1:G)
         APER <- APER + sum(misc[g,-g])/sum(misc[g,]) * prior[g] 
        APER # 0.51, nononononononononono
 
        
  
# Splines logistic regression:
# • Training and validation set:
#   - Number of rows in the dataset
    n <- nrow(reduced_dataset)

#   - Create placeholders for predictions and actual values
    predictions <- rep(NA, n)
    actual_values <- rep(NA, n)
    predicted_classes <- rep(NA, n)
    
#   - Perform leave-one-out cross-validation
      for (i in 1:n) {
        # Exclude the i-th observation
        train_data <- reduced_dataset[-i, ]
        test_data <- reduced_dataset[i, ]
        
        spline_DPI_m <- ns(train_data$DPI_m, df = 6, Boundary.knots = range(train_data$DPI_m))  
        spline_DUS_m <- ns(train_data$DUS_m, df = 6, Boundary.knots = range(train_data$DUS_m))  
        model_matrix <- cbind(spline_DPI_m, spline_DUS_m)
        
        # Fit a local logistic regression model on the training data
        fit <- glm(as.factor(Category) ~ model_matrix, family = "binomial",  data = train_data)
        
        # Predict on the left-out data point
        prediction <- predict(fit, newdata = test_data, type = "response")
        predictions[i] <- prediction
        
        predicted_classes[i] <- ifelse(predictions[i] > 0.7, 1, 0)
        # predicted_classes 
        
        # Store the actual value
        actual_values[i] <- ifelse(test_data[25] == "Control", 0, 1)
      } 
    
    
#  - Evaluate the performance (e.g., accuracy, ROC curve, etc.) using predictions and actual values
#     For example, calculate accuracy
    correct_predictions <- sum(predicted_classes == actual_values)
    accuracy <- correct_predictions / n
    
#  - Output the accuracy
    print(paste("Accuracy:", accuracy)) # 0.65
        
#  -  Retrain the model: create natural cubic splines for specific predictors
    spline_DPI_m <- ns(reduced_dataset$DPI_m, df = 6, Boundary.knots = range(reduced_dataset$DPI_m))  
    spline_DUS_m <- ns(reduced_dataset$DUS_m, df = 6, Boundary.knots = range(reduced_dataset$DUS_m))  
    model_matrix <- cbind(spline_DPI_m, spline_DUS_m)
        
# -  Fit logistic regression model with splines
    fit <- glm(as.factor(Category) ~ model_matrix, family = "binomial", data = reduced_dataset)
    summary(fit)
    

#  - Predicting probabilities for new data (test set):
    predictions <- predict(fit, newdata = RBD, type = "response")
    predictions
    
#  - Predicting classes (0 or 1) based on a threshold 
    predicted_classes <- ifelse(predictions > 0.7, 1, 0)
    predicted_classes 
    
#  - check the error (APER):
    prior <- c(0.6515152, 0.3484848)
    G <- 2 # number of the groups
    misc <- table(class.true=cluster.ew, class.assigned=predicted_classes[1:50]) # cluster.ew and predicted_classes do not the same length. Why?
    APER <- 0
    for(g in 1:G)
      APER <- APER + sum(misc[g,-g])/sum(misc[g,]) * prior[g] 
    APER # 0.36
    
    
# GAM:
    # • Training and validation set:
    #   - Number of rows in the dataset
    n <- nrow(reduced_dataset)
    
    #   - Create placeholders for predictions and actual values
    predictions <- rep(NA, n)
    actual_values <- rep(NA, n)
    predicted_classes <- rep(NA, n)
    
    #   - Perform leave-one-out cross-validation
    for (i in 1:n) {
      # Exclude the i-th observation
      train_data <- reduced_dataset[-i, ]
      test_data <- reduced_dataset[i, ]
      
     
      
      # Fit a local logistic regression model on the training data
      model_gam=gam(as.factor(Category) ~ ns(DPI_m,df=8) + ns(DUS_m,df=8) , family = "binomial", data = train_data)  
      
      # Predict on the left-out data point
      prediction <- predict(model_gam, newdata = test_data, type = "response")
      predictions[i] <- prediction
      
      predicted_classes[i] <- ifelse(predictions[i] > 0.7, 1, 0)
      # predicted_classes 
      
      # Store the actual value
      actual_values[i] <- ifelse(test_data[25] == "Control", 0, 1)
    } 
    
    #  - Evaluate the performance (e.g., accuracy, ROC curve, etc.) using predictions and actual values
    #     For example, calculate accuracy
    correct_predictions <- sum(predicted_classes == actual_values)
    accuracy <- correct_predictions / n
    
    #  - Output the accuracy
    print(paste("Accuracy:", accuracy)) # 0.67
    
    
   model_gam=gam(as.factor(Category) ~ ns(DPI_m,df=8) + ns(DUS_m,df=8), family = "binomial", data = reduced_dataset)   
   summary(model_gam)
   predictions <- predict(model_gam,newdata=RBD,  type = "response")
   predicted_classes<- ifelse(predictions > 0.7, 1, 0) 
   predicted_classes
   
   prior <- c(0.6515152, 0.3484848)
   G <- 2 # number of the groups
   misc <- table(class.true=cluster.ew, class.assigned=predicted_classes)# cluster.ew and predicted_classes do not the same length. Why?
   misc
   APER <- 0
   for(g in 1:G)
     APER <- APER + sum(misc[g,-g])/sum(misc[g,]) * prior[g] 
   APER # 0.5323095
   
   







