
seed <- 20012024

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
  feature_data <- filtered_data[, feature]
  group_labels <- filtered_data$Category
  
  obs_diff <- median(feature_data[group_labels == "PD"]) - median(feature_data[group_labels == "Control"])
  
  perm_diffs <- replicate(n_permutations, {
    shuffled_labels <- sample(group_labels)
    median(feature_data[shuffled_labels == "PD"]) - median(feature_data[shuffled_labels == "Control"])
  })
  
  p_value <- mean(abs(perm_diffs) >= abs(obs_diff))
  
  if (p_value < 0.05) {
    significant_p_values[feature] <- p_value
  }
}

print(significant_p_values)

# Among the non-normally distributed speech variables, "RST_r" "DPI_r" "DUS_r" "DPI_m" 
# have significant p-values



### Boxplots of significant features (HC vs PD)

library(ggplot2)
data <- speech_dataset[speech_dataset$Category != "RBD",]
features <- c("RST_r", "DPI_r", "DUS_r", "DPI_m", "RST_m")

for (feature in features) {
  plot <- ggplot(data, aes(x=Category, y=data[,feature], fill=Category)) +
    scale_fill_manual(values=c("lightskyblue","lightcyan4")) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = paste("Box Plot of", feature), x = "Category", y = feature)
  
  print(plot)
}



### Modeling step: semiparametric model

library(mgcv)
library(splines)


## Building a model with: Gender RST_r DPI_r DUS_r RST_m DPI_m 

df_model_g <- filtered_data[,c(25,2,4,7,14,16)]
df_model_g$Category <- ifelse(filtered_data$Category == 'PD',1,0)

df_model_g <- data.frame(Gender=parkinson[c(pd, control),2], df_model_g)
df_model_g$Gender <- ifelse(df_model_g$Gender == 'F',1,0)

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
updrs_RBD <- ifelse(updrsIII.new > 3, 1, 0)

table(true_label=updrs_RBD, prediction=model_RBD)


cols <- ifelse(updrs_RBD == 1, 'red', 'black')
plot(predictions$fit, col=cols, pch=16)
abline(0.5, 0, col='red', lty=2, lwd=2)

#           prediction
# true_label   0  1
#          0  18  9
#          1  11 12



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

table(true_label=updrs_RBD, prediction=model_RBD)


