
parkinson <- read.csv("data/parkinson.csv")

colnames(parkinson) <- parkinson[1,]
parkinson <- parkinson[-1,]
rownames(parkinson) <- 1:dim(parkinson)[1]

sani <- which(startsWith(parkinson$`Participant code`, "HC"))
malati <- which(startsWith(parkinson$`Participant code`, "PD"))
rischio <- which(startsWith(parkinson$`Participant code`, "RBD"))

parkinson <- parkinson[,-1]

parkinson <- type.convert(parkinson, as.is=T)

parkinson$Category <- as.factor(c( rep("PD", 30) , rep("RBD", 50), rep("Control",50) ))


# Extracting 'UPDRS III Total' column and removing 'Action tremor'
updrsIII <- as.numeric(parkinson[rischio, 13])
updrsIII.new <- updrsIII - (as.numeric(parkinson[rischio, 21]) + as.numeric(parkinson[rischio, 22]))


# Dataset with only speech variables:
speech_dataset <- parkinson[,c(41:65)]
colnames(speech_dataset) <- c('EST_r','RST_r','AST_r','DPI_r','DVI_r','GVI_r','DUS_r','DUF_r','RLR_r',
                             'PIR_r','RSR_r','LRE_r','EST_m','RST_m','AST_m','DPI_m','DVI_m','GVI_m','DUS_m',
                             'DUF_m','RLR_m','PIR_m','RSR_m','LRE_m','Category')


### Boxplots

library(ggplot2)
data <- speech_dataset
features <- names(data)[names(data) != "Category"]

for (feature in features) {
  plot <- ggplot(data, aes_string(x = "Category", y = feature, color = "Category")) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = paste("Box Plot of", feature), x = "Category", y = feature)
  
  print(plot)
}


### Checking for speech variables with Gaussian distribution

filtered_data <- subset(data, Category %in% c("PD", "Control"))
features <- names(filtered_data)[names(filtered_data) != "Category"]

normality_test_results <- list()
normal_variables <- c()

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


### Testing the normally distributed variables for the difference in mean between HC and PD

significancy_test_results <- list()
significant_normal_variables <- c()


for (feature in normal_variables) {
  pd_data <- filtered_data[filtered_data$Category == "PD", feature]
  hc_data <- filtered_data[filtered_data$Category == "Control", feature]
  
  test <- t.test(pd_data, hc_data, alternative="two.sided", mu=0)
  
  normality_test_results[[feature]] <- list(t_test_p_value = test$p.value)
  
  if (test$p.value < 0.05) {
    significant_normal_variables <- c(significant_normal_variables, feature)
  }
}

print(significant_normal_variables)

# Among the normally distributed speech variables, only "RST_m" seems to be significant



### Testing the non-normally distributed variables for the difference in mean between HC and PD
### (using a permutation test)

not_normal_variables <- setdiff(features, normal_variables)

n_permutations <- 1000
significant_p_values <- c()

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

# Among the not-normally distributed speech variables, "RST_r" "DPI_r" "DUS_r" "DPI_m" 
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


### Semiparametric model

library(car)
library(mgcv)
library(splines)

## Building a model with these variables: RST_r DPI_r DUS_r DPI_m RST_m

df_model <- filtered_data[,c(25,2,4,7,14,16)]
df_model$Category <- ifelse(df_model$Category == 'PD',1,0)

with(df_model, scatterplotMatrix(data.frame(Category, RST_r, DPI_r, DUS_r, DPI_m, RST_m)))

# linear relationship: RST_r, RST_m
# non-linear relationship: DPI_r, DPI_m, DUS_r

model_semipar <- gam(Category ~ RST_r + RST_m + s(DPI_r) + s(DPI_m) + s(DUS_r), 
                     family=binomial, data=df_model)
summary(model_semipar)

test <- speech_dataset[rischio, c(2,14,4,16,7)]

predictions <- predict(model_semipar, test, se=TRUE, type="response")
# predictions <- predict(model_semipar, test, se=TRUE)
# pfit <- exp(predictions$fit)/(1 + exp(predictions$fit))


model_RBD <- ifelse(predictions$fit > .5, 1, 0)
updrs_RBD <- ifelse(updrsIII.new > 3, 1, 0)

table(true_label=updrs_RBD, prediction=model_RBD)

cols <- ifelse(updrs_RBD == 1, 'red', 'black')
plot(predictions$fit, col=cols, pch=16)
abline(0.5, 0, col='red', lty=2, lwd=2)
# The threshold itself is okay. It is the model that doesn't perform too well



## Adding 'Gender' to the previous model (Gender RST_r DPI_r DUS_r DPI_m RST_m)

df_model_g <- filtered_data[,c(25,2,4,7,14,16)]
df_model_g$Category <- ifelse(df_model_g$Category == 'PD',1,0)

df_model_g <- data.frame(Gender=parkinson[c(1:30, 81:130),2], df_model_g)
df_model_g$Gender <- ifelse(df_model_g$Gender == 'F',1,0)

with(df_model_g, scatterplotMatrix(data.frame(Category, Gender, RST_r, DPI_r, DUS_r, DPI_m, RST_m)))

# linear relationship: Gender, RST_r, RST_m
# non-linear relationship: DPI_r, DPI_m, DUS_r

model_semipar_g <- gam(Category ~ Gender + RST_r + RST_m + s(DPI_r) + s(DPI_m) + s(DUS_r), 
                       family=binomial, data=df_model_g)
summary(model_semipar_g)

test <- speech_dataset[rischio, c(2,14,4,16,7)]
test <- data.frame(Gender=ifelse(parkinson[31:80, 2] == 'F',1,0),
                   test)

predictions <- predict(model_semipar_g, test, se=TRUE, type="response")
# predictions <- predict(model_semipar_g, test, se=TRUE)
# pfit <- exp(predictions$fit)/(1 + exp(predictions$fit))


model_RBD <- ifelse(predictions$fit > .5, 1, 0)
updrs_RBD <- ifelse(updrsIII.new > 3, 1, 0)

table(true_label=updrs_RBD, prediction=model_RBD)

cols <- ifelse(updrs_RBD == 1, 'red', 'black')
plot(predictions$fit, col=cols, pch=16)
abline(0.5, 0, col='red', lty=2, lwd=2)



### Conformal prediction

library(conformalInference)

train_gam <- function(x, y, out=NULL){
  colnames(x) <- c('Gender','RST_r','RST_m','DPI_r','DPI_m','DUS_r')
  train_data <- data.frame(y, x)
  model_gam <- gam(as.factor(y) ~ Gender + RST_r + RST_m + s(DPI_r) + s(DPI_m) + s(DUS_r),
                   family=binomial, data=train_data)
}

predict_gam <- function(obj, new_x){
  new_x <- data.frame(new_x)
  colnames(new_x) <- c('Gender','RST_r','RST_m','DPI_r','DPI_m','DUS_r')
  predict.gam(obj, new_x, type="response")
}

model_gam <- gam(Category ~ Gender + RST_r + RST_m + s(DPI_r) + s(DPI_m) + s(DUS_r), 
                 family=binomial, data=df_model_g)

target <- df_model_g$Category
covariates <- df_model_g[,c(1,3,6,4,7,5)]

c_preds <- conformal.pred(covariates, target, as.matrix(test), alpha=0.05, verbose=TRUE, 
                          train.fun=train_gam, predict.fun=predict_gam, num.grid.pts=200)

data.frame(lower=c_preds$lo, 
           prediction=c_preds$pred, 
           upper=c_preds$up)


model_RBD <- ifelse(c_preds$pred > .5, 1, 0)
updrs_RBD <- ifelse(updrsIII.new > 3, 1, 0)

table(true_label=updrs_RBD, prediction=model_RBD)

cols <- ifelse(updrs_RBD == 1, 'red', 'black')
plot(c_preds$pred , col=cols, pch=16)
abline(0.5, 0, col='red', lty=2, lwd=2)

