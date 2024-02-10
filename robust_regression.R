library(robustbase)
library(robust)
library(psych)
library(MASS)
library(ellipse)
library(here)
library(DescTools)
library(knitr)
library(RobStatTM)

parkinson <- read.csv("data/parkinson.csv")

colnames(parkinson) <- parkinson[1,]
parkinson <- parkinson[-1,]
rownames(parkinson) <- 1:dim(parkinson)[1]

parkinson <- parkinson[,-1]

parkinson <- type.convert(parkinson, as.is=T)

# Extracting 'UPDRS III Total' column and removing 'Action tremor'
updrsIII <- as.numeric(parkinson[31:80, 13])
updrsIII.new <- updrsIII - (as.numeric(parkinson[31:80, 21]) + as.numeric(parkinson[31:80, 22]))


# Dataset with only speech variables:
speech <- read.csv("data/speech.csv")

norbd <- speech[-which(speech$Category == "RBD"), -26]
y <- as.numeric(norbd$Category == "PD")
norbd$Category <- NULL

### robust analysis

p <- dim(norbd)[2]

fit_MCD <- covMcd(x = norbd, alpha = .75, nsamp = "best")

ind_best_subset <- fit_MCD$best
ind_best_subset

ind_rew_obs <- which(
    mahalanobis(
      x = norbd,
      center = fit_MCD$raw.center,
      cov = fit_MCD$raw.cov
    ) <= qchisq(p = .975, df = p)
  )

names(ind_rew_obs) <- NULL

ind_rew_obs == ind_best_subset   #true for every member

cov.mcd <- fit_MCD$cov

sample_mean <- colMeans(norbd)
sample_cov <- cov(norbd)

classical_dist <- sqrt(mahalanobis(x = norbd, center = sample_mean, cov = sample_cov))
robust_dist <- sqrt(fit_MCD$mah)

plot(classical_dist, robust_dist)
sapply(c(1:80)[-ind_best_subset], 
       function(i) text(classical_dist[i], robust_dist[i], labels = i, pos = 3))
abline(a=0, b=1, col='black', lty = 3)

#we can label as outliers those points that don't fall on this line

no.outliers.df <- read.csv("data/speech_no_outliers.csv")

#These row indexes are those we removed in the outliers script:
#3  4  6  8 11 23 28 41 50 52 59 72 74 75

n <- dim(norbd)[1]
(1:n)[-ind_best_subset] #these are the row indexes of the outliers computed by MCD
# 3  4  6  8 11 18 23 28 33 34 41 50 52 72

#there are only a few differences


#########ROBUST REGRESSION

rbd <- speech[speech$Category == "RBD",-27]
rbd$Gender <- as.factor(rbd$Gender)

norbd$Gender <- speech$Gender[speech$Category != "RBD"]

y <- as.numeric(speech$Category[speech$Category != "RBD"] == "PD")


#the most influential variables were 
#c("RST_r", "DPI_r" ,"RST_m", "DPI_m", "DUS_m", "Age", "Gender")



fit1 <- glmrob(y ~  RST_r + DPI_r + RST_m + DPI_m + DUS_m + Age + Gender  , 
                  family = binomial, method="BY", data=norbd)
summary(fit1)$coefficients  #only gender is significant



fit2 <- glmrob(y ~  RST_r + DPI_r + RST_m + DPI_m + DUS_m + Gender  , 
               family = binomial, method="BY", data=norbd)
summary(fit2)$coefficients  #only gender is significant


fit3 <- glmrob(y ~  RST_r + DPI_r + RST_m + DPI_m + Gender  , 
               family = binomial, method="BY", data=norbd)
summary(fit3)$coefficients  #only gender is significant


fit4 <- glmrob(y ~ DPI_r + RST_m + DPI_m + Gender  , 
               family = binomial, method="BY", data=norbd)
summary(fit4)$coefficients  #DPI_r and gender are significant


fit5 <- glmrob(y ~ DPI_r + DPI_m + Gender  , 
               family = binomial, method="Mqle", data=norbd)
summary(fit5)$coefficients #DPI_m is not significant

preds <- fit5$coefficients[1] + fit5$coefficients[2] * rbd$DPI_r + 
  fit5$coefficients[3]*rbd$DPI_m + fit5$coefficients[4]*(rbd$Gender == "M")
preds.fit5 <- exp(preds)/(1 + exp(preds))

table(true_label = as.numeric(updrsIII.new>3), 
      prediction = as.numeric(preds.fit5>0.5))
#             prediction
# true_label   0  1
#           0 20  7
#           1 12 11

############# DA RIVEDERE
##### plot of the decision boundaries of model 5

DPI_r.grid <- seq(range(norbd$DPI_r)[1] - 0.05 * diff(range(norbd$DPI_r)),
                  range(norbd$DPI_r)[2] + 0.05 * diff(range(norbd$DPI_r)),
                  length.out=100)
DPI_m.grid <- seq(range(norbd$DPI_m)[1] - 0.05 * diff(range(norbd$DPI_m)),
                  range(norbd$DPI_m)[2] + 0.05 * diff(range(norbd$DPI_m)),
                  length.out=100)

grid <- expand.grid(DPI_r.grid, DPI_m.grid)


pfitM <- fit5$coefficients[1] + fit5$coefficients[2] * norbd$DPI_r + 
  fit5$coefficients[3]*norbd$DPI_m + fit5$coefficients[4]*(norbd$Gender == "M")

grid$Male <- exp(pfitM)/ (1 + exp(pfitM))

pfitF <- fit5$coefficients[1] + fit5$coefficients[2] * norbd$DPI_r + 
  fit5$coefficients[3]*norbd$DPI_m + fit5$coefficients[4]*(norbd$Gender == "F")

grid$Female <- exp(pfitF)/ (1 + exp(pfitF))



par(mfrow=c(1,2))
with(norbd, 
     plot(DPI_m, DPI_r,
          col = ifelse(Gender == "M", "grey", 0), pch=16, cex=0.7))
contour(DPI_m.grid, DPI_r.grid, 
        matrix(grid$Male, ncol = length(DPI_r.grid), nrow = length(DPI_m.grid)),
        levels = 0.5, add = TRUE, col = "blue", lwd = 2)

with(norbd, 
     plot(DPI_m, DPI_r, ylim=c(200, 400), xlim=c(90, 250),
          col = ifelse(Gender == "F", "grey", 0), pch=16, cex=0.7))
contour(DPI_m.grid, DPI_r.grid, 
        matrix(grid$Female, ncol = length(DPI_r.grid), nrow = length(DPI_m.grid)),
        levels = 0.5, add = TRUE, col = "pink", lwd = 2)

#####


fit6 <- glmrob(y ~ DPI_r + Gender  , 
               family = binomial, method="BY", data=norbd)
summary(fit6)$coefficients #all terms are significant

preds <- fit6$coefficients[1] + fit6$coefficients[2] * rbd$DPI_r + 
  fit6$coefficients[3]*(rbd$Gender == "M")
preds.fit6 <- exp(preds)/(1 + exp(preds))

table(true_label = as.numeric(updrsIII.new>3), 
      prediction = as.numeric(preds.fit6>0.5))
#             prediction
# true_label   0  1
#           0 21  6
#           1 15  8



## Overall the robust regression doesn't perform well as the non robust one,
#this is probably due to the fact that the GAMs or GLMs have been fitted on
#the dataset without the outliers. 


