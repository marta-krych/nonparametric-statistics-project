
### Quick pre-processing part

parkinson <- read.csv(file="data/parkinson_preprocessed.csv")
parkinson <- parkinson[,-c(1)]
colnames(parkinson) <- c("Age", "Gender", "Duration of disease from first symptoms (years)",
                         "Antidepressant therapy", "Benzodiazepine medication",
                         "Clonazepam (mg/day)", "Hoehn-Yahr scale", "UPDRS III total",
                         "Speech", "Facial Expression", "Total tremor", 
                         "Total rigidity", "Total spasm", "Total stability",
                         "Entropy of speech timing (RP)", "Rate of speech timing (RP)",
                         "Acceleration of speech timing (RP)", "Gaping in between voiced intervals (RP)",
                         "Duration of unvoiced stops (RP)", "Decay of unvoiced fricatives (RP)",
                         "Relative loudness of respiration (RP)", "Pause intervals per respiration (RP)",
                         "Latency of respiratory exchange (RP)", "Entropy of speech timing (M)",
                         "Rate of speech timing (M)", "Acceleration of speech timing (M)", 
                         "Decay of unvoiced fricatives (M)", "Latency of respiratory exchange (M)")
# RP: reading passage
# M: monologue

cols.int <- c(3,7:14)
parkinson[cols.int] <- sapply(parkinson[cols.int], as.integer)

# 1-30: PD (2) // 31-80: RBD (1) // 81-130: HEALTHY (0)

category <- c(rep(2,30), rep(1,50), rep(0,50))
parkinson <- cbind(category, parkinson)
colnames(parkinson)[1] <- "Group"

parkinson$Gender <- as.factor(parkinson$Gender)
parkinson$Group <- as.factor(parkinson$Group)



####  SPEECH  ####

B <- 1e3
seed <- 23112023
set.seed(seed)

speech <- parkinson[,c(1,16:29)]


### One-way ANOVA

g <- nlevels(speech$Group)
n <- dim(speech)[1]

idx <- 12   # <---- put here the right index (between 2 and 15)
plot(speech$Group, speech[,idx], xlab='group', col=heat.colors(g), main=colnames(speech)[idx]) # group-wise boxplot


## Parametric approach:
# (we first check whether the assumption of parametric ANOVA hold, if that's
#  not the case, we will use a nonparametric approach instead)

# 1) univariate normality in each group
Ps <- c(shapiro.test(speech[speech$Group == 0, idx])$p,
        shapiro.test(speech[speech$Group == 1, idx])$p,
        shapiro.test(speech[speech$Group == 2, idx])$p)
Ps

# 2) Bartlett's test for homogeneity of variances
bartlett.test(speech[,idx], speech$Group)

fit <- aov(speech[,idx] ~ speech$Group)
summary(fit)


## Nonparametric approach:

T0 <- summary(fit)[[1]][1,4]
T_stat <- numeric(B) 

for(perm in 1:B){
  permutation <- sample(1:n)
  variable_perm <- speech[permutation,idx]
  fit_perm <- aov(variable_perm ~ speech$Group)
  
  T_stat[perm] <- summary(fit_perm)[[1]][1,4]
}

p_val <- sum(T_stat>=T0)/B
p_val

# significant idx: 3, 10, 12
colnames(speech)[c(3,10,12)]



### One-way MANOVA

fit <- manova(as.matrix(speech[,2:15]) ~ speech$Group)
summary.manova(fit, test="Wilks") 

T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
T_stat <- numeric(B)

for(perm in 1:B){
  permutation <- sample(1:n)
  group.perm <- speech$Group[permutation]
  fit.perm <- manova(as.matrix(speech[,2:15]) ~ group.perm)
  T_stat[perm] <- -summary.manova(fit.perm,test="Wilks")$stats[1,2]
}

p_val <- sum(T_stat>=T0)/B
p_val  # ---> SIGNIFICANT!



### BOOTSTRAP CIs 
# Bootstrap confidence intervals for the difference of the mean of same speech variable in different groups 
# ==> mainly to compare RBD and PD groups (on indexes found using ANOVA: 3,10,12)

idx <- 3   # index of speech variable (between 2 and 15)

x1 <- speech[speech$Group == 2, idx]
x2 <- speech[speech$Group == 1, idx]

# Plot data
par(mfrow=c(1,2))
boxplot(x1, ylim=range(c(x1,x2)), main='PD', col='lightblue4')
boxplot(x2, ylim=range(c(x1,x2)), main='RBD', col='lightblue')

# Computing the bootstrap distribution of the difference of the two sample means

x1.obs <- x1
x2.obs <- x2
diff.obs <- mean(x1) - mean(x2)

T.boot.diff <- numeric(B)

for(b in 1:B)
{
  x1.b <- sample(x1.obs, replace = T)
  x2.b <- sample(x2.obs, replace = T)
  T.boot.diff[b] <- mean(x1.b) - mean(x2.b)
}

plot(ecdf(T.boot.diff), main='Sample Mean PD - Sample Mean RBD')
abline(v = diff.obs, lty=2)

# RP intervals

alpha <- 0.05

right.quantile <- quantile(T.boot.diff, 1 - alpha/2)
left.quantile  <- quantile(T.boot.diff, alpha/2)

CI.RP <- c(diff.obs - (right.quantile - diff.obs), diff.obs - (left.quantile - diff.obs))
CI.RP

abline(v=CI.RP)


# idx: 3
# CI for difference PD-RBD: [-36.048333   6.446667] ==> no significant difference!
# CI for difference PD-healthy: [-56.91400 -11.81667] ==> SIGNIFICANT difference
# CI for difference RBD-healthy: [-37.8235  -3.8515] ==> SIGNIFICANT difference

# idx: 10
# CI for difference PD-RBD: [-7.685333 47.929000] ==> no significant difference!
# CI for difference PD-healthy: [-36.46000  22.99383] ==> no significant difference
# CI for difference RBD-healthy: [-50.5655  -4.9075] ==> SIGNIFICANT difference

# idx: 12
# CI for difference PD-RBD: [-23.81367  22.74767] ==> no significant difference!
# CI for difference PD-healthy: [-64.2615 -14.8435] ==> SIGNIFICANT difference
# CI for difference RBD-healthy: [-55.9875 -18.5590] ==> SIGNIFICANT difference



