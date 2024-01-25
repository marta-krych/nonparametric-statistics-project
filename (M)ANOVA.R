### Quick pre-processing part

parkinson <- read.csv(file="data/parkinson_preprocessed.csv")
colnames(parkinson)[1:35] <- c("Age (years)", "Gender", "Duration of disease from first symptoms (years)",
                               "Antidepressant therapy", "Benzodiazepine medication",
                               "Clonazepam (mg/day)", "Hoehn-Yahr scale", "UPDRS III total",
                               "Speech", "Facial Expression", "Tremor at Rest - Head",
                               "Tremor at Rest - RUE", "Tremor at Rest - LUE", "Tremor at Rest - RLE",
                               "Tremor at Rest - LLE", "Action or Postural Tremor - RUE",
                               "Action or Postural Tremor - LUE", "Rigidity - Neck", "Rigidity - RUE",
                               "Rigidity - LUE", "Rigidity - RLE", "Rigidity - LLE", "Finger Taps - RUE",
                               "Finger Taps - LUE", "Hand Movements - RUE", "Hand Movements - LUE",
                               "Rapid Alternating Movements - RUE", "Rapid Alternating Movements - LUE",
                               "Leg Agility - RLE", "Leg Agility - LLE", "Arising from Chair", "Posture",
                               "Gait", "Postural Stability", "Body Bradykinesia and Hypokinesia")

parkinson[,c(3,7)] <- sapply(parkinson[,c(3,7)], as.numeric)
parkinson[,8:35] <- sapply(parkinson[,8:35], as.integer)

parkinson$Category <- as.factor(parkinson$Category)
parkinson$Gender <- as.factor(parkinson$Gender)
parkinson$`Antidepressant therapy` <- as.factor(parkinson$`Antidepressant therapy`)
parkinson$`Benzodiazepine medication` <- as.factor(parkinson$`Benzodiazepine medication`)



####  SPEECH  ####

# Analysis to identify if there are speech variables for which there is a 
# significant difference in the mean of the groups (PD, RBD and Control)

B <- 1e3
seed <- 23112023
set.seed(seed)

speech <- parkinson[,36:60]



### One-way MANOVA

fit <- manova(as.matrix(speech[,1:24]) ~ speech$Category)
summary.manova(fit, test="Wilks") 

T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
T_stat <- numeric(B)

set.seed(seed)

for(perm in 1:B){
  permutation <- sample(1:n)
  group.perm <- speech$Category[permutation]
  fit.perm <- manova(as.matrix(speech[,2:15]) ~ group.perm)
  T_stat[perm] <- -summary.manova(fit.perm, test="Wilks")$stats[1,2]
}

p_val <- sum(T_stat>=T0)/B
p_val  # 0

# Indeed, there are some variables in which there is significant difference among
# the groups (PD, RBD and Control)
# To check specifically which are those variables, we will use a one-way permutation ANOVA



### One-way ANOVA

g <- nlevels(speech$Category)
n <- dim(speech)[1]

speech$Category <- factor(speech$Category, levels=c('Control','RBD','PD'))

idx <- 1   # <---- put here the right index (between 1 and 24)
plot(speech$Category, speech[,idx], xlab='group', 
     col=c("lightskyblue","lightsteelblue","lightcyan4"), main=colnames(speech)[idx]) # group-wise boxplot


## Parametric approach:
# (we first check whether the assumption of parametric ANOVA hold, if that's
#  not the case, we will use a nonparametric approach instead)

# 1) univariate normality in each group
Ps <- c(shapiro.test(speech[speech$Category == 'Control', idx])$p,
        shapiro.test(speech[speech$Category == 'RBD', idx])$p,
        shapiro.test(speech[speech$Category == 'PD', idx])$p)
Ps

# 2) Bartlett's test for homogeneity of variances
bartlett.test(speech[,idx], speech$Category)

fit <- aov(speech[,idx] ~ speech$Category)
summary(fit)


## Nonparametric approach:

T0 <- summary(fit)[[1]][1,4]
T_stat <- numeric(B) 

set.seed(seed)

for(perm in 1:B){
  permutation <- sample(1:n)
  variable_perm <- speech[permutation,idx]
  fit_perm <- aov(variable_perm ~ speech$Category)
  
  T_stat[perm] <- summary(fit_perm)[[1]][1,4]
}

p_val <- sum(T_stat>=T0)/B
p_val

# significant idx: 2 (0.005), 4 (0.001), 5 (0.039), 12 (0.06), 14 (0), 16 (0.002), 
#                  17 (0.07), 19 (0.001), 22 (0.016), 23 (0.04)

colnames(speech)[c(2,4,5,12,14,16,17,19,22,23)]
# "RST_r" "DPI_r" "DVI_r" "LRE_r" "RST_m" "DPI_m" "DVI_m" "DUS_m" "PIR_m" "RSR_m"





