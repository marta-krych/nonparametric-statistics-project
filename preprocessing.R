setwd("C:/Users/simon/Desktop/universita/MAGISTRALE/nonparametric statistics/project")

parkinson <- read.csv("parkinson.csv")

colnames(parkinson) <- parkinson[1,]
parkinson <- parkinson[-1,]
rownames(parkinson) <- 1:dim(parkinson)[1]

sani <- which(startsWith(parkinson$`Participant code`, "HC"))
malati <- which(startsWith(parkinson$`Participant code`, "PD"))
rischio <- which(startsWith(parkinson$`Participant code`, "RBD"))

parkinson <- parkinson[,-1]

parkinson <- type.convert(parkinson, as.is=T)

parkinson <- parkinson[,-c(7,8,10)]

parkinson$Category <- as.factor(c( rep("PD", 30) , rep("RBD", 50), rep("Control",50) ))
col <- c(rep("red",30), rep("blue",50), rep("green",50))


##################### categorical variables

categorical <- parkinson[1:80,1:37]
group <- as.factor( c(rep("malati",30), rep("rischio",50)))

# attach(categorical)
# 
# which(`Positive history of Parkinson disease in family` == "Yes")  #remove column
# 
# which(`Antidepressant therapy` != "No")   #do not distinguish between different medication
# `Antidepressant therapy`[which(`Antidepressant therapy` != "No")]
# 
# `Benzodiazepine medication`[which(`Benzodiazepine medication` != "No")]  
# 
# group[which(`Benzodiazepine medication` != "No")]
# 
# which(`Clonazepam (mg/day)`>0)
# which(`Benzodiazepine medication` == "Yes (Rivotril)")
# 
# # if Benzodiazepine medication == Yes (Rivotril), then Clonazepam (mg/day)>0
# # because Rivotril is in the Clonazepam family of medication,
# # so we redefine 
# 
# `Age (years)` - as.numeric(`Age of disease onset (years)`)
# `Duration of disease from first symptoms (years)`
# # remove Age of disease onset (years)
# 
# detach(categorical)

categorical <- categorical[,-c(3,4)]

categorical[which(categorical[,4] != "No"),4] = "Yes"

categorical[which(categorical[,5] != "No" & categorical[,5] != "Yes (Rivotril)" ), 5] = "Other"
categorical[which(categorical[,5] == "Yes (Rivotril)"),5] = "Rivotril"

categorical.new <- categorical[,1:10]

##### tremor 

tremor <- sapply(categorical[,11:17], as.numeric)
summary(tremor)

sum.tremor <- rowSums(tremor)
categorical.new <- cbind(categorical.new, sum.tremor)
cbind(sum.tremor, col[1:80])

boxplot(sum.tremor ~ group)

fit <- aov(sum.tremor ~ group)
T0 <- summary(fit)[[1]][1,4]

B <- 1000
set.seed(2023)

T_stat <- numeric(B) 
n <- length(sum.tremor)

for(perm in 1:B){
  # Permutation:
  permutation <- sample(1:n)
  sum.tremor_perm <- sum.tremor[permutation]
  fit_perm <- aov(sum.tremor_perm ~ group)
  
  # Test statistic:
  T_stat[perm] <- summary(fit_perm)[[1]][1,4]
}

hist(T_stat,xlim=range(c(T_stat,T0)))
abline(v=T0,col=3,lwd=2)
p.val <- sum(T_stat>=T0)/B
p.val  #the sum of the tremor variables is significant

cor.tremor <- cor(tremor)

pc.tremor <- princomp(tremor, scores=T)
summary(pc.tremor)

plot(cumsum(pc.tremor$sd^2)/sum(pc.tremor$sd^2), type='b', axes=F, xlab='number of components', 
     ylab='contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(tremor),labels=1:ncol(tremor),las=2)

tremor.reduced <- pc.tremor$scores[,1:3]
colnames(tremor.reduced) <- c("Tremor 1", "Tremor 2", "Tremor 3")
pairs(tremor.reduced, col=col)
pairs(tremor, col=col)

categorical.new <- cbind(categorical.new, tremor.reduced)


############## rigidity 

rigidity <- data.frame(sapply(categorical[,18:22], as.numeric))
colnames(rigidity) <- colnames(categorical[,18:22])

cormat.rigidity <- cor( rigidity)

sum.rigidity <- rowSums(rigidity)
categorical.new <- cbind(categorical.new, sum.rigidity)

cbind(sum.rigidity, col[1:80])

boxplot(sum.rigidity ~ group)

fit <- aov(sum.rigidity ~ group)
T0 <- summary(fit)[[1]][1,4]

B <- 1000
set.seed(2023)

T_stat <- numeric(B) 
n <- length(sum.rigidity)

for(perm in 1:B){
  # Permutation:
  permutation <- sample(1:n)
  sum.rigidity_perm <- sum.rigidity[permutation]
  fit_perm <- aov(sum.rigidity_perm ~ group)
  
  # Test statistic:
  T_stat[perm] <- summary(fit_perm)[[1]][1,4]
}

hist(T_stat,xlim=range(c(T_stat,T0)))
abline(v=T0,col=3,lwd=2)
p.val <- sum(T_stat>=T0)/B
p.val  #the sum of the rigidity variables is significant


pc.rigidity <- princomp(rigidity, scores=T)
summary(pc.rigidity) 

plot(cumsum(pc.rigidity$sd^2)/sum(pc.rigidity$sd^2), type='b', axes=F, xlab='number of components', 
     ylab='contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(rigidity),labels=1:ncol(rigidity),las=2)

rigidity.reduced <- pc.rigidity$scores[,1:2]
colnames(rigidity.reduced) <- c("Rigidity 1", "Rigidity 2")
plot(rigidity.reduced, col=col)
#pairs(rigidity, col=col)

categorical.new <- cbind(categorical.new, rigidity.reduced)


############### spasms 

spasm <- data.frame(sapply(categorical[,23:28], as.numeric))
colnames(spasm) <- colnames(categorical[,23:28])

cormat.spasm <- cor( spasm)


sum.spasm <- rowSums(spasm)
categorical.new <- cbind(categorical.new, sum.spasm)

cbind(sum.spasm, col[1:80])

boxplot(sum.spasm ~ group)

fit <- aov(sum.spasm ~ group)
T0 <- summary(fit)[[1]][1,4]

B <- 1000
set.seed(2023)

T_stat <- numeric(B) 
n <- length(sum.spasm)

for(perm in 1:B){
  # Permutation:
  permutation <- sample(1:n)
  sum.spasm_perm <- sum.spasm[permutation]
  fit_perm <- aov(sum.spasm_perm ~ group)
  
  # Test statistic:
  T_stat[perm] <- summary(fit_perm)[[1]][1,4]
}

hist(T_stat,xlim=range(c(T_stat,T0)))
abline(v=T0,col=3,lwd=2)
p.val <- sum(T_stat>=T0)/B
p.val  #the sum of the spasm variables is significant

pc.spasm <- princomp(spasm, scores=T)
summary(pc.spasm) 

plot(cumsum(pc.spasm$sd^2)/sum(pc.spasm$sd^2), type='b', axes=F, xlab='number of components', 
     ylab='contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(spasm),labels=1:ncol(spasm),las=2)

spasm.reduced <- pc.spasm$scores[,1:2]
colnames(spasm.reduced) <- c("Spasm 1", "Spasm 2")
plot(spasm.reduced, col=col)

categorical.new <- cbind(categorical.new, spasm.reduced)

##################### stability

stability <- data.frame(sapply(categorical[,29:35], as.numeric))
colnames(stability) <- colnames(categorical[,29:35])

cormat.stability <- cor( stability)


sum.stability <- rowSums(stability)
categorical.new <- cbind(categorical.new, sum.stability)
colnames(categorical.new) <- c( colnames(categorical.new)[1:10], "Total.Tremor",
                                "Total.Rigidity", "Total.Spasm", "Total.Stability")

cbind(sum.stability, col[1:80])

boxplot(sum.stability ~ group)

fit <- aov(sum.stability ~ group)
T0 <- summary(fit)[[1]][1,4]

B <- 1000
set.seed(2023)

T_stat <- numeric(B) 
n <- length(sum.stability)

for(perm in 1:B){
  # Permutation:
  permutation <- sample(1:n)
  sum.stability_perm <- sum.stability[permutation]
  fit_perm <- aov(sum.stability_perm ~ group)
  
  # Test statistic:
  T_stat[perm] <- summary(fit_perm)[[1]][1,4]
}

hist(T_stat,xlim=range(c(T_stat,T0)))
abline(v=T0,col=3,lwd=2)
p.val <- sum(T_stat>=T0)/B
p.val #the sum of the stability variables is significant

pc.stability <- princomp(stability, scores=T)
summary(pc.stability) 

plot(cumsum(pc.stability$sd^2)/sum(pc.stability$sd^2), type='b', axes=F, xlab='number of components', 
     ylab='contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(stability),labels=1:ncol(stability),las=2)

stability.reduced <- pc.stability$scores[,1:3]
colnames(stability.reduced) <- c("Stability 1", "Stability 2", "Stability 3")
pairs(stability.reduced, col=col)

categorical.new <- cbind(categorical.new, stability.reduced)


######################

categorical.hc <- parkinson[81:130,1:12 ] 
categorical.hc$Total.Tremor <- categorical.hc$Total.Rigidity <- categorical.hc$Total.Spasm <- 
  categorical.hc$Total.Stability <- "-" 
  
categorical.hc <- categorical.hc[,-c(3,4)]

categorical.new <- rbind(categorical.new, categorical.hc)

#####################  numerical variables

numerical <- data.frame( sapply(parkinson[,38:61], as.numeric) )
colnames(numerical) <- colnames(parkinson[,38:61])
cormat <- cor(numerical)

library(corrplot)
corrplot(cormat, method = "number", type = "upper", add = F, tl.col = "black", tl.pos = "n",is.corr=T,
         diag=T)

#pairs(numerical, col=col)

numerical.new <- NULL
names <- NULL
stop <- FALSE

while (!stop) {
  
  highly_correlated <- abs(cormat) > 0.6
  
  if (sum(highly_correlated)>ncol(highly_correlated)) 
    stop <- TRUE
  
  for (i in 1:ncol(highly_correlated)) {
    s <- sum(highly_correlated[1:i,i])
    if (s<=1) {
      names <- c(names, colnames(numerical)[i])
      numerical.new <- cbind(numerical.new, numerical[,i])
    }
  }
  cormat <- cor(numerical.new)
}

numerical.new <- as.data.frame(numerical.new)
colnames(numerical.new) <- names

pairs(numerical.new, col=col)


parkinson.clean <- as.data.frame( cbind(categorical.new, numerical.new))

write.csv(parkinson.clean, file="parkinson_preprocessed.csv")




