 library(readr)
 
## Import the dataset and fix it:
 
 dataset <- read.csv("dataset.csv")
 View(dataset)
 
 colnames(dataset) <- dataset[1,]
 dataset <- dataset[-1,]
 rownames(dataset) <- 1:dim(dataset)[1]
 
 sani <- which(startsWith(dataset$`Participant code`, "HC"))
 malati <- which(startsWith(dataset$`Participant code`, "PD"))
 rischio <- which(startsWith(dataset$`Participant code`, "RBD"))
 
 dataset <- dataset[,-1]
 dataset <- type.convert(dataset, as.is=T)
 dataset <- dataset[,-c(7,8,10)]
 
 
 
 dataset$Category <- as.factor(c( rep("Ill", 30) , rep("At Risk", 50), rep("Healthy",50) ))
 
 
 for (i in 38:61) 
   boxplot(dataset[,i] ~ dataset$Category , main=colnames(dataset)[i]) # How to improve it?
 
 for (i in 38:61) 
   boxplot(scale(x=dataset[,i], center=T, scale=T) ~ dataset$Category , main=colnames(dataset)[i])

## Exploratory analysis for the numerical variables (PCA):
 
#   • Let's go:
       pc.dataset <- princomp(dataset[,c(38:61)], scores=T)
       pc.dataset
       summary(pc.dataset)
       
#   • loadings
       load.dataset <- pc.dataset$loadings
       load.dataset
       
#   • summary plot:
       x11()
       layout(matrix(c(2,3,1,3),2,byrow=T))
       plot(pc.dataset, las=2, main='Principal components', ylim=c(0,4.5e4))
       barplot(sapply(dataset[,c(1,8,38:61)],sd)^2, las=2, main='Original Variables', ylim=c(0,4.5e4), ylab='Variances')
       plot(cumsum(pc.dataset$sd^2)/sum(pc.dataset$sd^2), type='b', axes=F, xlab='number of components', 
            ylab='contribution to the total variance', ylim=c(0,1), main = 'Componet vs Cumulative variance')
       abline(h=1, col='blue')
       abline(h=0.8 , lty=2, col='blue') 
       box()
       axis(2,at=0:10/10,labels=0:10/10)
       axis(1,at=1:ncol(dataset[,c(1,8,38:61)]),labels=1:ncol(dataset[,c(1,8,38:61)]),las=2)
       
#      It seems that keeping 4 PCs is enough to explain over the 90% of the total variability.
#      Looking the matrix of the loadings i can see that some of the features can be deleted.
#      The strategy if to look the correlation with the other features and then proceed.
#      Now find some other information looking the scores.       
       
#     •  scores:
       scores.dataset <- pc.dataset$scores
       scores.dataset
       
       x11()
       boxplot(dataset[,c(1,8,38:61)], las=2, col='gold', main='Original variables')
       
       x11()
       boxplot(scores.dataset, las=2, col='gold', main='Principal components')
       
       x11()
       plot(scores.dataset[,1],scores.dataset[,2],type="n",xlab="pc1",ylab="pc2", asp=1)
       text(scores.dataset[,1],scores.dataset[,2],dimnames(dataset[,c(1,8,38:61)])[[1]], cex=0.7)
       
# NB: the PCA lacks of interpretability, look the correlation matrix :
       cormat <- cor(dataset[,c(38:61)])
       library(corrplot)
       x11()
       corrplot(cormat, method = "number", type = "upper", add = F, tl.col = "black", tl.pos = "n",is.corr=T,diag=T) 
       
# the chosen numerical variables:
      pseudo_numerical <- dataset[,-c(38,40,43,45,46,50,52,55,57,58)]
      numerical <- pseudo_numerical[,c(38:51)]
#     numerical <- numerical[,-c(7,14)]

## I took the categorical variable preproccessing part by Simone Giacomello's script:
    
       categorical <- dataset[1:80,1:37]
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
       cbind(sum.tremor, col[1:80]) # here is the problem to create the dataframe
       
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
       
       categorical.hc <- dataset[81:130,1:12 ] 
       categorical.hc$Total.Tremor <- categorical.hc$Total.Rigidity <- categorical.hc$Total.Spasm <- 
         categorical.hc$Total.Stability <- "-" 
       
       categorical.hc <- categorical.hc[,-c(3,4)]
       colnames(categorical.hc[,11:14]) <- c("Total.Tremor", "Total.Rigidity", "Total.Spasm", "Total.Stability")
       
       
       categorical.new <- categorical.new[,-c(12,13,14,16,17,19,20,22,23,24)]
       colnames(categorical.new) <- colnames(categorical.hc)
       
       categorical.new <- rbind(categorical.new, categorical.hc)
 
      
## here the preprocessed dataset:
dataset.new <- as.data.frame( cbind(categorical.new, numerical))
View(dataset.new)

dataset.new <- dataset.new[,-6] # drop the column clonazepam (mg/day) because is different from 0 in and only if
                                # benzodiazepine medication is Revotril.
write.csv(dataset.new, file="dataset_preprocessed_1.csv")



