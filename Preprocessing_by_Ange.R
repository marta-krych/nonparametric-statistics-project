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

## Exploratory analysis:
#  • PCA
#    - take the numerical values:
       dataset_num <- dataset[,-c(2,3,6,7)]
       
#    - because with numerical value the command do not work i take the continous features. Let's go:
       pc.dataset <- princomp(dataset[,c(38:61)], scores=T)
       pc.dataset
       summary(pc.dataset)
       
#    - loadings
       load.dataset <- pc.dataset$loadings
       load.dataset
       
#    - summary plot:
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
       
#     - scores:
       scores.dataset <- pc.dataset$scores
       scores.dataset
       
       x11()
       boxplot(dataset[,c(1,8,38:61)], las=2, col='gold', main='Original variables')
       
       x11()
       boxplot(scores.dataset, las=2, col='gold', main='Principal components')
       
       x11()
       plot(scores.dataset[,1],scores.dataset[,2],type="n",xlab="pc1",ylab="pc2", asp=1)
       text(scores.dataset[,1],scores.dataset[,2],dimnames(dataset[,c(1,8,38:61)])[[1]], cex=0.7)
       
#       • look the correlation matrix:
       cormat <- cor(dataset[,c(38:61)])
       
       library(corrplot)
       x11()
       corrplot(cormat, method = "number", type = "upper", add = F, tl.col = "black", tl.pos = "n",is.corr=T,diag=T)

       
## here the obtained dataset:
       dataset <- dataset[,-c(38,40,43,45,46,50,52,55,57,58)]
       View(dataset)
       