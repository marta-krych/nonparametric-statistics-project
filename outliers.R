library(aplpack)

speech <- read.csv("data/speech.csv")

speech$Category <- factor(speech$Category)
speech$Gender <- factor(speech$Gender)

rbd <- speech[which(speech$Category == "RBD"),]
norbd <- speech[-which(speech$Category == "RBD"),-c(25:27)]


norbd.full <- speech[-which(speech$Category == "RBD"),]

####OUTLIERS DETECTION
#I look for outliers only in the categories PD and Control since these patients 
#will be used to build a regression model to assess the risk for RBD patients

OUTLIER.INDICES <- numeric(0)
p <- dim(norbd)[2]

for (i in 1:p) {
  
  for (j in i:p) {
    
    df.comb <- norbd[,c(i,j)]
    outliers <- bagplot(df.comb, create.plot = F)$pxy.outlier
    indices <- which(apply(df.comb,1,function(x) all(x %in% outliers)))
    
    OUTLIER.INDICES <- c(OUTLIER.INDICES, indices)
    
  }
}

unique(OUTLIER.INDICES)   #41 values, too many!!!
OUTLIER.INDICES

#number of repetitions of each index 
counts.outliers <- table(OUTLIER.INDICES)
#sort them in a descending order to identify which are outlying in many variables
counts.outliers <- sort(counts.outliers,decreasing = T)
counts.outliers

#I define outliers as the units that are outliers in more than 20 bivariate bagplots
outliers <- head(counts.outliers, 14)
MOST.OUTLYING <- as.numeric(names(outliers))
MOST.OUTLYING

#4  6  3 41 72 50 52 28 75  8 23 59 74 11

n <- dim(speech)[1]
col <- rep("blue", n)
col[MOST.OUTLYING] <- "red"
pairs(norbd[,-c(25:27)], col=col)

#no outliers dataset
no.rbd.speech.clean <- norbd.full[-MOST.OUTLYING,]

speech.clean <- rbind(no.rbd.speech.clean, rbd)
ids <- as.numeric(rownames(speech.clean))

speech.clean <- speech.clean[order(ids),]

write.csv(speech.clean, file="data/speech_no_outliers.csv", row.names = FALSE)







