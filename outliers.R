library(DepthProc)
library(aplpack)

parkinson <- read.csv("data/parkinson.csv")

colnames(parkinson) <- parkinson[1,]
parkinson <- parkinson[-1,]
rownames(parkinson) <- 1:dim(parkinson)[1]

parkinson <- parkinson[,-1]

parkinson <- type.convert(parkinson, as.is=T)

parkinson$Category <- as.factor(c( rep("PD", 30) , rep("RBD", 50), rep("Control",50) ))

#dataset with only speech variables:
speech_dataset <- parkinson[,c(41:65)]

#change the colnames of the speech variables with their acronyms
colnames(speech_dataset)
colnames(speech_dataset) <- c('EST_r','RST_r','AST_r','DPI_r','DVI_r','GVI_r','DUS_r','DUF_r','RLR_r',
                             'PIR_r','RSR_r','LRE_r','EST_m','RST_m','AST_m','DPI_m','DVI_m','GVI_m','DUS_m',
                             'DUF_m','RLR_m','PIR_m','RSR_m','LRE_m','Category')
colnames(speech_dataset)

#save the rbd patients into a separate table
rbd.speech_dataset <- speech_dataset[which(speech_dataset$Category == "RBD"),]

norbd.speech_dataset <- speech_dataset[-which(speech_dataset$Category == "RBD"),]

####OUTLIERS DETECTION

OUTLIER.INDICES <- numeric(0)
p <- dim(norbd.speech_dataset[,-25])[2]

for (i in 1:p) {
  
  for (j in (i+1):p) {
    
    df.comb <- norbd.speech_dataset[,c(i,j)]
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
outliers <- head(counts.outliers, 14)
MOST.OUTLYING <- as.numeric(names(outliers))
MOST.OUTLYING

n <- dim(rbd.speech_dataset)[1]
col <- rep("blue", n)
col[MOST.OUTLYING] <- "red"
pairs(norbd.speech_dataset[,-25], col=col)

#no outliers dataset
norbd.speech_dataset.clean <- norbd.speech_dataset[-MOST.OUTLYING,]

speech_dataset.clean <- rbind(norbd.speech_dataset.clean, rbd.speech_dataset)

write.csv(speech_dataset.clean, file="data/speech_no_outliers.csv")








