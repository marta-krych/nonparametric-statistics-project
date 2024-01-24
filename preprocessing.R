#setwd("C:/Users/simon/Desktop/universita/MAGISTRALE/nonparametric statistics/project")

parkinson <- read.csv("data/parkinson.csv")

colnames(parkinson) <- parkinson[1,]
parkinson <- parkinson[-1,]
rownames(parkinson) <- 1:dim(parkinson)[1]

parkinson <- parkinson[,-1]  #remove the ID column 

parkinson <- type.convert(parkinson, as.is=T)

parkinson <- parkinson[,-c(7,8,10)]  #all units are equal in these columns

parkinson$Category <- as.factor(c( rep("PD", 30) , rep("RBD", 50), rep("Control",50) ))
col <- c(rep("red",30), rep("blue",50), rep("green",50))


##################### categorical variables

categorical <- parkinson[1:80,1:37] #these variables are NA in HC patients

attach(categorical)

#only 3 Yes values, so this variable is useless
which(`Positive history of Parkinson disease in family` == "Yes")  

#Only 10 out of 80 patients have Antidepressant medication
which(`Antidepressant therapy` != "No")   
#9 of these 10 are different medications, so we only keep Yes or No as a variable
`Antidepressant therapy`[which(`Antidepressant therapy` != "No")]

`Benzodiazepine medication`[which(`Benzodiazepine medication` != "No")]

parkinson$Category[which(`Benzodiazepine medication` != "No")]

# if Benzodiazepine medication == Yes (Rivotril), then Clonazepam (mg/day)>0
# because Rivotril is in the Clonazepam family of medication
which(`Clonazepam (mg/day)`>0)
which(`Benzodiazepine medication` == "Yes (Rivotril)")


`Age (years)` - as.numeric(`Age of disease onset (years)`)
`Duration of disease from first symptoms (years)`
# remove Age of disease onset (years)

detach(categorical)

categorical <- categorical[,-c(3,4)]

categorical[which(categorical[,4] != "No"),4] = "Yes"

categorical[which(categorical[,5] != "No" & categorical[,5] != "Yes (Rivotril)" ), 5] = "Other"
categorical[which(categorical[,5] == "Yes (Rivotril)"),5] = "Rivotril"

categorical <- rbind(categorical, parkinson[81:130, colnames(categorical)])

#####################  speech variables

speech <- data.frame( sapply(parkinson[,38:61], as.numeric) )
colnames(speech) <- c("EST_r","RST_r","AST_r","DPI_r","DVI_r","GVI_r","DUS_r", 
                      "DUF_r","RLR_r","PIR_r","RSR_r" ,"LRE_r" ,"EST_m","RST_m",
                      "AST_m","DPI_m","DVI_m","GVI_m","DUS_m","DUF_m","RLR_m",
                      "PIR_m", "RSR_m" , "LRE_m" )
  
speech$Age <- parkinson$`Age (years)`
speech$Gender <- parkinson$Gender
speech$Category <- parkinson$Category

#speech dataset
write.csv(speech, file="data/speech.csv", row.names = FALSE)


#preprocessed dataset
parkinson.clean <- as.data.frame( cbind(categorical, speech))
parkinson.clean$Category <- parkinson$Category

write.csv(parkinson.clean, file="data/parkinson_preprocessed.csv", row.names = FALSE)




