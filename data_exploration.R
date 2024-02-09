
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

parkinsonPR <- parkinson[parkinson$Category != 'Control', ] # only RBD and PD patients




### Tables with age and gender 

# 1-30: PD -- 31-80: RBD

age_classes <- cut(parkinson$`Age (years)`, seq(30,90,10))
table(parkinson$Category, age_classes)

table(parkinson$Category, parkinson$Gender)



### Distribution of duration of disease for PD and RBD

library(ggplot2)
# library(devtools)
# devtools::install_github("kassambara/easyGgplot2")
library(easyGgplot2)

duration_disease <- data.frame(Category=parkinsonPR$Category,
                               Duration=parkinsonPR$`Duration of disease from first symptoms (years)`)

p <- ggplot2.histogram(data=duration_disease, xName='Duration',
                       groupName='Category', legendPosition="right",
                       faceting=TRUE, facetingVarNames="Category", 
                       binwidth=1, groupColors=c("lightcyan4","lightsteelblue"))

p <- ggplot2.customize(p, mainTitle="Duration of the disease from first symptoms",
                       xtitle="years", xtitleFont=c(12,"plain","black"), xTickLabelFont=c(10,"plain","black"), 
                       ytitle="count", ytitleFont=c(12,"plain","black"), yTickLabelFont=c(10,"plain","black"),
                       legendPosition="right")
plot(p)




### UPDRS III Total (for PD and RBD patients)

score_classes <- cut(parkinsonPR$`UPDRS III total`, c(-Inf,0,32,59,132))
table(parkinsonPR$Category, score_classes)

updrs_total <- data.frame(Category=parkinsonPR$Category,
                          Total=parkinsonPR[,8])

p <- ggplot2.histogram(data=updrs_total, xName='Total',
                       groupName='Category', legendPosition="right",
                       faceting=TRUE, facetingVarNames="Category", 
                       binwidth=1, groupColors=c("lightcyan4","lightsteelblue")) + 
  geom_vline(aes(xintercept=32.5), color="lightsteelblue4", linetype="dashed") + 
  geom_vline(aes(xintercept=58.5), color="lightsteelblue4", linetype="dashed") +
  geom_vline(aes(xintercept=0.5), color="lightsteelblue4", linetype="dashed")

p <- ggplot2.customize(p, xtitle="UPDRS III Total score", xtitleFont=c(12,"plain","black"), 
                       xTickLabelFont=c(10,"plain","black"), ytitle="count", ytitleFont=c(12,"plain","black"), 
                       yTickLabelFont=c(10,"plain","black"), legendPosition="right") + 
  annotate("text", x=14, y=6.6, label= "mild", size=4) + annotate("text", x=45, y=6.6, label= "moderate", size=4)

plot(p)




### UPDRS III Speech (for PD and RBD patients)

speech_updrs <- data.frame(Category=parkinsonPR$Category,
                           Speech=parkinsonPR$Speech)

p1 <- ggplot2.barplot(data=speech_updrs, xName='Speech',
                      groupName='Category', legendPosition="right",
                      faceting=TRUE, facetingVarNames="Category",
                      groupColors=c("lightcyan4","lightsteelblue"),
                      ylim=c(0,55)) + geom_text(aes(label = ..count..), stat="count", vjust=-0.8, colour="grey40")

p1 <- ggplot2.customize(p1, xtitle="UPDRS III Speech score", xtitleFont=c(12,"plain","black"), 
                        xTickLabelFont=c(10,"plain","black"), ytitle="count", ytitleFont=c(12,"plain","black"), 
                        yTickLabelFont=c(10,"plain","black"), legendPosition="right")
plot(p1)




### Hoehn-Yahr score (for PD patients)

hoehn_total <- data.frame(Total=parkinsonPR[1:30,7])

p1 <- ggplot2.barplot(data=hoehn_total, xName='Total', fill="lightcyan4") +
  geom_text(aes(label = ..count..), stat="count", vjust=-0.8, colour="grey40")

p1 <- ggplot2.customize(p1, xtitle="Hoehn-Yahr score", xtitleFont=c(12,"plain","black"), xTickLabelFont=c(10,"plain","black"), 
                       ytitle="count", ytitleFont=c(12,"plain","black"), yTickLabelFont=c(10,"plain","black"), ylim=c(0,20))
plot(p1)




