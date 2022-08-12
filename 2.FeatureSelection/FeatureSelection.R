setwd("D:/Final/FeatureSelection/")

library(GenomicRanges)
library(genomation)
library(randomForest)
library(ggplot2)
library(ggrepel)

load("../DataPreparation/DataSets_AllFeature.Rdata")

TrainSet <- rbind.data.frame(as.data.frame(mcols(Train_True.bed)),
                             as.data.frame(mcols(Train_False.bed)))

TestSet<- rbind.data.frame(as.data.frame(mcols(Test_True.bed)),
                           as.data.frame(mcols(Test_False.bed)))

PredictionSet <- rbind.data.frame(as.data.frame(mcols(chr22.TrueEhaner.bed)),
                                  as.data.frame(mcols(chr22.FalseEnhancer.bed)))

# 1. train Random Foreset model#######################################
# factor for classification
TrainSet$is_Enh <- factor(TrainSet$is_Enh, levels = c(0,1))

# Using random forest for variable selection
rfModel <- randomForest(is_Enh ~ ., data = TrainSet, ntree=500)


# 2. Select feature by importance (MeanDecreaseGini)####################
# Getting the list of important variables
FeatureImportance <- rfModel$importance[,1]
FeatureImportance <- FeatureImportance[order(FeatureImportance,decreasing = T)]

quantile(FeatureImportance)
hist(FeatureImportance,breaks = 50)
importance.cutoff <- 100
abline(v=importance.cutoff, col = "red")

SelectedFeature <- names(FeatureImportance[FeatureImportance > importance.cutoff])
SelectedFeature


pdf("FeatureSelection.pdf",10,10)
plot.data <- as.data.frame(FeatureImportance)
plot.data$name <- rownames(plot.data)

plot.data$selected <- "No"
plot.data[SelectedFeature,]$selected <- "Yes"

plot.data$class <- "Histone"
plot.data[c("P300Std","Hdac2a300705aStd","Hdac1sc6298Std",
            "Ezh239875Std","Hdac6a301341a","Suz12051317","Cbpsc369"),]$class <- "Enzyme"
plot.data["ATAC",]$class <- "Accessibility"
plot.data["CtcfStd",]$class <- "TF"
table(plot.data$class)

plot.data$x_cood <- rnorm(21,mean=1,sd=0.05)
ggplot(plot.data, aes(x=x_cood, y=FeatureImportance))+
  geom_violin( alpha=0.6, fill="aquamarine")+
  geom_point(aes(color = selected),size=4)+
  xlim(c(0.5,1.5))+
  geom_label_repel(aes(label = name, fill = class),
    fontface = 'bold', color = 'white',
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines"),
    segment.color = 'grey50'  )
dev.off()


# 3. output data for classifation########################################
TrainSet <- TrainSet[,c("is_Enh",SelectedFeature)]
TestSet <- TestSet[,c("is_Enh",SelectedFeature)]
PredictionSet <- PredictionSet[,c("is_Enh",SelectedFeature)]

save(TrainSet, TestSet, PredictionSet,SelectedFeature, file = "DataSets_SelectedFeature.Rdata")

dir.create("DataSets_SelectedFeature")
write.table(TrainSet, file = "DataSets_SelectedFeature/TrainSet.tab",
            sep="\t", quote = F, row.names = F)
write.table(TestSet, file = "DataSets_SelectedFeature/TestSet.tab",
            sep="\t", quote = F, row.names = F)
write.table(PredictionSet, file = "DataSets_SelectedFeature/PredictionSet.tab",
            sep="\t", quote = F, row.names = F)
