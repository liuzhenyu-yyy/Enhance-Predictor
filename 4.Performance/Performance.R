setwd("D:/Final/Performance//")

load("../FeatureSelection/DataSets_SelectedFeature.Rdata")
load("../ModelTraining/model1.RData")

library(genomation)
library(GenomicRanges)
library(e1071)
library(stringr)

#1. predict Enhancer on K562 cell line Chr 22##############################
Predict_x <- PredictionSet[,2:13]
Predict_y <- PredictionSet[,1]

pred.model1 <- predict(object = model1, newdata = Predict_x)
Freq.model1 <- table(Predict_y, pred.model1)
Freq.model1
accuracy.model1 <- sum(diag(Freq.model1))/sum(Freq.model1)
accuracy.model1

PredictionSet$prediction <- pred.model1
#2. predict Enhancer for H1 cell############################################
SelectedFeature
SelectedFeature <- SelectedFeature[!SelectedFeature=="ATAC"]
FeatureFileName <- paste("wgEncodeBroadHistoneK562",SelectedFeature,"Sig.bigWig",
                         sep="")
AllFile <- dir("../data/K562_BroadHistone/")
FeatureFileName%in%AllFile

#read bed file
ChromHMM.H1.bed <- readGeneric("../data/ChromHMM.H1.bed", 
                                 header = F,zero.based = T,
                                 meta.cols = list(annotation=4))
#select enhancers
is_enhancer <- str_detect(ChromHMM.H1.bed$annotation, "Strong_Enhancer")
table(is_enhancer)
ChromHMM.Enh.H1.bed <- ChromHMM.H1.bed[is_enhancer]

write.table(as.data.frame(ChromHMM.Enh.H1.bed)[,1:3], 
            file = "ChromHMM.Enh.H1.bed", quote = F,
            col.names = F, sep = "\t", row.names = F)

#False enhancers
FalseEnhancer.H1.bed <- readGeneric("FalseEnhancer.H1.bed", 
                                    header = F,zero.based = T)

#merge bed
ChromHMM.Enh.H1.bed$is_Enh <- 1
FalseEnhancer.H1.bed$is_Enh <- 0

Merge.bed <- c(ChromHMM.Enh.H1.bed, FalseEnhancer.H1.bed)
table(Merge.bed$is.Enh)

#calculate feature
for(i in 1:length(SelectedFeature)){
  temp <- ScoreMatrixBin(target = paste("../data/K562_BroadHistone/",
                                        FeatureFileName[i], sep=''),
                         windows = Merge.bed,
                         bin.num = 1)
  mcols(Merge.bed)[,SelectedFeature[i]] <- temp@.Data
  print(i)
}

#normalize to control(input) data
control.td <- ScoreMatrixBin(target = "../data/K562_BroadHistone/wgEncodeBroadHistoneK562ControlStdSig.bigWig",
                             windows = Merge.bed,
                             bin.num = 1)

for(feature in SelectedFeature){
  mcols(Merge.bed)[,feature] <- mcols(Merge.bed)[,feature]-control.td
}

#Addition feature: chromatin accsibility (ATAC-Seq signal)
#From GEO: GSE70482
ATAC_rep1 <- ScoreMatrixBin(target = "../data/K562_ATAC-Seq/GSM1782764_K562_50K_ATAC-seq_nan_nan_nan_1_1_hg19.bigWig",
                            windows = Merge.bed,
                            bin.num = 1)
ATAC_rep2 <- ScoreMatrixBin(target = "../data/K562_ATAC-Seq/GSM1782765_K562_50K_ATAC-seq_nan_nan_nan_2_1_hg19.bigWig",
                            windows = Merge.bed,
                            bin.num = 1)
#mean value of two repilcates
Merge.bed$ATAC <- (ATAC_rep1@.Data + ATAC_rep2@.Data)/ 2

#prediction
H1.data <- as.data.frame(Merge.bed)[,names(TrainSet)]

Predict.H1_x <- H1.data[,2:13]
Predict.H1_y <- H1.data[,1]

pred.H1.model1 <- predict(object = model1, newdata = Predict.H1_x)
Freq.H1.model1 <- table(Predict.H1_y, pred.H1.model1)
Freq.H1.model1
accuracy.H1.model1 <- sum(diag(Freq.H1.model1))/sum(Freq.H1.model1)
accuracy.H1.model1

H1.data$prediction <- pred.H1.model1
Merge.bed$prediction <- pred.H1.model1
save(Predict.H1_x, Predict.H1_y, Freq.H1.model1, pred.H1.model1, accuracy.H1.model1,
     Predict_x, Predict_y, Freq.model1, pred.model1, accuracy.model1,
     file="ModelPerformance.RData")

#find false negative
FN.id <- which(H1.data$prediction==0 & H1.data$is_Enh==1)
length(FN.id)
write.table(as.data.frame(Merge.bed)[FN.id,1:3], 
            file = "H1.FalseNegative.bed", quote = F,
            col.names = F, sep = "\t", row.names = F)

length(ChromHMM.Enh.H1.bed)
sum(TrainSet$is_Enh==1)+sum(TestSet$is_Enh==1)

#3. output predict result
K562.Chr22.prediction <- data.frame("name"=rownames(PredictionSet))
K562.Chr22.prediction <- cbind.data.frame(K562.Chr22.prediction, PredictionSet[,2:13])
K562.Chr22.prediction$is_Enh <- PredictionSet$is_Enh
K562.Chr22.prediction$prediction <- PredictionSet$prediction
write.table(K562.Chr22.prediction, 
            file = "K562.Chr22.prediction.tab", quote = F,
            col.names = T, sep = "\t", row.names = T)
write.table(H1.data, 
            file = "H1.prediction.tab", quote = F,
            col.names = F, sep = "\t", row.names = F)

#4. output Config
#prepare for final software 
InputData <- names(TrainSet)[2:13]
SelectedFeature
FeatureFileName

save(InputData, SelectedFeature, FeatureFileName, model1,
     file = "../Config.RData")
