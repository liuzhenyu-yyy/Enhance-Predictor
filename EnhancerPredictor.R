args=commandArgs(T)

load("Config.RData")
library(genomation)
library(GenomicRanges)
library(e1071)
warnings('off')
#args <- c("Performance/ChromHMM.Enh.H1.bed","Hi.out.bed")

#read bed file
Input.bed <- readGeneric(args[1], 
                         header = F,zero.based = T)[1:2000]

#calculate feature value
for(i in 1:length(SelectedFeature)){
  print(paste("Calculating feature value:",SelectedFeature[i],'...'))
  temp <- ScoreMatrixBin(target = paste("./data/K562_BroadHistone/",
                                        FeatureFileName[i], sep=''),
                         windows = Input.bed,
                         bin.num = 1)
  mcols(Input.bed)[,SelectedFeature[i]] <- temp@.Data
}

#calculate ATAC signal
print(paste("Calculating feature value:","ATAC","..."))
ATAC_rep1 <- ScoreMatrixBin(target = "./data/K562_ATAC-Seq/GSM1782764_K562_50K_ATAC-seq_nan_nan_nan_1_1_hg19.bigWig",
                            windows = Input.bed,
                            bin.num = 1)
ATAC_rep2 <- ScoreMatrixBin(target = "./data/K562_ATAC-Seq/GSM1782765_K562_50K_ATAC-seq_nan_nan_nan_2_1_hg19.bigWig",
                            windows = Input.bed,
                            bin.num = 1)
#mean value of two repilcates
Input.bed$ATAC <- (ATAC_rep1@.Data + ATAC_rep2@.Data)/ 2

#normalize to input
print("Normalizing Signal Density...")
control.td <- ScoreMatrixBin(target = "./data/K562_BroadHistone/wgEncodeBroadHistoneK562ControlStdSig.bigWig",
                             windows = Input.bed,
                             bin.num = 1)

for(feature in SelectedFeature){
  mcols(Input.bed)[,feature] <- mcols(Input.bed)[,feature]-control.td
}

#prediction
print("Classification...")
Input.data <- as.data.frame(Input.bed)[,InputData]

pred.model1 <- predict(object = model1, newdata = Input.data)

Input.data$prediction <- pred.model1

#output results
print("Generating output file...")
Output.data <- as.data.frame(Input.bed)
Output.data$prediction <-  pred.model1
write.table(Output.data, 
            file = args[2], quote = F,
            col.names = T, sep = "\t", row.names = F)

print("Done.")