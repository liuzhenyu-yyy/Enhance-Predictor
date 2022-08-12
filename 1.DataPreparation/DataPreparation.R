setwd("D:/Final/DataPreparation/")

library(stringr)
library(GenomicRanges)
library(genomation)
library(ChIPseeker)
library(data.table)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

TxDb.hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene

#1.prepare positive data (true enhancers)###########

#downloaded from UCSC Vista enhancers and ChromHMM data sets
VistaEnhancer.bed <- readGeneric("../data/VistaEnhancer.bed", 
                                 header = F,zero.based = T)  
ChromHMM.K562.bed <- readGeneric("../data/ChromHMM.K562.bed", 
                                 header = F,zero.based = T,
                                 meta.cols = list(annotation=4))  

#select enhancer from ChromHMM data set
is_enhancer <- str_detect(ChromHMM.K562.bed$annotation, "Strong_Enhancer")
table(is_enhancer)
ChromHMM.Enh.bed <- ChromHMM.K562.bed[is_enhancer]

#select 23 Chromosomes for analysis
SelectedChrom <- names(table(VistaEnhancer.bed@seqnames))
SelectedChrom

#compare two different data sets
length(ChromHMM.Enh.bed)
table(ChromHMM.Enh.bed@seqnames)

length(VistaEnhancer.bed)
table(VistaEnhancer.bed@seqnames) #too small train set

quantile(width(VistaEnhancer.bed))
quantile(width(ChromHMM.Enh.bed)) #some region are too short
  
rm(VistaEnhancer.bed, ChromHMM.K562.bed) #use ChromHMM as golden standard

#quality control of enhancers
hist(log2(width(ChromHMM.Enh.bed)), breaks = 50)
width.cutoff <- 8
abline(v = width.cutoff, col="red")

ChromHMM.Enh.bed <- ChromHMM.Enh.bed[log2(width(ChromHMM.Enh.bed))>width.cutoff]
table(ChromHMM.Enh.bed@seqnames)  
quantile(width(ChromHMM.Enh.bed)) #reasonable width

#annotate peak to TSSs
ChromHMM.Enh.bed.ann <- annotatePeak(ChromHMM.Enh.bed, tssRegion = c(-3000, 3000),
                                TxDb = TxDb.hg19)

#select distal region (>3000kb) as true enhaner 
is.distal <- abs(ChromHMM.Enh.bed.ann@anno$distanceToTSS)>3000
table(is.distal)

TrueEnhancer.bed <- ChromHMM.Enh.bed[is.distal]

table(TrueEnhancer.bed@seqnames)
quantile(width(ChromHMM.Enh.bed))


#2.prepare negative data (false enhancers)##########
#output TrueEnhancer.bed
write.table(as.data.frame(TrueEnhancer.bed)[,1:3], 
            file = "TrueEnhancer.bed", quote = F,
            col.names = F, sep = "\t", row.names = F)

#generate tss 3k region (promoters)
tss.3k.bed <- promoters(TxDb.hg19, upstream = 3000, downstream = 3000)
start(tss.3k.bed)[start(tss.3k.bed)<0] <- 0
tss.3k.bed <- tss.3k.bed[tss.3k.bed@seqnames %in% SelectedChrom]

#exclude tss region and True enhancer region
exclude.bed <- reduce(c(TrueEnhancer.bed, tss.3k.bed)) #regions excluded in shuffling
exclude.bed <- exclude.bed[exclude.bed@seqnames%in%SelectedChrom]
write.table(as.data.frame(exclude.bed)[,1:3], 
            file = "exclude.bed", quote = F,
            col.names = F, sep = "\t", row.names = F)


#shuffle bed using bedtools
#bedtools shuffle -i TrueEnhancer.bed -g hg19.genome.size -excl exclude.bed -chrom > FalseEnhancer.bed
FalseEnhancer.bed <- readGeneric("./FalseEnhancer.bed", header = F,zero.based = T)  

#check properties of false enhancers
#same genomic distribution
table(FalseEnhancer.bed@seqnames)
table(TrueEnhancer.bed@seqnames) 
plot.data <- data.frame(table(TrueEnhancer.bed@seqnames))
names(plot.data) <- c("Chromosome","Positive")
plot.data$Negative <- table(FalseEnhancer.bed@seqnames)
pdf("ChromosomeDistrubution.pdf",4.3,3)
ggplot(plot.data, aes(x = Positive, y =Negative , color = Chromosome))+
  geom_point() +
  xlab("peak number: positive")+
  ylab("peak number: negative")
dev.off()

#same width distribution
quantile(width(FalseEnhancer.bed))
quantile(width(TrueEnhancer.bed))
end(FalseEnhancer.bed) <- end(FalseEnhancer.bed) + 1
pdf("Width.pdf")
par(mfrow=c(2,1))
hist(width(FalseEnhancer.bed))
hist(width(TrueEnhancer.bed))
dev.off()
#no overlap with true enhancers or promoters
table(countOverlaps(FalseEnhancer.bed,TrueEnhancer.bed))
table(countOverlaps(FalseEnhancer.bed,tss.3k.bed))
list1 <-GenomicRanges::GRangesList("Posotive"=TrueEnhancer.bed,
                                   "Negative"=FalseEnhancer.bed)
pdf("cover.pdf")
covplot(FalseEnhancer.bed,chrs = "chr22",xlim = c(22000,32000))
covplot(TrueEnhancer.bed,chrs = "chr22",xlim = c(22000,32000))
covplot(tss.3k.bed,chrs = "chr22",xlim = c(22000,32000))
dev.off()

#3.calculate feature data###########################################

FalseEnhancer.bed$is_Enh <- 0
TrueEnhancer.bed$is_Enh <- 1

Merge.bed <- c(FalseEnhancer.bed, TrueEnhancer.bed)
table(Merge.bed$is_Enh)

#downloaded from USCS Broad Histone 
# http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeBroadHistone

AllFiles <- dir("../data/K562_BroadHistone/")
AllFiles

#get feature name
AllFeatures <- AllFiles
for(i in 1:length(AllFeatures)){
  AllFeatures[i] <- substr(AllFeatures[i],start = 25, stop = str_length(AllFeatures[i])-10)
}
AllFeatures

#calculate feature value
for(i in 1:length(AllFeatures)){
  temp <- ScoreMatrixBin(target = paste("../data/K562_BroadHistone/",
                                        AllFiles[i], sep=''),
                         windows = Merge.bed,
                         bin.num = 1)
  mcols(Merge.bed)[,AllFeatures[i]] <- temp@.Data
  print(i)
}

#normalize to control(input) data
for(feature in AllFeatures){
  mcols(Merge.bed)[,feature] <- mcols(Merge.bed)[,feature]-mcols(Merge.bed)[,"ControlStd"]
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

AllFeatures <- c(AllFeatures, "ATAC")
AllFeatures <- AllFeatures[AllFeatures!="ControlStd"]

mcols(Merge.bed) <- mcols(Merge.bed)[,c("is_Enh",AllFeatures)]

#4. Split data set for training and testing 
TrueEnhancer.bed <- Merge.bed[Merge.bed$is_Enh==1]
FalseEnhancer.bed <- Merge.bed[Merge.bed$is_Enh==0]
names(TrueEnhancer.bed) <- paste("TE", 1:length(TrueEnhancer.bed), sep='')
names(FalseEnhancer.bed) <- paste("FE", 1:length(FalseEnhancer.bed), sep='')

#chr22 for final prediction
chr22.TrueEhaner.bed <- TrueEnhancer.bed[TrueEnhancer.bed@seqnames == "chr22"]
chr22.FalseEnhancer.bed <- FalseEnhancer.bed[FalseEnhancer.bed@seqnames == "chr22"]


# sampling train set and test set
Non_Chr22.TrueEhaner <- names(TrueEnhancer.bed[TrueEnhancer.bed@seqnames!="chr22"])
Non_Chr22.FlaseEhaner <- names(FalseEnhancer.bed[FalseEnhancer.bed@seqnames!="chr22"])

Train_True <- sample(Non_Chr22.TrueEhaner, length(Non_Chr22.TrueEhaner)/2)
Test_True <- Non_Chr22.TrueEhaner[!Non_Chr22.TrueEhaner%in%Train_True]

Train_False <- sample(Non_Chr22.FlaseEhaner, length(Non_Chr22.FlaseEhaner)/2)
Test_False <- Non_Chr22.FlaseEhaner[!Non_Chr22.FlaseEhaner%in%Train_False]

Train_True.bed <- TrueEnhancer.bed[Train_True]
Test_True.bed <- TrueEnhancer.bed[Test_True]
Train_False.bed <- FalseEnhancer.bed[Train_False]
Test_False.bed <- FalseEnhancer.bed[Test_False]

#4. output for downstream analysis
save(Train_True.bed, Test_True.bed,
     Train_False.bed, Test_False.bed,
     chr22.TrueEhaner.bed, chr22.FalseEnhancer.bed,
     AllFeatures, 
     file = "DataSets_AllFeature.Rdata")

dir.create("DataSets_AllFeature")
write.table(mcols(Train_True.bed), file = "DataSets_AllFeature/Train_True.tab",
            sep="\t", quote = F, row.names = F)
write.table(mcols(Train_False.bed), file = "DataSets_AllFeature/Train_False.tab",
            sep="\t", quote = F, row.names = F)
write.table(mcols(Test_True.bed), file = "DataSets_AllFeature/Test_True.tab",
            sep="\t", quote = F, row.names = F)
write.table(mcols(Test_False.bed), file = "DataSets_AllFeature/Test_False.tab",
            sep="\t", quote = F, row.names = F)