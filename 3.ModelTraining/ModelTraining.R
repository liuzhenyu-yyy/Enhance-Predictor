setwd("D:/Final/ModelTraining/")

load("../FeatureSelection/DataSets_SelectedFeature.Rdata")

library(e1071)

#1. simple SVM classification###########################################
#Select the best kernel function and classification methon
SVM_test <- function(x,y){
  type <- c('C-classification','nu-classification')
  kernel <- c('linear','polynomial','radial','sigmoid')
  pred <- array(0, dim=c(nrow(x),2,4))
  errors <- matrix(0,2,4)
  dimnames(errors) <- list(type, kernel)
  for(i in 1:2){
    for(j in 1:4){
      print("Test SVM parameters:")
      print(paste("Classification:",type[i]))
      print(paste("kernal:",kernel[j]))
      pred[,i,j] <- predict(object = svm(x, y, type = type[i], kernel = kernel[j]), 
                            newdata = x)
      errors[i,j] <- sum(pred[,i,j] != as.integer(y))
      print("done......")
      print("########################")
    }
  }
  return(errors)
}

train_x <- TrainSet[,2:13]
train_y <- TrainSet[,1]
test_x <- TestSet[,2:13]
test_y <- TestSet[,1]

SVMtest <- SVM_test(x=train_x,y=train_y)
SVMtest

model1 <- svm(x=train_x, y=train_y, type = 'C-classification', kernel = 'radial')

pred.model1 <- predict(object = model1, newdata = test_x)
Freq.model1 <- table(test_y, pred.model1)
Freq.model1
accuracy.model1 <- sum(diag(Freq.model1))/sum(Freq.model1)
accuracy.model1

TrainSet$prediction <- pred.model1
write.table(TrainSet, file = "TrainSet_Prediction.tab",
            sep="\t", quote = F, row.names = F)

#adjust hyperparameter
#This step takes a lot of time.
#downsampling to speed up calculation
DownSample <- sample(rownames(TrainSet), 500)
tune.svm(is_Enh ~., data =TrainSet[DownSample,], gamma = 10^(-100:0), cost = 10^(0:3))

model2 <- svm(x=train_x, y=train_y, type = 'C-classification', kernel = 'radial',
              gamma = 0.01, cost = 100)

pred.model2 <- predict(object = model2, newdata = test_x)
Freq.model2 <- table(test_y, pred.model2)
Freq.model2
accuracy.model2 <- sum(diag(Freq.model2))/sum(Freq.model2)
accuracy.model2

#2.perform 10 fold cross validation#########################################
library(caret)

# Define train control for k fold cross validation
AllData <- rbind.data.frame(TestSet,TrainSet)
AllData$is_Enh <- factor(AllData$is_Enh)
train_control <- trainControl(method="cv", number=10)

#10 fold cross test
#This command takes a lot of time
model_CV <- train(is_Enh~., data=AllData, trControl=train_control, method="svmRadialSigma")

# Summarise Results
print(model_CV)


#3. Save the model for performance evaluation.
save(model1, SVMtest, Freq.model1, accuracy.model1,
     file = "model1.RData")

