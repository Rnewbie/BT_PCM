library(RCurl)
library(caret)
library(C50)
library(dplyr)
library(ROCR)
library(ggplot2)
library(randomForest)
index = c(20,44,26,31,11,42,17,34,29,39,9,23,38,22,24,3,33,8,18,21,5,37,28,14,15,35,30,
          41,1,43,177,211,262,58,151,218,88,231,180,274,256,104,84,205,65,175,279,49,59,
          145,174,242,196,283,51,219,273,260,206,74)
### protein alone
p <- getURL("https://raw.githubusercontent.com/Rnewbie/BT_PCM/master/data/request/protein.csv")
protein_alone_data <- read.csv(text = p, header = TRUE)
protein_alone_train <- protein_alone_data[index,]
protein_alone_test <- protein_alone_data[-index,]
c <- getURL("https://raw.githubusercontent.com/Rnewbie/BT_PCM/master/data/request/compound.csv")
compound_alone_data <- read.csv(text = c, header = TRUE)
compound_alone_train <- compound_alone_data[index,]
compound_alone_test <- compound_alone_data[-index,]
cxp <- getURL("https://raw.githubusercontent.com/Rnewbie/BT_PCM/master/data/request/CxP.csv")
CxP_data <- read.csv(text = cxp, header = TRUE)
CxP_train <- CxP_data[index,]
CxP_test <- CxP_data[-index,]
cxc <- getURL("https://raw.githubusercontent.com/Rnewbie/BT_PCM/master/data/request/CxC.csv")
CxC_data <- read.csv(text = cxc, header = TRUE)
CxC_train <- CxC_data[index,]
CxC_test <- CxC_data[-index,]
pxp <- getURL("https://raw.githubusercontent.com/Rnewbie/BT_PCM/master/data/request/PxP.csv")
PxP_data <- read.csv(text = pxp, header = TRUE)
PxP_train <- PxP_data[index,]
PxP_test <- PxP_data[-index,]


### Boss email Morning RF, ANN, SVM (three kernel)
training_randomforest <- function(x) {
  fit <- randomForest(activity~., data = x)
  actual <- x$activity
  prediction <- predict(fit, x)
  pred <- prediction(predict(fit, newdata = x, type = "prob")[, 2], x$activity)
  perf.auc <- performance(pred, measure = "auc")
  AUC <- unlist(perf.auc@y.values)
  confusionmatrix  <- confusionMatrix(prediction, actual)
  table <- confusionmatrix$table
  results <- as.numeric(table)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  return(round(c(ACC, SENS, SPEC, AUC, MCC), digits = 2))
}

### 10 fold cross validation RF

CV_randomforest <- function(x) {
  myData <- x
  myData$activity <- factor(as.character(myData$activity))
  myRes <- data.frame()
  k = 10
  index <- sample(1:k, nrow(myData), replace = TRUE)
  folds <- 1:k
  for (j in 1:k)
    training <- subset(myData, index %in% folds[-j])
  testing <- subset(myData, index %in% c(j))
  fit <- randomForest(activity~., data = training)
  actual <- testing$activity
  prediction <- predict(fit, testing)
  pred <- prediction(predict(fit, newdata = testing, type = "prob")[, 2], testing$activity)
  perf.auc <- performance(pred, measure = "auc")
  AUC <- unlist(perf.auc@y.values)
  confusionmatrix  <- confusionMatrix(prediction, actual)
  table <- confusionmatrix$table
  results <- as.numeric(table)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  return(round(c(ACC, SENS, SPEC, AUC, MCC), digits = 2))
}



#### randomForest Testing

testing_randomforest <- function(x) {
  fit <- randomForest(activity~., data = x)
  actual <- CxC_test$activity
  prediction <- predict(fit, CxC_test)
  pred <- prediction(predict(fit, newdata = CxC_test, type = "prob")[, 2], CxC_test$activity)
  perf.auc <- performance(pred, measure = "auc")
  AUC <- unlist(perf.auc@y.values)
  confusionmatrix  <- confusionMatrix(prediction, actual)
  table <- confusionmatrix$table
  results <- as.numeric(table)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  return(round(c(ACC, SENS, SPEC, AUC, MCC), digits = 2))
}



#### Artificial neural network

training_neural <- function(x) {
  fit <- nnet(activity~., data = x, size = 5, MaxNWts = 4000, na.action = na.omit)
  actual <- x$activity
  prediction <- predict(fit, x, type = "class")
  pred <- prediction(predict(fit, newdata = x, type = "raw")[, 1], x$activity)
  perf.auc <- performance(pred, measure = "auc")
  AUC <- unlist(perf.auc@y.values)
  confusionmatrix  <- confusionMatrix(prediction, actual)
  table <- confusionmatrix$table
  results <- as.numeric(table)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  return(round(c(ACC, SENS, SPEC, AUC, MCC), digits = 2))
}





#### CV ANN

CV_neural <- function(x) {
  myData <- x
  myData$activity <- factor(as.character(myData$activity))
  myRes <- data.frame()
  k = 10
  index <- sample(1:k, nrow(myData), replace = TRUE)
  folds <- 1:k
  for (j in 1:k)
    training <- subset(myData, index %in% folds[-j])
  testing <- subset(myData, index %in% c(j))
  fit <- nnet(activity~., data = x, size = 5, MaxNWts = 4000, na.action = na.omit)
  actual <- testing$activity
  prediction <- predict(fit, testing, type = "class")
  pred <- prediction(predict(fit, newdata = testing, type = "raw")[, 1], factor(testing$activity))
  perf.auc <- performance(pred, measure = "auc")
  AUC <- unlist(perf.auc@y.values)
  confusionmatrix  <- confusionMatrix(prediction, actual)
  table <- confusionmatrix$table
  results <- as.numeric(table)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  return(round(c(ACC, SENS, SPEC, AUC, MCC), digits = 2))
}


#### testing for artificial neural network

testing_neural <- function(x) {
  fit <- nnet(activity~., data = x, size = 5, MaxNWts = 4000, na.action = na.omit)
  actual <- CxC_test$activity
  prediction <- predict(fit, CxC_test, type = "class")
  pred <- prediction(predict(fit, newdata = CxC_test, type = "raw")[, 1], CxC_test$activity)
  perf.auc <- performance(pred, measure = "auc")
  AUC <- unlist(perf.auc@y.values)
  confusionmatrix  <- confusionMatrix(prediction, actual)
  table <- confusionmatrix$table
  results <- as.numeric(table)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  return(round(c(ACC, SENS, SPEC, AUC, MCC), digits = 2))
}



#### SVM radial

training_SVM_radial <- function(x) {
  fit <- svm(activity~., data = x, na.action = na.omit, scale = FALSE, kernel = "radial", probability = TRUE)
  actual <- x$activity
  prediction <- predict(fit, x, type = "class")
  prediction_auc <- predict(fit, newdata = x, probability = TRUE)
  attributes_auc <- data.frame(attributes(prediction_auc))[, 5]
  pred <- prediction(predictions = attributes_auc, x$activity)
  perf.auc <- performance(pred, measure = "auc")
  AUC <- unlist(perf.auc@y.values)
  confusionmatrix  <- confusionMatrix(prediction, actual)
  table <- confusionmatrix$table
  results <- as.numeric(table)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  return(round(c(ACC, SENS, SPEC, AUC, MCC), digits = 2))
}



### CV Radial

CV_radial <- function(x) {
  myData <- x
  myData$activity <- factor(as.character(myData$activity))
  myRes <- data.frame()
  k = 10
  index <- sample(1:k, nrow(myData), replace = TRUE)
  folds <- 1:k
  for (j in 1:k)
    training <- subset(myData, index %in% folds[-j])
  testing <- subset(myData, index %in% c(j))
  fit <- svm(activity~., data = training, na.action = na.omit, scale = FALSE, kernel = "radial", probability = TRUE)
  actual <- testing$activity
  prediction <- predict(fit, testing, type = "class")
  prediction_auc <- predict(fit, newdata = testing, probability = TRUE)
  attributes_auc <- attributes(prediction_auc)$probabilities[, 2]
  pred <- prediction(predictions = attributes_auc, testing$activity)
  perf.auc <- performance(pred, measure = "auc")
  AUC <- unlist(perf.auc@y.values)
  confusionmatrix  <- confusionMatrix(prediction, actual)
  table <- confusionmatrix$table
  results <- as.numeric(table)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  return(round(c(ACC, SENS, SPEC, AUC, MCC), digits = 2))
}


#### Testing SVM radial

testing_SVM_radial <- function(x) {
  fit <- svm(activity~., data = x, na.action = na.omit, scale = FALSE, kernel = "radial", probability = TRUE)
  actual <- CxC_test$activity
  prediction <- predict(fit, CxC_test, type = "class")
  prediction_auc <- predict(fit, newdata = CxC_test, probability = TRUE)
  attributes_auc <- data.frame(attributes(prediction_auc))[, 5]
  pred <- prediction(predictions = attributes_auc, CxC_test$activity)
  perf.auc <- performance(pred, measure = "auc")
  AUC <- unlist(perf.auc@y.values)
  confusionmatrix  <- confusionMatrix(prediction, actual)
  table <- confusionmatrix$table
  results <- as.numeric(table)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  return(round(c(ACC, SENS, SPEC, AUC, MCC), digits = 2))
}


### Training SVM linear 

training_SVM_linear <- function(x) {
  fit <- svm(activity~., data = x, na.action = na.omit, scale = FALSE, kernel = "linear", probability = TRUE)
  actual <- x$activity
  prediction <- predict(fit, x, type = "class")
  prediction_auc <- predict(fit, newdata = x, probability = TRUE)
  attributes_auc <- data.frame(attributes(prediction_auc))[, 5]
  pred <- prediction(predictions = attributes_auc, x$activity)
  perf.auc <- performance(pred, measure = "auc")
  AUC <- unlist(perf.auc@y.values)
  confusionmatrix  <- confusionMatrix(prediction, actual)
  table <- confusionmatrix$table
  results <- as.numeric(table)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  return(round(c(ACC, SENS, SPEC, AUC, MCC), digits = 2))
}


### 10 fold CV SVM linear


CV_linear <- function(x) {
  myData <- x
  myData$activity <- factor(as.character(myData$activity))
  myRes <- data.frame()
  k = 10
  index <- sample(1:k, nrow(myData), replace = TRUE)
  folds <- 1:k
  for (j in 1:k)
    training <- subset(myData, index %in% folds[-j])
  testing <- subset(myData, index %in% c(j))
  fit <- svm(activity~., data = training, na.action = na.omit, scale = FALSE, kernel = "linear", probability = TRUE)
  actual <- testing$activity
  prediction <- predict(fit, testing, type = "class")
  prediction_auc <- predict(fit, newdata = testing, probability = TRUE)
  attributes_auc <- attributes(prediction_auc)$probabilities[, 2]
  pred <- prediction(predictions = attributes_auc, testing$activity)
  perf.auc <- performance(pred, measure = "auc")
  AUC <- unlist(perf.auc@y.values)
  confusionmatrix  <- confusionMatrix(prediction, actual)
  table <- confusionmatrix$table
  results <- as.numeric(table)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  return(round(c(ACC, SENS, SPEC, AUC, MCC), digits = 2))
}


### Testing SVM linear

testing_SVM_linear <- function(x) {
  fit <- svm(activity~., data = x, na.action = na.omit, scale = FALSE, kernel = "linear", probability = TRUE)
  actual <- CxC_test$activity
  prediction <- predict(fit, CxC_test, type = "class")
  prediction_auc <- predict(fit, newdata = CxC_test, probability = TRUE)
  attributes_auc <- data.frame(attributes(prediction_auc))[, 5]
  pred <- prediction(predictions = attributes_auc, CxC_test$activity)
  perf.auc <- performance(pred, measure = "auc")
  AUC <- unlist(perf.auc@y.values)
  confusionmatrix  <- confusionMatrix(prediction, actual)
  table <- confusionmatrix$table
  results <- as.numeric(table)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  return(round(c(ACC, SENS, SPEC, AUC, MCC), digits = 2))
}


### Training SVM polynomial

training_SVM_polynomial <- function(x) {
  fit <- svm(activity~., data = x, na.action = na.omit, scale = FALSE, kernel = "polynomial", probability = TRUE)
  actual <- x$activity
  prediction <- predict(fit, x, type = "class")
  prediction_auc <- predict(fit, newdata = x, probability = TRUE)
  attributes_auc <- data.frame(attributes(prediction_auc))[, 5]
  pred <- prediction(predictions = attributes_auc, x$activity)
  perf.auc <- performance(pred, measure = "auc")
  AUC <- unlist(perf.auc@y.values)
  confusionmatrix  <- confusionMatrix(prediction, actual)
  table <- confusionmatrix$table
  results <- as.numeric(table)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  return(round(c(ACC, SENS, SPEC, AUC, MCC), digits = 2))
}


### SVM polynomial 10-CV

CV_polynomial <- function(x) {
  myData <- x
  myData$activity <- factor(as.character(myData$activity))
  myRes <- data.frame()
  k = 10
  index <- sample(1:k, nrow(myData), replace = TRUE)
  folds <- 1:k
  for (j in 1:k)
    training <- subset(myData, index %in% folds[-j])
  testing <- subset(myData, index %in% c(j))
  fit <- svm(activity~., data = training, na.action = na.omit, scale = FALSE, kernel = "polynomial", probability = TRUE)
  actual <- testing$activity
  prediction <- predict(fit, testing, type = "class")
  prediction_auc <- predict(fit, newdata = testing, probability = TRUE)
  attributes_auc <- attributes(prediction_auc)$probabilities[, 2]
  pred <- prediction(predictions = attributes_auc, testing$activity)
  perf.auc <- performance(pred, measure = "auc")
  AUC <- unlist(perf.auc@y.values)
  confusionmatrix  <- confusionMatrix(prediction, actual)
  table <- confusionmatrix$table
  results <- as.numeric(table)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  return(round(c(ACC, SENS, SPEC, AUC, MCC), digits = 2))
}


### SVM testing Polynomial

testing_SVM_polynomial <- function(x) {
  fit <- svm(activity~., data = x, na.action = na.omit, scale = FALSE, kernel = "polynomial", probability = TRUE)
  actual <- CxC_test$activity
  prediction <- predict(fit, CxC_test, type = "class")
  prediction_auc <- predict(fit, newdata = CxC_test, probability = TRUE)
  attributes_auc <- data.frame(attributes(prediction_auc))[, 5]
  pred <- prediction(predictions = attributes_auc, CxC_test$activity)
  perf.auc <- performance(pred, measure = "auc")
  AUC <- unlist(perf.auc@y.values)
  confusionmatrix  <- confusionMatrix(prediction, actual)
  table <- confusionmatrix$table
  results <- as.numeric(table)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  return(round(c(ACC, SENS, SPEC, AUC, MCC), digits = 2))
}


































































































































































































































































View(compound_alone_data)
boss_question <- function(x) {
  Active <- subset(x, activity = "Active")
  Inactive <- subset(x, activity = "Inactive")
  results <- list(100)
  for (i in 1:100) {
    train_active <- sample_n(Active, size = 35)
    test_active <- sample_n(Active, size = 9)
    train_inactive <- sample_n(Inactive, size = 195)
    test_inactive <- sample_n(Inactive, size = 49)
    train <- rbind(train_active, train_inactive)
    test <- rbind(test_active, test_inactive)
    model_train <- C5.0(activity~., data = train)
    actual <- train$activity
    prediction <- predict(model_train, train)
    confusionmatrix  <- confusionMatrix(prediction, actual)
    results[[i]] <- confusionmatrix$table
  }
  return(results)
}

#### ROC curve for the training data
### For compound
results_compound <- C5.0(activity~., data = compound_alone_train)
actual <- compound_alone_train$activity
prediction_prob <- predict(results_compound, newdata = compound_alone_train, type = "prob")[, 1]
prediction_class <- predict(results_compound, compound_alone_train)
confusionmatrix <- confusionMatrix(prediction_class, actual)
table <- confusionmatrix$table
results <- as.numeric(table)
rocCurve <- pROC::roc(actual, prediction_prob)
plot.roc(rocCurve, col="blue", legacy.axes = TRUE, xlim = c(1,0), ylim = c(0,1))### Blue for compound
### For protein
results_protein <- C5.0(activity~., data = protein_alone_train)
actual <- protein_alone_train$activity
prediction_prob <- predict(results_protein, protein_alone_train, type = "prob")[, 1]
prediction_class <- predict(results_protein, protein_alone_train)
confusionmatrix <- confusionMatrix(prediction_class, actual)
table <- confusionmatrix$table
results <- as.numeric(table)

rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
plot.roc(rocCurve, col = "red", add = TRUE, legacy.axes = TRUE)

## For crossterms
results_CxP <- C5.0(activity~., data = CxP_train)
actual <- CxP_train$activity
prediction_prob <- predict(results_CxP, CxP_train, type = "prob")[, 1]
prediction_class <- predict(results_CxP, CxP_train)
confusionmatrix <- confusionMatrix(prediction_class, actual)
table <- confusionmatrix$table
results <- as.numeric(table)
rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
plot.roc(rocCurve, col = "green", add = TRUE, legacy.axes = TRUE)

### for Protein-Protein
results_PxP <- C5.0(activity~., data = PxP_train)
actual <- PxP_train$activity
prediction_prob <- predict(results_PxP, PxP_train, type = "prob")[, 1]
prediction_class <- predict(results_PxP, PxP_train)
confusionmatrix <- confusionMatrix(prediction_class, actual)
table <- confusionmatrix$table
results <- as.numeric(table)
rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
plot.roc(rocCurve, col = "black", add = TRUE, legacy.axes = TRUE)

### For compound-compound
results_CxC <- C5.0(activity~., data = CxC_train)
actual <- CxC_train$activity
prediction_prob <- predict(results_CxC, CxC_train, type = "prob")[, 1]
prediction_class <- predict(results_CxC, CxC_train)
confusionmatrix <- confusionMatrix(prediction_class, actual)
table <- confusionmatrix$table
results <- as.numeric(table)

rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
plot.roc(rocCurve, col = "yellow", add = TRUE, legacy.axes = TRUE)

legend("bottomright", c("C", "P", "CxP", "PxP", "CxC"), fill = c("blue","red", "green", "black", "yellow"))


### ROC curve for the 10 fold cross validation

### compound
myData <- compound_alone_train
myData$activity <- factor(as.character(myData$activity))
k = 10
index <- sample(1:k, nrow(myData), replace = TRUE)
folds <- 1:k
myRes <- data.frame()
for (i in 1:k) {
  training <- subset(myData, index %in% folds[-i])
  testing <- subset(myData, index %in% c(i))
  mymodel <- C5.0(activity~., data = training)
  actual <- testing$activity
  prediction <- predict(mymodel, testing, type = "prob")[, 1]
  results <- data.frame(Predicted = prediction, Actual = actual)
  myRes = rbind(myRes, results)
}

actual <- myRes$Actual
prediction_prob <- myRes$Predicted
rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
plot.roc(rocCurve, col="blue", legacy.axes = TRUE, xlim = c(1,0), ylim = c(0,1)) ### BLUE

### protein
myData <- protein_alone_train
myData$activity <- factor(as.character(myData$activity))
k = 10
index <- sample(1:k, nrow(myData), replace = TRUE)
folds <- 1:k
myRes <- data.frame()
for (i in 1:k) {
  training <- subset(myData, index %in% folds[-i])
  testing <- subset(myData, index %in% c(i))
  mymodel <- C5.0(activity~., data = training)
  actual <- testing$activity
  prediction <- predict(mymodel, testing, type = "prob")[, 1]
  results <- data.frame(Predicted = prediction, Actual = actual)
  myRes = rbind(myRes, results)
}

actual <- myRes$Actual
prediction_prob <- myRes$Predicted
rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
plot.roc(rocCurve, col = "red", add = TRUE, legacy.axes = TRUE)

### Cross-terms
myData <- CxP_train
myData$activity <- factor(as.character(myData$activity))
k = 10
index <- sample(1:k, nrow(myData), replace = TRUE)
folds <- 1:k
myRes <- data.frame()
for (i in 1:k) {
  training <- subset(myData, index %in% folds[-i])
  testing <- subset(myData, index %in% c(i))
  mymodel <- C5.0(activity~., data = training)
  actual <- testing$activity
  prediction <- predict(mymodel, testing, type = "prob")[, 1]
  results <- data.frame(Predicted = prediction, Actual = actual)
  myRes = rbind(myRes, results)
}

actual <- myRes$Actual
prediction_prob <- myRes$Predicted
rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
plot.roc(rocCurve, col = "green", add = TRUE, legacy.axes = TRUE)

### cross_terms Protein Protein
myData <- PxP_train
myData$activity <- factor(as.character(myData$activity))
k = 10
index <- sample(1:k, nrow(myData), replace = TRUE)
folds <- 1:k
myRes <- data.frame()
for (i in 1:k) {
  training <- subset(myData, index %in% folds[-i])
  testing <- subset(myData, index %in% c(i))
  mymodel <- C5.0(activity~., data = training)
  actual <- testing$activity
  prediction <- predict(mymodel, testing, type = "prob")[, 1]
  results <- data.frame(Predicted = prediction, Actual = actual)
  myRes = rbind(myRes, results)
}

actual <- myRes$Actual
prediction_prob <- myRes$Predicted
rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
plot.roc(rocCurve, col = "black", add = TRUE, legacy.axes = TRUE)

### compound compound cross terms 
myData <- CxC_train
myData$activity <- factor(as.character(myData$activity))
k = 10
index <- sample(1:k, nrow(myData), replace = TRUE)
folds <- 1:k
myRes <- data.frame()
for (i in 1:k) {
  training <- subset(myData, index %in% folds[-i])
  testing <- subset(myData, index %in% c(i))
  mymodel <- C5.0(activity~., data = training)
  actual <- testing$activity
  prediction <- predict(mymodel, testing, type = "prob")[, 1]
  results <- data.frame(Predicted = prediction, Actual = actual)
  myRes = rbind(myRes, results)
}

actual <- myRes$Actual
prediction_prob <- myRes$Predicted
rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
plot.roc(rocCurve, col = "yellow", add = TRUE, legacy.axes = TRUE)

legend("bottomright", c("C", "P", "CxP", "PxP", "CxC"), fill = c("blue","red", "green", "black", "yellow"))


### ROC curve for testing data set
### For compound
results_compound <- C5.0(activity~., data = compound_alone_train)
actual <- compound_alone_test$activity
prediction_prob <- predict(results_compound, compound_alone_test, type = "prob")[, 1]
prediction_class <- predict(results_compound, compound_alone_test)
rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
plot.roc(rocCurve, col="blue", legacy.axes = TRUE, xlim = c(1,0), ylim = c(0,1))### Blue for compound
### For protein
results_protein <- C5.0(activity~., data = protein_alone_train)
actual <- protein_alone_test$activity
prediction_prob <- predict(results_protein, protein_alone_test, type = "prob")[, 1]
prediction_class <- predict(results_protein, protein_alone_test)
rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
plot.roc(rocCurve, col = "red", add = TRUE, legacy.axes = TRUE)

## For crossterms
results_CxP <- C5.0(activity~., data = CxP_train)
actual <- CxP_test$activity
prediction_prob <- predict(results_CxP, CxP_test, type = "prob")[, 1]
prediction_class <- predict(results_CxP, CxP_test)
rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
plot.roc(rocCurve, col = "green", add = TRUE, legacy.axes = TRUE)

### for Protein-Protein
results_PxP <- C5.0(activity~., data = PxP_train)
actual <- PxP_test$activity
prediction_prob <- predict(results_PxP, PxP_test, type = "prob")[, 1]
prediction_class <- predict(results_PxP, PxP_test)
rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
plot.roc(rocCurve, col = "black", add = TRUE, legacy.axes = TRUE)

### For compound-compound
results_CxC <- C5.0(activity~., data = CxC_train)
actual <- CxC_test$activity
prediction_prob <- predict(results_CxC, CxC_test, type = "prob")[, 1]
prediction_class <- predict(results_CxC, CxC_test)
rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
plot.roc(rocCurve, col = "yellow", add = TRUE, legacy.axes = TRUE)


### trying with ROC Curve







data <- data.frame(results)
m = ncol(data)
ACC  <- matrix(nrow = m, ncol = 1)
SENS  <- matrix(nrow = m, ncol = 1)
SPEC  <-matrix(nrow = m, ncol = 1)
MCC <- matrix(nrow = m, ncol = 1)

for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4  =  sqrt(MCC2)*sqrt(MCC3)
  
  
  MCC[i,1]  = MCC1/MCC4
}

### ggplot ROC CURVE

#### ROC curve for the training data
### For compound
all <- list()
results_compound <- C5.0(activity~., data = compound_alone_train)
actual <- compound_alone_train$activity
pred <- prediction(predict(results_compound, newdata = compound_alone_train, type = "prob")[, 2], compound_alone_train$activity)
perform <- performance(pred, "tpr", "fpr")
toPlot <- data.frame(tpr = unlist(slot(perform, "y.values")), fpr = unlist(slot(perform, "x.values")))
#res <- ggplot(toPlot) + geom_line(aes(x=fpr, y=tpr)) + geom_abline(intercept=0, slope=1, colour="orange") + ylab("Sensitivity") + xlab("1 - Specificity") + ggtitle("")
#print(res)
all <- c(all, list(toPlot))
#rocCurve <- pROC::roc(actual, prediction_prob)
#plot.roc(rocCurve, col="blue", legacy.axes = TRUE, xlim = c(1,0), ylim = c(0,1))### Blue for compound
### For protein
results_protein <- C5.0(activity~., data = protein_alone_train)
actual <- protein_alone_train$activity
pred <- prediction(predict(results_protein, newdata = protein_alone_train, type = "prob")[, 2], protein_alone_train$activity)
perform <- performance(pred, "tpr", "fpr")
toPlot <- data.frame(tpr = unlist(slot(perform, "y.values")), fpr = unlist(slot(perform, "x.values")))
all <- c(all, list(toPlot))





#rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
#plot.roc(rocCurve, col = "red", add = TRUE, legacy.axes = TRUE)

## For crossterms
results_CxP <- C5.0(activity~., data = CxP_train)
actual <- CxP_train$activity
pred <- prediction(predict(results_CxP, newdata = CxP_train, type = "prob")[, 2], CxP_train$activity)
perform <- performance(pred, "tpr", "fpr")
toPlot <- data.frame(tpr = unlist(slot(perform, "y.values")), fpr = unlist(slot(perform, "x.values")))
all <- c(all, list(toPlot))
#rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
#plot.roc(rocCurve, col = "green", add = TRUE, legacy.axes = TRUE)

### for Protein-Protein
results_PxP <- C5.0(activity~., data = PxP_train)
actual <- PxP_train$activity
pred <- prediction(predict(results_PxP, newdata = PxP_train, type = "prob")[, 2], PxP_train$activity)
perform <- performance(pred, "tpr", "fpr")
toPlot <- data.frame(tpr = unlist(slot(perform, "y.values")), fpr = unlist(slot(perform, "x.values")))
all <- c(all, list(toPlot))
#rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
#plot.roc(rocCurve, col = "black", add = TRUE, legacy.axes = TRUE)

### For compound-compound
results_CxC <- C5.0(activity~., data = CxC_train)
actual <- CxC_train$activity
pred <- prediction(predict(results_CxC, newdata = CxC_train, type = "prob")[, 2], CxC_train$activity)
perform <- performance(pred, "tpr", "fpr")
toPlot <- data.frame(tpr = unlist(slot(perform, "y.values")), fpr = unlist(slot(perform, "x.values")))
all <- c(all, list(toPlot))


#rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
#plot.roc(rocCurve, col = "yellow", add = TRUE, legacy.axes = TRUE)

### Getting a Complete plot for training

names(all) <- c("C", "P", "CxP", "PxP", "CxC")
for(i in 1:5) {
  all[[i]] <- cbind(all[[i]], rep(names(all)[i], nrow(all[[i]])))
  colnames(all[[i]])[5] <- "Models"
}
for(i in 1:5) {
  all[[i]] <- cbind(all[[i]], rep(names(all)[i], nrow(all[[i]])))
  colnames(all[[i]])[5] <- "Models"
}

all <- do.call(rbind, all)
names(all) <- NULL
df <- data.frame(all)

C <- all$C
C <- data.frame(C, Model = c("C", "C", "C", "C", "C", "C"))
P <- all$P
P <- data.frame(P, Model = c("P", "P", "P", "P", "P", "P"))
CxP <- all$CxP
CxP <- data.frame(CxP, Model = c("CxP", "CxP", "CxP", "CxP"))
PxP <- all$PxP
PxP <- data.frame(PxP, Model = c("PxP", "PxP", "PxP", "PxP"))
CxC <- all$CxC
CxC <- data.frame(CxC, Model = c("CxC", "CxC", "CxC", "CxC", "CxC", "CxC"))

training_ROC <- ggplot(ROC_train) + geom_line(aes(x=fpr, y=tpr, colour=Model)) +
  geom_abline(intercept=0, slope=1, colour="orange") + 
  ylab("Sensitivity") + xlab("1 - Specificity") + ggtitle("")


### 10-fold Cross Validation
### Compound 
all <- list()
myData <- compound_alone_train
myData$activity <- factor(as.character(myData$activity))
folds <- createFolds(myData$activity, k = 10)
cv_results <- lapply(folds, function(x) {
  training <- myData[x, ]
  testing <- myData[-x, ]
  mymodel <- C5.0(activity~., data = testing)
  pred <- prediction(predict(mymodel, newdata = training, type = "prob")[, 2], training$activity)
  return(pred)
})

perf.auc <- performance(cv_results$Fold01, measure = "auc")
perf.auc
unlist(perf.auc@y.values)

perform <- performance(cv_results$Fold9, measure = "tpr", x.measure = "fpr")
plot_CV <- data.frame(tpr = unlist(slot(perform, "y.values")), fpr = unlist(slot(perform, "x.values")))

compound_ROC_data 
C <- data.frame(compound_ROC_data, Model = c("C", "C", "C"))
P_ROC_data 
P <- data.frame(P_ROC_data, Model = c("P", "P", "P"))
CxP_ROC_data
CxP <- data.frame(CxP_ROC_data, Model = c("CxP", "CxP", "CxP"))
PxP_ROC_data 
PxP <- data.frame(PxP_ROC_data, Model = c("PxP", "PxP", "PxP", "PxP"))
CxC_ROC_data 
CxC <- data.frame(CxC_ROC_data, Model = c("CxC", "CxC", "CxC"))






C <- all$C
C <- data.frame(C, Model = c("C", "C", "C", "C", "C", "C"))
P <- all$P
P <- data.frame(P, Model = c("P", "P", "P", "P", "P", "P"))
CxP <- all$CxP
CxP <- data.frame(CxP, Model = c("CxP", "CxP", "CxP", "CxP"))
PxP <- all$PxP
PxP <- data.frame(PxP, Model = c("PxP", "PxP", "PxP", "PxP"))
CxC <- all$CxC
CxC <- data.frame(CxC, Model = c("CxC", "CxC", "CxC", "CxC", "CxC", "CxC"))



ROC_10_CV <- rbind(C, P, CxP, PxP, CxC)
write.csv(ROC_10_CV, file = "ROC_10_CV.csv", row.names = FALSE)

ROC_train <- read.csv("ROC_train.csv", header = TRUE)
training_ROC <- ggplot(ROC_train) + geom_line(aes(x=fpr, y=tpr, colour=Model)) +
  geom_abline(intercept=0, slope=1, colour="orange", linetype = "dotted") + 
  ylab("Sensitivity") + xlab("1 - Specificity") + ggtitle("Training Set")
ROC_10_CV <- read.csv("ROC_10_CV.csv")
CV_ROC <- ggplot(ROC_10_CV) + geom_line(aes(x=fpr, y=tpr, colour=Model)) +
  geom_abline(intercept=0, slope=1, colour="orange", linetype = "dotted") + 
  ylab("Sensitivity") + xlab("1 - Specificity") + ggtitle("10-fold CV")

ROC_test <- read.csv("ROC_test.csv", header = TRUE)
testing_ROC <- ggplot(ROC_test) + geom_line(aes(x=fpr, y=tpr, colour=Model)) +
  geom_abline(intercept=0, slope=1, colour="orange", linetype = "dotted") + 
  ylab("Sensitivity") + xlab("1 - Specificity") + ggtitle("External Validation")

plot <- plot_grid(training_ROC, CV_ROC, testing_ROC, labels = c("A", "B", "C"), nrow = 3, ncol = 1)

grid.arrange(training_ROC, CV_ROC, testing_ROC, ncol=1)




#rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
#plot.roc(rocCurve, col="blue", legacy.axes = TRUE, xlim = c(1,0), ylim = c(0,1))
### Protein Alone


myData <- CxP_train
myData$activity <- factor(as.character(myData$activity))
folds <- createFolds(myData$activity, k = 9)
cv_results <- lapply(folds, function(x) {
  training <- myData[x, ]
  testing <- myData[-x, ]
  mymodel <- C5.0(activity~., data = testing)
  pred <- prediction(predict(mymodel, newdata = testing, type = "prob")[, 2], testing$activity)
  return(pred)
})

perf.auc <- performance(cv_results$Fold1, measure = "auc")
unlist(perf.auc@y.values)

perform <- performance(pred, "tpr", "fpr")
plot_CV <- data.frame(tpr = unlist(slot(perform, "y.values")), fpr = unlist(slot(perform, "x.values")))


ggplot(cv_results$Fold03) + geom_line(aes(x=fpr, y=tpr)) + 
  geom_abline(intercept=0, slope=1, colour="orange") + ylab("Sensitivity") + xlab("1 - Specificity") + 
  ggtitle("")


all_CV <- c(all_CV, list(myRes))

### 10-fold CV Cross-Terms
myData <- CxP_train
myData$activity <- factor(as.character(myData$activity))
myData <- compound_alone_train
myData$activity <- factor(as.character(myData$activity))
folds <- createFolds(myData$activity, k = 10)
cv_results <- lapply(folds, function(x) {
  training <- myData[x, ]
  testing <- myData[-x, ]
  mymodel <- C5.0(activity~., data = training)
  pred <- prediction(predict(mymodel, newdata = training, type = "prob")[, 2], training$activity)
  return(pred)
})

CxP_ROC_data <- cv_results$Fold08
ggplot(cv_results$Fold08) + geom_line(aes(x=fpr, y=tpr)) + 
  geom_abline(intercept=0, slope=1, colour="orange") + ylab("Sensitivity") + xlab("1 - Specificity") + 
  ggtitle("")

all <- c(all, list(myRes))

### 10-fold PxP
myData <- PxP_train
myData$activity <- factor(as.character(myData$activity))
myData <- compound_alone_train
myData$activity <- factor(as.character(myData$activity))
folds <- createFolds(myData$activity, k = 10)
cv_results <- lapply(folds, function(x) {
  training <- myData[x, ]
  testing <- myData[-x, ]
  mymodel <- C5.0(activity~., data = training)
  pred <- prediction(predict(mymodel, newdata = training, type = "prob")[, 2], training$activity)
  return(pred)
})


### 10-fold CxC
myData <- CxC_train
myData$activity <- factor(as.character(myData$activity))
myData <- compound_alone_train
myData$activity <- factor(as.character(myData$activity))
folds <- createFolds(myData$activity, k = 10)
cv_results <- lapply(folds, function(x) {
  training <- myData[x, ]
  testing <- myData[-x, ]
  mymodel <- C5.0(activity~., data = training)
  pred <- prediction(predict(mymodel, newdata = testing, type = "prob")[, 2], testing$activity)
  return(pred)
})


all_CV <- c(all_CV, list(myRes))

##
names(all_CV) <- c("C", "P", "CxP", "PxP", "CxC")


### testing Data
### For compound
all <- list()
results_compound <- C5.0(activity~., data = compound_alone_train)
actual <- compound_alone_test$activity
pred <- prediction(predict(results_compound, newdata = compound_alone_test, type = "prob")[, 2], compound_alone_test$activity)
perform <- performance(pred, "tpr", "fpr")
toPlot <- data.frame(tpr = unlist(slot(perform, "y.values")), fpr = unlist(slot(perform, "x.values")))
#res <- ggplot(toPlot) + geom_line(aes(x=fpr, y=tpr)) + geom_abline(intercept=0, slope=1, colour="orange") + ylab("Sensitivity") + xlab("1 - Specificity") + ggtitle("")
#print(res)
all <- c(all, list(toPlot))
#rocCurve <- pROC::roc(actual, prediction_prob)
#plot.roc(rocCurve, col="blue", legacy.axes = TRUE, xlim = c(1,0), ylim = c(0,1))### Blue for compound
### For protein
results_protein <- C5.0(activity~., data = protein_alone_train)
actual <- protein_alone_test$activity
pred <- prediction(predict(results_protein, newdata = protein_alone_test, type = "prob")[, 2], protein_alone_test$activity)
perform <- performance(pred, "tpr", "fpr")
toPlot <- data.frame(tpr = unlist(slot(perform, "y.values")), fpr = unlist(slot(perform, "x.values")))
all <- c(all, list(toPlot))





#rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
#plot.roc(rocCurve, col = "red", add = TRUE, legacy.axes = TRUE)

## For crossterms
results_CxP <- C5.0(activity~., data = CxP_train)
actual <- CxP_test$activity
pred <- prediction(predict(results_CxP, newdata = CxP_test, type = "prob")[, 2], CxP_test$activity)
perform <- performance(pred, "tpr", "fpr")
toPlot <- data.frame(tpr = unlist(slot(perform, "y.values")), fpr = unlist(slot(perform, "x.values")))
all <- c(all, list(toPlot))
#rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
#plot.roc(rocCurve, col = "green", add = TRUE, legacy.axes = TRUE)

### for Protein-Protein
results_PxP <- C5.0(activity~., data = PxP_train)
actual <- PxP_test$activity
pred <- prediction(predict(results_PxP, newdata = PxP_test, type = "prob")[, 2], PxP_test$activity)
perform <- performance(pred, "tpr", "fpr")
toPlot <- data.frame(tpr = unlist(slot(perform, "y.values")), fpr = unlist(slot(perform, "x.values")))
all <- c(all, list(toPlot))
#rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
#plot.roc(rocCurve, col = "black", add = TRUE, legacy.axes = TRUE)

### For compound-compound
results_CxC <- C5.0(activity~., data = CxC_train)
actual <- CxC_test$activity
pred <- prediction(predict(results_CxC, newdata = CxC_test, type = "prob")[, 2], CxC_test$activity)
perform <- performance(pred, "tpr", "fpr")
toPlot <- data.frame(tpr = unlist(slot(perform, "y.values")), fpr = unlist(slot(perform, "x.values")))
all <- c(all, list(toPlot))

names(all) <- c("C", "P", "CxP", "PxP", "CxC")

C <- all$C
C <- data.frame(C, Model = c("C", "C", "C", "C", "C", "C"))
P <- all$P
P <- data.frame(P, Model = c("P", "P", "P", "P", "P", "P"))
CxP <- all$CxP
CxP <- data.frame(CxP, Model = c("CxP", "CxP", "CxP", "CxP"))
PxP <- all$PxP
PxP <- data.frame(PxP, Model = c("PxP", "PxP", "PxP", "PxP"))
CxC <- all$CxC
CxC <- data.frame(CxC, Model = c("CxC", "CxC", "CxC", "CxC", "CxC", "CxC"))

ROC_test <- rbind(C, P, CxP, PxP, CxC)
write.csv(ROC_test, file = "ROC_test.csv", row.names = FALSE)

testing_ROC <- ggplot(ROC_test) + geom_line(aes(x=fpr, y=tpr, colour=Model)) +
  geom_abline(intercept=0, slope=1, colour="orange") + 
  ylab("Sensitivity") + xlab("1 - Specificity") + ggtitle("")

### ten fold Clasisfier

all <- list()
results_compound <- C5.0(activity~., data = compound_alone_train)
actual <- compound_alone_test$activity

pred <- prediction(predict(results_compound, newdata = compound_alone_test, type = "prob")[, 2], compound_alone_test$activity)
perform <- performance(pred, "tpr", "fpr")
toPlot <- data.frame(tpr = unlist(slot(perform, "y.values")), fpr = unlist(slot(perform, "x.values")))
#res <- ggplot(toPlot) + geom_line(aes(x=fpr, y=tpr)) + geom_abline(intercept=0, slope=1, colour="orange") + ylab("Sensitivity") + xlab("1 - Specificity") + ggtitle("")
#print(res)
all <- c(all, list(toPlot))
#rocCurve <- pROC::roc(actual, prediction_prob)
#plot.roc(rocCurve, col="blue", legacy.axes = TRUE, xlim = c(1,0), ylim = c(0,1))### Blue for compound
### For protein
results_protein <- C5.0(activity~., data = protein_alone_train)
actual <- protein_alone_test$activity
prediction <- predict(results_protein, newdata = protein_alone_test, type  = "prob")[, 2]
pred <- prediction(prediction, protein_alone_test$activity)
pred <- prediction(predict(results_protein, newdata = protein_alone_test, type = "prob")[, 2], protein_alone_test$activity)
perform <- performance(pred, "tpr", "fpr")
toPlot <- data.frame(tpr = unlist(slot(perform, "y.values")), fpr = unlist(slot(perform, "x.values")))
all <- c(all, list(toPlot))

perform <- performance(myRes, "tpr", "fpr")
actual <- myRes$Actual
prediction_prob <- myRes$Predicted
rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))




#rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
#plot.roc(rocCurve, col = "red", add = TRUE, legacy.axes = TRUE)

## For crossterms
results_CxP <- J48(activity~., data = CxP_train)
results_CxP <- evaluate_Weka_classifier(results_CxP, newdata = CxP_Train, numFolds = 10,
                                        class = TRUE, complexity = FALSE, seed = 1)
actual <- CxP_test$activity
pred <- prediction(predict(results_CxP, newdata = CxP_test, type = "prob")[, 2], CxP_test$activity)
perform <- performance(pred, "tpr", "fpr")
toPlot <- data.frame(tpr = unlist(slot(perform, "y.values")), fpr = unlist(slot(perform, "x.values")))
all <- c(all, list(toPlot))
#rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
#plot.roc(rocCurve, col = "green", add = TRUE, legacy.axes = TRUE)

### for Protein-Protein
results_PxP <- C5.0(activity~., data = PxP_train)
actual <- PxP_test$activity
pred <- prediction(predict(results_PxP, newdata = PxP_test, type = "prob")[, 2], PxP_test$activity)
perform <- performance(pred, "tpr", "fpr")
toPlot <- data.frame(tpr = unlist(slot(perform, "y.values")), fpr = unlist(slot(perform, "x.values")))
all <- c(all, list(toPlot))
#rocCurve <- pROC::roc(actual, prediction_prob, levels = rev(levels(actual)))
#plot.roc(rocCurve, col = "black", add = TRUE, legacy.axes = TRUE)

### For compound-compound
results_CxC <- C5.0(activity~., data = CxC_train)
actual <- CxC_test$activity
pred <- prediction(predict(results_CxC, newdata = CxC_test, type = "prob")[, 2], CxC_test$activity)
perform <- performance(pred, "tpr", "fpr")
toPlot <- data.frame(tpr = unlist(slot(perform, "y.values")), fpr = unlist(slot(perform, "x.values")))
all <- c(all, list(toPlot))

names(all) <- c("C", "P", "CxP", "PxP", "CxC")

C <- all$C
C <- data.frame(C, Model = c("C", "C", "C", "C", "C", "C"))
P <- all$P
P <- data.frame(P, Model = c("P", "P", "P", "P", "P", "P"))
CxP <- all$CxP
CxP <- data.frame(CxP, Model = c("CxP", "CxP", "CxP", "CxP"))
PxP <- all$PxP
PxP <- data.frame(PxP, Model = c("PxP", "PxP", "PxP", "PxP"))
CxC <- all$CxC
CxC <- data.frame(CxC, Model = c("CxC", "CxC", "CxC", "CxC", "CxC", "CxC"))

ROC_10_Fold <- rbind(C, P, CxP, PxP, CxC)
write.csv(ROC_test, file = "ROC_test.csv", row.names = FALSE)

testing_ROC <- ggplot(toPlot) + geom_line(aes(x=fpr, y=tpr, colour=Model)) +
  geom_abline(intercept=0, slope=1, colour="orange") + 
  ylab("Sensitivity") + xlab("1 - Specificity") + ggtitle("")


ggplot(toPlot) + geom_line(aes(x = fpr, y = tpr))








