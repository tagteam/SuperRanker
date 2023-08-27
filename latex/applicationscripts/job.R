library(randomForest)
library(caret)
library(glmnet)
library(pROC)

#Performs random 5-folds B-repeated CV and returns optimal lambda
glmnetWrapper <- function(x, y, alpha, B, type) {
  cvOut <- lapply(1:B, function(b) {
    #cat(b, "\n")
    cv.glmnet(x, y, family="binomial", type.measure=type, nfolds=5, alpha=alpha, nlambda=500)
  })
  #Do the hack again
  #sapply(cvOut, function(a) length(a$lambda))
  lambdaSeq <- Reduce(intersect, lapply(cvOut, function(a) a$lambda))
  cvMat <- t(sapply(cvOut, function(a) a$cvm[a$lambda %in% lambdaSeq]))
  lambdaSeq[which.min(colMeans(cvMat))]
}


doglmnet <- function(x, y, xValid=xOle, yValid=yOle, alpha, type) {
  lambda <- glmnetWrapper(x, y, alpha=alpha, B=20, type=type)
  mod <- glmnet(x, y, family="binomial", alpha=alpha, lambda=lambda)
  
  risk <- predict(mod, xValid, type="response")
  riskOrder <- order(risk, decreasing=TRUE)
  AUC <- as.numeric(auc(yValid, risk[,1]))

  coef <- coef(mod)[-1]
  if(alpha == 1) { #censor when lasso
    order <- order(abs(coef), decreasing=TRUE)[1:sum(coef != 0)] #censor when coef = 0
  } else { #no censoring when ridge
    order <- order(abs(coef), decreasing=TRUE)
  }

  list(riskOrder=riskOrder, importance=order, risk=risk[,1], coef=as.numeric(coef), auc=AUC, lambda=lambda)  
}


tctl <- trainControl(method="repeatedcv", number=5, repeats=20, verboseIter=TRUE)

dopls <- function(x, y, xValid=xOle, yValid=yOle, ncomp=10) {
  datTrain <- data.frame(y=factor(c("benign", "malign")[y+1]), x=x)
  datValid <- data.frame(x=xValid)

  plsMod <- train(y ~ ., data=datTrain, method="pls", trControl=tctl, tuneGrid=data.frame(ncomp=1:ncomp))
  plsncomp <- as.numeric(plsMod$bestTune)
  
  plsRisk <- predict(plsMod, newdata=datValid, type="prob")[, "malign"]
  plsRiskOrder <- order(plsRisk, decreasing=TRUE)
  plsAUC <- as.numeric(auc(yValid, plsRisk))
  
  plsImp <- varImp(plsMod)$importance[,"Overall"]
  plsOrder <- order(plsImp, decreasing=TRUE)

  list(riskOrder=plsRiskOrder, importance=plsOrder, risk=plsRisk, imp=plsImp, auc=plsAUC, ncomp=plsncomp)	
}


dorf <- function(x, y, xValid=xOle, yValid=yOle, ntree=5000) {
  rfTune <- tuneRF(x, factor(y), ntreeTry=1000, trace=FALSE, plot=FALSE)
  rfmtry <- rfTune[which.min(rfTune[,2]), 1]
  rfMod <- randomForest(x, factor(y), ntree=ntree, mtry=rfmtry, importance=TRUE)

  rfRisk <- predict(rfMod, newdata=xValid, type="prob")[, "1"]
  rfRiskOrder <- order(rfRisk, decreasing=TRUE)
  rfAUC <- as.numeric(auc(yValid, rfRisk))

  rfImp <- importance(rfMod, type=2)
  rfOrder <- order(rfImp, decreasing=TRUE)[1:sum(rfImp != 0)] #censor when importance = 0

  list(riskOrder=rfRiskOrder, importance=rfOrder, risk=rfRisk, imp=as.vector(rfImp), auc=rfAUC, mtry=rfmtry)
}



#################
currentIndex <- as.numeric(Sys.getenv("SGE_TASK_ID")) # Index of current process
processTotal <- as.numeric(Sys.getenv("SGE_TASK_LAST")) # Total number of processes
jobname <- Sys.getenv("JOB_NAME")

seed.in <- 1565745
set.seed(seed.in)
seeed <- sample(1:1000000,size=processTotal,replace=FALSE) # Generate array of seeds
seed <- seeed[currentIndex]


#################
load("preprocessed.RData")
xMal <- scale(xMal)
xOle <- scale(xOle)

cat("Seed = ", seed, "\n")
set.seed(seed)

subI <- sample(1:nrow(xMal), round(0.9*nrow(xMal))) #subsample random 90% for training
xTrain <- xMal[subI,]
yTrain <- yMal[subI]

  
#################
rfRes <- dorf(xTrain, yTrain)

plsRes <- dopls(xTrain, yTrain)

lassoDev <- doglmnet(xTrain, yTrain, alpha=1, type="deviance")

lassoClass <- doglmnet(xTrain, yTrain, alpha=1, type="class")

lassoAUC <- doglmnet(xTrain, yTrain, alpha=1, type="auc")

ridgeDev <- doglmnet(xTrain, yTrain, alpha=0, type="deviance")

ridgeClass <- doglmnet(xTrain, yTrain, alpha=0, type="class")

ridgeAUC <- doglmnet(xTrain, yTrain, alpha=0, type="auc")



#################
out <- list(rf=rfRes, pls=plsRes,
            lassoDev=lassoDev, lassoClass=lassoClass, lassoAUC=lassoAUC,
            ridgeDev=ridgeDev, ridgeClass=ridgeClass, ridgeAUC=ridgeAUC)

filename <- paste('output/', jobname,'-', currentIndex,'.RData',sep="")
save(out,file=filename)
