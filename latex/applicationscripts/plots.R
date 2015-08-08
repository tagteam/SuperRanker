library(SuperRanker)

rm(list=ls())
load("joinedOutput.RData")
load("DACOVAlabels.RData")

B <- 2 #number of randomization for censored lists


########## PREPARE SRA OBJECTS ##########
rfRisk <- sra(lapply(res, function(a) a$rf$riskOrder))
plsRisk <- sra(lapply(res, function(a) a$pls$riskOrder))

lassoDevRisk <- sra(lapply(res, function(a) a$lassoDev$riskOrder))
lassoClassRisk <- sra(lapply(res, function(a) a$lassoClass$riskOrder))
lassoAucRisk <- sra(lapply(res, function(a) a$lassoAUC$riskOrder))

ridgeDevRisk <- sra(lapply(res, function(a) a$ridgeDev$riskOrder))
ridgeClassRisk <- sra(lapply(res, function(a) a$ridgeClass$riskOrder))
ridgeAucRisk <- sra(lapply(res, function(a) a$ridgeAUC$riskOrder))


rfImpSRA <- sra(lapply(res, function(a) a$rf$importance), nitems=5000, B=B)
plsImpSRA <- sra(lapply(res, function(a) a$pls$importance))

lassoDevImpSRA <- sra(lapply(res, function(a) a$lassoDev$importance), nitems=5000, B=B)
lassoClassImpSRA <- sra(lapply(res, function(a) a$lassoClass$importance), nitems=5000, B=B)
lassoAucImpSRA <- sra(lapply(res, function(a) a$lassoAUC$importance), nitems=5000, B=B)

ridgeDevImpSRA <- sra(lapply(res, function(a) a$ridgeDev$importance))
ridgeClassImpSRA <- sra(lapply(res, function(a) a$ridgeClass$importance))
ridgeAucImpSRA <- sra(lapply(res, function(a) a$ridgeAUC$importance))


########## RISK PLOT ##########
plot(rfRisk, ylim=c(0,30), xlim=c(0,120))
title("Predicted risks of malignant tumor in 113 patients (DACOVA)")

plot(lassoDevRisk, add=TRUE, lines.col="red")
plot(lassoClassRisk, add=TRUE, lines.col="red", lines.lty=2)
plot(lassoAucRisk, add=TRUE, lines.col="red", lines.lty=3)

plot(plsRisk, add=TRUE, lines.col="blue")

plot(ridgeDevRisk, add=TRUE, lines.col="forestgreen")
plot(ridgeClassRisk, add=TRUE, lines.col="forestgreen", lines.lty=2)
plot(ridgeAucRisk, add=TRUE, lines.col="forestgreen", lines.lty=3)

legend("topright", c("Random Forest", "PLS-DA", "Lasso deviance", "Lasso class", "Lasso AUC", "Ridge deviance", "Ridge class", "Ridge AUC"),
       col=c("black", "blue", "red", "red", "red", "forestgreen", "forestgreen", "forestgreen"), 
       lwd=3, lty=c(1,1,1,2,3,1,2,3), bty="n")


########## IMPORTANCE PLOT ##########
plot(rfImpSRA, ylim=c(0,2000), xlim=c(0,600))
title("Discrimination importance of 5000 predictors")

plot(lassoDevImpSRA, add=TRUE, lines.col="red")
plot(lassoClassImpSRA, add=TRUE, lines.col="red", lines.lty=2)
plot(lassoAucImpSRA, add=TRUE, lines.col="red", lines.lty=3)

plot(plsImpSRA, add=TRUE, lines.col="blue")

plot(ridgeDevImpSRA, add=TRUE, lines.col="forestgreen")
plot(ridgeClassImpSRA, add=TRUE, lines.col="forestgreen", lines.lty=2)
plot(ridgeAucImpSRA, add=TRUE, lines.col="forestgreen", lines.lty=3)


########## RIDGE CENSORING ##########
cutOff <- quantile(abs(unlist(out)), 0.001)
outOrderCens <- lapply(out, function(a) {
  b <- abs(a)
  b <- b[b > cutOff]
  order(b, decreasing=TRUE)	
})

median(sapply(outOrderCens, function(a) 5000 - length(a))) #median number of censorings
range(sapply(outOrderCens, function(a) 5000 - length(a)))



########## AUC + BRIER PLOT ##########
aucMat <- cbind(sapply(res, function(a) a$rf$auc),
                sapply(res, function(a) a$pls$auc),
                sapply(res, function(a) a$lassoDev$auc),
                sapply(res, function(a) a$ridgeDev$auc),
                sapply(res, function(a) a$ridgeAUC$auc))
colnames(aucMat) <- c("Random Forest", "PLS-DA", "Lasso deviance", "Ridge deviance", "Ridge AUC")
boxplot(aucMat, ylab="AUC", ylim=c(0.4,0.8), axes=FALSE, main="AUC distribution of DACOVA predictions")
axis(2, at=seq(0.4,0.8,0.1))
axis(1, at=1:5, labels=colnames(aucMat), cex.axis=0.8)
abline(h=0.5, lty=3)



brierScaled <- function(label, pred) {
  brier <- mean((label - pred)^2)
  brier
}


########## LASSO SUMMARY ##########
median(sapply(res, function(a) length(a$lassoDev$importance)))
range(sapply(res, function(a) length(a$lassoDev$importance)))

median(sapply(res, function(a) length(a$lassoClass$importance)))
range(sapply(res, function(a) length(a$lassoClass$importance)))

median(sapply(res, function(a) length(a$lassoAUC$importance)))
range(sapply(res, function(a) length(a$lassoAUC$importance)))


save(list=ls()[-which(ls() == "res")], file="plots.RData")