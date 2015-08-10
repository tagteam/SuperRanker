library(SuperRanker)

rm(list=ls())
load("joinedOutput.RData")
load("DACOVAlabels.RData")
res <- res[1:1000]

B <- 1000 #number of randomization for censored lists

########## PREPARE OBJECTS ##########
rfRisk <- sra(lapply(res, function(a) a$rf$riskOrder))
plsRisk <- sra(lapply(res, function(a) a$pls$riskOrder))

lassoDevRisk <- sra(lapply(res, function(a) a$lassoDev$riskOrder))
lassoClassRisk <- sra(lapply(res, function(a) a$lassoClass$riskOrder))
lassoAucRisk <- sra(lapply(res, function(a) a$lassoAUC$riskOrder))

ridgeDevRisk <- sra(lapply(res, function(a) a$ridgeDev$riskOrder))
ridgeClassRisk <- sra(lapply(res, function(a) a$ridgeClass$riskOrder))
ridgeAucRisk <- sra(lapply(res, function(a) a$ridgeAUC$riskOrder))

##
rfImpSRA <- sra(lapply(res, function(a) a$rf$importance), nitems=5000, B=B)
plsImpSRA <- sra(lapply(res, function(a) a$pls$importance))

lassoDevImpSRA <- sra(lapply(res, function(a) a$lassoDev$importance), nitems=5000, B=B)
lassoClassImpSRA <- sra(lapply(res, function(a) a$lassoClass$importance), nitems=5000, B=B)
lassoAucImpSRA <- sra(lapply(res, function(a) a$lassoAUC$importance), nitems=5000, B=B)

ridgeDevImpSRA <- sra(lapply(res, function(a) a$ridgeDev$importance))
ridgeClassImpSRA <- sra(lapply(res, function(a) a$ridgeClass$importance))
ridgeAucImpSRA <- sra(lapply(res, function(a) a$ridgeAUC$importance))

##
aucMat <- cbind(sapply(res, function(a) a$rf$auc),
                sapply(res, function(a) a$pls$auc),
                sapply(res, function(a) a$lassoDev$auc),
                sapply(res, function(a) a$ridgeDev$auc),
                sapply(res, function(a) a$ridgeAUC$auc))
colnames(aucMat) <- c("Random Forest", "PLS-DA", "Lasso deviance", "Ridge deviance", "Ridge AUC")

brierMat <- cbind(sapply(res, function(a) mean((a$rf$risk - dacovaLabels)^2)),
                  sapply(res, function(a) mean((a$pls$risk - dacovaLabels)^2)),
                  sapply(res, function(a) mean((a$lassoDev$risk - dacovaLabels)^2)),
                  sapply(res, function(a) mean((a$ridgeDev$risk - dacovaLabels)^2)),
                  sapply(res, function(a) mean((a$ridgeAUC$risk - dacovaLabels)^2)))
colnames(brierMat) <- c("Random Forest", "PLS-DA", "Lasso deviance", "Ridge deviance", "Ridge AUC")

##
save(list=ls()[-which(ls() == "res")], file="plots.RData")


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
plot(rfImpSRA, ylim=c(0,1500), xlim=c(0,600))
title("Discrimination importance of 5000 predictors")

plot(lassoDevImpSRA, add=TRUE, lines.col="red")
plot(lassoClassImpSRA, add=TRUE, lines.col="red", lines.lty=2)
plot(lassoAucImpSRA, add=TRUE, lines.col="red", lines.lty=3)

plot(plsImpSRA, add=TRUE, lines.col="blue")

plot(ridgeDevImpSRA, add=TRUE, lines.col="forestgreen")
plot(ridgeClassImpSRA, add=TRUE, lines.col="forestgreen", lines.lty=2)
plot(ridgeAucImpSRA, add=TRUE, lines.col="forestgreen", lines.lty=3)


########## AUC PLOT ##########
boxplot(aucMat, ylab="AUC", ylim=c(0.4,0.8), axes=FALSE, main="AUC distribution of DACOVA predictions")
axis(2, at=seq(0.4,0.8,0.1))
axis(1, at=1:5, labels=colnames(aucMat), cex.axis=0.8)
abline(h=0.5, lty=3)


########## RIDGE CENSORING ##########
ridgeDevCens <- lapply(res, function(a) {
  b <- abs(a$ridgeDev$coef)
  b <- b[b > quantile(b, 0.001)]
  order(b, decreasing=TRUE)	
})
ridgeDevCencRisk <- sra(ridgeDevCens, nitems=5000, B=B)

plot(ridgeDevImpSRA, ylim=c(0,1500), xlim=c(0,600))
plot(ridgeDevCencRisk, add=TRUE, lines.col="red")
title("Discrimination importance for Ridge with thresholding")
legend("topright", c("No censoring", "Censored at 0.1% quantile"), col=c(1,2), lwd=3, lty=1, bty="n")



########## LASSO SUMMARY ##########
median(sapply(res, function(a) length(a$lassoDev$importance)))
range(sapply(res, function(a) length(a$lassoDev$importance)))

median(sapply(res, function(a) length(a$lassoClass$importance)))
range(sapply(res, function(a) length(a$lassoClass$importance)))

median(sapply(res, function(a) length(a$lassoAUC$importance)))
range(sapply(res, function(a) length(a$lassoAUC$importance)))


save(list=ls()[-which(ls() == "res")], file="plots.RData")


