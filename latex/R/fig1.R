library(SuperRanker)
library(glmnet)
library(minerva)
library(changepoint)
library(randomForestSRC)

# Read data
library(multtest)
data(golub)


y <- golub.cl
x <- golub 

producelists <- function(x, y) {
    index <- seq(1, nrow(golub))
                                        # d <- data.frame(y, x)

                                        # Marginale t-tests
    mt.p <- sapply(index, function(i) { t.test(x[i,] ~ y)$p.value } )
    list1 <- order(mt.p)
    
                                        # Marginale logreg-tests
    mlogreg.p <- sapply(index, function(i) { drop1(glm(y ~ x[i,], family=binomial), test="Chisq")[2,5] } )
    list2 <- order(mlogreg.p)
    
                                        # Elastic net
    X <- scale(t(x))
    enet <- glmnet(X, y, family="binomial", alpha=.8)
    nyres <- cv.glmnet(X, y, family="binomial", alpha=.8)
    coefficients <- coef(enet, s=nyres$lambda.1se)[-1]
    nonzeros <- sum(coefficients!=0)
    list3 <- order(abs(coefficients), decreasing=TRUE)
    
                                        # MIC
    MIC <- sapply(index, function(i) { mine(x[i,], y)$MIC})
    list4 <- order(MIC, decreasing=TRUE)
    
    
                                        # Random Forest
###    dd <- data.frame(y=factor(y), t(x))
                                        # dd <- data.frame(y=y, t(x))
###    f1 <- rfsrc(y ~ ., data=dd, ntree=100)
###    variables <- abs(f1$importance[,1])
###    num.undecided <- sum(variables==0)
###    list5 <- order(variables, decreasing=TRUE)
###    list5[(length(variables)-num.undecided):length(variables)] <- 0
    
    cbind(list1,list2,list3,list4)
}

inputmatrix <- producelists(x, y)
colnames(inputmatrix) <- c("T", "LogReg", "ElasticNet", "MIC", "RF")[1:ncol(inputmatrix)]

res <- sqrt(sra(inputmatrix))


ysize <- 1000
plot(res[1:ysize], lwd=3, col="red", type="l", ylim=c(0,1300), ylab="Sequential rank agreement")

## 
# col1 <- makeTransparent()
null <- sqrt(random_list_sra(inputmatrix, B=1, n=100))
matlines(null[1:ysize,], col="gray", lty=1)

## 

null2 <- sapply(1:400, function(i) {
    sqrt(sra(producelists(x, sample(y)), B=1))
})
# matlines(null2[1:ysize,], col="lightblue", lty=1)

makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

bcolor <- makeTransparent("lightblue", alpha=80)


www <- smooth_sra(null2)
polygon(c(1:length(www$lower), rev(1:length(www$lower))),
        c(www$lower, rev(www$upper)),
        col=bcolor, border=NA)

xxx <- as.data.frame(inputmatrix[1:10,])
xxx

library(xtable)
print(xtable(xxx))
#print(xtable(xxx, caption = "Top 10 list of results", label = "tab:one",
#digits = c(0, 0, 2, 0, 2, 3, 3)), table.placement = "tbp",
#caption.placement = "top")
