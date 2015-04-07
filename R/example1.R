#
#
# Read data
#
#
library(glmnet)
library(minerva)
library(changepoint)
library(randomForestSRC)

source("sra.R")

# Read data
library(multtest)
data(golub)


y <- golub.cl
x <- golub 


vs.vimp.md <- function(obj, which.outcome, conservative = TRUE) {
   msub <- max.subtree(obj, conservative = conservative)
   vimp <- obj$importance
   if (is.null(vimp)) {
       vimp <- matrix(1, nrow(msub$order), ncol(obj$xvar))
   }
   names <- rownames(msub$order)
   names[vimp[, which.outcome] > 0 & msub$order[, 1] <= msub$threshold]    
} 

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
    dd <- data.frame(y=factor(y), t(x))
                                        # dd <- data.frame(y=y, t(x))
    f1 <- rfsrc(y ~ ., data=dd, ntree=100)
    variables <- abs(f1$importance[,1])
    num.undecided <- sum(variables==0)
    list5 <- order(variables, decreasing=TRUE)
    list5[(length(variables)-num.undecided):length(variables)] <- 0
    
    cbind(list1,list2,list3,list4,list5)
}

inputmatrix <- producelists(x, y)


res <- sra(inputmatrix)

plot(sqrt(res$agreement)[1:200], lwd=3, col="red", type="l")


####
####
####

makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

random.list <- function(x, nsim=100) {
    npreds <- nrow(x)

    res <- sapply(1:nsim, function(i) {
        
        xx <- sapply(1:ncol(x), function(ii) { sample(1:npreds) } )
        sra(xx, B=1, listlength=nrow(x))$agreement
    } )
    res
}

newcol <- makeTransparent

null1 <- random.list(inputmatrix)

matlines(null1, type="l", col="lightgray")


nsims <- 20

nulleval <- sapply(1:nsims, function(i) {
    py <- sample(y)

    inputmatrix <- producelists(x, py)
    
    res <- sra(inputmatrix)$agreement
    res
})


    
