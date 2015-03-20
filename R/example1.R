#
#
# Read data
#
#


source("sra.R")


# Read data
snps <- read.table("http://people.umass.edu/foulkes/asg/data/snpSimCVDNum.txt", header=TRUE)
traits <- read.table("http://people.umass.edu/foulkes/asg/data/ySimCVD.txt", header=TRUE)


y <- traits$yA  
x <- snps
index <- seq(1, NCOL(snps))
d <- data.frame(y, x)
#md <- mice(d)

# Marginal t-tests

mt.p <- sapply(index, function(i) { anova(lm(y ~ x[,i]))[1,5] } )
list1 <- order(mt.p)

# Elastic net

library(glmnet)
#library(mice)
#res <- coef(glmnet(as.matrix(x), y, alpha=0), s=lambda)[-1]

# MIC

library(minerva)
MIC <- sapply(index, function(i) { cc <- complete.cases(snps[,i],y) ; mine(x[cc,i], y[cc])$MIC})
list3 <- order(MIC, decreasing=TRUE)

## Forest
library(randomForestSRC)
f1 <- rfsrc(y~., data=d, ntree=1000)
variables <- abs(f1$importance)
num.undecided <- sum(variables==0)
list4 <- order(variables, decreasing=TRUE)
list4[(length(variables)-num.undecided):length(variables)] <- 0


res <- sra(list(list1, list3, list4))
plot(sqrt(res$agreement),xlim=c(0,20))




####
####
####


nulleval <- sapply(1:2, function(i) {
    py <- sample(y)

    #1
    mt.p <- sapply(index, function(i) { anova(lm(py ~ x[,i]))[1,5] } )
    list1 <- order(mt.p)

    #3
    MIC <- sapply(index, function(i) { cc <- complete.cases(snps[,i],py) ; mine(x[cc,i], py[cc])$MIC})
    list3 <- order(MIC, decreasing=TRUE)

    #4
    dd <- d
    dd$y <-py
    f1 <- rfsrc(y~., data=dd, ntree=1000)
    variables <- abs(f1$importance)
    num.undecided <- sum(variables==0)
    list4 <- order(variables, decreasing=TRUE)
    list4[(length(variables)-num.undecided):length(variables)] <- 0

    res <- sra(list(list1, list3, list4))$agreement
    res
})

nulleval1 <- sapply(1:2, function(i) {
    nobs <- length(list1)
    res <- sra(list(rep(0, nobs), rep(0, nobs), rep(0, nobs)), listlength=nobs)$agreement
})

