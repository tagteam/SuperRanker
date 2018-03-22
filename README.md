# SuperRanker

Development version of the R package SuperRanker

To install the development version of SuperRanker run the following commands
from within R
```{r}
library(devtools)
install_github('TagTeam/SuperRanker')
```

## Example


```{r}
mlist <- matrix(cbind(1:8,c(1,2,3,5,6,7,4,8),c(1,5,3,4,2,8,7,6)),ncol=3)
sra(mlist)
```


### Comparisons between sra and the topK method

See [this script](https://www.github.com) for a comparison between sequential rank agreement the topK method. This corresponds to section 5 of the manuscript


### Long example

The following code produces the results for the Table and Figure 1 in the manuscript. They do take some time to run. This reproduces the results in section 3 of the paper.


```{r}
library(SuperRanker)  # For the SuperRanker package
library(glmnet)       # For elastic net
library(minerva)      # For MIC
library(parallel)     # To speed up computations
library(multtest)     # For the golub data
library(MESS)         # For transparent colours

# Read data
data(golub)

y <- golub.cl
x <- golub


## This function runs the 4 analyses
## For a given design matrix (predictor set) X, and an outcome Y
producelists <- function(x, y) {
    nitems <- nrow(x)
    index <- seq(1, nrow(golub))
    ## d <- data.frame(y, x)

    ## Marginale t-tests
    mt.p <- sapply(index, function(i) { t.test(x[i,] ~ y)$p.value } )
    list1 <- order(mt.p)

    ## Marginale logreg-tests
    mlogreg.p <- sapply(index, function(i) { drop1(glm(y ~ x[i,], family=binomial), test="Chisq")[2,5] } )
    list2 <- order(mlogreg.p)

    ## Elastic net
    X <- scale(t(x))
    enet <- glmnet(X, y, family="binomial", alpha=.3)
    nyres <- cv.glmnet(X, y, family="binomial", alpha=.3)
    coefficients <- coef(enet, s=nyres$lambda.min)[-1]
    nonzeros <- sum(coefficients!=0)
    list3 <- order(abs(coefficients), decreasing=TRUE)
    if (nonzeros<nitems) {
        list3[(nonzeros+1):nitems] <- NA
    }
    ## MIC
    MIC <- sapply(index, function(i) { mine(x[i,], y)$MIC})
    list4 <- order(MIC, decreasing=TRUE)

    cbind(list1,list2,list3,list4)
}

## Original run
set.seed(1)
epsilon <- 0
B <- 200  # Number of resamples to use for censored lists. 
          # This comes into play because of elastic net


## inputmatrix corresponds to the table 
inputmatrix <- producelists(x, y)
colnames(inputmatrix) <- c("T", "LogReg", "ElasticNet", "MIC", "RF")[1:ncol(inputmatrix)]
res <- sra(inputmatrix, B=B, epsilon=epsilon)   ## Result from original analyses

## Now make the same run but where we only have top 20 lists
topk <- 20
inputmatrix2 <- inputmatrix
inputmatrix2[(topk+1):nrow(inputmatrix),] <- NA
res2 <- sra(inputmatrix2, B=B, epsilon=epsilon)

#pdf(file="newfig1a.pdf")
par(mar=c(5,5,0,0)+.1)
ysize <- length(res)
## Plot the line
plot(res[1:ysize], col="black", lwd=3, type="l", 
     ylim=c(1,1300), ylab="Sequential rank agreement", 
     xlab="List depth", bty="n", log="xy", cex.axis=1.7, cex.lab=1.7)

## Compute random list backgrounds
nreference <- 100  ## Number used
null <- random_list_sra(inputmatrix, B=B, n=nreference, epsilon=0.25)
# save(null, file="R/fig1null.rda")
# load("R/fig1null.rda")

bcolor <- MESS::col.alpha("red", alpha=120/255)
bcolor2 <- MESS::col.alpha("blue", alpha=120/255)

## Plot the random lists
www <- smooth_sra(null)
polygon(c(1:length(www$lower), rev(1:length(www$lower))),
        c(www$lower, rev(www$upper)),
        col=bcolor2, border=NA)

## Load the data for the methods
makereference <- function(x, y, n=100, B=B) {
    mclapply(1:n, function(i) { sra(producelists(x, sample(y)), B=B, epsilon=epsilon) }, 
             mc.cores=3)
}
## Uncomment the two lines below to redo the analysis
null2 <- makereference(x, y, n=100, B=B)
#save(null2, file="R/fig1null2.rda")
#load("~/ku/R/SuperRanker/latex/R/fig1null2.rda")
ttt <- do.call("cbind", null2)

## Plot the reference sras
www2 <- smooth_sra(ttt)
polygon(c(1:length(www2$lower), rev(1:length(www2$lower))),
        c(www2$lower, rev(www2$upper)),
        col=bcolor, border=NA)

# dev.off()


#pdf(file="newfig1b.pdf")
#par(mar=c(5,5,0,0)+.1)

ysize <- length(res)
nreference <- 20
## Plot the line
plot(res2[1:ysize], col="black", lwd=3, type="l", ylim=c(1,1300), 
     ylab="Sequential rank agreement", xlab="List depth", bty="n", log="xy",
     cex.axis=1.7, cex.lab=1.7)

## Compute random list backgrounds
null <- random_list_sra(inputmatrix2, B=B, n=nreference)
#save(null, file="R/fig1nullc.rda")
#load("R/fig1nullc.rda")


## Plot the random lists
www <- smooth_sra(null)
polygon(c(1:length(www$lower), rev(1:length(www$lower))),
        c(www$lower, rev(www$upper)),
        col=bcolor2, border=NA)

## Load the data for the methods
makereference <- function(x, y, n=100, B=20) {
    mclapply(1:n, function(i) { mymat <- producelists(x, sample(y)) ; mymat[(topk+1):nrow(mymat),] <- NA ; sra(mymat, B=B, epsilon=epsilon) }, mc.cores=3  )
}
## Uncomment the two lines below to redo the analysis
null2 <- makereference(x, y, n=100, B=B)
#save(null2, file="R/fig1null2c.rda")
#load("R/fig1null2c.rda")


ttt <- do.call("cbind", null2)

## Plot the reference sras
www2 <- smooth_sra(ttt)
polygon(c(1:length(www2$lower), rev(1:length(www2$lower))),
        c(www2$lower, rev(www2$upper)),
        col=bcolor, border=NA)

```
