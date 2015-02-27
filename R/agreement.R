##' Quantify the agreement regarding common elements of a
##' list of lists.
##' 
##' @title Measure agreement between a set of lists
##' @param object Either matrix where each column is a list or a list of lists. Elements are integers between 1 and the length of the lists. The lists need not be the same length, but the shorter lists need to be padded with NAs
##' @param metric Function used to evaluate the agreement among the ranks
##' @param depth Interest margin, i.e., measurement of agreement is limited to the first \code{depth} elements of each list. Useful for long lists where only the top is of interest.
##' @param B Number of resamples to use for missing items
##' @return Vector of lengths \code{depth} where element e contains the result
##' of metric applied to the ranks of the lists at position e.
##' @examples
##' agreement(list(c(1,2,3,4,5,7,6,8), c(1,2,3,5,6,4,7,8)))
##' agreement(list(c(1,2,3,4,5,7,6,8,9,10,11,12,13), c(5,11,4,7,8,3,12,13,6,10,9,2,1)),depth=13,metric=sd)
##' agreement(list(c(1,2,3,4,5,7,6,8,9,10,11,12,13), c(5,11,4,7,8,3,12,13,6,10,9,2,1)),depth=4,metric=sd)
##' 
##' ## Simulate some data
##' library(lava)
##' m <- lvm()
##' nvar <- 100
##' addvar(m) <- c("Y",paste("X",1:nvar,sep=""))
##' regression(m,Y~X1+X2+X3+X4+X5) <- rnorm(5,0.5,0.5)
##' nobs <- 100
##' d <- sim(m,nobs)
##' 
##' ## Run analyses
##' p1 <- unlist(lapply(names(d[,-1]), function(xi){summary(lm(formula(paste("Y ~ ",xi)),d))$coefficients[2,4]}))
##' l1 <- order(p1)
##' lambda <- .3
##' xx <- scale(d[,-1])
##' y <- d[,1]
##' library(glmnet)
##' # Lasso
##' res  <- coef(glmnet(xx,y), s=lambda)[-1]
##' res2 <- coef(glmnet(xx,y, alpha=0), s=lambda)[-1]
##' l2 <- order(abs(res), decreasing=TRUE)
##' l3 <- order(abs(res2), decreasing=TRUE)
##' l3 <- order(abs(coef(glmnet(xx,y), s=.1)[-1]), decreasing=TRUE)
##' z1 <- agreement(list(l1, l2), depth=20)
##' plot(z1)
##' 
##' # Forest
##' library(randomForestSRC)
##' f1 <- rfsrc(Y~.,data=d,ntree=1000)
##' z2 <- agreement(list(order(abs(f1$importance),decreasing=TRUE),l2))
##' plot(z2,xlim=c(0,20))
##' 
##' ## f1 <- rfsrc(Y~X1+X2+X77+X79.,data=d,ntree=10)
##' l4 <- order(max.subtree(f1, conservative = TRUE,max.order=1)$order)
##' L4 <- max.subtree(f1,maxorder=, conservative = TRUE,max.order=0)$order
##' L4a <- apply(jitter(L4),2,order)
##' agreement(L4a)
##' agreement(L4a[,1:4])
##'
##' boot.list <- lapply(1:10,function(b){
##' d.b <- d[sample(1:NROW(d),replace=TRUE),]
##' xx.b <- scale(d.b[,-1])
##' y.b <- d.b[,1]
##' res.b <- coef(glmnet(xx.b,y.b), s=lambda)[-1]
##' order(abs(res.b), decreasing=TRUE)
##' })
##' a.re <- agreement(boot.list)
##' sub.list <- lapply(1:10,function(b){
##' d.b <- d[sample(1:NROW(d),replace=FALSE,size=80),]
##' xx.b <- scale(d.b[,-1])
##' y.b <- d.b[,1]
##' res.b <- coef(glmnet(xx.b,y.b), s=lambda)[-1]
##' order(abs(res.b), decreasing=TRUE)
##' })
##' a.sub <- agreement(sub.list)
##' plot(a.re)
##' plot(a.sub,add=TRUE,lines.col=2)
##' boot.forest <- lapply(1:10,function(b){
##' d.b <- d[sample(1:NROW(d),replace=TRUE),]
##' f1 <- rfsrc(Y~.,data=d.b,ntree=100,nsplit=10)
##' order(abs(f1$importance),decreasing=TRUE)
##' })
##' forest.re <- agreement(boot.forest)
##' sub.forest <- lapply(1:10,function(b){
##' d.b <- d[sample(1:NROW(d),replace=FALSE,size=80),]
##' f1 <- rfsrc(Y~.,data=d.b,ntree=100,nsplit=10)
##' order(abs(f1$importance),decreasing=TRUE)
##' })
##' forest.sub <- agreement(sub.forest)
##' plot(forest.re,xlim=c(0,20))
##' plot(forest.sub,add=TRUE,lines.col=2)
##'
##' 
##' @export 
##' @author Claus Ekstrøm <ekstrom@@sund.ku.dk>, Thomas A. Gerds <tag@@biostat.ku.dk>
agreement <- function(object, metric=var, depth=NULL, listlength=NULL, B=10) {
    metric <- match.fun(metric)
    
    # Make sure that the input object ends up as a matrix with integer columns all
    # consisting of elements from 1 and up to listlength
    if (is.matrix(object))
        rankmat <- object
    else
        rankmat <- as.matrix(do.call("cbind",object))
    
    if (missing(depth)) depth <- NROW(rankmat)
    
    # Now do sanity check of input
    # Check if there are dupliates in the columns
    if (any(apply(rankmat, 2, duplicated, incomparables=c(0, NA)))) {
      stop("an element must only occur exactly once in each list")
    }
    # Check and compute the maximum list length
    ll <- max(c(NROW(rankmat), max(rankmat, na.rm=TRUE)))
    if (is.null(listlength)) {
      listlength <- ll
    }    
    if (listlength<ll) {
      stop("You are requesting a list length that is smaller than the largest element found")
    }    
    # Should probably also check integers and such
        
    # Compute the full metric for each item in a full sized list
    # The speed here could be increased (for example by converting everything to itemlist scale first)
    full.list.metric <- function(m, metric) {
      itemmatrix <- apply(m, 2, order)
      apply(itemmatrix, 1, metric)
    }
    resample <- function(x, ...) x[sample.int(length(x), ...)]

    
    # General approach
    # 1) Fill in missing items in each list
    # 2) Compute the metric for the full lists
    # 3) Compute agreement 
    # 4) Finally average everything

    nseq <- seq(listlength)
    nlists <- NCOL(rankmat)
  

    # Compute a list of missing items for each list
    missing.items <- lapply(as.data.frame(rankmat), function(x) { nseq[-x] })    
    
    fullrankmat <- matrix(rep(0, nlists*listlength), ncol=nlists)
    fullrankmat[1:NROW(rankmat),] <- rankmat
      
    tmpres <- sapply(1:B, function(i) {
      # Ad 1)
      for (j in 1:nlists) {
        if (length(missing.items[[j]])>0) {      
          fullrankmat[(listlength-length(missing.items[[j]])+1):listlength,j] <- resample(missing.items[[j]])
        }
      }
      res <- full.list.metric(fullrankmat, metric)
      res
    })
    
    rankmat.l <- lapply(seq_len(nrow(rankmat)), function(i) rankmat[i,])  
    maxdepth <- min(nrow(rankmat), ifelse(is.null(depth), nrow(rankmat), depth))
    uniq.l <- lapply(rankmat.l, function(x) {sort(unique(x))})
    
    agreement <- sapply(1:maxdepth, function(x) {
      myvar <- unique(unlist(uniq.l[1:x]))
      mean(tmpres[myvar])
    })

    out <- list(agreement=agreement, metric=quote(metric), depth=maxdepth, B=B)
    class(out) <- "agreement"
    out
}
