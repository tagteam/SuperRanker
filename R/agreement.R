##' Quantify the agreement regarding common elements of a
##' list of lists.
##' 
##' @title Measure agreement between a set of lists
##' @param object Either matrix where each column is a list or a list of lists.
##' @param metric Function used to be applied to the ranks
##' @param depth Interest margin, i.e., agreement is limited to the first \code{depth} elements of each list.
##' @return Vector of lengths \code{depth} where element e contains the result
##' of metric applied to the ranks of the lists at position e.
##' @seealso 
##' @examples
##' agreement(list(c(1,2,3,4,5,7,6,8), c(1,2,3,5,6,4,7,8)))
##' agreement(list(c(1,2,3,4,5,7,6,8,9,10,11,12,13), c(5,11,4,7,8,3,12,13,6,10,9,2,1)),depth=4)
##' 
##' ## Simulate some data
##' library(lava)
##' m <- lvm()
##' nvar <- 1000
##' addvar(m) <- c("Y",paste("X",1:nvar,sep=""))
##' ## x <- matrix(rnorm(nobs*nvar), ncol=nvar)
##' ## nsig <- 1:5
##' regression(m,Y~X1+X2+X3+X4+X5) <- rnorm(5,0.5,0.5)
##' ## y <- rnorm(nobs, mean = 0*rowSums(x[,nsig]))
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
##' f1 <- rfsrc(Y~.,data=d,ntree=10)
##' l4 <- order(max.subtree(f1, conservative = TRUE,max.order=1)$order)
##' L4 <- max.subtree(f1,maxorder=, conservative = TRUE,max.order=0)$order
##' agreement(L4)
##' @export 
##' @author Claus Ekstrøm <ekstrom@@sund.ku.dk>, Thomas A. Gerds <tag@@biostat.ku.dk>
agreement <- function(object, metric=var, depth=NULL) {
    metric <- match.fun(metric)
    if (missing(depth)) depth <- length(list[[1]])
    if (is.matrix(object))
        rankmat <- object
    else
        rankmat <- as.matrix(do.call("cbind",lapply(object,function(l)l[1:depth])))
    # Compute rank for each item per list
    itemranks <- apply(rankmat, 2, order)
    rankmat.l <- lapply(seq_len(nrow(rankmat)), function(i) rankmat[i,])  
    maxdepth <- min(nrow(rankmat), ifelse(is.null(depth), nrow(rankmat), depth))
    uniq.l <- lapply(rankmat.l, function(x) {sort(unique(x))})
    res <- apply(itemranks, 1, metric)  
    agreement <- sapply(1:maxdepth, function(x) {
        myvar <- unique(unlist(uniq.l[1:x]))
        mean(res[myvar])
    })
    out <- list(agreement=agreement,metric=quote(metric),depth=maxdepth)
    class(out) <- "agreement"
    out
}
