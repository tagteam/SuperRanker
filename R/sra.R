#' Compute the sequential rank agreement
#'
#' @param object Either matrix where each column is a ranked list of
#' items or a list of ranked lists of items. Elements are integers
#' between 1 and the length of the lists. The lists should have the
#' same length but censoring can be used by setting the list to zero
#' from a point onwards. See details for more information.
#' @param na.strings Code for missing items 
#' @param B An integer giving the number of randomization to sample
#' over in the case of censored observations
#' @return A vector of the sequential rank agreement
##' @examples
##' 
##' alist <- list(a=1:8,b=sample(1:8),c=sample(1:8))
##' sra(alist)
##' 
##' blist <- list(x1=letters,x2=sample(letters),x3=sample(letters))
##' sra(blist)
##' 
##' ## censored lists are either too short
##' clist <- list(x1=c("a","b","c","d","e","f","g","h"),
##'               x2=c("h","c","f","g","b"),
##'               x3=c("d","e","a"))
##' set.seed(17)
##' sra(clist,na.strings="z",B=10)
##' 
##' ## or use a special code for missing elements
##' Clist <- list(x1=c("a","b","c","d","e","f","g","h"),
##'               x2=c("h","c","f","g","b","z","z","z"),
##'               x3=c("d","e","a","z","z","z","z","z"))
##' set.seed(17)
##' sra(Clist,na.strings="z",B=10)
##' 
#' @author Claus Ekstrøm <ekstrom@@sund.ku.dk> and Thomas A Gerds <tag@@biostat.ku.dk>
#' 
#' @export
sra <- function(object, na.strings=NA, B=1) {
    # Make sure that the input object ends up as a matrix with integer columns all
    # consisting of elements from 1 and up to listlength
    if (is.matrix(object))
        object <- apply(object,2,function(x)x)
    else
        stopifnot(is.list(object)) # data.frame is a list

    nlists <- length(object)

    ## sanity checks
    object <- lapply(object,function(x){
                         out <- which(is.na(x))
                         na.strings <- na.strings[!is.na(na.strings)]
                         if (length(na.strings)>0)
                             out <- c(out,grep(paste0("^",na.strings,"$"),x))
                         ## remove censored items with side effect:
                         ## in case where all lists have trailing censored information
                         ## this is pruned
                         if (length(out)>0)
                             x <- x[-out]
                         # stop at duplicated items
                         if (any(duplicated(x)))
                             stop(paste0("Duplicated items found in list ",n))
                         print(x)
                         x
                     })
    ## check class of elements, then coerce to integer
    cc <- sapply(object,class)
    if (length(cc <- unique(cc))>1) stop(paste("All elements of object must have the same class. Found:",paste(cc,collapse=", ")))
    if (match(cc,c("integer","character","numeric","factor"),nomatch=0)==0)
        stop("Class of lists in object should be one of 'integer', 'character', 'numeric' or 'factor'.")
    labels <- unique(unlist(object,recursive=TRUE,use.names=FALSE))
    nitems <- length(labels)
    if (match(c("sraNULL"),labels,nomatch=0)>0) stop("Item name sraNULL is reserved for missing items")
    object <- lapply(object,function(x){as.integer(factor(x,levels=c("sraNULL",labels)))-1})

    ## items are coded as 1, 2, 3, ...
    ## missing items (sraNULL) are coded as 0
    items <- seq(nitems) 

    ## Compute a list of missing items for each list
    missing.items <- lapply(object, function(x) {items[match(items,x,nomatch=0)==0]})
    missing.items1 <- lapply(object, function(x) {items[-x]})
    nmiss <- sapply(missing.items,length)
    
    ## fill too short lists with 0 (code for missing)
    ll <- sapply(object,length)
    listlength <- max(ll)
    tooshort <- any(ll<max(listlength))
    if (tooshort)
        object <- lapply(object,function(x){c(x,rep(0,listlength-length(x)))})

    ## set B to 1, if there is no censoring
    ##             or if only one element is censored
    iscensored <- any(nmiss!=0)
    if (B!=1 && (!iscensored || (max(nmiss)==1))) {B <- 1}

    # Special version of sample needed
    resample <- function(x, ...) x[sample.int(length(x), ...)]
    tmpres <- sapply(1:B, function(b) {
                         obj.b <- lapply(1:nlists,function(j){
                                             list <- object[[j]]
                                             if (nmiss[[j]]>0){
                                                 if (nmiss[[j]]==1){ ## fill in missing item
                                                     list[list==0] <- missing.items[[j]]
                                                 } else{ ## replace censored items with random order
                                                       list[list==0] <- resample(missing.items[[j]])
                                                   }
                                             }
                                             list
                                         })
                         ## bind lists
                         rankmat <- do.call("cbind",obj.b)
                         res <- sracppfull(rankmat)
                         res
                     })
    agreement <- rowMeans(tmpres)
    names(agreement) <- labels
    class(agreement) <- "sra"
    attr(agreement, "B") <- B
    agreement
}



#' Simulate sequential rank agreement for randomized unrelated lists
#'
#' @param obj Either a vector or matrix
#' @param B Either a vector or matrix
#' @param n the number of sequential rank agreement curves to produce
#' @return A matrix with n columns each representing the sequential rank agreement obtained from 
#' @author Claus Ekstrøm <ekstrom@@sund.ku.dk>
#' @export
random_list_sra <- function(object, B=1, n=1) {

    ## Make sure that the input object ends up as a matrix with integer columns all
    ## consisting of elements from 1 and up to listlength
    if (!is.matrix(object))
        object <- as.matrix(do.call("cbind",object))

    nitems <- nrow(object)
    notmiss <- apply(object, 2, function(x) {sum(x>0)} )
    res <- sapply(1:n, function(i) {
        ## Do a permutation with the same number of missing
        for (j in 1:ncol(object)) {
            object[,j] <- c(sample(nitems, size=notmiss[j]), rep(0, nitems-notmiss[j]))
        }
        sra(object, B=B) 
    })
    res
    
}



#' Smooth quantiles of a matrix of sequential ranked agreements
#'
#' @param obj A matrix
#' @param confidence the limits to compute
#' @return A list containing two vectors for the smoothed lower and upper limits
#' @author Claus Ekstrøm <ekstrom@@sund.ku.dk>
#' @export
smooth_sra <- function(object, confidence=0.95) {

    alpha <- (1-confidence)/2
    limits <- apply(object, 1, function(x) { quantile(x, probs=c(alpha, 1-alpha)) })
#    res <- apply(limits, 1, function(x) { loess(x)} )
#    list(lower=res[[1]]$y, upper=res[[2]]$y)
    list(lower=limits[1,], upper=limits[2,])
}
