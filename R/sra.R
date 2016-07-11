#' Compute the sequential rank agreement
#'
#' @param object Either matrix where each column is a ranked list of
#'     items or a list of ranked lists of items. Elements are integers
#'     between 1 and the length of the lists. The lists should have
#'     the same length but censoring can be used by setting the list
#'     to zero from a point onwards. See details for more information.
#' @param B An integer giving the number of randomization to sample
#'     over in the case of censored observations
#' @param na.strings A vector of strings/values that represent missing
#'     values in addition to NA. Defaults to NULL which means only NA
#'     are censored values.
#' @param nitems Claus: please describe this argument
#' @param type The type of measure to use. Either sd (standard
#'     deviation - the default) or mad (median absolute deviance)
#' @param ... Arguments passed to methods.
#' @return A vector of the sequential rank agreement
##' @examples
##'
##' mlist <- matrix(cbind(1:8,c(1,2,3,5,6,7,4,8),c(1,5,3,4,2,8,7,6)),ncol=3)
##' sra(mlist)
##'
##' mlist <- matrix(cbind(1:8,c(1,2,3,5,6,7,4,8),c(1,5,3,4,2,8,7,6)),ncol=3)
##' sra(mlist, nitems=20, B=10)
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
#' @rdname sra
#' @export
sra <- function(object,B,na.strings,nitems,type,...) {
  UseMethod("sra")
}

#' @rdname sra
#' @export
sra.default <- function(object,B,na.strings,nitems,type,...) {
    stop("Input must be either a matrix, a data.frame or a list.")
}

#' @rdname sra
#' @export
sra.matrix <- function(object, B=1, na.strings=NULL, nitems=nrow(object), type=c("sd", "mad"),...) {
    if (!is.matrix(object))
        stop("Input object must be a matrix")

    ## Convert all missing types to NAs
    if (!is.null(na.strings)) {
        object[object %in% na.strings] <- NA
    }

    unique.items <- length(unique(object[!is.na(object)]))

    ## Check that we dont have more unique values than rows
    if (unique.items > nitems) {
        stop("Found more unique items in the matrix than rows/nitems. Increase nitems to match")
    }

    ## Expand the columns in the matrix to have length unique.items
    if (nitems>nrow(object)) {
        glue <- matrix(rep(NA, NCOL(object)*(nitems - nrow(object))), ncol=NCOL(object))
        object <- rbind(object, glue)
    }
    object <- lapply(1:NCOL(object),function(j)object[,j]) # Convert matrix to list
    sra.list(object, B=B, nitems=nitems, type=type)
}



#' @rdname sra
#' @export
sra.list <- function(object, B=1, na.strings=NULL, nitems=max(sapply(object, length)), type=c("sd", "mad"),...) {
    # Make sure that the input object ends up as a matrix with integer columns all
    # consisting of elements from 1 and up to listlength

    stopifnot(is.list(object))

    nlists <- length(object)
    nitems <- nitems ## force evaluation here

    type <- match.arg(type)

    ## Sanity checks
    object <- lapply(1:length(object),function(j){
                         x <- object[[j]]
                         ## Add the NA items to be removed
                         out <- which(is.na(x))
                         if (!is.null(na.strings)) {
                             na.strings <- na.strings[!is.na(na.strings)]
                             out <- c(out,grep(paste0("^",na.strings,"$"),x))
                         }
                         ## remove censored items with side effect:
                         ## in case where all lists have trailing censored information
                         ## this is pruned
                         if (length(out)>0)
                             x <- x[-out]
                         ## stop at duplicated items
                         if (any(duplicated(x)))
                             stop(paste0("Duplicated items found in list ",j))
                         x
                     })
    ## check class of elements, then coerce to integer

    cc <- sapply(object, class)
    if (length(cc <- unique(cc))>1)
        stop(paste("All elements of object must have the same class. Found:",paste(cc,collapse=", ")))

    if (match(cc,c("integer","character","numeric","factor"),nomatch=0)==0)
        stop("Class of lists in object should be one of 'integer', 'character', 'numeric' or 'factor'.")

    labels <- unique(unlist(object,recursive=TRUE,use.names=FALSE))

    nitems <- max(nitems, length(labels))

    object <- lapply(object,function(x){as.integer(factor(x,levels=labels))})

    ## items are coded as 1, 2, 3, ...
    ## missing items (sraNULL) are coded as 0
    items <- seq(nitems)

    ## Compute a list of missing items for each list
    missing.items <- lapply(object, function(x) {items[match(items,x,nomatch=0)==0]})
    nmiss <- sapply(missing.items,length)

    ## fill too short lists with 0 (code for missing)
    ll <- sapply(object,length)
    listlength <- max(ll)
    tooshort <- any(ll<nitems)
    if (tooshort)
        object <- lapply(object,function(x){c(x,rep(0,nitems-length(x)))})

    ## set B to 1, if there is no censoring
    ##             or if only one element is censored
    iscensored <- any(nmiss!=0)
    if (B!=1 && (!iscensored || (max(nmiss)==1))) {B <- 1}

    itype <- 0
    if (type == "mad")
        itype <- 1

    # Special version of sample needed
    resample <- function(x, ...) x[sample.int(length(x), ...)]
    tmpres <- sapply(1:B, function(b) {
        obj.b <- lapply(1:nlists,function(j){
            list <- object[[j]]
            if (nmiss[[j]]>0){
                list[list==0] <- resample(missing.items[[j]])
            }
            list
        })
        ## bind lists
        rankmat <- do.call("cbind",obj.b)
        res <- sracppfull(rankmat, type=itype)
        res
    })
    if (itype==0) {
        agreement <- sqrt(rowMeans(tmpres))
    } else {
        agreement <- rowMeans(tmpres)
    }
    names(agreement) <- items
    class(agreement) <- "sra"
    attr(agreement, "B") <- B
    attr(agreement, "type") <- type
    agreement

}


#' Simulate sequential rank agreement for randomized unrelated lists
#'
#' Simulate sequential rank agreement from completely uninformative lists (ie., raw permutations of items) and compute the corresponding sequential rank agreement curves.
#' The following attributes are copied from the input object: number of lists, number of items and amount of censoring.
#'
#' @param object A matrix or list of vectors representing ranked lists.
#' @param B An integer giving the number of randomizations to sample
#'     over in the case of censored observations
#' @param n Integer: the number of permutation runs. For each permutation run we permute each of the lists in object
#' and compute corresponding the sequential rank agreement curves 
#' @param na.strings A vector of character values that represent
#'     censored observations
#' @param type The type of measure to use. Either sd (standard
#'     deviation - the default) or mad (median absolute deviance)
#' @return A matrix with n columns and the same number of rows as for the input object. Each column contains one
#' simulated sequential rank agreement curve from one permutation run.
#' @author Claus Ekstrøm <ekstrom@@sund.ku.dk>
#' @examples
#' # setting with 3 lists
#' mlist <- matrix(cbind(1:8,c(1,2,3,5,6,7,4,8),c(1,5,3,4,2,8,7,6)),ncol=3)
#' # compute sequential rank agreement of lists
#' sra(mlist)
#' # compute sequential rank agreement of 5 random permutations
#' random_list_sra(mlist, n=5)
#'
#' @export
random_list_sra <- function(object, B=1, n=1, na.strings=NULL, type=c("sd", "mad")) {

    type <- match.arg(type)

    ## Make sure that the input object ends up as a matrix with integer columns all
    ## consisting of elements from 1 and up to listlength
    if (!is.matrix(object))
        object <- as.matrix(do.call("cbind",object))

    ## Convert all missing types to NAs
    if (!is.null(na.strings)) {
        object[object %in% na.strings] <- NA
    }

    nitems <- nrow(object)
    notmiss <- apply(object, 2, function(x) {sum(!is.na(x))} )
    res <- sapply(1:n, function(i) {
        ## Do a permutation with the same number of missing
        for (j in 1:ncol(object)) {
            object[,j] <- c(sample(nitems, size=notmiss[j]), rep(NA, nitems-notmiss[j]))
        }
        sra(object, B=B, type=type)
    })
    res

}



#' Smooth quantiles of a matrix of sequential ranked agreements. 
#'
#' @param object A matrix
#' @param confidence the limits to compute
#' @return A list containing two vectors for the smoothed lower and upper limits
#' @author Claus Ekstrøm <ekstrom@@sund.ku.dk> 
#' @examples
#' # setting with 3 lists
#' mlist <- matrix(cbind(1:8,c(1,2,3,5,6,7,4,8),c(1,5,3,4,2,8,7,6)),ncol=3)
#' # compute rank agreement of 5 random permutations
#' null=random_list_sra(mlist,n=15)
#' # now extract point-wise quantiles according to confidence level
#' smooth_sra(null)
#' @export
smooth_sra <- function(object, confidence=0.95) {

    alpha <- (1-confidence)/2
    limits <- apply(object, 1, function(x) {stats::quantile(x, probs=c(alpha, 1-alpha)) })
    list(lower=limits[1,], upper=limits[2,])
}



#' Compute a Kolmogorov-Smirnoff-like test for Smooth quantiles of a matrix of sequential rank agreements
#'
#' @param object An object created with \code{sra}.
#' @param nullobject An object created with \code{random_list_sra}.
#' @param weights Either a single value or a vector of the same length as the number of item with the weight that should be given to specific depths.
#' @return A single value corresponding to the p-value
#' @author Claus Ekstrøm <ekstrom@@sund.ku.dk>
#' @examples
#' # setting with 3 lists
#' mlist <- matrix(cbind(1:8,c(1,2,3,5,6,7,4,8),c(1,5,3,4,2,8,7,6)),ncol=3)
#' # compute sequential rank agreements
#' x=sra(mlist)
#' # compute rank agreement of 5 random permutations
#' null=random_list_sra(mlist,n=15)
#' # now extract point-wise quantiles according to confidence level
#' test_sra(x,null)
#' # compare to when we use the result of the first permutation run
#' test_sra(null[,1],null[,-1])
#' 
#' @export
test_sra <- function(object, nullobject, weights=1) {
    ## Sanity checks
    if (! (length(weights) %in% c(1, length(object))))
        stop("the vector of weights must have the same length as the number of items")
    
    ## Test statistic
    T <- max(weights*abs(object - apply(nullobject, 1, mean)))

    ## Now compute the individual jackknife variations from the null object
    B <- ncol(nullobject)
    nullres <- sapply(1:B, function(i) {
        max(weights*abs(nullobject[,i] - apply(nullobject[,-i], 1, mean)))
    })

    res <- sum(nullres>=T)/(B+1)
    attr(res, "B") <- B
    return(res)

}
