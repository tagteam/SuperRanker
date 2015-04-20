#' Compute the sequential rank agreement
#'
#' @param obj Either matrix where each column is a ranked list of items or a list of ranked lists of items. Elements are integers between 1 and the length of the lists. The lists should have the same length but censoring can be used by setting the list to zero from a point onwards. See details for more information.
#' @param B An integer giving the number of randomization to sample over in the case of censored observations
#' @return A vector of the sequential rank agreement
#' @author Claus Ekstrøm <ekstrom@@sund.ku.dk>
#' @export
sra <- function(object, B=1) {

    # Make sure that the input object ends up as a matrix with integer columns all
    # consisting of elements from 1 and up to listlength
    if (is.matrix(object))
        rankmat <- object
    else
        rankmat <- as.matrix(do.call("cbind",object))

    listlength <- nrow(rankmat)
    nlists <- NCOL(rankmat)
    nseq <- seq(listlength)

    # Special version of sample needed
    resample <- function(x, ...) x[sample.int(length(x), ...)]
    
    # Compute a list of missing items for each list
    missing.items <- lapply(as.data.frame(rankmat), function(x) { nseq[-x] })

    ## Should make a sanity check that zeros are from a point onwards
    if (!all(sapply(1:nlists, function(x) {  res <- TRUE
                                             if (length(missing.items[[x]])>0) {
                                                 if (any(rankmat[(nitems-missing.items[[x]]):nitems,x]) )
                                                     { res <- FALSE }
                                             }
                                             res
                                         }))) {
        stop("Censored ranked lists should be coded 0 from a rank onwards")
    }

    ## If there is no censoring then we should set B to 1
    if (max(sapply(missing.items, function(i) length(i) ))==0) {
        B <- 1
    }
    
    tmpres <- sapply(1:B, function(i) {
        for (j in 1:nlists) {
            if (length(missing.items[[j]])>0) {      
                rankmat[(listlength-length(missing.items[[j]])+1):listlength,j] <- resample(missing.items[[j]])
            }
        }
        res <- sracppfull(rankmat)
        res
    })

    agreement <- apply(tmpres, 1, mean)

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
