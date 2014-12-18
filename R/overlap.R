##' Quantify the average overlap of two lists
##' in the head of the lists and until given depth
##'
##' 
##' @title Quantify overlap of two lists
##' @param l1 First list 
##' @param l2 Second list
##' @param depth Interest margin, i.e., agreement is limited to the first \code{depth} elements of each list.
##' @return A vector of length \code{depth} whose elements are the average overlap.
##' @seealso list.agreement
##' @examples
##' overlap(c(1,2,3,4,5,7,6,8), c(1,2,3,5,6,4,7,8))
##' @export 
##' @author Claus Ekstrøm <ekstrom@@sund.ku.dk>, Thomas A. Gerds <tag@@biostat.ku.dk>
overlap <- function(l1, l2, depth=NULL) {
    rankmat <- as.matrix(cbind(l1, l2))
    nlists <- ncol(rankmat)
    maxdepth <- min(nrow(rankmat), ifelse(is.null(depth), nrow(rankmat), depth))
    overlap <- sapply(1:maxdepth, function(x) { sum(duplicated(as.vector(rankmat[1:x,]))) / x})
    cumsum(overlap)/seq_len(maxdepth)
}
