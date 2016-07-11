#' Compute the average overlap
#'
#' @param obj Either a vector or matrix
#' @return A vector of the average overlap
#' @examples
#' # setting with 3 lists
#' mlist <- matrix(cbind(1:8,c(1,2,3,5,6,7,4,8),c(1,5,3,4,2,8,7,6)),ncol=3)
#' average_overlap(mlist)
#' @export
average_overlap <- function(obj) {
    if (is.matrix(obj)) {
        obj <- overlap(obj)
    } else {
        stop("average_overlap only works for matrices right now")
    }
    cumsum(obj)/seq_along(obj)
}
