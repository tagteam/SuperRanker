#' Compute the average overlap
#'
#'
#' @param obj Either a vector or matrix
#' @return A vector of the average overlap
#' @export
average_overlap <- function(obj) {
    if (is.matrix(obj)) {
        obj <- overlap(obj)
    } else {
        stop("average_overlap only works for matrices right now")
    }
    cumsum(obj)/seq_along(obj)
}
