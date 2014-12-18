##' Plot the agreement between lists as a function of the list depth
##' 
##' @title Plot agreement 
##' @param x Agreement object
##' @param xlim x-axis limits
##' @param ylim y-axis limits
##' @param xlab x-axis lab
##' @param ylab y-axis lab
##' @param ... Processed by \code{\link{prodlim::SmartControl}}.
##' @return Graph
##' @seealso 
##' @examples
##' a <- agreement(list(c(1,2,3,4,5,7,6,8,9,10,11,12,13), c(5,11,4,7,8,3,12,13,6,10,9,2,1)))
##' plot(a)
#' @method plot agreement
#' @S3method plot agreement
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
plot.agreement <- function(x,xlim,ylim,xlab="List depth",ylab="Agreement",add=FALSE,...){
    if (missing(xlim)) xlim=c(0,x$depth)
    if (missing(ylim)) ylim=range(x$agreement)
    axis1.DefaultArgs <- list()
    axis2.DefaultArgs <- list(side=2)
    plot.DefaultArgs <- list(x=0,y=0,type = "n",ylim = ylim,xlim = xlim,xlab = xlab,ylab = ylab)
    lines.DefaultArgs <- list(type="l",lwd=3,x=1:x$depth,y=x$agreement)
    smartA <- prodlim::SmartControl(call=list(...),
                                    keys=c("plot","lines","axis1","axis2"),
                                    ignore=c("x","ylim","xlim","xlab","ylab"),
                                    defaults=list("plot"=plot.DefaultArgs,
                                        "lines"=lines.DefaultArgs,
                                        "axis1"=axis1.DefaultArgs,
                                        "axis2"=axis2.DefaultArgs),
                                    forced=list("plot"=list(axes=FALSE),"axis1"=list(side=1)),
                                    ignore.case=TRUE,
                                    replaceDefaults=FALSE,
                                    verbose=TRUE)
    if (add==FALSE){
        do.call("plot",smartA$plot)
        ## plot(0,0,type="n",axes=FALSE,xlab="Depth",ylab="")
        do.call("axis",smartA$axis1)
        do.call("axis",smartA$axis2)
    }
    do.call("lines",smartA$lines)
}
