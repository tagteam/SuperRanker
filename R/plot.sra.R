##' Plot the agreement between lists as a function of the list depth
##'
##' @title Plot sequential rank agreement
##' @param x Agreement object
##' @param xlim x-axis limits
##' @param ylim y-axis limits
##' @param xlab x-axis lab
##' @param ylab y-axis lab
##' @param add Logical. If \code{TRUE} add graph to existing plot.
##' @param ... Processed by function \code{prodlim::SmartControl}.
##' @return Graph
##' @examples
##' R1=c(1,2,3,4,5,7,6,8,9,10,11,12,13)
##' R2=c(5,11,4,7,8,3,12,13,6,10,9,2,1)
##' a <- sra(list(R1,R2))
##' plot(a)
##' arand = colMeans(do.call("rbind",lapply(1:20,function(b){
##'     sra(list(sample(R1),sample(R1)))
##' })))
##' lines(1:length(R1),arand,col=2,lwd=3)
##'
##' l <- c(1,2,3,4,5,7,6,8,9,10,11,12,13)
##' l <- 1:100
##' aa <- sapply(1:20,function(i){
##'     sra(list(sample(l),sample(l),sample(l)))[i]
##' })
##' c(mean(aa),sd(aa))
##'
##' @importFrom graphics plot axis lines
#' @export
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
plot.sra <- function(x, xlim, ylim, xlab="List depth", ylab="Sequential rank agreement", add=FALSE, ...){
    if (missing(xlim)) xlim=c(1,length(x))
    if (missing(ylim)) ylim=c(0,max(x))
    axis1.DefaultArgs <- list()
    axis2.DefaultArgs <- list(side=2)
    plot.DefaultArgs <- list(x=0,y=0,type = "n",ylim = ylim,xlim = xlim,xlab = xlab,ylab = ylab)
    lines.DefaultArgs <- list(type="l",lwd=3,x=1:length(x),y=x)
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
