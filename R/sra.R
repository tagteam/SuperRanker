sra <- function(object, metric=var, depth=NULL, listlength=NULL, B=10) {
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
    maxdepth <- min(nrow(fullrankmat), ifelse(is.null(depth), nrow(fullrankmat), depth))

      
    tmpres <- sapply(1:B, function(i) {
      # Ad 1)
      for (j in 1:nlists) {
          if (length(missing.items[[j]])>0) {
          fullrankmat[(listlength-length(missing.items[[j]])+1):listlength,j] <- sample(missing.items[[j]])
        }
      }
      res <- full.list.metric(fullrankmat, metric)

      rankmat.l <- lapply(seq_len(nrow(fullrankmat)), function(i) fullrankmat[i,])
      uniq.l <- lapply(rankmat.l, function(x) {sort(unique(x))})

      agreement <- sapply(1:maxdepth, function(x) {
          myvar <- unique(unlist(uniq.l[1:x]))
          mean(res[myvar])
      })

      agreement
      
  })
    
    out <- list(agreement=apply(tmpres, 1, mean), metric=quote(metric), depth=maxdepth, B=B)
    class(out) <- "agreement"
    out
}
