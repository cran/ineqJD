inequalityCurves.decomposition <-
function(x, l = 1:dim(x$decomposition)[2], k = 1:dim(x$decomposition)[4], ...) {

    if (anyDuplicated(l) | anyDuplicated(k)) {
      stop("multiple selections are not allowed")
    }
    
    r <- dim(x$decomposition)[3]
    
    p <- rowSums(x$dataProcessed$Phl)
    p <- p / p[r]

    h <- x$decomposition[, l, , k]
    h <- array(h, c(dim(x$decomposition)[1], length(l), r, length(k)))
    h <- apply(h, 3, sum)

    res <- stepfun(p, c(h, h[r]), right = TRUE)

    g_names <- dimnames(x$decomposition)[[2]]
    names(g_names) <- g_names
    s_names <- dimnames(x$decomposition)[[4]]
    names(s_names) <- s_names
    attributes(res)$call <- deparse(sys.calls()[[sys.nframe() - 1]])
    attributes(res)$index <- x$index
    attributes(res)$min <- min(h)
    attributes(res)$max <- max(h)
    attributes(res)$groups <- g_names
    attributes(res)$sources <- s_names
    attributes(res)$selected_groups <- g_names[l]
    attributes(res)$selected_sources <- s_names[k]
    class(res) <- c("inequality_curves", "stepfun", "function")
    return(res)
}
