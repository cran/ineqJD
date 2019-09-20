summary.decomposition <-
function(object, ...) {
    

    tmp <- dim(object$decomposition)
    g <- tmp[2]
    r <- tmp[3]
    s <- tmp[4]

    g_names <- dimnames(object$decomposition)[[2]]
    s_names <- dimnames(object$decomposition)[[4]]

    t <- 2:r
    nhl <- object$dataProcessed$Phl
    nhl[t, ] <- NA
    nhl[t, ] <- object$dataProcessed$Phl[t, ] - object$dataProcessed$Phl[t - 1, ]
    nh <- rowSums(nhl)
    N <- sum(nh)
  
    tmp <- object$decomposition
    tmp[] <- NA
    for(k in 1:s){
        for(h in 1:r){
            tmp[, , h, k] <- object$decomposition[, , h, k] * nh[h] / N
        }
    }
    
    joint <- sapply(1:s, function(k) apply(array(tmp[ , , , k], c(g, g, r)), c(1, 2), sum))
    joint <- array(joint, c(g , g, s), dimnames = list(l1 = g_names, l2 = g_names, s = s_names))

    sources <- apply(joint, 3, sum)
    groups_sources <- apply(joint, c(2, 3), sum)
    names(dimnames(groups_sources)) <- c("l", "s")

    pairs <- apply(joint, c(1, 2), sum)
    groups <- colSums(pairs)
    within <- diag(pairs)
    between <- groups - within
    synthetic <- sum(pairs)

    res <- list(
        index = object$index,
        joint = joint,
        pairs = pairs,
        within = within,
        between = between,
        groups = groups,
        groups_sources = groups_sources,
        sources = sources,
        synthetic = synthetic
    )

    rm(tmp)
    

    class(res) <- "summary.decomposition"
    return(res)
}
