bonferroni.dataProcessed <-
function(x) {
  r <- length(x$yh)
  t <- 2:r
  g <- dim(x$Phl)[2]
  s <- dim(x$Qhlk)[3]

  g_names <- colnames(x$Phl)

  nhl <- x$Phl
  nhl[t, ] <- NA
  nhl[t, ] <- x$Phl[t, ] - x$Phl[t - 1, ]

  nh <- rowSums(nhl)
  nl <- colSums(nhl)
  N <- sum(nl)
  fl <- nl / N

  ol <- apply(nhl, 2, function (l) min(which(l > 0)))

  pl_h <- x$Phl / rowSums(x$Phl)

  wllh <- array(sapply(1:r, function(h) outer(fl, pl_h[h, ])), c(g, g, r))
  dimnames(wllh) <- list(l1 = g_names, l2 = g_names, h = 1:r)
  
  Shlk <- x$Qhlk
  Shlk[t, , ] <- NA
  Shlk[t, , ] <- x$Qhlk[t, , ] - x$Qhlk[t-1, , ]

  tmp <- Minfhlk <- x$Qhlk
  tmp[] <- Minfhlk[] <- NA
  for (k in 1:s) {
    for (h in 1:r) {
      for (l in 1:g) {
        tmp[h, l, k] <- x$Qhlk[h, l, k] / x$Phl[h, l]
        Minfhlk[h, l, k] <- ifelse(
          is.nan(tmp[h, l, k]) == TRUE,
          Shlk[ol[l], l, k] / nhl[ol[l], l],
          tmp[h, l, k]
        )
      }
    }
  }
  
  Minfhl <- apply(Minfhlk, c(1, 2), sum)
  MY <- sum(Minfhl[r, ] %*% x$Phl[r, ]) / sum(x$Phl[r, ])
  Mlk <- array(Minfhlk[r, , ], c(g, s))
  dimnames(Mlk) <- dimnames(Minfhlk)[-1]
  tmp <- array(
    sapply(1:s,
    function(k)
      sapply(1:r,
        function(h) outer(Mlk[, k], Minfhlk[h, , k], "-") / MY
      )
    ),
    c(g, g, r, s)
  )
  Vllhk <- array(
    sapply(1:s,
      function(k) tmp[, , , k] * wllh),
    c(g, g, r, s)
  )
  dimnames(Vllhk) <- list(
    l1 = g_names,
    l2 = g_names,
    h = 1:r,
    k = dimnames(x$Qhlk)[[3]]
  )

  rm(tmp)
  res <- list(
    index = "Bonferroni",
    decomposition = Vllhk,
    dataProcessed = x
  )

  class(res) <- "decomposition"
  
  return(res)
}
