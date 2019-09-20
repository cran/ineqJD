zenga.dataProcessed <-
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

  ol <- apply(nhl, 2, function (l) min(which(l > 0)))
  ul <- apply(nhl, 2, function (l) max(which(l > 0)))

  pl_h <- x$Phl / rowSums(x$Phl)

  Shlk <- x$Qhlk
  Shlk[t, , ] <- NA
  Shlk[t, , ] <- x$Qhlk[t, , ] - x$Qhlk[t - 1, , ]

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
  Mlk <- Minfhlk[r, , ]
  
  # upper group
  ag_h <- sapply(1:g, function(l) (nl[l] - x$Phl[, l]) / (N - rowSums(x$Phl)))
  ag_h[r, ] <- nhl[r, ] / nh[r]
  
  wllh <- array(sapply(1:r, function(h) outer(ag_h[h,], pl_h[h, ])), c(g, g, r))
  dimnames(wllh) <- list(l1 = g_names, l2 = g_names, h = 1:r)
  
  # Dif_w.g_Phg , matrice delle differenze [w.g (num. tot nel gruppo g)]-[Phg (num cumulata nel gruppo g)]
  Dif_nl_Phl <- array(sapply(1:g, function(l) nl[l ] - x$Phl[, l]), c(r, g))
  
  # Dif_TXjwig_QhgXjwi , array delle differenze TXjwig-QhgXjwi: r righe, g colonne, c matrici
  Dif_T_Qhlk <- array(sapply(1:s, function(k) sapply(1:g, function(l) x$Qhlk[r, l, k]-x$Qhlk[, l, k])),c(r, g, s))
  
  tmp <- Msuphlk <- x$Qhlk
  tmp[] <- Msuphlk[] <- NA
  for (k in 1:s) {
    for (h in 1:r) {
      for (l in 1:g) {
        tmp[h, l, k] <- Dif_T_Qhlk[h, l, k]/Dif_nl_Phl[h, l]
        Msuphlk[h, l, k] <- ifelse(is.nan(tmp[h, l, k]) == TRUE, Shlk[ul[l], l, k] / nhl[ul[l], l], tmp[h, l, k])
      }
    }
  }
  
  tmp <- array(sapply(1:s, function(k) Msuphlk[ , , k] * ag_h), c(r, g, s))
  tmp <- apply(tmp, c(1, 2), sum)
  Msuph <- rowSums(tmp)
  
  tmp <- array(
    sapply(1:s, 
           function(k) sapply(1:r, 
            function(h) outer(Msuphlk[h, , k], t(Minfhlk[h, , k]), "-") / Msuph[h])),
              c(g, g, r, s)
    )
  Bllhk <- array(sapply(1:s, function(k) tmp[ , , , k] * wllh), c(g, g, r, s))
  
  dimnames(Bllhk) <- list(
    l1 = g_names,
    l2 = g_names,
    h = 1:r,
    k = dimnames(x$Qhlk)[[3]]
  )
  
  rm(tmp)

  
  res <- list(
    index = "Zenga",
    decomposition = Bllhk,
    dataProcessed = x
  )
  
  class(res) <- "decomposition"
  
  return(res)
}
