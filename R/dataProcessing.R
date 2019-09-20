dataProcessing <-
function(
  units, # sources values for each unit, vector if income, matrix if sources
  groups = rep("G1", nrow(as.matrix(units))), # vector, preferibly factor
  weights = rep(1, nrow(as.matrix(units))) # vector
) {
  ###############
  # input check #
  ###############
  units <- as.matrix(units)
  groups <- as.factor(groups)
  N <- nrow(units)
  if (length(groups) != N)
    stop("the length of groups vector differs from the number of units")
  if (length(weights) != N)
    stop("the length of weights vector differs from the number of units")
  if (any(weights < 0))
    stop("weights cannot be negative")
  ####################
  # data preparation #
  ####################
  yi <- rowSums(units)
  if (any(yi < 0))
    warning("data contain negative sums by sources")
  # sorting
  ord <- order(yi)
  yi <- yi[ord]
  units <- as.matrix(units[ord, ])
  groups <- groups[ord]
  # grouping
  g_names <- levels(groups)
  g <- length(g_names)
  s <- ncol(units)
  if (is.null(colnames(units)))
    colnames(units) <- paste0("X", 1:s)
  s_names <- colnames(units)
  wil <- outer(1:N, g_names, function(i, l) weights[i] * (groups == l))
  colnames(wil) <- g_names
  Pil <- apply(wil, 2, cumsum)
  Qilk <- sapply(s_names, function(k) apply(units[, k] * wil, 2, cumsum))
  Qilk <- array(Qilk, c(N, g, s))
  dimnames(Qilk) <- list(i = 1:N, l = g_names, k = s_names)
  sel_h <- !duplicated(yi, fromLast = TRUE)
  r <- sum(sel_h)
  res <- list(
    yh = yi[sel_h],
    Phl = array(Pil[sel_h, ], c(r, g)),
    Qhlk = array(Qilk[sel_h, , ], c(r, g, s))
  )
  dimnames(res$Phl) <- list(h = 1:r, l = g_names)
  dimnames(res$Qhlk) <- list(h = 1:r, l = g_names, k = s_names)
  class(res) <- "dataProcessed"
  return(res)
}
