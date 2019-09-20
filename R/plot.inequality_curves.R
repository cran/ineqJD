plot.inequality_curves <-
function(
    x,
    pch = 16,
    from = 0,
    to = 1,
    xlim = NULL,
    ylim = NULL, # the y limits of the plot
    xaxs = "i",
    yaxs = "i",
    xlab = "p",
    ylab = NULL, # NULL see details
    main = attributes(x)$index,
    sub = paste0(
        "grp: ",
        paste(attributes(x)$groups, collapse = ", "),
        "; src: ",
        paste(attributes(x)$sources, collapse = ", ")
    ),
    ...) {
    if (is.null(ylab)) {
        ylab <- switch(
            attributes(x)$index,
            Zenga = "I(",
            Gini = "G(",
            Bonferroni = "V("
        )
        ylab <- paste0(ylab, xlab)
        if (length(attributes(x)$selected_groups) < length(attributes(x)$groups)) {
            ylab <- paste0(ylab, ", l = {", paste(attributes(x)$selected_groups, collapse = ", "), "}")
        }
        if (length(attributes(x)$selected_sources) < length(attributes(x)$sources)) {
            ylab <- paste0(ylab, ", k = {", paste(attributes(x)$selected_sources, collapse = ", "), "}")
        }
        ylab <- paste0(ylab, ")")
    }
    if (is.null(xlim)) {
        xlim <- c(from, to)
    }
    if (is.null(ylim)) {
        ylim <- c(min(0, attributes(x)$min), max(1, attributes(x)$max))
    }
    plot.stepfun(
        x,
        pch = pch,
        xlim = xlim,
        ylim = ylim,
        xaxs = xaxs,
        yaxs = yaxs,
        xlab = xlab,
        ylab = ylab,
        main = main,
        sub = sub,
        ...)
}
