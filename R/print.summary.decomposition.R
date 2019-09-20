print.summary.decomposition <-
function(x, ...) {
    app <- 4
    cat("\nSynthetic", x$index, "inequality index:", round(x$synthetic, app), "\n")
    cat("\nJoint contributions:\n")
    print(round(x$joint, app))
    cat("### GROUPS ###\n")
    print(round(x$pairs, app))
    cat("\nWithin:\n")
    print(round(x$within, app))
    cat("\nBetween:\n")
    print(round(x$between, app))
    cat("\nMarginal:\n")
    print(round(x$groups, app))
    cat("\n### GROUPS AND SOURCES ###\n")
    print(round(x$groups_sources, app))
    cat("\n### SOURCES ###\n")
    print(round(x$sources, app))
    cat("\n")
}
