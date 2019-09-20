print.dataProcessed <-
function(x, max = 1000, ...) {
    d <- dim(x$Qhlk)
    cat(
        "\nr = ", d[1], " different values of Y\n",
        "\ng = ", d[2], " group(s)\n\t",
        paste(
            dimnames(x$Qhlk)[[2]], collapse = "\n\t"
        ),
        "\n\ns = ", d[3], " source(s)\n\t",
        paste(
            dimnames(x$Qhlk)[[3]], collapse = "\n\t"
        ),
        "\n",
        sep = ""
    )

    cat("\n$yh values of Y, h = 1, ..., r\n")
    print(x$yh, max = max)
    cat("\n$Phl cumlated weights, l = 1, ..., g\n")
    print(x$Phl, max = max)
    cat("\n$Qhlk cumulated sources, k = 1, ..., s\n")
    print(x$Qhlk, max = max)
}
