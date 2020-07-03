library('Ternary')

par(mar=rep(0.2, 4))
TernaryPlot(alab = 'w_1', blab = 'w_2', clab = 'w_*', axis.labels = seq(0, 1, by = 0.1), grid.minor.lines=0)

FunctionToContour <- function (a, b, c) {
  2*a*b
}

values <- TernaryPointValues(FunctionToContour, resolution=100L)
ColourTernary(values)
TernaryContour(FunctionToContour, resolution=200L)
