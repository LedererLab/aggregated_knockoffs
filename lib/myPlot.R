#! R
# Filename: myPlot.R

# A plot function: display the relationship between x and y 
# Args: 
#    x: a vector
#    y: a vector or a matrix whose length or rows should be equal to the length of x
#    main: main title.
#    xlab: a character. The default is "Target FDR".
#    ylab: a character. It can be one of "Actual FDR", "Power" and "selection accuracy".
#    legend.x: the x-axis position of legend.
#    legend.y: the y-axis position of legend.
#    legendnames: should be a character vector with length = ncol(x)
#
# Output:
#    a plot displaying the relationship between x and y (a plot with grid)

myPlot <- function(x, y, main, ylab = c("Actual FDR", "Power", "selection accuracy"), 
                   xlab = "Target FDR", legend.x = NULL, legend.y = NULL,
                   legendnames = NULL){
  pch <- c(1, 19)
  col <- c("darkorchid", "orange")
  plot(NULL, NULL, xlim = c(min(0, min(x)), max(1, max(x))), cex = 1, xlab = "", ylab = "", 
       cex.axis = 1.2, type = "n", ylim = c(min(0, min(y)), max(1, max(y))), las = 1)
  points(x, y[ , 1], pch = pch[1], cex = 0.8, col = col[1])
  points(x, y[ , 2], pch = pch[2], cex = 0.8, col = col[2])
  if(ylab == "Actual FDR"){
    abline(a = 0, b = 1)
  }
  title(ylab = ylab, xlab = xlab, cex.lab = 1.5, line = 3)
  title(main = main, cex.main = 2, line = 1, font.main = 2)
  grid(col = "darkgray")
  if(is.null(legend.x) & is.null(legend.y)){
    if(ylab == "Actual FDR" | ylab == "selection accuracy"){
      legend.x = - 0.05
      legend.y = 1.05
    }
    if(ylab == "Power"){
      legend.x = 0.45
      legend.y = 0.3
    }
  }
  legend(x = legend.x, y = legend.y, legend = legendnames, cex = 1.5, col = col,              
         pch = pch, bty = "n", y.intersp = 0.7)
}









