#----------------------------------------------------------------
# A plot function: display the relationship between x and y 
# Args: 
#    x: a vector
#    y: a vector or a matrix whose length or rows should be equal to the length of x
#    main: main title, such as "lasso".
#    xlab: a character. The default is "Target FDR".
#    ylab: a character. The default is "selection accuracy".
#    legendnames: should be a character vector with length = ncol(x)
#
# Output:
#    a plot displaying the relationship between x and y (a plot with grid)

myPlot.SA <- function(x, y, main, xlab = "Target FDR", ylab = "selection accuracy",
                      legendnames = NULL){
  pch <- c(1, 19)
  col <- c("darkorchid","orange")
  plot(NULL, NULL, xlim = c(min(0, min(x)), max(1, max(x))), cex = 1, xlab = "", ylab = "", 
       cex.axis = 1.2, type = "n", ylim = c(min(0, min(y)), max(1, max(y))), las = 1)
  points(x, y[ , 1], pch = pch[1], cex = 0.8, col = col[1])
  points(x, y[ , 2], pch = pch[2], cex = 0.8, col = col[2])
  title(ylab = ylab, xlab = xlab, cex.lab = 1.5, line = 3)
  grid(col = "darkgray")
  title(main = main, cex.main = 2, line = 1, font.main = 2)
  legend(x = 0.5, y = 1.05, legend = legendnames, cex = 1.5, col = col,              
         pch = pch, bty = "n", y.intersp = 0.7)
}

#--------------------------------------------------------------------------
# A plot function: display the relationship between x and y 
# Args: 
#    x: a vector
#    y: a vector or a matrix whose length or rows should be equal to the length of x
#    main: main title, such as "lasso".
#    xlab: a character. The default is "Target FDR".
#    ylab: a character. The default is "Actual FDR".
#    legendnames: should be a character vector with length = ncol(x)
#
# Output:
#    a plot displaying the relationship between x and y (a plot with diagonal line)

myPlot.FDR <- function(x, y, main, ylab = "Actual FDR", xlab = "Target FDR",
                       legendnames = NULL){
  pch <- c(1, 19)
  col <- c("darkorchid","orange")
  plot(NULL, NULL, xlim = c(min(0, min(x)), max(1, max(x))), cex = 1, xlab = "", ylab = "", 
       cex.axis = 1.2, type = "n", ylim = c(min(0, min(y)), max(1, max(y))), las = 1)
  points(x, y[ , 1], pch = pch[1], cex = 0.8, col = col[1])
  points(x, y[ , 2], pch = pch[2], cex = 0.8, col = col[2])
  abline(a = 0, b = 1)
  title(ylab = ylab, xlab = xlab, cex.lab = 1.5, line = 3)
  title(main = main, cex.main = 2, line = 1, font.main = 2)
  legend(x = - 0.05, y = 1.05, legend = legendnames, cex = 1.5, col = col,              
         pch = pch, bty = "n", y.intersp = 0.7)
}

#---plot Power---#
myPlot.Power <- function(x, y, main, ylab = "Power", xlab = "Target FDR",
                       legendnames = NULL){
  pch <- c(1, 19)
  col <- c("darkorchid", "orange")
  plot(NULL, NULL, xlim = c(min(0, min(x)), max(1, max(x))), cex = 1, xlab = "", ylab = "", 
       cex.axis = 1.2, type = "n", ylim = c(min(0, min(y)), max(1, max(y))), las = 1)
  points(x, y[ , 1], pch = pch[1], cex = 0.8, col = col[1])
  points(x, y[ , 2], pch = pch[2], cex = 0.8, col = col[2])
  title(ylab = ylab, xlab = xlab, cex.lab = 1.5, line = 3)
  title(main = main, cex.main = 2, line = 1, font.main = 2)
  grid(col = "darkgray")
  legend(x = 0.45, y = 0.3, legend = legendnames, cex = 1.5, col = col,              
         pch = pch, bty = "n", y.intersp = 0.7)
}

#--my plot function: to plot the relationship between FDR, Power and target FDR--#
myPlot <- function(x, y, main, ylab = c("Actual FDR", "Power", "selection accuracy"), 
                   xlab = "Target FDR", legend.x = NULL, legend.y = NULL,
                   legendnames = NULL){
  pch <- c(1, 19, 17)
  col <- c("skyblue", "darkorchid", "orange")
  plot(NULL, NULL, xlim = c(min(0, min(x)), max(1, max(x))), cex = 1, xlab = "", ylab = "", 
       cex.axis = 1.2, type = "n", ylim = c(min(0, min(y)), max(1, max(y))), las = 1)
  points(x, y[ , 1], pch = pch[1], cex = 0.8, col = col[1])
  points(x, y[ , 2], pch = pch[2], cex = 0.8, col = col[2])
  points(x, y[ , 3], pch = pch[3], cex = 0.8, col = col[3])  
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









