# A plot function: plot the relationship between k and FDR or Power or selection accuracy
#
# Args:
#    y: a vector or a matrix (numeric)
#    xlab: default is "k"; it can be changed.
#    ylab: default is "FDR"; it can be changed, such as "Power" and "selection accuracy".
#    main: main title, such as "(n,p,s0)=(200,100,10)".
#    legendnames: should be a character vector, such as c("fdr = 0.1", "fdr = 0.2", "fdr = 0.5").
#    pch: default is NULL
#    col: default is NULL
#    ablines: should be a target FDR vector
#
# Outputs:
#    a plot

myPlot.kcurves <- function(y, xlab = "k", ylab = "FDR", main = main, legend.x = 0, legend.y = 0.8,
                           legendnames = NULL, pch = NULL, col = NULL, ablines = NULL){
  if(is.null(legendnames)){
    stop("legendnames should be a character vector with the length of nrow(y)")
  }
  if(is.null(ablines)){
    stop("ablines should be the target FDR vector with the length of nrow(y)")
  }
  nrows.y <- nrow(y)
  ncols.y <- ncol(y)
  if(is.null(pch)){
    pch <- c(1 : nrows.y)
  }
  if(is.null(col)){
    col <- c(1 : nrows.y)
  }
  plot(NULL, NULL, xlim = c(0, ncol(y) + 0.5), cex.main = 1.2, main = main, cex = 0.5, 
       xlab = xlab, ylab = ylab, cex.axis = .80, type = "b", ylim = c(0, max(1, max(y))), las = 1)
  for (j in 1 : nrows.y) {
    lines(x = 1 : ncols.y, y = y[j, ], type = "b", lty = 1, pch = pch[j], col = col[j])
  }
  if(ylab == "FDR"){
    legend(legend.x, legend.y, legend = legendnames, xjust = 1, yjust = 1,
           cex = 0.5, lty = rep(1, nrows.y), pch = pch, col = col)
    for (j in 1 : nrows.y) {
      abline(h = ablines[j], lty = 2)
    }
  }
  if(ylab == "Power"){
    legend(legend.x, legend.y, legend = legendnames, xjust = 1, yjust = 0, 
           cex = 0.5, lty = rep(1, nrows.y), pch = pch, col = col)
  }
  if(ylab == "selection accuracy"){
    legend(legend.x, legend.y, legend = legendnames, xjust = 1, yjust = 1, 
           cex = 0.5, lty = rep(1, nrows.y), pch = pch, col = col)
  }
}

#----------------------------------------------------------------
# A plot function: display the relationship between x and y
# Args: 
#    x: a vector or a matrix
#    y: a vector or a matrix should have the same dimension with x
#    main: main title, such as "(n,p,s0)=(200,100,10)".
#    xlab: a character, such as "FDR"
#    ylab: a character, such as "Power"
#    legendnames: should be a character vector with length = ncol(x); such as c("BCknockoff","BCknockoff+","geometric")
#    pch: default is NULL
#    col: default is NULL; col could be c("darkorchid","orange","deepskyblue","deeppink","green3","firebrick1","blue")
#
# Output:
#    a plot displaying the relationship between x and y

myPlot.FP <- function(x, y, main, xlab = "FDR", ylab = "Power", 
                      legend.x = 0, legend.y = 0.8, legendnames = NULL, 
                      pch = NULL, col = NULL){
  if(nrow(x) != nrow(y) || ncol(x) != ncol(y)){
    stop("x and y should have the same dimension")
  }
  if(is.null(legendnames)){
    stop("legendnames should be a character vector with the length = ncol(x)")
  }
  nrows.x <- nrow(x)
  ncols.x <- ncol(x)
  if(is.null(col)){
    if(ncols.x > 7){
      col <- rainbow(ncols.x)  # col could be c("darkorchid","orange","deepskyblue","deeppink","green3","firebrick1","blue")
    }else{
      col <- c("darkorchid","orange","deepskyblue","deeppink","green3","firebrick1","blue")[1 : ncols.x]
    }
  }
  if(is.null(pch)){
    pch <- c(15 : (15 + ncols.x))
  }
  if(ylab == "Power"){
    plot(NULL, NULL, xlim = c(min(0, min(x)), max(1, max(x))), cex = 1, xlab = "", ylab = "", 
         cex.axis = 1.2, type = "n", ylim = c(min(0, min(y)), max(1, max(y))), las = 1)
  }
  if(ylab == "selection accuracy"){
    plot(NULL, NULL, xlim = c(min(0, min(x)), max(1, max(x))), cex = 1, xlab = "", ylab = "", 
         cex.axis = 1.2, type = "n", ylim = c(min(0, min(y)), max(1, max(y))), las = 1)
  }
  for (j in 1 : ncols.x) {
    points(x[ , j], y[ , j], pch = pch[j], cex = 0.5, col = col[j])
  }
  title(ylab = ylab, xlab = xlab, cex.lab = 1.5, line = 3)
  grid(col = "darkgray")
  title(main = main, cex.main = 2, line = 1, font.main = 2)
  legend(legend.x, legend.y, legend = legendnames, cex = 1.5, col = col,              
         pch = pch, bty = "n", y.intersp = 0.7)
}

#----------------------------------------------------------------
# A plot function: display the relationship between each columns of x and y
# Args: 
#    x: a vector or a matrix
#    y: a vector with length = nrow(x)
#    main: main title, such as "(n,p,s0)=(200,100,10)".
#    xlab: a character, such as "Actual FDR"
#    ylab: a character, such as "Target FDR"
#    legendnames: should be a character vector with length = ncol(x); such as c("BCknockoff","BCknockoff+","geometric")
#    pch: default is NULL
#    col: default is NULL; col could be c("darkorchid","orange","deepskyblue","deeppink","green3","firebrick1","blue")
#
# Output:
#    a plot displaying the relationship between each columns of x and y

myPlot.x.Target <- function(x, y, main, ylab = "Actual FDR", xlab = "Target FDR",
                            legend.x = 0, legend.y = 0.8, legendnames = NULL, 
                            pch = NULL, col = NULL){
  if(is.null(legendnames)){
    stop("legendnames should be a character vector with the length = ncol(x)")
  }
  nrows.y <- nrow(y)
  ncols.y <- ncol(y)
  if(is.null(col)){
    if(ncols.y > 7){
      col <- rainbow(ncols.y)  # col could be c("darkorchid","orange","deepskyblue","deeppink","green3","firebrick1","blue")
    }else{
      col <- c("darkorchid","orange","deepskyblue","deeppink","green3","firebrick1","blue")[1 : ncols.y]
    }
  }
  if(is.null(pch)){
    pch <- c(15 : (15 + ncols.y - 1))
  }
  plot(NULL, NULL, xlim = c(min(0, min(x)), max(1, max(x))), cex = 1, xlab = "", ylab = "",
      cex.axis = 1.2, type = "n", ylim = c(min(0, min(y)), max(1, max(y))), las = 1)
  for (j in 1 : ncols.y) {
    points(x, y[ , j], pch = pch[j], cex = 0.5, col = col[j])
  }
  if(ylab == "Actual FDR"){
    abline(a = 0, b = 1)
  }
  title(ylab = ylab, xlab = xlab, line = 3, cex.lab = 1.5)
  title(main = main, cex.main = 2, line = 1, font.main = 2)
  legend(legend.x, legend.y, legend = legendnames, cex = 1.5, col = col,           
         pch = pch, bty = "n", y.intersp = 0.7)
}












