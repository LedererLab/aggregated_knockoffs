# Code written by: Lun Li and Johannes Lederer
TransformData <- function(rawdata, method = 'LiLederer', g = 1) {
  # This function accepts a set of data and perform request transformation.
  # There are two options: our method of transformation (Lilederer), or Kurtz's
  # method (Kurtzetal).
  #
  # Args:
  #   rawdata: A matrix of data.
  #   method : Transformation with 2 options: 'LiLederer' and 'Kurtzetal'.
  #   g      : Integer indicates the number of clusters used for the GMM model.
  #            Only used for 'LiLederer'.
  #             
  # Returns  : A matrix of transformed data using the requested method.
  if (class(rawdata) != "matrix") {
    stop("invalid data type")
  }
  if (method == 'LiLederer') {
    return(centerByGMM(getClr(rawdata), g)) # our method
  } else if(method == 'Kurtzetal') {
    return(getClr(rawdata, FALSE))
  } else{
    stop("unknown method")
  }
}