#! R
# Filename: ImportData_genus.R
ImportData_genus <- function(filename, threshold = 0, replace.zero = TRUE) {
  # This function reads a raw file and returns relevant information from the 
  # data file depend on the need of the user. 
  #
  # Args:
  #   filename: Character string specifying the file name (with or without path) 
  #             of the raw dataset.
  #   threshold: A numeric value indicates the threshold for genus data. If the
  #              percentage of non-zero entry for a genus is lower than this
  #              value, it will not be included in the final returned value. 
  #   replace.zero: a logical value. If TURE, the zero terms will be replaced 
  #                 by the 0.5 times the minimal abundance of all the genus.
  #             
  # Returns: A list containing:
  #   genus:   a data.frame of numeirc genus data with correct names.
  
  data.raw <- read.table(filename, sep = "\t", header = T,
                         fill = T, comment.char="", quote = "")
  genus.indx <- which(startsWith(colnames(data.raw), "k__Bacteria")) 
  if (length(genus.indx) == 0) {
    stop("data does not contain 'k__Bacteria' with default names")
  }
  genus <- apply(data.raw[,genus.indx], 2, as.numeric)
  # collect only the genus level data (it should modified for different data)
  genus <- genus[, 270 : 2290]
  
  # select the phyla that its samples have at most 'threshold' percentage zero items
  if(threshold != 0) {
    phy.nz <- apply(genus != 0, 2, sum) / nrow(genus)
    genus <- genus[ ,phy.nz > threshold] 
    genus <- t(apply(genus, 1, function(x) x/sum(x)))
  }
  # replace the zero terms by the 0.5 times the minimal abundance of all the genus
  if(replace.zero){
    for(row in 1:nrow(genus)){
      row.zero.indx <- which(genus[row, ] == 0)
      if(length(row.zero.indx) != 0){
        genus[row, row.zero.indx] <- 0.5 * min(genus[row, -row.zero.indx])
      }
    }
    genus <- t(apply(genus, 1, function(x) x/sum(x)))
  }
  # extract the data of variables "COUNTRY_OF_BIRTH", "AGE_YEARS", "VIOSCREEN_GENDER" and "BMI"
  Names <- c(colnames(genus), "COUNTRY_OF_BIRTH", "AGE_YEARS", "VIOSCREEN_GENDER", "BMI")             ### added by Fang
  cob.indx <- c(which(colnames(data.raw) == "COUNTRY_OF_BIRTH"),
                which(colnames(data.raw) == "AGE_YEARS"),
                which(colnames(data.raw) == "VIOSCREEN_GENDER"),
                which(colnames(data.raw) == "BMI"))             
  genus <- cbind(genus, data.raw[ , cob.indx])                  
  names(genus) <- Names                                         
  return(genus)
}




