#! R
# Filename: ImportData.R
ImportData <- function(filename, threshold = 0, replace.zero = TRUE) {
  # This function reads a raw file and returns relevant information from the 
  # data file depend on the need of the user. 
  #
  # Args:
  #   filename: Character string specifying the file name (with or without path) 
  #             of the raw dataset.
  #   threshold: A numeric value indicates the threshold for phylum data. If the
  #              percentage of non-zero entry for a phylum is lower than this
  #              value, it will not be included in the final returned value. 
  #   replace.zero: a logical value. If TURE, the zero terms will be replaced 
  #                 by the 0.5 times the minimal abundance of all the phylum.
  #             
  # Returns: A list containing:
  #   phylum:   a data.frame of numeirc phylum data with correct names.

  data.raw <- read.table(filename, sep = "\t", header = T,
                         fill = T, comment.char="", quote = "")
  phylum.indx <- which(startsWith(colnames(data.raw), "k__Bacteria")) 
  if (length(phylum.indx) == 0) {
    stop("data does not contain phylum with default names")
  }
  phylum <- apply(data.raw[,phylum.indx], 2, as.numeric)
  names.temp <- colnames(phylum)
  for (i in 1:length(names.temp)) {
    foo <- strsplit(names.temp[i], "_")[[1]]
    foo <- foo[length(foo)]
    # Note: there might be cases that the name of phylum is not specified.
    if (startsWith(foo, ".")) {
      colnames(phylum)[i] <- substr(foo, 2, nchar(foo) - 1)
    } else {
      colnames(phylum)[i] <- foo
    }
  }
  # change the name of bacteria which has the same with another
  repeatname.indx = which(colnames(phylum) == "Caldithrix")[2]
  colnames(phylum)[repeatname.indx] <- "Caldithrix.1"
  # select the phyla that its samples have at most 'threshold' percentage zero items
  if(threshold != 0) {
    phy.nz <- apply(phylum != 0, 2, sum) / nrow(phylum)
    phylum <- phylum[ ,phy.nz > threshold] 
    phylum <- t(apply(phylum, 1, function(x) x/sum(x)))
  }
  # replace the zero terms by the 0.5 times the minimal abundance of all the phylum
  if(replace.zero){
    for(row in 1:nrow(phylum)){
      row.zero.indx <- which(phylum[row, ] == 0)
      if(length(row.zero.indx) != 0){
        phylum[row, row.zero.indx] <- 0.5 * min(phylum[row, -row.zero.indx])
      }
    }
    phylum <- t(apply(phylum, 1, function(x) x/sum(x)))
  }
  # extract the data of variables "COUNTRY_OF_BIRTH", "AGE_YEARS", "VIOSCREEN_GENDER" and "BMI"
  Names <- c(colnames(phylum), "COUNTRY_OF_BIRTH", "AGE_YEARS", "VIOSCREEN_GENDER", "BMI")             ### added by Fang
  cob.indx <- c(which(colnames(data.raw) == "COUNTRY_OF_BIRTH"),
                which(colnames(data.raw) == "AGE_YEARS"),
                which(colnames(data.raw) == "VIOSCREEN_GENDER"),
                which(colnames(data.raw) == "BMI"))             
  phylum <- cbind(phylum, data.raw[ , cob.indx])                  
  names(phylum) <- Names                                         
  return(phylum)
}




