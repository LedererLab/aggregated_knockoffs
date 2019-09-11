# Code written by: Lun Li and Johannes Lederer
ImportData <- function(filename, codebook, all.data = FALSE,
                       clean = TRUE, threshold = 0) {
  # This function reads a raw file and returns relevant information from the 
  # data file depend on the need of the user. The user can chose to obtain:
  # only phyla data; phyla data and uncleaned characteristics of the subjects,
  # or phyla data and cleaned characteristics. Note that the cleaned data will
  # have less subject than original data, as the subjects with unreasonable
  # BMI's are removed from the returned stet. 
  #
  # Args:
  #   filename: Character string specifying the file name (with or without path) 
  #             of the raw dataset.
  #   codebook: Character  string specifying the file name (with or without 
  #             path) of the sorted code book.  
  #   all.data: Logical. If = TRUE, return all variables in section 'Returns';
  #             otherwise only return the phylum data. Default = FALSE.
  #   clean:    Logical. If = TRUE, the returned phylum and features will be 
  #             cleaned with funcion CleanData.
  #   threshold:A numeric value indicates the threshold for phylum data. If the
  #             percentage of non-zero entry for a phylum is lower than this
  #             value, it will not be included in the final returned value. 
  #             
  # Returns: A list containing:
  #   phylum:   a data.frame of numeirc phylum data with correct names.
  #   features: data.frame that includes features of each subject.
  #   entry:    data.frame of entry names from the codebook, data type, and 
  #             index of each term.
  #   descript: a descriptive list for each entry. Can be accessed by calling 
  #             'descript$...'.
  data.raw <- read.table(filename, sep = "\t", header = T,
                         fill = T, comment.char="", quote = "")
  phylum.indx <- c(which(startsWith(colnames(data.raw), "k__Bacteria")), which(colnames(data.raw) == "BMI"))  ### edited by Fang
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
  phylum <- as.data.frame(phylum)
  if(threshold != 0) {
    phy.nz <- apply(phylum != 0, 2, sum) / nrow(phylum)
    phylum <- phylum[,phy.nz > threshold] 
    phylum <- t(apply(phylum, 1, function(x) x/sum(x)))
  }
  if(!all.data) {
    Names <- c(colnames(phylum), "COUNTRY_OF_BIRTH", "AGE_YEARS", "VIOSCREEN_GENDER")             ### added by Fang
    cob.indx <- c(which(colnames(data.raw) == "COUNTRY_OF_BIRTH"),
                  which(colnames(data.raw) == "AGE_YEARS"),
                  which(colnames(data.raw) == "VIOSCREEN_GENDER"))  ### added by Fang
    phylum <- cbind(phylum, data.raw[ , cob.indx])                  ### added by Fang
    names(phylum) <- Names                                          ### added by Fang
    return(list(phylum = phylum))
  } else {
    # import entry names from the codebook
    entry <- read.csv(codebook, sep = ",", fill = T)
    descript <- as.list(as.character(entry[, 3]))
    names(descript) <- as.character(entry[, 1])
    entry <- cbind(entry[, 1:2], 1:nrow(entry))
    entry[1:2] <- apply(entry[1:2], 2, as.character)
    entry[,2] <- tolower(entry[, 2])
    names(entry) <- c(names(entry)[1:2], "index")
    entry <- entry
    descript <- descript
    # get required data
    matched <- match(entry[,1],names(data.raw))
    features <- data.raw[, matched]
    rm(data.raw)
  }
  if(clean) {
    cleanHelper <- function(x) {
      # This helper function combines the same levels from a factor and return
      # a variable with regulated levels  
      #
      # Args:
      #   x: The passed in variable that needs fix
      #
      # Returns: A variable with regulated levels
      x <- tolower(as.character(x))
      x[x == "yes"] <- "true"
      x[x == "no"] <- "false"
      x[x == "unknown" | x == "unspecified" | x == "no_data" |
          x == "not sure" | x == "other"] <- NA
      return(factor(x))
    }
    
    num.indx <- which(entry[,2] == "int"|entry[,2] == "float")
    for(i in num.indx) {
      features[,i] <- suppressWarnings(as.numeric(as.character(features[,i])))
    }
    
    # Only use the data with reasonable weight/height
    use.indx <- !is.na(features$BMI) # only take data that has a BMI value
    use.indx <- use.indx & features$WEIGHT_KG >= 40 & features$WEIGHT_KG <= 150 & 
      features$HEIGHT_CM >= 140 & features$HEIGHT_CM <= 220 &
      !is.na(features$WEIGHT_KG) & !is.na(features$HEIGHT_CM)
    
    # fix all date type
    features$COLLECTION_DATE <- as.Date(as.character(features$COLLECTION_DATE), 
                                        "%m/%d/%Y")
    ## fix all boolean type boolean
    bool.indx <- which(entry[,2] == "bool")
    for(i in bool.indx) {
      features[,i] <- cleanHelper(features[,i])
    }
    
    ## fix all strings
    str.indx <- which(entry[,2] == "str")
    str.indx <- str.indx[-length(str.indx)]
    for(i in str.indx) {
      features[,i] <- cleanHelper(features[,i])
    }
    temp <- tolower(as.character(features$VIOSCREEN_CALCIUM_FREQ))
    temp[temp == "unknown" | temp == "unspecified" | temp == "no_data"] <- NA
    temp[temp == "less than once per week"] <- 0.5
    features$FLU_VACCINE_DATE <- cleanHelper(features$FLU_VACCINE_DATE)
    features$VIOSCREEN_CALCIUM_FREQ <- as.numeric(temp)
    features$VIOSCREEN_CALCIUM_DOSE <- as.factor(features$VIOSCREEN_CALCIUM_DOSE)
    features$VIOSCREEN_MULTIVITAMIN_FREQ <- 
      as.factor(features$VIOSCREEN_MULTIVITAMIN_FREQ)
    features <- features[use.indx,]
    phylum <- as.data.frame(phylum[use.indx,])
  }
  return(list(phylum = phylum, features = features, entry = entry,
              descript = descript))
}

cleanHelper <- function(x) {
  # This helper function combines the same levels from a factor and return
  # a variable with regulated levels  
  #
  # Args:
  #   x: The passed in variable that needs fix
  #
  # Returns: A variable with regulated levels
  x <- tolower(as.character(x))
  x[x == "yes"] <- "true"
  x[x == "no"] <- "false"
  x[x == "unknown" | x == "unspecified" | x == "no_data" |
      x == "not sure" | x == "other"] <- NA
  return(factor(x))
}
