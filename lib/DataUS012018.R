# Import Phylum level data before January 2018
phylum.raw <- ImportData("RealData/ag-cleaned_L2.txt", threshold = 0, 
                         all.data = F)[[1]]
#---------------------------------------------------------------------
## "COUNTRY_OF_BIRTH"---US
indx.cob <- which(colnames(phylum.raw) == "COUNTRY_OF_BIRTH")
sam.US <- which(phylum.raw[ , indx.cob] == "United States")
length(sam.US)
#------------Study on US samples---------------
phylum.US <- phylum.raw[sam.US, -indx.cob]
## "BMI"
indx.BMI <- which(colnames(phylum.US) == "BMI")
ext.indx <- which(phylum.US[,indx.BMI] <= 60 & phylum.US[,indx.BMI] > 15)     # 15 <= BMI <= 60
phylum.US <- phylum.US[ext.indx, ]
## "AGE_YEARS"
indx.age <- which(colnames(phylum.US) == "AGE_YEARS")
del.indx <- c(which(phylum.US[, indx.age] == "Unknown"), which(phylum.US[, indx.age] == "Unspecified"))
phylum.US <- phylum.US[-del.indx, ]

US.age <- phylum.US[, indx.age]
US.age.total <- as.numeric(levels(US.age)[US.age])       
ext.age.20_69 <- which(US.age.total >= 20 & US.age.total <= 69) # 20 <= AGE <= 69
phylum.US <- phylum.US[ext.age.20_69, ]
## "VIOSCREEN_GENDER"
indx.gender <- which(colnames(phylum.US) == "VIOSCREEN_GENDER")

# "BMI"
ext.indx.1 <- which(phylum.US[,indx.BMI] < 18.5 & phylum.US[,indx.BMI] >= 15)
ext.indx.2 <- which(phylum.US[,indx.BMI] <= 25 & phylum.US[,indx.BMI] >= 18.5)
ext.indx.3 <- which(phylum.US[,indx.BMI] <= 30 & phylum.US[,indx.BMI] > 25)
ext.indx.4 <- which(phylum.US[,indx.BMI] <= 60 & phylum.US[,indx.BMI] > 30)

# X.raw and Y.BMI for US sample
X.raw.US <- phylum.US[ , -c(indx.BMI, indx.age, indx.gender)]
Y.BMI.US <- phylum.US[ , indx.BMI]
