#! R 
# Filename: DataUS012018.R

# Import Phylum level data before January 2018
phylum.raw <- ImportData("RealData/ag-cleaned_L2.txt", threshold = 0, replace.zero = T)
phylum.raw <- as.matrix(phylum.raw)

#--------------------------
#--"COUNTRY_OF_BIRTH"--US
indx.cob <- which(colnames(phylum.raw) == "COUNTRY_OF_BIRTH")
sam.US <- which(phylum.raw[ , indx.cob] == "United States")
length(sam.US)

#---Study on US samples---
phylum.US <- phylum.raw[sam.US, -indx.cob]

#--"BMI"--#
indx.BMI <- which(colnames(phylum.US) == "BMI")
del.indx.bmi <- c(which(phylum.US[, indx.BMI] == "Unknown"), which(phylum.US[, indx.BMI] == "Unspecified"))
if(length(del.indx.bmi) != 0){
  phylum.US <- phylum.US[-del.indx.bmi, ]
}
bmi = as.numeric(phylum.US[,indx.BMI])
ext.indx <- which(bmi <= 60 & bmi > 15)     # 15 <= BMI <= 60
phylum.US <- phylum.US[ext.indx, ]

#--"AGE_YEARS"--#
indx.age <- which(colnames(phylum.US) == "AGE_YEARS")
del.indx.age <- c(which(phylum.US[, indx.age] == "Unknown"), which(phylum.US[, indx.age] == "Unspecified"))
if(length(del.indx.age) != 0){
  phylum.US <- phylum.US[-del.indx.age, ]
}
US.age <- as.numeric(phylum.US[, indx.age])      
ext.age.20_69 <- which(US.age >= 20 & US.age <= 69) # 20 <= AGE <= 69
phylum.US <- phylum.US[ext.age.20_69, ]

#--"VIOSCREEN_GENDER"--#
indx.gender <- which(colnames(phylum.US) == "VIOSCREEN_GENDER")

#--"BMI" index--#
ext.indx.1 <- which(phylum.US[,indx.BMI] < 18.5 & phylum.US[,indx.BMI] >= 15)
ext.indx.2 <- which(phylum.US[,indx.BMI] <= 25 & phylum.US[,indx.BMI] >= 18.5)
ext.indx.3 <- which(phylum.US[,indx.BMI] <= 30 & phylum.US[,indx.BMI] > 25)
ext.indx.4 <- which(phylum.US[,indx.BMI] <= 60 & phylum.US[,indx.BMI] > 30)

#--X.raw and Y.BMI for US sample--#
X.raw.US <- phylum.US[ , -c(indx.BMI, indx.age, indx.gender)]
X.raw.US <- apply(X.raw.US, 2, function(x) as.numeric(x))
Y.BMI.US <- as.numeric(phylum.US[, indx.BMI])





