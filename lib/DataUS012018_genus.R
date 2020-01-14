#! R 
# Filename: DataUS012018_genus.R

# Import Phylum level data before January 2018
genus.raw <- ImportData_genus("RealData/ag-cleaned_L6.txt", threshold = 0, replace.zero = T)
genus.raw <- as.matrix(genus.raw)

#--------------------------
#--"COUNTRY_OF_BIRTH"--US
indx.cob <- which(colnames(genus.raw) == "COUNTRY_OF_BIRTH")
sam.US <- which(genus.raw[ , indx.cob] == "United States")
length(sam.US)

#---Study on US samples---
genus.US <- genus.raw[sam.US, -indx.cob]

#--"BMI"--#
indx.BMI <- which(colnames(genus.US) == "BMI")
del.indx.bmi <- c(which(genus.US[, indx.BMI] == "Unknown"), which(genus.US[, indx.BMI] == "Unspecified"))
if(length(del.indx.bmi) != 0){
  genus.US <- genus.US[-del.indx.bmi, ]
}
bmi = as.numeric(genus.US[,indx.BMI])
ext.indx <- which(bmi <= 60 & bmi > 15)     # 15 <= BMI <= 60
genus.US <- genus.US[ext.indx, ]

#--"AGE_YEARS"--#
indx.age <- which(colnames(genus.US) == "AGE_YEARS")
del.indx.age <- c(which(genus.US[, indx.age] == "Unknown"), which(genus.US[, indx.age] == "Unspecified"))
if(length(del.indx.age) != 0){
  genus.US <- genus.US[-del.indx.age, ]
}
US.age <- as.numeric(genus.US[, indx.age])      
ext.age.20_69 <- which(US.age >= 20 & US.age <= 69) # 20 <= AGE <= 69
genus.US <- genus.US[ext.age.20_69, ]

#--"VIOSCREEN_GENDER"--#
indx.gender <- which(colnames(genus.US) == "VIOSCREEN_GENDER")

#--"BMI" index--#
ext.indx.1 <- which(genus.US[,indx.BMI] < 18.5 & genus.US[,indx.BMI] >= 15)
ext.indx.2 <- which(genus.US[,indx.BMI] <= 25 & genus.US[,indx.BMI] >= 18.5)
ext.indx.3 <- which(genus.US[,indx.BMI] <= 30 & genus.US[,indx.BMI] > 25)
ext.indx.4 <- which(genus.US[,indx.BMI] <= 60 & genus.US[,indx.BMI] > 30)

#--X.raw and Y.BMI for US sample--#
X.raw.US <- genus.US[ , -c(indx.BMI, indx.age, indx.gender)]
X.raw.US <- apply(X.raw.US, 2, function(x) as.numeric(x))
Y.BMI.US <- as.numeric(genus.US[, indx.BMI])





