#! R 
# Filename: plotlinear.R

#-------------------------------------------------------#
#---plot the curves of FDR and Power w.r.t Target FDR---#
#-------------------------------------------------------#
#----Lasso----#
#-------------#

#------plots of KO and modified AKO-------#
FDR2 <- cbind(FDR, FDR.AKO.m)
Pwr2 <- cbind(Pwr, Pwr.AKO.m)
legendnames.m = c("KO", "AKO")
#--define file names and legend names--#  
# filename1 = "Results/linear2/FDRm.pdf"  
# filename2 = "Results/linear2/Powerm.pdf" 
filename1 = "Results/NIPS_response/linear/FDR15.pdf"  
filename2 = "Results/NIPS_response/linear/Power15.pdf" 

# display the relationship between FDR and TargetFDR
pdf(file = filename1, width = 5, height = 5)
myPlot(x = fdr, y = FDR2, main = "linear", ylab = "Actual FDR", legendnames = legendnames.m)
dev.off()

# display the relationship between Power and TargetFDR 
pdf(file = filename2, width = 5, height = 5)
myPlot(x = fdr, y = Pwr2, main = "linear", ylab = "Power", legendnames = legendnames.m,
         legend.x = 0.7, legend.y = 0.2)
dev.off()



#-------------------------------------------------------------#
#--plots of KO, modified AKO, average AKO, and modified average AKO--#
FDR2 <- cbind(FDR, FDR.AKO.m, FDR.AKO.aver.m)
Pwr2 <- cbind(Pwr, Pwr.AKO.m, Pwr.AKO.aver.m)
legendnames.m = c("KO", "AKO", "AKO.ave")
#--define file names and legend names--#  
filename3 = "Results/NIPS_response/linear/FDR17.pdf"  
filename4 = "Results/NIPS_response/linear/Power17.pdf" 

# display the relationship between FDR and TargetFDR
pdf(file = filename3, width = 5, height = 5)
myPlot.m(x = fdr, y = FDR2, main = "linear", ylab = "Actual FDR", legendnames = legendnames.m)
dev.off()

# display the relationship between Power and TargetFDR 
pdf(file = filename4, width = 5, height = 5)
myPlot.m(x = fdr, y = Pwr2, main = "linear", ylab = "Power", legendnames = legendnames.m,
       legend.x = 0.55, legend.y = 0.3)
dev.off()



#--------------------------------------#
#--plots of KO, modified AKO, and MKO--#
FDR2 <- cbind(FDR, FDR.AKO.m, FDR.MKO)
Pwr2 <- cbind(Pwr, Pwr.AKO.m, Pwr.MKO)
legendnames.m = c("KO", "AKO", "MKO")
#--define file names and legend names--#  
filename5 = "Results/AISTATS/linear/FDR2.pdf"  
filename6 = "Results/AISTATS/linear/Power2.pdf" 

# display the relationship between FDR and TargetFDR
pdf(file = filename5, width = 5, height = 5)
myPlot.m(x = fdr, y = FDR2, main = "linear", ylab = "Actual FDR", legendnames = legendnames.m)
dev.off()

# display the relationship between Power and TargetFDR 
pdf(file = filename6, width = 5, height = 5)
myPlot.m(x = fdr, y = Pwr2, main = "linear", ylab = "Power", legendnames = legendnames.m,
         legend.x = -0.05, legend.y = 1.05)
dev.off()

