#! R 
# Filename: plotlinear.R

#-------------------------------------------------------#
#---plot the curves of FDR and Power w.r.t Target FDR---#
#-------------------------------------------------------#
#----Lasso----#
#-------------#
#--define file names and legend names--#  
filename1 = "Results/linear2/FDRm.pdf"  
filename2 = "Results/linear2/Powerm.pdf" 


#------plots of KO and modified AKO-------#
FDR2 <- cbind(FDR, FDR.AKO.m)
Pwr2 <- cbind(Pwr, Pwr.AKO.m)
legendnames.m = c("KO", "AKO")

# display the relationship between FDR and TargetFDR
pdf(file = filename1, width = 5, height = 5)
myPlot(x = fdr, y = FDR2, main = "linear", ylab = "Actual FDR", legendnames = legendnames.m)
dev.off()

# display the relationship between Power and TargetFDR 
pdf(file = filename2, width = 5, height = 5)
myPlot(x = fdr, y = Pwr2, main = "linear", ylab = "Power", legendnames = legendnames.m,
       legend.x = 0.7, legend.y = 0.2)
dev.off()

