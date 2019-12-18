#! R 
# Filename: plotlogistic.R

#--------------------------------------------------------------------#
#---plot the curves of FDR and Power w.r.t Target FDR---#
#--------------------------------------------------------------------#
#----logistic----#
#----------------#
# define file names and legend names
filename1 = "Results/logistic2/FDR.pdf"  
filename2 = "Results/logistic2/Power.pdf"     
filename3 = "Results/logistic2/FDRm.pdf"  
filename4 = "Results/logistic2/Powerm.pdf"      
legendnames = c("KO", "AKO")
legendnames.m = c("KO", "modified AKO")

#------AKO------#
FDR1 <- cbind(FDR, FDR.AKO)
Pwr1 <- cbind(Pwr, Pwr.AKO)
SA1 <- cbind(SA, SA.AKO)

# display the relationship between FDR and TargetFDR
pdf(file = filename1, width = 5, height = 5)
myPlot.FDR(x = fdr, y = FDR1, main = "logistic",legendnames = legendnames)
dev.off()

pdf(file = filename2, width = 5, height = 5)
myPlot.Power(x = fdr, y = Pwr1, main = "logistic", legendnames = legendnames)
dev.off()


#------modified AKO-------#
FDR2 <- cbind(FDR, FDR.AKO.m)
Pwr2 <- cbind(Pwr, Pwr.AKO.m)
SA2 <- cbind(SA, SA.AKO.m)

# display the relationship between FDR and TargetFDR
pdf(file = filename3, width = 5, height = 5)
myPlot.FDR(x = fdr, y = FDR2, main = "logistic", legendnames = legendnames.m)
dev.off()

# display the relationship between Power and TargetFDR 
pdf(file = filename4, width = 5, height = 5)
myPlot.Power(x = fdr, y = Pwr2, main = "logistic", legendnames = legendnames.m)
dev.off()





