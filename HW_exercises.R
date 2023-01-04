library("HardyWeinberg")
####Silene data
Seedlings <- c(AA = 79, AB = 43, BB  = 21)
HW.test.Seedlings <- HWChisq(Seedlings, cc = 0, verbose = TRUE)
HW.test.Seedlings
Adults<-  c(AA = 70, AB = 60, BB  = 13)
HW.test.Adults <- HWChisq(Adults, cc = 0, verbose = TRUE)
HW.test.Adults
#### SNP rs16891982 data 
AFR <- c("CC" = 617, "CG" = 41, "GG" = 3)
HW.test.AFR <- HWChisq(AFR, cc = 0, verbose = TRUE)
HW.test.AFR
EUR  <- c("CC" = 4, "CG" =  54, "GG" = 445)
HW.test.EUR <- HWChisq(EUR, cc = 0, verbose = TRUE)
HW.test.EUR
###Human HapMap data
setwd("~/exercises/HardyWeinberg")
SNP_data <-as.matrix(read.table("CEU_500.hw", header=TRUE))
SNP_sim_data <-as.matrix(read.table("CEU_500.sim", header=TRUE))
#How does the data look like?
head(SNP_data)
#Plot the genotype frequencies for each of the 500 populations
#in a deFinetti diagram and indicate populations that differ 
#signicantly from Hardy-Weinberg proportions.
par(mfrow=c(1,2))
HWTernaryPlot(SNP_data,   region = 1, 
              curvecols=c("black", "red","green","black","purple"), 
              vbounds = FALSE, main ="Original CEU data", 
              cex = 0.5,  cex.main=1.5, font.main = 1)
HWTernaryPlot(SNP_sim_data, region = 1, 
              curvecols=c("black", "red","green","black","purple"), 
              vbounds = FALSE, main ="Simulated data", 
              cex = 0.5, cex.main=1.5, font.main = 1)
#### count the number of significant tests
#First make a vector with chi test values 
#using the function HWChisqStats
chitest<-HWChisqStats(SNP_data,pvalues=FALSE)
#Then make a vector with significant tests
sigchitest <-chitest[chitest>3.84]
length(sigchitest)
#####QQ plots
par(mfrow=c(1,2))
HWQqplot(SNP_data,   
         main ="Q-Q plot for original data",  
         cex.main=1.5, font.main = 1)
HWQqplot(SNP_sim_data, 
         main ="Q-Q plot for simulated data", 
         cex.main=1.5, font.main = 1)

