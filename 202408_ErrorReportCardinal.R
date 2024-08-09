library(BiocManager)
library(Cardinal)

# Minimal Preprocessing ---------------------------------------------------

zeroA <-readMSIData("0A.imzML")
zeroA <- zeroA |>
  normalize(method="tic") |>
  peakProcess(SNR=3,filterFreq=FALSE) # 885 ref peaks w/ tol = 234 ppm

zeroA_ssc <- spatialShrunkenCentroids(zeroA, r = 1, s=0,k=2, weights = "gaussian")
image(zeroA_ssc) # 2 regions - one for tissue and another for tape/background

leaf<-subsetPixels(zeroA, zeroA_ssc@model[["class"]]==2) 
bkgd<-subsetPixels(zeroA, zeroA_ssc@model[["class"]]==1)
leaf_ssc <- spatialShrunkenCentroids(leaf,r=1, s=8, 
                                     k=7, weights = "gaussian") #Segments tissue into various regions as expected
bkgd_ssc <- spatialShrunkenCentroids(bkgd, r=1, s=8, 
                                     k=7, weights = "gaussian") #Err 
#Error in h(simpleError(msg, call)) : 
# error in evaluating the argument 'x' in selecting a method for function 't':
# 'x' must be an array of at least two dimensions

zeroA$regions <- makeFactor(leaf = (zeroA_ssc@model[["class"]]==2),
                            bkgd = (zeroA_ssc@model[["class"]]==1))
leaf2<-subsetPixels(zeroA, regions == "leaf")
bkgd2<-subsetPixels(zeroA, regions == "bkgd")
leaf2_ssc <- spatialShrunkenCentroids(leaf2,r=1, s=8, 
                                     k=7, weights = "gaussian") #Segments tissue into various regions as expected
bkgd2_ssc <- spatialShrunkenCentroids(bkgd2, r=1, s=8, 
                                     k=7, weights = "gaussian") #Err
#Error in h(simpleError(msg, call)) : 
# error in evaluating the argument 'x' in selecting a method for function 't':
# 'x' must be an array of at least two dimensions

# Specified Tolerance -----------------------------------------------------

zeroA2 <-readMSIData("0A.imzML")
zeroA2 <- zeroA2 |>
  normalize(method="tic") |>
  peakProcess(SNR=3,filterFreq=FALSE, tolerance = 5, units="ppm") #1389 ref peaks w/ tol = 5 ppm

zeroA2_ssc <- spatialShrunkenCentroids(zeroA2, r = 1, s=0,k=2, weights = "gaussian")
image(zeroA2_ssc) # 2 regions - one for tissue and another for tape/background 

leaf3<-subsetPixels(zeroA2, zeroA2_ssc@model[["class"]]==2)
bkgd3<-subsetPixels(zeroA2, zeroA2_ssc@model[["class"]]==1)
leaf3_ssc <- spatialShrunkenCentroids(leaf3,r=1, s=8, 
                                     k=7, weights = "gaussian") #Err
#Error in h(simpleError(msg, call)) : 
# error in evaluating the argument 'x' in selecting a method for function 't':
# 'x' must be an array of at least two dimensions
bkgd3_ssc <- spatialShrunkenCentroids(bkgd3, r=1, s=8, 
                                     k=7, weights = "gaussian") #Err
#Error in h(simpleError(msg, call)) : 
# error in evaluating the argument 'x' in selecting a method for function 't':
# 'x' must be an array of at least two dimensions
