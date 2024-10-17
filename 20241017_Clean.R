#Cardinal v3.6.5 w/ R/4.4.0 RStudio/2024.04.1+748
.libPaths("/work/LAS/yjlee-lab/buckm065/Rlibs44")
setwd("/work/LAS/yjlee-lab/buckm065/CardinalProject/Replicates")
library(BiocManager)
library(matter)
library(Cardinal)
library(gplots)
library(ggplot2)
library(Rtsne)

# col=matter::cpal("Cividis")

# Unlabeled ---------------------------------------------------------------

zeroA <-readMSIData("0A.imzML") #Read Unlabeled DW
zeroA <- zeroA |>
  normalize(method="tic") |> #TIC Normalization
  peakProcess(SNR=3,tolerance = 7, units = "ppm") #Peak Picking and Alignment

zeroA_ssc <- spatialShrunkenCentroids(zeroA, r = 1, s=0,k=2, 
                                      weights = "gaussian")
image(zeroA_ssc)

zeroA$regions <- makeFactor(tissue = (zeroA_ssc@model[["class"]]==1),
                            bkgd = (zeroA_ssc@model[["class"]]==2))
image(zeroA, "regions")
cntrl_tissue <- subsetPixels(zeroA, regions == "tissue")
cntrl_bkgd <- subsetPixels(zeroA, regions == "bkgd")

# Background
cntrl_bkgd <- cntrl_bkgd |>
  summarizePixels(c(TIC="sum")) |>
  summarizeFeatures(c(Mean = "mean"))
image(cntrl_bkgd, "TIC", smooth="adaptive", enhance="histogram",
      col=matter::cpal("Cividis"))
plot(cntrl_bkgd, "Mean", annPeaks=5, xaxt="n")
axis(1, at = seq(100,1200, by=200), las=1)

# Tissue
cntrl_tissue <- cntrl_tissue |>
  summarizePixels(c(TIC= "sum"))|>
  summarizeFeatures(c(Mean = "mean"))
image(cntrl_tissue, "TIC", smooth="adaptive", enhance="histogram",
      col=matter::cpal("Cividis"))
plot(cntrl_tissue, "Mean", annPeaks=6, xaxt="n")
axis(1, at = seq(100,1200, by=200), las=1)

# Spatial Shrunken Centroids
cntrl_tissue_ssc <- spatialShrunkenCentroids(cntrl_tissue, r=1, 
                                             k = 7, s= 2^(0:5), 
                                             weights="gaussian")
cntrl_tissue_ssc
image(cntrl_tissue_ssc, i = 1:6)
cntrl_tissue_mdl4 <- cntrl_tissue_ssc[[4]]
image(cntrl_tissue_mdl4)
image(cntrl_tissue_mdl4,superpose=FALSE)
plot(cntrl_tissue_mdl4, type= "statistic", lwd=2, select = 1)
plot(cntrl_tissue_mdl4, type= "statistic", lwd=2, superpose=FALSE, 
     key=FALSE, xlim = c(100,1200), annPeaks= 2)
cntrl_mdl4_top <- topFeatures(cntrl_tissue_mdl4)

cntrl_tissue$anatomy <- makeFactor(
  Background = (cntrl_tissue_mdl4@model[["class"]] == 4),
  Daughter = (cntrl_tissue_mdl4@model[["class"]] == 5),
  innParent = (cntrl_tissue_mdl4@model[["class"]] == 2),
  outParent = (cntrl_tissue_mdl4@model[["class"]] == 3),
  edgeFrond = (cntrl_tissue_mdl4@model[["class"]] == 6),
  Background2 = cntrl_tissue_mdl4@model[["class"]] == 1)
image(cntrl_tissue, "anatomy")
cntrl_tissue$anatomy <- factor(cntrl_tissue$anatomy, 
                               levels = c("Background", "Daughter",
                                          "innParent","outParent",
                                          "edgeFrond","Background2"))

# Principle Component Analysis
cntrl_tissue_pca <- PCA(cntrl_tissue, ncomp=6)
cntrl_tissue_pca
image(cntrl_tissue_pca, superpose=FALSE, layout= c(2,3), smooth="adaptive", 
      enhance="hist")

plot(cntrl_tissue_pca, type="x", groups=cntrl_tissue$anatomy, shape=20)

plot(cntrl_tissue_pca@model[["x"]][,1],cntrl_tissue_pca@model[["x"]][,6],  
     col=cntrl_tissue$anatomy, pch=20, xlab= "PC1", ylab="PC6")
legend("topright", legend=unique(cntrl_tissue$anatomy), pch=16, 
       col=unique(cntrl_tissue$anatomy))

plot(cntrl_tissue_pca@model[["x"]][,1],cntrl_tissue_pca@model[["x"]][,5],  
     col=cntrl_tissue$anatomy, pch=20, xlab= "PC1", ylab="PC5")
legend("topright", legend=levels(cntrl_tissue$anatomy), pch=16, 
       col=unique(cntrl_tissue$anatomy))

plot(cntrl_tissue_pca@model[["x"]][,1],cntrl_tissue_pca@model[["x"]][,4],  
     col=cntrl_tissue$anatomy, pch=20, xlab= "PC1", ylab="PC4")
legend("topright", legend=unique(cntrl_tissue$anatomy), pch=16, 
       col=unique(cntrl_tissue$anatomy)) 

plot(cntrl_tissue_pca@model[["x"]][,1],cntrl_tissue_pca@model[["x"]][,3],  
     col=cntrl_tissue$anatomy, pch=20, xlab= "PC1", ylab="PC3")
legend("topright", legend=unique(cntrl_tissue$anatomy), pch=16, 
       col=unique(cntrl_tissue$anatomy)) 

plot(cntrl_tissue_pca@model[["x"]][,1],cntrl_tissue_pca@model[["x"]][,2],  
     col=cntrl_tissue$anatomy, pch=20, xlab= "PC1", ylab="PC2")
legend("topright", legend=unique(cntrl_tissue$anatomy), pch=16,
       col=unique(cntrl_tissue$anatomy))

plot(cntrl_tissue_pca@model[["x"]][,2],cntrl_tissue_pca@model[["x"]][,3],  
     col=cntrl_tissue$anatomy, pch=20, xlab= "PC2", ylab="PC3")
legend("topright", legend=unique(cntrl_tissue$anatomy), pch=16, 
       col=unique(cntrl_tissue$anatomy))

plot(cntrl_tissue_pca@model[["x"]][,2],cntrl_tissue_pca@model[["x"]][,4],  
     col=cntrl_tissue$anatomy, pch=20, xlab= "PC2", ylab="PC4")
legend("topright", legend=unique(cntrl_tissue$anatomy), pch=16, 
       col=unique(cntrl_tissue$anatomy)) 

plot(cntrl_tissue_pca@model[["x"]][,2],cntrl_tissue_pca@model[["x"]][,5],  
     col=cntrl_tissue$anatomy, pch=20, xlab= "PC2", ylab="PC5")
legend("topright", legend=levels(cntrl_tissue$anatomy), pch=16, 
       col=unique(cntrl_tissue$anatomy)) 

plot(cntrl_tissue_pca@model[["x"]][,2],cntrl_tissue_pca@model[["x"]][,6],  
     col=cntrl_tissue$anatomy, pch=20, xlab= "PC2", ylab="PC6")
legend("topright", legend=levels(cntrl_tissue$anatomy), pch=16, 
       col=unique(cntrl_tissue$anatomy)) 

plot(cntrl_tissue_pca@model[["x"]][,3],cntrl_tissue_pca@model[["x"]][,6],  
     col=cntrl_tissue$anatomy, pch=20, xlab= "PC3", ylab="PC6")
legend("topright", legend=levels(cntrl_tissue$anatomy), pch=16, 
       col=unique(cntrl_tissue$anatomy)) 

plot(cntrl_tissue_pca@model[["x"]][,3],cntrl_tissue_pca@model[["x"]][,5],  
     col=cntrl_tissue$anatomy, pch=20, xlab= "PC3", ylab="PC5")
legend("topright", legend=levels(cntrl_tissue$anatomy), pch=16, 
       col=unique(cntrl_tissue$anatomy))

plot(cntrl_tissue_pca@model[["x"]][,3],cntrl_tissue_pca@model[["x"]][,4],  
     col=cntrl_tissue$anatomy, pch=20, xlab= "PC3", ylab="PC4")
legend("topright", legend=levels(cntrl_tissue$anatomy), pch=16, 
       col=unique(cntrl_tissue$anatomy))

plot(cntrl_tissue_pca@model[["x"]][,4],cntrl_tissue_pca@model[["x"]][,5],  
     col=cntrl_tissue$anatomy, pch=20, xlab= "PC4", ylab="PC5")
legend("topright", legend=levels(cntrl_tissue$anatomy), pch=16, 
       col=unique(cntrl_tissue$anatomy))

plot(cntrl_tissue_pca@model[["x"]][,4],cntrl_tissue_pca@model[["x"]][,6],  
     col=cntrl_tissue$anatomy, pch=20, xlab= "PC4", ylab="PC6")
legend("topright", legend=levels(cntrl_tissue$anatomy), pch=16, 
       col=unique(cntrl_tissue$anatomy))

plot(cntrl_tissue_pca@model[["x"]][,5],cntrl_tissue_pca@model[["x"]][,6],  
     col=cntrl_tissue$anatomy, pch=20, xlab= "PC5", ylab="PC6")
legend("topright", legend=levels(cntrl_tissue$anatomy), pch=16, 
       col=unique(cntrl_tissue$anatomy))

# Exporting images in mass
# plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE) 
# plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
# file.copy(from=plots.png.paths, to="/work/LAS/yjlee-lab/buckm065/CardinalProject/Replicates/20240911_images")
# plots.png.details <- file.info(plots.png.paths)
# plots.png.details <- plots.png.details[order(plots.png.details$mtime),]
# sorted.png.names <- gsub(plots.dir.path, "/work/LAS/yjlee-lab/buckm065/CardinalProject/Replicates/20240911_images", row.names(plots.png.details), fixed=TRUE)
# numbered.png.names <- paste0("/work/LAS/yjlee-lab/buckm065/CardinalProject/Replicates/20240911_images", 1:length(sorted.png.names), ".png")
# file.rename(from=sorted.png.names, to=numbered.png.names)

plot(cntrl_tissue_pca, type="rotation", select=2, groups=cntrl_tissue$anatomy, lwd=2,
     annPeaks = 10, key = FALSE)
plot(cntrl_tissue_pca, type = "scree")

image(cntrl_tissue, mz=381.0798, enhance="hist")

# Nonnegative Matrix Factorization

cntrl_tissue_nmf <- NMF(cntrl_tissue, ncomp = 6, niter=30)
cntrl_tissue_nmf
image(cntrl_tissue_nmf, smooth="adaptive", enhance="hist")
image(cntrl_tissue_nmf, smooth="adaptive", enhance="hist", superpose=FALSE, 
      layout= c(2,3), key= FALSE)
plot(cntrl_tissue_nmf, lwd = 2, superpose= FALSE, layout=c(2,3))
plot(cntrl_tissue_nmf, type="x", groups=cntrl_tissue$anatomy)

# Projection to Latent Structures

cntrl_tissue_pls <- PLS(cntrl_tissue, y=cntrl_tissue$anatomy, ncomp= 6)
image(cntrl_tissue_pls, type="response", layout = c(2,3), scale = TRUE)
plot(cntrl_tissue_pls, type="coefficients", lwd = 2, annPeaks = "circle",
     superpose=FALSE)
plot(cntrl_tissue_pls, type="scores", groups = cntrl_tissue$anatomy, lwd=2)

# Spatially-aware Dirichlet Gaussian Mixture Model
zeroA <- zeroA |>
  summarizePixels(c(TIC= "sum"))|>
  summarizeFeatures(c(Mean = "mean"))
zeroA$tissue <- ifelse(zeroA$regions=="tissue",zeroA$TIC, 0)
zeroA$bkgd <- ifelse(zeroA$regions=="bkgd", zeroA$TIC, 0)
image(zeroA, c("tissue", "bkgd"),superpose=TRUE, free="xy", 
      col=dpal("Dark2"),enhance="histogram", scale=TRUE)
zeroA$samples <- interaction(run(zeroA), zeroA$regions)

cntrl_sdgmm <- spatialDGMM(zeroA, r = 10, k = 2, groups=zeroA$samples)
zeroA_means <- meansTest(zeroA, ~regions, samples=zeroA$samples)
zeroA_means
# Two Days 13C ------------------------------------------------------------

twoA <-readMSIData("2dayA_13C.imzML") #Read Unlabeled DW
twoA <- twoA |>
  normalize(method="tic") |> #TIC Normalization
  peakProcess(SNR=3,tolerance = 7, units = "ppm", 
              filterFreq = 0.01) #Peak Picking and Alignment

twoA_ssc <- spatialShrunkenCentroids(twoA, r = 1, s=0,k=2, 
                                      weights = "gaussian")
image(twoA_ssc)

twoA$regions <- makeFactor(tissue = (twoA_ssc@model[["class"]]==2),
                            bkgd = (twoA_ssc@model[["class"]]==1))
image(twoA, "regions")
twoA_tissue <- subsetPixels(twoA, regions == "tissue")
twoA_bkgd <- subsetPixels(twoA, regions == "bkgd")

# Background
twoA_bkgd <- twoA_bkgd |>
  summarizePixels(c(TIC="sum")) |>
  summarizeFeatures(c(Mean = "mean"))
image(twoA_bkgd, "TIC", smooth="adaptive", enhance="histogram",
      col=matter::cpal("Cividis"))
plot(twoA_bkgd, "Mean", annPeaks=5, xaxt="n")
axis(1, at = seq(800,1100, by=200), las=1)

# Tissue
twoA_tissue <- twoA_tissue |>
  summarizePixels(c(TIC= "sum"))|>
  summarizeFeatures(c(Mean = "mean"))
image(twoA_tissue, "TIC", smooth="adaptive", enhance="histogram",
      col=matter::cpal("Cividis"))
plot(twoA_tissue, "Mean", annPeaks=6, xaxt="n")
axis(1, at = seq(800,1100, by=200), las=1)

# Spatial Shrunken Centroids
twoA_tissue_ssc <- spatialShrunkenCentroids(twoA_tissue, r=1, 
                                             k = 7, s= 2^(0:5), 
                                             weights="gaussian")
twoA_tissue_ssc
image(twoA_tissue_ssc, i = 1:6)
twoA_tissue_mdl4 <- twoA_tissue_ssc[[4]]
image(twoA_tissue_mdl4, smooth="adaptive", enhance="hist",select=5)
image(twoA_tissue_mdl4,superpose=FALSE)
plot(twoA_tissue_mdl4, type= "statistic", lwd=2)
plot(twoA_tissue_mdl4, type= "statistic", lwd=2, superpose=FALSE,
     annPeaks = 5)
twoA_tissue_top <- topFeatures(twoA_tissue_mdl4)

twoA_tissue$anatomy <- makeFactor(
  Background = (twoA_tissue_mdl4@model[["class"]] == 4),
  Background2= (twoA_tissue_mdl4@model[["class"]] == 5),
  Daughter = (twoA_tissue_mdl4@model[["class"]] == 7),
  Parent = (twoA_tissue_mdl4@model[["class"]] == 2),
  Intermediate = (twoA_tissue_mdl4@model[["class"]] == 3),
  edgeFrond = (twoA_tissue_mdl4@model[["class"]] == 6),
  edgeFrond2 = (twoA_tissue_mdl4@model[["class"]]==1) )

image(twoA_tissue, "anatomy", smoot="adaptive", enhance="hist")


# Colocalization of Segment 5 Peaks
image(twoA, mz=806.89, enhance="hist", smooth = "adapt")
image(zeroA, mz=806.8946,, enhance="hist", smooth = "adapt")

twoA_coloc <- colocalized(twoA, mz=806.89)
twoA_coloc
image(twoA, mz=twoA_coloc$mz[1:5], smooth="adapt", enhance="hist", scale=TRUE)


segmentationTest(twoA_tissue, ~anatomy, samples=run(twoA_tissue))

# Principle Component Analysis
twoA_tissue_pca <- PCA(twoA_tissue, ncomp=6)
twoA_tissue_pca
image(twoA_tissue_pca, superpose=FALSE, layout= c(2,3), smooth="adaptive", 
      enhance="hist")
plot(twoA_tissue_pca, type="x", groups=twoA_tissue$anatomy, shape=20)
plot(twoA_tissue_pca, type="rotation", groups=twoA_tissue$anatomy, shape=20,
     superpose=FALSE, layout = c(2,3), scale = TRUE, annPeaks = 10)
plot(twoA_tissue_pca, type = "scree")

# Nonnegative Matrix Factorization

twoA_tissue_nmf <- NMF(twoA_tissue, ncomp = 6, niter=30)
twoA_tissue_nmf
image(twoA_tissue_nmf, smooth="adaptive", enhance="hist")
image(twoA_tissue_nmf, smooth="adaptive", enhance="hist", superpose=FALSE, 
      layout= c(2,3), key= FALSE)
plot(twoA_tissue_nmf, lwd = 2, annPeaks = 5, superpose= FALSE, layout=c(2,3))
plot(twoA_tissue_nmf, type="x", groups=twoA_tissue$anatomy, shape=20)

# Projection to Latent Structures

twoA_tissue_pls <- PLS(twoA_tissue, y=twoA_tissue$anatomy, ncomp= 6)
image(twoA_tissue_pls, type="response", scale = TRUE)
# image(twoA_tissue_pls, type="response", layout = c(2,3), superpose=FALSE)
plot(twoA_tissue_pls, type="coefficients", lwd = 2, annPeaks = 5,
     superpose=FALSE)
plot(twoA_tissue_pls, type="scores", groups = twoA_tissue$anatomy, shape=20)

# Spatially-aware Dirichlet Gaussian Mixture Model

twoA_sdgmm <- spatialDGMM(twoA, r = 1, k = 2)

# Three Days 13C ----------------------------------------------------------

threeB <-readMSIData("3dayB_13C.imzML") #Read Unlabeled DW
threeB <- threeB |>
  normalize(method="tic") |> #TIC Normalization
  peakProcess(SNR=3,tolerance = 7, units = "ppm", 
              filterFreq = 0.01) #Peak Picking and Alignment

threeB_ssc <- spatialShrunkenCentroids(threeB, r = 1, s=0,k=2, 
                                     weights = "gaussian")
image(threeB_ssc)

threeB$regions <- makeFactor(tissue = (threeB_ssc@model[["class"]]==2),
                           bkgd = (threeB_ssc@model[["class"]]==1))
image(threeB, "regions")
threeB_tissue <- subsetPixels(threeB, regions == "tissue")
threeB_bkgd <- subsetPixels(threeB, regions == "bkgd")

# Background
threeB_bkgd <- threeB_bkgd |>
  summarizePixels(c(TIC="sum")) |>
  summarizeFeatures(c(Mean = "mean"))
image(threeB_bkgd, "TIC", smooth="adaptive", enhance="histogram",
      col=matter::cpal("Cividis"))
plot(threeB_bkgd, "Mean", annPeaks=5, xaxt="n")
axis(1, at = seq(800,1100, by=200), las=1)

# Tissue
threeB_tissue <- threeB_tissue |>
  summarizePixels(c(TIC= "sum"))|>
  summarizeFeatures(c(Mean = "mean"))
image(threeB_tissue, "TIC", smooth="adaptive", enhance="histogram",
      col=matter::cpal("Cividis"))
plot(threeB_tissue, "Mean", annPeaks=6, xaxt="n")
axis(1, at = seq(800,1100, by=200), las=1)

# Spatial Shrunken Centroids
threeB_tissue_ssc <- spatialShrunkenCentroids(threeB_tissue, r=1, 
                                            k = 7, s= 2^(0:5), 
                                            weights="gaussian")
threeB_tissue_ssc
image(threeB_tissue_ssc, i = 1:6)
threeB_tissue_mdl4 <- threeB_tissue_ssc[[4]]
image(threeB_tissue_mdl4)
image(threeB_tissue_mdl4,superpose=FALSE)
plot(threeB_tissue_mdl4, type= "statistic", lwd=2)
plot(threeB_tissue_mdl4, type= "statistic", lwd=2, superpose=FALSE,
     annPeaks = 5)
threeB_tissue_top <- topFeatures(threeB_tissue_mdl4)

threeB_tissue$anatomy <- makeFactor(
  Background = threeB_tissue_mdl4@model[["class"]] == 1,
  Background2= threeB_tissue_mdl4@model[["class"]] == 6,
  Daughter = threeB_tissue_mdl4@model[["class"]] == 4,
  Parent = threeB_tissue_mdl4@model[["class"]] == 2,
  Intermediate = threeB_tissue_mdl4@model[["class"]] == 3,
  edgeFrond = threeB_tissue_mdl4@model[["class"]] == 5 )
image(threeB_tissue, "anatomy")

# Principle Component Analysis
threeB_tissue_pca <- PCA(threeB_tissue, ncomp=6)
threeB_tissue_pca
image(threeB_tissue_pca, superpose=FALSE, layout= c(2,3), smooth="adaptive", 
      enhance="hist")
plot(threeB_tissue_pca, type="x", groups=threeB_tissue$anatomy, shape=20)
plot(threeB_tissue_pca, type="rotation", groups=threeB_tissue$anatomy, shape=20,
     superpose=FALSE, layout = c(2,3), scale = TRUE, annPeaks = 10, key=FALSE)
plot(threeB_tissue_pca, type = "scree")

# Nonnegative Matrix Factorization

threeB_tissue_nmf <- NMF(threeB_tissue, ncomp = 6, niter=30)
threeB_tissue_nmf
image(threeB_tissue_nmf, smooth="adaptive", enhance="hist")
image(threeB_tissue_nmf, smooth="adaptive", enhance="hist", superpose=FALSE, 
      layout= c(2,3), key= FALSE)
plot(threeB_tissue_nmf, lwd = 2, annPeaks = 5, superpose= FALSE, layout=c(2,3))
plot(threeB_tissue_nmf, type="x", groups=threeB_tissue$anatomy, shape = 20)

# Projection to Latent Structures

threeB_tissue_pls <- PLS(threeB_tissue, y=threeB_tissue$anatomy, ncomp= 6)
image(threeB_tissue_pls, type="response", layout = c(2,3), scale = TRUE)
plot(threeB_tissue_pls, type="coefficients", lwd = 2, annPeaks = "circle",
     superpose=FALSE)
plot(threeB_tissue_pls, type="scores", groups = threeB_tissue$anatomy, lwd=2)

# Spatially-aware Dirichlet Gaussian Mixture Model

threeB_sdgmm <- spatialDGMM(threeB, r = 1, k = 2)



# Five Days D -------------------------------------------------------------

fiveB <-readMSIData("5dayB_2H.imzML") #Read Unlabeled DW
fiveB <- fiveB |>
  normalize(method="tic") |> #TIC Normalization
  peakProcess(SNR=2,tolerance = 7, units = "ppm", 
              filterFreq = 25) #Peak Picking and Alignment

fiveB_ssc <- spatialShrunkenCentroids(fiveB, r = 1, s=0,k=2, 
                                       weights = "gaussian")
image(fiveB_ssc)

fiveB$regions <- makeFactor(tissue = (fiveB_ssc@model[["class"]]==1),
                             bkgd = (fiveB_ssc@model[["class"]]==2))
image(fiveB, "regions")
fiveB_tissue <- subsetPixels(fiveB, regions == "tissue")
fiveB_bkgd <- subsetPixels(fiveB, regions == "bkgd")

# Background
fiveB_bkgd <- fiveB_bkgd |>
  summarizePixels(c(TIC="sum")) |>
  summarizeFeatures(c(Mean = "mean"))
image(fiveB_bkgd, "TIC", smooth="adaptive", enhance="histogram",
      col=matter::cpal("Cividis"))
plot(fiveB_bkgd, "Mean", annPeaks=5, xaxt="n")
axis(1, at = seq(100,1200, by=200), las=1)

# Tissue
fiveB_tissue <- fiveB_tissue |>
  summarizePixels(c(TIC= "sum"))|>
  summarizeFeatures(c(Mean = "mean"))
image(fiveB_tissue, "TIC", smooth="adaptive", enhance="histogram",
      col=matter::cpal("Cividis"))
plot(fiveB_tissue, "Mean", annPeaks=6, xaxt="n")
axis(1, at = seq(100,1200, by=200), las=1)

# Spatial Shrunken Centroids
fiveB_tissue_ssc <- spatialShrunkenCentroids(fiveB_tissue, r=1, 
                                              k = 7, s= 2^(0:5), 
                                              weights="gaussian")
fiveB_tissue_ssc
image(fiveB_tissue_ssc, i = 1:6)
fiveB_tissue_mdl4 <- fiveB_tissue_ssc[[3]]
image(fiveB_tissue_mdl4)
image(fiveB_tissue_mdl4,superpose=FALSE)
plot(fiveB_tissue_mdl4, type= "statistic", lwd=2)
plot(fiveB_tissue_mdl4, type= "statistic", lwd=2, superpose=TRUE, xlim=c(322,400), key = TRUE)
plot(fiveB_tissue_mdl4, type= "centers", lwd=2)

fiveB_tissue_top <- topFeatures(fiveB_tissue_mdl4)


fiveB_tissue$anatomy <- makeFactor(
  Background = fiveB_tissue_mdl4@model[["class"]] == 2,
  Background2= fiveB_tissue_mdl4@model[["class"]] == 4,
  Daughter = fiveB_tissue_mdl4@model[["class"]] == 7,
  Parent = fiveB_tissue_mdl4@model[["class"]] == 1,
  Intermediate = fiveB_tissue_mdl4@model[["class"]] == 6,
  edgeFrond= fiveB_tissue_mdl4@model[["class"]] == 5,
  DaughterEdge = fiveB_tissue_mdl4@model[["class"]] == 3 )
image(fiveB_tissue, "anatomy")

# Principle Component Analysis
fiveB_tissue_pca <- PCA(fiveB_tissue, ncomp=6)
fiveB_tissue_pca
image(fiveB_tissue_pca, superpose=FALSE, layout= c(2,3), smooth="adaptive", 
      enhance="hist", key = FALSE)
plot(fiveB_tissue_pca, type = "x", groups = fiveB_tissue$anatomy, shape = 20)
plot(fiveB_tissue_pca@model[["x"]][,1],fiveB_tissue_pca@model[["x"]][,6],  
     col=fiveB_tissue$anatomy, pch=20)
plot(fiveB_tissue_pca, type="rotation", groups=fiveB_tissue$anatomy, lwd=2,
     superpose=FALSE, layout = c(2,3), scale = TRUE, annPeaks = 10, key = FALSE)
plot(fiveB_tissue_pca, type = "scree")
# cls = c("#4F7AA8","#76B7B2","#EDC947","#F28E27","#59A14E","#B07AA1","#E05759")

plot(fiveB_tissue_pca@model[["x"]][,1],fiveB_tissue_pca@model[["x"]][,2],  
     col=fiveB_tissue$anatomy,pch=20, xlab= "PC1", ylab="PC2")
legend("topright", legend=levels(fiveB_tissue$anatomy), pch=16, 
       col=unique(fiveB_tissue$anatomy))

plot(fiveB_tissue_pca@model[["x"]][,1],fiveB_tissue_pca@model[["x"]][,3],  
     col=fiveB_tissue$anatomy,pch=20, xlab= "PC1", ylab="PC3")
legend("topright", legend=levels(fiveB_tissue$anatomy), pch=16, 
       col=unique(fiveB_tissue$anatomy))

plot(fiveB_tissue_pca@model[["x"]][,1],fiveB_tissue_pca@model[["x"]][,4],  
     col=fiveB_tissue$anatomy,pch=20, xlab= "PC1", ylab="PC4")
legend("topright", legend=levels(fiveB_tissue$anatomy), pch=16, 
       col=unique(fiveB_tissue$anatomy))

plot(fiveB_tissue_pca@model[["x"]][,1],fiveB_tissue_pca@model[["x"]][,5],  
     col=fiveB_tissue$anatomy,pch=20, xlab= "PC1", ylab="PC5")
legend("topright", legend=levels(fiveB_tissue$anatomy), pch=16, 
       col=unique(fiveB_tissue$anatomy))

plot(fiveB_tissue_pca@model[["x"]][,2],fiveB_tissue_pca@model[["x"]][,6],  
     col=fiveB_tissue$anatomy,pch=20, xlab= "PC2", ylab="PC6")
legend("topright", legend=levels(fiveB_tissue$anatomy), pch=16, 
       col=unique(fiveB_tissue$anatomy))

plot(fiveB_tissue_pca@model[["x"]][,2],fiveB_tissue_pca@model[["x"]][,5],  
     col=fiveB_tissue$anatomy,pch=20, xlab= "PC2", ylab="PC5")
legend("topright", legend=levels(fiveB_tissue$anatomy), pch=16, 
       col=unique(fiveB_tissue$anatomy))

plot(fiveB_tissue_pca@model[["x"]][,2],fiveB_tissue_pca@model[["x"]][,4],  
     col=fiveB_tissue$anatomy,pch=20, xlab= "PC2", ylab="PC4")
legend("topright", legend=levels(fiveB_tissue$anatomy), pch=16, 
       col=unique(fiveB_tissue$anatomy))

plot(fiveB_tissue_pca@model[["x"]][,2],fiveB_tissue_pca@model[["x"]][,3],  
     col=fiveB_tissue$anatomy,pch=20, xlab= "PC2", ylab="PC3")
legend("topright", legend=levels(fiveB_tissue$anatomy), pch=16, 
       col=unique(fiveB_tissue$anatomy))

plot(fiveB_tissue_pca@model[["x"]][,3],fiveB_tissue_pca@model[["x"]][,4],  
     col=fiveB_tissue$anatomy,pch=20, xlab= "PC3", ylab="PC4")
legend("topright", legend=levels(fiveB_tissue$anatomy), pch=16, 
       col=unique(fiveB_tissue$anatomy))

plot(fiveB_tissue_pca@model[["x"]][,3],fiveB_tissue_pca@model[["x"]][,5],  
     col=fiveB_tissue$anatomy,pch=20, xlab= "PC3", ylab="PC5")
legend("topright", legend=levels(fiveB_tissue$anatomy), pch=16, 
       col=unique(fiveB_tissue$anatomy))

plot(fiveB_tissue_pca@model[["x"]][,3],fiveB_tissue_pca@model[["x"]][,6],  
     col=fiveB_tissue$anatomy,pch=20, xlab= "PC3", ylab="PC6")
legend("topright", legend=levels(fiveB_tissue$anatomy), pch=16, 
       col=unique(fiveB_tissue$anatomy))

plot(fiveB_tissue_pca@model[["x"]][,4],fiveB_tissue_pca@model[["x"]][,6],  
     col=fiveB_tissue$anatomy,pch=20, xlab= "PC4", ylab="PC6")
legend("topright", legend=levels(fiveB_tissue$anatomy), pch=16, 
       col=unique(fiveB_tissue$anatomy))

plot(fiveB_tissue_pca@model[["x"]][,4],fiveB_tissue_pca@model[["x"]][,5],  
     col=fiveB_tissue$anatomy,pch=20, xlab= "PC4", ylab="PC5")
legend("topright", legend=levels(fiveB_tissue$anatomy), pch=16, 
       col=unique(fiveB_tissue$anatomy))

plot(fiveB_tissue_pca@model[["x"]][,5],fiveB_tissue_pca@model[["x"]][,6],  
     col=fiveB_tissue$anatomy,pch=20, xlab= "PC5", ylab="PC6")
legend("topright", legend=levels(fiveB_tissue$anatomy), pch=16, 
       col=unique(fiveB_tissue$anatomy))

# Nonnegative Matrix Factorization

fiveB_tissue_nmf <- NMF(fiveB_tissue, ncomp = 6, niter=30)
fiveB_tissue_nmf
image(fiveB_tissue_nmf, smooth="adaptive", enhance="hist")
image(fiveB_tissue_nmf, smooth="adaptive", enhance="hist", superpose=FALSE, 
      layout= c(2,3), key= FALSE)
plot(fiveB_tissue_nmf, lwd = 2, annPeaks = 5, superpose= FALSE, 
     layout=c(2,3))
plot(fiveB_tissue_nmf@model[["x"]][,3],fiveB_tissue_nmf@model[["x"]][,2],  col=fiveB_tissue$anatomy, pch=20, 
     log ="xy", xlim = c(1E-06, 1E+01),  ylim=c(1E-06, 1E+01))

# Projection to Latent Structures

fiveB_tissue_pls <- PLS(fiveB_tissue, y=fiveB_tissue$anatomy, ncomp= 6)
image(fiveB_tissue_pls, type="response", layout = c(2,3), scale = TRUE)
plot(fiveB_tissue_pls, type="coefficients", lwd = 2, annPeaks = "circle",
     superpose=FALSE)
plot(fiveB_tissue_pls, type="scores", groups = fiveB_tissue$anatomy, shape=20)


# Spatially-aware Dirichlet Gaussian Mixture Model

fiveB_sdgmm <- spatialDGMM(fiveB, r = 1, k = 2)

# ROI
fiveB2 <-readMSIData("5dayB_2H.imzML") #Read Unlabeled DW
fiveB2 <- fiveB2 |>
  normalize(method="tic") |> #TIC Normalization
  peakProcess(SNR=2,tolerance = 10, units = "ppm", 
              filterFreq = 25) #Peak Picking and Alignment
fiveB2_ROI <- selectROI(fiveB2, mz = 104.1075,mode= "region")

fiveB2$regions <- makeFactor(tissue = (fiveB2_ROI==TRUE),
                     background = (fiveB2_ROI ==FALSE))
image(fiveB2, "regions")

fiveB2_tissue <- subsetPixels(fiveB2, regions == "tissue")
fiveB2_bkgd <- subsetPixels(fiveB2, regions == "background")

fiveB2_tissue_ssc <- spatialShrunkenCentroids(fiveB2_tissue, r=1, 
                                             k = 7, s= 2^(0:5), 
                                             weights="gaussian")
fiveB2_tissue_ssc
image(fiveB2_tissue_ssc, i = 1:6)
fiveB2_tissue_mdl3 <- fiveB2_tissue_ssc[[3]]
image(fiveB2_tissue_mdl3)
image(fiveB2_tissue_mdl3,superpose=FALSE)
plot(fiveB2_tissue_mdl3, type= "statistic", lwd=2)
plot(fiveB2_tissue_mdl3, type= "statistic", lwd=2, superpose=FALSE,
     annPeaks = 5, key = FALSE)
plot(fiveB2_tissue_mdl3, type= "statistic", xlim = c(322,400), lwd=2)
fiveB2_tissue_top <- topFeatures(fiveB2_tissue_mdl3)

# Background
fiveB2_bkgd <- fiveB2_bkgd |>
  summarizePixels(c(TIC="sum")) |>
  summarizeFeatures(c(Mean = "mean"))
image(fiveB2_bkgd, "TIC", smooth="adaptive", enhance="histogram",
      col=matter::cpal("Cividis"))
plot(fiveB2_bkgd, "Mean", annPeaks=5, xaxt="n")
axis(1, at = seq(100,1200, by=200), las=1)

# Tissue
fiveB2_tissue <- fiveB_tissue |>
  summarizePixels(c(TIC= "sum"))|>
  summarizeFeatures(c(Mean = "mean"))
image(fiveB2_tissue, "TIC", smooth="adaptive", enhance="histogram",
      col=matter::cpal("Cividis"))
plot(fiveB2_tissue, "Mean", annPeaks=6, xaxt="n")
axis(1, at = seq(100,1200, by=200), las=1)

# setup.layout(c(1,1))
fiveB_mean <- summarizeFeatures(fiveB_tissue, "mean")
threeB_mean <- summarizeFeatures(threeB_tissue, "mean")
twoA_mean <- summarizeFeatures(twoA_tissue, "mean")
fiveB_bkgd_mean <- summarizeFeatures(fiveB_bkgd, "mean")
plot(subsetFeatures(fiveB_mean, mz>=810, mz <= 860),"mean", annPeaks = 7)
plot(subsetFeatures(threeB_mean, mz>=810, mz <= 860),"mean", annPeaks = 20)
plot(subsetFeatures(twoA_mean, mz>=810, mz <= 860),"mean", annPeaks = 20)
plot(subsetFeatures(fiveB_bkgd_mean, mz>=322, mz <= 400),"mean", annPeaks = 7)
image(fiveB_tissue, mz=c(813.4914, 849.7225),
      main = "Labeling Localizations",
      enhance = "histogram",
      smooth = "adaptive",
      normalize = "linear",
      superpose=TRUE)

image(threeB_tissue, mz=c(813.4919,
                          821.519,
                          822.522,
                          823.525,
                          825.532,
                          835.566,
                          836.569,
                          837.572,
                          838.576,
                          839.579,
                          845.599,
                          846.603,
                          847.606,
                          848.609,
                          856.636,
                          857.640,
                          858.643),
      enhance = "histogram",
      smooth = "adaptive",
      normalize = "linear",
      superpose=TRUE)
image(twoA_tissue, mz=c(813.4919,
                        821.519,
                        822.522,
                        823.525,
                        825.532,
                        835.566,
                        836.569,
                        837.572,
                        838.576,
                        839.579,
                        845.599,
                        846.603,
                        847.606,
                        848.609,
                        856.636,
                        857.640),
      enhance = "histogram",
      smooth = "adaptive",
      normalize = "linear",
      superpose=TRUE)

image(fiveB_tissue, mz=c(104.1075,
                         106.120,
                         107.126,
                         108.1326,
                         109.139,
                         110.145,
                         111.151,
                         112.158),
      enhance = "histogram",
      smooth = "adaptive",
      normalize = "linear",
      superpose=TRUE)

image(fiveB_tissue, mz=c(813.4919,
                         818.5235,
                         834.6264,
                         835.6335,
                         836.6397,
                         843.6812,
                         844.6884,
                         851.736,
                         852.7424,
                         853.7492),
      enhance = "histogram",
      smooth = "adaptive",
      normalize = "linear",
      superpose=TRUE)
# Export ------------------------------------------------------------------

write.csv(fiveB_tissue_top, "/work/LAS/yjlee-lab/buckm065/CardinalProject/Replicates/fiveB_top.csv")
write.csv(fiveB2_tissue_top, "/work/LAS/yjlee-lab/buckm065/CardinalProject/Replicates/fiveB2_top.csv")
write.csv(threeB_tissue_top, "/work/LAS/yjlee-lab/buckm065/CardinalProject/Replicates/threeB_top.csv")
write.csv(twoA_tissue_top, "/work/LAS/yjlee-lab/buckm065/CardinalProject/Replicates/twoA_top.csv")
write.csv(cntrl_mdl4_top, "/work/LAS/yjlee-lab/buckm065/CardinalProject/Replicates/cntrl_top.csv")


cntrl_bkgd_mean <- cbind(cntrl_bkgd@featureData@listData[["mz"]],cntrl_bkgd@featureData@listData[["Mean"]])
write.csv(cntrl_bkgd_mean, "/work/LAS/yjlee-lab/buckm065/CardinalProject/Replicates/cntrl_bkgd_mean.csv")

cntrl_tissue_mean <- cbind(cntrl_tissue@featureData@listData[["mz"]],cntrl_tissue@featureData@listData[["Mean"]])
write.csv(cntrl_tissue_mean, "/work/LAS/yjlee-lab/buckm065/CardinalProject/Replicates/cntrl_tissue_mean.csv")

twoA_tissue_mean <- cbind(twoA_tissue@featureData@listData[["mz"]],twoA_tissue@featureData@listData[["Mean"]])
write.csv(twoA_tissue_mean, "/work/LAS/yjlee-lab/buckm065/CardinalProject/Replicates/twoA_tissue_mean.csv")

twoA_bkgd_mean <- cbind(twoA_bkgd@featureData@listData[["mz"]],twoA_bkgd@featureData@listData[["Mean"]])
write.csv(twoA_bkgd_mean, "/work/LAS/yjlee-lab/buckm065/CardinalProject/Replicates/twoA_bkgd_mean.csv")

threeB_tissue_mean <- cbind(threeB_tissue@featureData@listData[["mz"]],threeB_tissue@featureData@listData[["Mean"]])
write.csv(threeB_tissue_mean, "/work/LAS/yjlee-lab/buckm065/CardinalProject/Replicates/threeB_tissue_mean.csv")

threeB_bkgd_mean <- cbind(threeB_bkgd@featureData@listData[["mz"]],threeB_bkgd@featureData@listData[["Mean"]])
write.csv(threeB_bkgd_mean, "/work/LAS/yjlee-lab/buckm065/CardinalProject/Replicates/threeB_bkgd_mean.csv")

fiveB_bkgd_mean <- cbind(fiveB_bkgd@featureData@listData[["mz"]],fiveB_bkgd@featureData@listData[["Mean"]])
write.csv(fiveB_bkgd_mean, "/work/LAS/yjlee-lab/buckm065/CardinalProject/Replicates/fiveB_bkgd_mean.csv")

fiveB_tissue_mean <- cbind(fiveB_tissue@featureData@listData[["mz"]],fiveB_tissue@featureData@listData[["Mean"]])
write.csv(fiveB_tissue_mean, "/work/LAS/yjlee-lab/buckm065/CardinalProject/Replicates/fiveB_tissue_mean.csv")

fiveB2_tissue_mean <- cbind(fiveB2_tissue@featureData@listData[["mz"]],fiveB2_tissue@featureData@listData[["Mean"]])
write.csv(fiveB2_tissue_mean, "/work/LAS/yjlee-lab/buckm065/CardinalProject/Replicates/fiveB2_tissue_mean.csv")

fiveB2_bkgd_mean <- cbind(fiveB2_bkgd@featureData@listData[["mz"]],fiveB2_bkgd@featureData@listData[["Mean"]])
write.csv(fiveB2_bkgd_mean, "/work/LAS/yjlee-lab/buckm065/CardinalProject/Replicates/fiveB2_bkgd_mean.csv")


# LABELING

image(fiveB_tissue, mz= c(208.97278,209.97909,210.98524,211.99149,212.9834,213.9896),
      enhance = "histogram",
      superpose=TRUE,
      smooth = "adaptive")
plot(fiveB_tissue, "Mean", xlim=c(205,215), ylim=c(0,10), annPeaks=60)
image(zeroA, mz= 208.97289, enhance="histogram")
image(zeroA, mz= 813.4919, enhance="histogram")

image(fiveB_tissue, mz= 381, enhance="histogram")
image(fiveB_tissue, mz=c(382.0012,
                         383.0074,
                         384.0989,
                         385.1051,
                         386.1111,
                         387.1177,
                         388.1237,
                         389.1303,
                         390.1367,
                         391.1429),
      enhance = "histogram",
      smooth = "adaptive",
      normalize = "linear",
      superpose=TRUE)
image(fiveB_tissue, mz=c(208.97278,
                         209.97909,
                         210.98524,
                         211.99149,
                         212.9834,
                         213.9896),
      enhance = "histogram",
      smooth = "adaptive",
      normalize = "linear",
      superpose=TRUE)
image(fiveB_tissue, mz=c(313.0323,
                         314.0357,
                         315.0479,
                         316.0507),
      enhance = "histogram",
      smooth = "adaptive",
      normalize = "linear",
      superpose=TRUE)
# 206.82968?
plot(subsetFeatures(fiveB_tissue, mz>=620, mz <= 680),"mean", annPeaks = 5)
# tSNE --------------------------------------------------------------------
fiveB_tsne_5 <- Rtsne(as.matrix(fiveB_pData), perplexity = 5)
plot(fiveB_tsne_5$Y, col = fiveB_tissue$anatomy, pch = 20)
legend("topright", legend=unique(fiveB_tissue$anatomy), pch=16, 
       col=unique(fiveB_tissue$anatomy))

fiveB_tsne_30 <- Rtsne(as.matrix(fiveB_pData), perplexity = 30)
plot(fiveB_tsne_30$Y, col = fiveB_tissue$anatomy, pch = 20)
legend("topright", legend=unique(fiveB_tissue$anatomy), pch=16, 
       col=unique(fiveB_tissue$anatomy))

fiveB_tsne_30_500 <- Rtsne(as.matrix(fiveB_pData), perplexity = 30,
                         step = 500)
plot(fiveB_tsne_30_500$Y, col = fiveB_tissue$anatomy, pch = 20)
legend("topright", legend=unique(fiveB_tissue$anatomy), pch=16, 
       col=unique(fiveB_tissue$anatomy))

fiveB_tsne_50_200 <- Rtsne(as.matrix(fiveB_pData), perplexity = 50,
                           eta = 1000, theta = 0.25, max_iter = 5000)
plot(fiveB_tsne_50_200$Y, col = fiveB_tissue$anatomy, pch = 20)
legend("topright", legend=unique(fiveB_tissue$anatomy), pch=16, 
       col=unique(fiveB_tissue$anatomy))

fiveB_tsne_100 <- Rtsne(as.matrix(fiveB_pData), perplexity = 100)
plot(fiveB_tsne_100$Y, col = fiveB_tissue$anatomy, pch = 20)
legend("topright", legend=unique(fiveB_tissue$anatomy), pch=16, 
       col=unique(fiveB_tissue$anatomy))


# fiveB_tsne_ssc  <- Rtsne(t(as.matrix(spectra(fiveB_tissue))), 
#                          check_duplicates = FALSE, perplexity = 350, dims = 3, 
#                          max_iter = 3000, step = 50,  eta = 200)
# plot(fiveB_tsne_ssc$Y[,1],fiveB_tsne_ssc$Y[,2], col = fiveB_tissue$anatomy, pch = 20)
# legend("topright", legend=levels(fiveB_tissue$anatomy), pch=16, 
#        col=unique(fiveB_tissue$anatomy))


# Means Test --------------------------------------------------------------
# fiveB_sdgmm <- spatialDGMM(fiveB_tissue, r = 50, k = 5, weights = "adaptive")
# fiveB_means <- meansTest(fiveB_sdgmm,~class)

# sessionInfo -------------------------------------------------------------

sessionInfo()
