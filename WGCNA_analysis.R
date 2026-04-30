rm(list=ls())


# Load required libraries

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(gridExtra)
allowWGCNAThreads() 

# Gene expression matrix and format to WGCNA program -LC data (Diet restricted)

#setwd('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/GitHub_folder/DR-PeripheralClock-Longevity/Data')

# Load your batch corrected TPM values

TPM<-read.csv(".../Data/Batch_corrected_TPM.csv", row.names = 1)

# TPM<-read.csv("Batch_corrected_TPM.csv", row.names = 1)


# subset DR data only
DR_only=TPM[,25:48]
norm.counts<-t(DR_only)

# 1. Network Construction  9Weighted networks : aij=(sij)^beta, beta is the soft power. we need to select the best beta

# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function fo select best beta
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

# visualization to pick power

sft.data <- sft$fitIndices

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)

# Set soft threshold as 18 as that retains the highest mean connectivity while reaching an R2 value above 0.80
soft_power <- 18

# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)


# use following steps to not get erros
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 8000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.1,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3,
                          deepSplit = 4, # sensitivity control (0-4))
                          minModuleSize = 20,
                          TOMDenom = 'mean',
                          reassignThreshold = 0.05)


# Convert labels to colors for plotting
mergedColors = labels2colors(bwnet$colors)

plotDendroAndColors(
  bwnet$dendrograms[[1]],
  mergedColors[bwnet$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )


# assign back original correlation
cor <- temp_cor



# Number of genes in each module
gnes_per_modules<-as.data.frame(table(bwnet$colors))
#write.csv(gnes_per_modules, 'gnes_per_modules_LC25_Batch_corrected_RUGv_K_2.csv')

D1<-as.data.frame(bwnet$colors)

#write.csv(D1, "WGCNA_module_LC25_batch_corrected_RUGv_k_2.csv")

