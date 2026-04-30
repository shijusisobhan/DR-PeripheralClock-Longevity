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


#------------------------------------------------------------------------------------------------
  # MDC analysis
#----------------------------------------------------------------------------------------------

HC25_expr<-TPM[,c(1:24)]
LC25_expr<-TPM[,c(25:48)]


# -------------------------------
# Step 2: Compute adjacency matrices
# -------------------------------
# Use same soft-threshold power as LC25 WGCNA
softPower <- 18  # replace with your chosen value

adj_LC25 <- adjacency(t(LC25_expr), power = softPower, type = "signed")
adj_HC25 <- adjacency(t(HC25_expr), power = softPower, type = "signed")

# -------------------------------
# Step 3: Compute intramodular connectivity for each gene
# -------------------------------


moduleColors<-D1$`bwnet$colors`
names(moduleColors)<-rownames(D1)

k_LC25 <- intramodularConnectivity(adj_LC25, moduleColors)
k_HC25 <- intramodularConnectivity(adj_HC25, moduleColors)

# Extract "kWithin" (connectivity within module)
k_LC25_within <- k_LC25$kWithin
k_HC25_within <- k_HC25$kWithin
names(k_LC25_within) <- rownames(k_LC25)
names(k_HC25_within) <- rownames(k_HC25)



# -------------------------------
# Step 5: Compute module-level MDC (median kDiff)
# -------------------------------
moduleList <- unique(moduleColors)

moduleMDC <- sapply(moduleList, function(mod){
  # Get genes in this module
  genes <- names(moduleColors)[moduleColors == mod]
  
  # Subset connectivity vectors using exact gene names
  mean(k_LC25_within[genes]) / mean(k_HC25_within[genes])
})

# moduleMedC now has your MedC per module
# Interpretation:
#  MedC > 0 : module stronger in LC25
#  MedC < 0 : module stronger in HC25


# -------------------------------
# Step 4: Permutation for empirical null
# -------------------------------
nPerm <- 1000
set.seed(123)

permMDC <- matrix(NA, nrow = nPerm, ncol = length(moduleList))
colnames(permMDC) <- moduleList

all_expr <- cbind(LC25_expr, HC25_expr)  # combine all samples

for (i in 1:nPerm){
  # Shuffle sample labels
  perm_labels <- sample(ncol(all_expr))
  perm_LC <- all_expr[, perm_labels[1:24]]
  perm_HC <- all_expr[, perm_labels[25:48]]
  
  # Compute adjacency
  perm_adj_LC <- adjacency(t(perm_LC), power = softPower, type = "signed")
  perm_adj_HC <- adjacency(t(perm_HC), power = softPower, type = "signed")
  
  # Intramodular connectivity
  perm_k_LC <- intramodularConnectivity(perm_adj_LC, colors = moduleColors)$kWithin
  perm_k_HC <- intramodularConnectivity(perm_adj_HC, colors = moduleColors)$kWithin
  names(perm_k_LC) <- rownames(perm_adj_LC)
  names(perm_k_HC) <- rownames(perm_adj_HC)
  
  # Compute MDC per module (ratio)
  for (mod in moduleList){
    genes <- names(moduleColors)[moduleColors == mod]
    permMDC[i, mod] <- mean(perm_k_LC[genes]) / mean(perm_k_HC[genes])
  }
}

# -------------------------------
# Step 5: Compute FDR (q-values)
# -------------------------------
empiricalP <- sapply(moduleList, function(mod){
  if (moduleMDC[mod] > 1){
    # Gain of connectivity
    mean(permMDC[,mod] >= moduleMDC[mod])
  } else {
    # Loss of connectivity
    mean(permMDC[,mod] <= moduleMDC[mod])
  }
})

# FDR correction
qValues <- p.adjust(empiricalP, method = "fdr")

# -------------------------------
# Step 6: Compile results
# -------------------------------
MDC_results <- data.frame(
  Module = moduleList,
  MDC = moduleMDC,
  pValue = empiricalP,
  qValue = qValues
)

# Interpretation:
# MDC > 1 & q < 0.05 -> module gains connectivity in LC25
# MDC < 1 & q < 0.05 -> module loses connectivity in LC25 (stronger in HC25)
#write.csv(MDC_results, 'MDC_result_BatchCorrected_TPM_K2.csv')

