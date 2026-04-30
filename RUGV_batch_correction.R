
rm(list=ls())
# Libraries
library(EDASeq)
library(RUVSeq)
library(DESeq2)

setwd('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/Drosophila/Fat body data/FB_old/FB_old_corrected/FB_old_paper/Kallisto_TPM')



df<-read.csv("TPM_FB_HC25_LC25.csv")



# Assume your data frame is called df
# First column is gene names, the remaining 24 columns are samples

#-------------------------------------------------------------------------------
# Subset only the numeric sample columns
# Filter rows
threshold <- 1
est_HC <- df[, 1:25]

#Keep the column which have count >5 in more than 60% of sample (exclude cont <5 in more than 40% of sample)
HC_filter <- est_HC[rowSums(est_HC < threshold) <= 0.4 * ncol(est_HC), ]


# Subset only the numeric sample columns
est_LC <- df[, c(1,26:49)]

LC_filter <- est_LC[rowSums(est_LC < threshold) <= 0.4 * ncol(est_LC), ]


#--------------------------------------------------------------------------------------------

# find the union of genes
#Common_gene<-intersect(HC_filter$ext_gene, LC_filter$ext_gene) # Intial filtering is over

Common_gene<-union(HC_filter$ext_gene, LC_filter$ext_gene) # Intial filtering is over
#--------------------------------------------------------------------------------------------

# Filtering and exploratory data analysis

df<-read.csv("Estimated_count_FB_HC25_LC25.csv")
counts<-df[which(df$ext_gene %in% Common_gene),]
counts<-counts[!duplicated(counts$ext_gene), ]
rownames(counts)<-counts[,1]
counts<-counts[,-1]
  
counts <- round(as.matrix(counts))

#--------------------------------------------------
# 2. Define sample info
#--------------------------------------------------
condition <- factor(rep(c("Control", "DR"), each=24))
batch     <- rep(rep(c("Batch1", "Batch2"), each=12), 2)
timepoints <- rep(seq(2, 24, by=2), times=4)   # 12 time points repeated 4 times

pheno <- data.frame(condition=condition,time=factor(timepoints), batch=batch)
rownames(pheno) <- colnames(counts)

set <- newSeqExpressionSet(counts=counts, phenoData=pheno)


#--------------------------------------------------
# 3. Use edgeR to run first-pass DE analysis
#    Design includes condition + time
#--------------------------------------------------
design <- model.matrix(~condition + time, data=pheno)

y <- DGEList(counts=counts, group=condition)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit)

# Rank genes by p-value (least DE = most stable)
res <- topTags(lrt, n=nrow(counts))$table
neg_controls <- rownames(res)[tail(order(res$PValue), 2000)]  # pick 1000 least DE genes
#--------------------------------------------------
# 4. RUVg within-condition batch correction
#--------------------------------------------------
set_ruv <- RUVg(set, cIdx=neg_controls, k=2) 
# (k = # factors of unwanted variation, can try 1–3)

#--------------------------------------------------
# 5. EDASeq UQ normalization between conditions
#--------------------------------------------------
set_uq <- betweenLaneNormalization(set_ruv, which="upper")

#--------------------------------------------------
# 6. Extract final normalized counts
#--------------------------------------------------
corrected_counts <- normCounts(set_uq)

#write.csv(corrected_counts, "Batch_corrected_est_Count_RUVg_shiju_K_2.csv")

gene_lengths_df<-read.csv('Average_gene_length_Fatbody.csv', row.names = 1)
gene_lengths_filtered<-gene_lengths_df[rownames(gene_lengths_df) %in% rownames(corrected_counts),,drop=F]
#o guarantee matching order with corrected_counts
gene_lengths_filtered <- gene_lengths_df[rownames(corrected_counts), , drop = FALSE]

# Divide counts by gene lengths (in kb or base pairs — consistent units)
rate <- sweep(corrected_counts, 1, gene_lengths_filtered[,1], FUN = "/")

# Normalize per sample to sum to 1 million (TPM)
Corrected_tpm <- sweep(rate, 2, colSums(rate), FUN = "/") * 1e6

write.csv(Corrected_tpm, "Batch_corrected_TPM_RUVg_shiju_K_2_neg_1000.csv")



