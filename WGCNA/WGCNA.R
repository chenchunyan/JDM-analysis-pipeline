# WGCNA Analysis Script 
#
# This script performs a Weighted Gene Co-expression Network Analysis (WGCNA).

# Load necessary libraries
library(WGCNA)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(openxlsx)

# Enable multi-threading
enableWGCNAThreads()

# Set working directory
setwd("/ifs1/User//chenchunyan/WGCNA/")

# =====================================================================================
# 1. Data Loading, Filtering, and Cleaning
# =====================================================================================

# Load expression data
data = read.table("./GSE11083_PBMC_express.txt", header=T, sep="\t", row.names=1, check.names=F)
datExpr0 <- as.data.frame(t(data))

# Check for genes and samples with too many missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Filter genes based on expression level and coefficient of variation
# As per the paper: average expression level > 1 and coefficient of variation (variance/mean) > 0.3
mean_expr <- colMeans(datExpr0)
cv <- apply(datExpr0, 2, function(x) var(x) / mean(x))

genes_to_keep <- mean_expr > 1 & cv > 0.3
datExpr0 <- datExpr0[, genes_to_keep]
print(paste("Number of genes after filtering:", ncol(datExpr0)))

# =====================================================================================
# 2. Outlier Detection and Removal
# =====================================================================================

# Cluster samples to detect outliers
sampleTree = hclust(dist(datExpr0), method = "average")

pdf(file = "1.sampleClustering.pdf", width = 12, height = 9)
par(cex = 0.6, mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# Define a height cut to identify outliers.
# This cut should be chosen based on the dendrogram. Let's set a cut-off height.
# For this dataset, a height of 80 seems reasonable to cut the main outliers.
abline(h = 80, col = "red")
dev.off()

# Determine clusters and remove outlier samples
clust = cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
table(clust)
# Cluster 0 contains the outliers
keepSamples = (clust == 1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Plot the dendrogram of the remaining samples
sampleTree2 = hclust(dist(datExpr), method = "average")
pdf(file = "1.sampleClustering_after_removal.pdf", width = 12, height = 9)
par(cex = 0.6, mar = c(0, 4, 2, 0))
plot(sampleTree2, main = "Sample clustering after removing outliers", sub="", xlab="",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

print(paste("Number of samples kept:", nSamples))

# =====================================================================================
# 3. Soft-thresholding Power Selection
# =====================================================================================

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf(file="3.soft_threshold_selection.pdf", width=9, height=5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red") # R^2 cut-off of 0.9
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# From the plot, we select the power
softPower = sft$powerEstimate
if(is.na(softPower)){
  softPower = 4 # Manually set power if it's not automatically estimated
}
print(paste("Selected soft power:", softPower))

# =====================================================================================
# 4. Network Construction and Module Detection
# =====================================================================================

adjacency = adjacency(datExpr, power = softPower)
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

# Gene clustering on TOM-based dissimilarity
geneTree = hclust(as.dist(dissTOM), method = "average")

# Module identification using dynamic tree cut
minModuleSize = 100
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 4, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

pdf(file="4.gene_dendrogram_and_module_colors.pdf", width=8, height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# Merge modules whose expression profiles are very similar
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")

pdf(file="5.module_eigengene_clustering.pdf", width=7, height=6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
MEDissThres = 0.2 # Corresponds to a correlation of 0.8
abline(h=MEDissThres, col = "red")
dev.off()

merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
moduleColors = mergedColors # Use merged modules for further analysis

pdf(file="4.gene_dendrogram_and_merged_module_colors.pdf", width=12, height=9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


# =====================================================================================
# 5.1 Module Stability Analysis (minModuleSize = 100,deepSplit = 4,cutHeight = 0.9,sampleFrac  <- 0.96,bootstrap=100)
# =====================================================================================
### set Bootstrap parameter
set.seed(12345)
nBoot       <- 100          # ⬆️ increase Bootstrap 
sampleFrac  <- 0.96         # ⬆️ Increase the sampling ratio
nSamples    <- ncol(datExpr)
geneNames   <- rownames(datExpr)
refColors   <- dynamicColors
refLevels   <- unique(refColors)

moduleFreq  <- matrix(0,
                      nrow = length(geneNames),
                      ncol = length(refLevels),
                      dimnames = list(geneNames, refLevels))

for (b in 1:nBoot) {
  cat("Bootstrap iteration", b, "\n")
  bootSamp <- sample(nSamples, size = ceiling(sampleFrac * nSamples), replace = FALSE)
  exprBoot <- t(datExpr)[bootSamp, , drop = FALSE]

  adjBoot   <- adjacency(exprBoot, power = softPower)
  tomBoot   <- TOMsimilarity(adjBoot)
  dissBoot  <- 1 - tomBoot
  treeBoot  <- hclust(as.dist(dissBoot), method = "average")

  modBoot   <- cutreeDynamic(treeBoot,
                             distM = dissBoot,
                             deepSplit = 4,               #Reduce the cutting depth
                             pamRespectsDendro = FALSE,
                             cutHeight = 0.9,
                             minClusterSize = minModuleSize)

  colorBoot <- labels2colors(modBoot)
  colorBoot <- factor(colorBoot, levels = refLevels)

  for (col in refLevels) {
    genesIn <- geneNames[colorBoot == col & !is.na(colorBoot)]
    if (length(genesIn) > 0) {
      moduleFreq[genesIn, col] <- moduleFreq[genesIn, col] + 1
    }
  }
}

# Calculate the stability score
stabScore <- apply(moduleFreq, 1, max) / nBoot
names(stabScore) <- geneNames

#Output Stability Summary Table
stabTable <- data.frame(Gene = names(stabScore),
                        Stability = stabScore,
                        Module = refColors)
write.table(stabTable, "WGCNA_ModuleStability_Scores_50_96.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Output module stability statistics table
threshold <- 0.7
module_stats <- as.data.frame(table(Module = refColors))
module_stats$StableGeneCount <- tapply(stabScore > threshold, refColors, sum)
module_stats$StableRatio <- module_stats$StableGeneCount / module_stats$Freq
write.table(module_stats, "WGCNA_ModuleStability_Summary_50_96.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Drawing: Stability Distribution
library(ggplot2)
stabPlot <- data.frame(Stability = stabScore, Module = refColors)
p <- ggplot(stabPlot, aes(x = Stability, fill = Module)) +
  geom_histogram(binwidth = 0.05, alpha = 0.7, position = "identity") +
  theme_bw() +
  xlab("Stability Score") +
  ylab("Gene Count") +
  ggtitle("Distribution of Gene Stability Scores by Module") +
  scale_fill_brewer(palette = "Set1")

ggsave("Stability_Histogram_By_Module.pdf", p, width = 8, height = 6)

# You can choose to print the highly stable genes in the brown module.
brown_genes <- geneNames[refColors == "brown"]
brown_high_stable <- stabTable[stabTable$Gene %in% brown_genes & stabTable$Stability > threshold, ]
write.table(brown_high_stable, "Brown_Module_High_Stability_Genes.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# =====================================================================================
# 5.2 Module Stability Analysis (minModuleSize = 150,deepSplit = 4,cutHeight = 0.9,sampleFrac  <- 0.96,bootstrap=100)
# =====================================================================================
# ----  Network construction parameter settings ----
softPower <- 4
minModuleSize <- 150
deepSplit <- 4

adjacency <- adjacency(t(datExpr), power = softPower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Dynamic Pruning Recognition Module
moduleLabels <- cutreeDynamic(dendro = geneTree,
                               distM = dissTOM,
                               deepSplit = deepSplit,
                               pamRespectsDendro = FALSE,
                               cutHeight = 0.9,
                               minClusterSize = minModuleSize)

moduleColors <- labels2colors(moduleLabels)

# ---- Correlation between modules and phenotypes----
trait <- read.table("group2.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

MEs <- moduleEigengenes(t(datExpr), colors = moduleColors)$eigengenes
MEs <- orderMEs(MEs)
trait <- trait[match(rownames(MEs), rownames(trait)), , drop = FALSE]

moduleTraitCor <- cor(MEs, trait, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

write.table(moduleTraitCor, file = "ModuleTraitCor.txt", sep = "\t")
write.table(moduleTraitPvalue, file = "ModuleTraitPvalue.txt", sep = "\t")

# ---- Bootstrap module stability assessment ----
set.seed(12345)
nBoot <- 100
sampleFrac <- 0.96
nSamples <- ncol(datExpr)
geneNames <- rownames(datExpr)
refColors <- moduleColors
refLevels <- unique(refColors)

moduleFreq <- matrix(0, nrow = length(geneNames), ncol = length(refLevels), dimnames = list(geneNames, refLevels))

for (b in 1:nBoot) {
  cat("Bootstrap ", b, "\n")
  bootSamp <- sample(nSamples, size = ceiling(sampleFrac * nSamples), replace = FALSE)
  exprBoot <- t(datExpr)[bootSamp, , drop = FALSE]

  adjBoot <- adjacency(exprBoot, power = softPower)
  tomBoot <- TOMsimilarity(adjBoot)
  dissBoot <- 1 - tomBoot
  treeBoot <- hclust(as.dist(dissBoot), method = "average")

  modBoot <- cutreeDynamic(treeBoot,
                           distM = dissBoot,
                           deepSplit = deepSplit,
                           pamRespectsDendro = FALSE,
                           cutHeight = 0.9,
                           minClusterSize = minModuleSize)
  colorBoot <- labels2colors(modBoot)
  colorBoot <- factor(colorBoot, levels = refLevels)

  for (col in refLevels) {
    genesIn <- geneNames[colorBoot == col & !is.na(colorBoot)]
    if (length(genesIn) > 0) {
      moduleFreq[genesIn, col] <- moduleFreq[genesIn, col] + 1
    }
  }
}

stabScore <- apply(moduleFreq, 1, max) / nBoot

# ---- Output stability results ----
stability_table <- data.frame(Gene = names(stabScore),
                              Module = refColors,
                              Stability = stabScore)
write.table(stability_table, "WGCNA_ModuleGene_Stability.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#The proportion of stable genes in each module
stable_threshold <- 0.7
module_stats <- data.frame(Module = refLevels,
                           Freq = sapply(refLevels, function(col) sum(refColors == col)),
                           StableGeneCount = sapply(refLevels, function(col) sum(stability_table$Module == col & stability_table$Stability >= stable_threshold)))
module_stats$StableRatio <- module_stats$StableGeneCount / module_stats$Freq
write.table(module_stats, file = "WGCNA_ModuleStability_Summary.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# ----visualization ----
library(ggplot2)
ggplot(module_stats, aes(x = Module, y = StableRatio, fill = Module)) +
  geom_bar(stat = "identity") +
  labs(title = "Module stability score", y = ">70% Stable gene ratio", x = "Module") +
  theme_minimal()
ggsave("ModuleStabilityBarplot.pdf", width = 6, height = 4)

# =====================================================================================
# 6. Module-Trait Relationship Analysis
# =====================================================================================

# Load trait data
traitData = read.table("group2.txt", header=T, sep="\t", row.names=1, check.names=F)
allTraits = as.data.frame(traitData)

# Match samples in expression data and trait data
fpkmSamples = rownames(datExpr)
traitSamples = rownames(allTraits)
traitRows = match(fpkmSamples, traitSamples)
datTraits = allTraits[traitRows, ]

# Recalculate MEs with merged colors
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes

# Correlate MEs with traits
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Visualize module-trait relationships
pdf(file="6.module-trait_relationships.pdf", width=10, height=8)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# =====================================================================================
# 7. Gene Enrichment Analysis for Significant Modules
# =====================================================================================

# Identify modules significantly associated with the trait (e.g., JDM)
# Let's assume 'JDM' is the trait of interest from the 'group2.txt' file.
# We select the module with the highest absolute correlation with JDM.
jdm_column_index <- which(colnames(datTraits) == "JDM")
if (length(jdm_column_index) == 0) {
    stop("Trait 'JDM' not found in the trait data file.")
}
significant_module_name <- names(MEs)[which.max(abs(moduleTraitCor[, jdm_column_index]))]
print(paste("Most significant module for JDM is:", significant_module_name))


# Perform GO and KEGG enrichment analysis for all modules
modules <- unique(moduleColors)
wb <- createWorkbook()

for (module in modules) {
  if (module == "grey") next # Skip the grey module (ungrouped genes)
  
  module_genes <- names(datExpr)[moduleColors == module]
  
  # Convert gene symbols to Entrez IDs
  entrez_ids <- bitr(module_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  if (nrow(entrez_ids) > 0) {
    # GO Enrichment
    go_enrich <- enrichGO(gene          = entrez_ids$ENTREZID,
                          OrgDb         = org.Hs.eg.db,
                          ont           = "BP", # Biological Process
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.2,
                          readable      = TRUE)
    
    if (!is.null(go_enrich) && nrow(go_enrich@result) > 0) {
      addWorksheet(wb, paste0("GO_BP_", module))
      writeData(wb, sheet = paste0("GO_BP_", module), as.data.frame(go_enrich@result))
    }
    
    # KEGG Enrichment
    kegg_enrich <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                              organism     = 'hsa',
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05)
    
    if (!is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
      addWorksheet(wb, paste0("KEGG_", module))
      writeData(wb, sheet = paste0("KEGG_", module), as.data.frame(kegg_enrich@result))
    }
  }
}
saveWorkbook(wb, "module_enrichment_analysis.xlsx", overwrite = TRUE)


# =====================================================================================
# 8. Output Gene and Module Information
# =====================================================================================

geneModuleInfo <- data.frame(Gene = names(datExpr), ModuleColor = moduleColors)
write.table(geneModuleInfo, file = "Gene_Module_Info_revised.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

print("WGCNA analysis complete. Check the output files.")
