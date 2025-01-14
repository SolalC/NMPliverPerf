library(data.table)
library(tidyverse)
library(sva)
library(edgeR)
library(variancePartition)
library(ggpattern)
library(gprofiler2)
#####################
# Hautz:
hautzCount <- fread('/QRISdata/Q1144/Data/hautz/GSE263614_raw_counts.txt.gz')
hautzInfo <- fread('/QRISdata/Q1144/Data/hautz/hautzSampleInfo.csv')
hautzInfo.f <- hautzInfo %>% mutate(time = ifelse(time == 'PRE', 'CIT', time),
         viability = ifelse(perf == 'TP', 'Viable', 'Non-Viable'),
         viability = ifelse(time == 'CIT', 'CIT', viability),
         batch = 'hautz') %>%
  select(colID, sampleID, batch, time, viability)

colID <- colnames(hautzCount) %in% hautzInfo.f$colID
hautzCount.f <- hautzCount[,..colID]
hautzCount.f <- as.data.frame(hautzCount.f)
rownames(hautzCount.f) <- hautzCount$GENE
# library Size Parameters:
hautzParam <- data.frame(colID = colnames(hautzCount.f),
                         LS = colSums(hautzCount.f),
                         nZero = apply(hautzCount.f, 2, function(x) sum(x == 0)),
                         nNotZero = apply(hautzCount.f, 2, function(x) sum(x > 0))) %>% left_join(., hautzInfo.f)
#####################
# Raigani dataset:
files <- list.files("/QRISdata/Q1144/Data/raigani", full.names = T)
raiganiInfo <- fread(files[which(str_detect(files, "Information"))])
files <- files[str_detect(files, ".txt.gz")]
countRaigan <- matrix(nrow = 60675)
for (i in files) {
  temp <- read.table(i)
  colnames(temp) <- c('gene', i)
  countRaigan <- cbind(countRaigan, temp)
}
countRaigan <- countRaigan[, -1]
genes <- countRaigan[, seq(1, ncol(countRaigan), by = 2)]
countRaigan <- countRaigan[, seq(2, ncol(countRaigan), by = 2)]
rownames(countRaigan) <- genes[, 1]
colnames(countRaigan) <- str_remove(basename(colnames(countRaigan)), '_CountTable.txt.gz')
# RaiganiInfo:
raiganiInfo.f <- raiganiInfo %>% dplyr::filter(condition %in%
                                               c("Control", "Excluded Control"))
raiganiInfo.f <- raiganiInfo.f %>% mutate(
  time = ifelse(perf == "0 hours perfusion", "CIT", perf),
  time = ifelse(time == "3 hours perfusion", "3H", time),
  time = ifelse(time == "6 hours perfusion", "6H", time),
  viability = ifelse(perf_success == 'Fatty Adequate Hepatic Function Liver', 'Viable', perf),
  viability = ifelse(perf_success == 'Fatty Inadequate Hepatic Function Liver', 'Non-Viable', viability),
  viability = ifelse(time == 'CIT', 'CIT', viability),
  sampleID = str_replace(sampleID, '-', '_'),
  batch = 'raigani'
) %>% select(colID, sampleID, batch, time, viability)
# Get the dataset:
raiganicolID <- raiganiInfo.f$colID
countRaigan.f <- countRaigan[, raiganicolID]
# Convert to gene symbol:
ensembl.genes <- rownames(countRaigan.f)
gene.symbols <- gconvert(ensembl.genes,
                         organism = "hsapiens",
                         target = "ENTREZGENE", filter_na = F
)
gene.symbols <- gene.symbols[!duplicated(gene.symbols$input), ]
countRaigan.f$gene <- gene.symbols$target
countRaigan.f <- countRaigan.f[!is.na(countRaigan.f$gene), ]
countRaigan.f <- countRaigan.f[!duplicated(countRaigan.f$gene), ]
rownames(countRaigan.f) <- countRaigan.f$gene
countRaigan.f <- countRaigan.f[,-which(colnames(countRaigan.f) == 'gene')]
# library Size Parameters:
raiganParam <- data.frame(colID = colnames(countRaigan.f),
                         LS = colSums(countRaigan.f),
                         nZero = apply(countRaigan.f, 2, function(x) sum(x == 0)),
                         nNotZero = apply(countRaigan.f, 2, function(x) sum(x > 0))) %>% left_join(., raiganiInfo.f)
#####################
# in-house:
inHouseCount <- fread('/QRISdata/Q1144/Data/inHouse/NMPcountData.csv')
inHouseSampleInfo <- fread('/QRISdata/Q1144/Data/inHouse/NMPsampleCovariates.csv') 
inHouseSampleInfo.f <- inHouseSampleInfo %>%
  mutate(viability = ifelse(perf_success == 'Perf_Successful', 'Viable', perf_success),
         viability = ifelse(viability == 'Perf_Unsuccessful', 'Non-Viable', viability),
         time = ifelse(perf == 'CIT_start', 'CIT', perf),
         time = ifelse(time == 'CIT_end', 'CIT', time),
         time = ifelse(time == '3h_perf', '3H', time),
         time = ifelse(time == '4h_perf', '3H', time),
         time = ifelse(time == '6h_perf', '6H', time),
         sampleID = sample,
         batch = 'inHouse') %>%
  select(colID, sampleID, batch, time, viability)

colID <- inHouseSampleInfo.f$colID
inHouseCount.f <- as.data.frame(inHouseCount[, ..colID])
rownames(inHouseCount.f) <- inHouseCount$V1
# Get library size:
inHouseLS <- colSums(inHouseCount.f)
# Get the number of genes expressing zero:
inHouseParam <- data.frame(colID = colnames(inHouseCount.f),
                          LS = colSums(inHouseCount.f),
                          nZero = apply(inHouseCount.f, 2, function(x) sum(x == 0)),
                          nNotZero = apply(inHouseCount.f, 2, function(x) sum(x > 0))) %>% left_join(., inHouseSampleInfo.f)
#####################
# Merge the three datasets:
hautzGenes <- rownames(hautzCount.f)
raiganiGenes <- rownames(countRaigan.f)
inHouseGenes <- rownames(inHouseCount.f)
intersectGenes <- Reduce(intersect, list(hautzGenes, raiganiGenes, inHouseGenes))
# Filter to keep only the common genes:
hautzCount.ff <- hautzCount.f[intersectGenes,]
countRaigan.ff <- countRaigan.f[intersectGenes,]
inHouseCount.ff <- inHouseCount.f[intersectGenes,]
# Merge the three datasests:
mergedCount <- cbind(hautzCount.ff, countRaigan.ff, inHouseCount.ff)
# Merge the sample information:
mergedInfo <- rbind(hautzInfo.f, raiganiInfo.f, inHouseSampleInfo.f)
# Save the data:
write.csv(mergedCount, '/QRISdata/Q1144/Data/merged/mergedCountData.csv')
write.csv(mergedInfo, '/QRISdata/Q1144/Data/merged/mergedSampleInfo.csv')
# Compare the library sizes and capture gene expression from each of the study:
RNAparam <- rbind(hautzParam, raiganParam, inHouseParam) %>% mutate(nGene = nNotZero + nZero)
# Library Size:
ggplot(RNAparam) + 
  geom_point(aes(x = viability, y = LS, color = viability)) +
  facet_wrap(~batch) +
  theme_minimal() -> librarySize
# Number of zero:
ggplot(RNAparam) + 
  geom_point(aes(x = viability, y = nZero, color = viability)) +
  facet_wrap(~batch) +
  theme_minimal() -> nZeroPlot
# Ratio:
ggplot(RNAparam) + 
  geom_point(aes(x = viability, y = nZero/nGene, color = viability)) +
  facet_wrap(~batch) +
  theme_minimal() -> percentZero
################### End of the file.
# Create the DGElist object
dge <- DGEList(
  counts = hautz.f, 
  samples = hautzInfo.f,
  genes = rownames(hautz.f)
)
# Double check everything is ordered correctly
# all(colnames(dge$counts) == rownames(dge$samples))
# # Set null and alternative models (ignore batch)
# # full model matrix - including both the adjustment variables and the
# # variable of interest
# mod1 <- model.matrix(~ perf, data = dge$samples)
# # Null model matrix contains only the adjustment variables
# mod0 <- model.matrix(~ 1, data = dge$samples)
# # filter low number genes:
# filt <- apply(dge$counts, 1, function(x) length(x[x > 5]) >= 2)
# filtered <- as.matrix(dge$counts[filt, ])
# # Run SVA for sequencing data - restrict to only 2 SVs
# svseq <- svaseq(filtered, mod1, mod0)
# colnames(svseq$sv) <- paste0("SV", seq(1, ncol(svseq$sv)))
# # Add the SVA to the dataframe:
# hautzInfo.f <- cbind(hautzInfo.f, svseq$sv)
# rownames(hautzInfo.f) <- hautzInfo.f$colID
# # DEG:
# design <- ~ (1 | sampleID) + SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + 
#   SV9 + SV10 + SV11 + SV12 +  perf
# # filter genes by number of counts
# isexpr <- rowSums(cpm(hautz.f) > 0.1) >= 5
# # Standard usage of limma/voom
# dge <- DGEList(hautz.f[isexpr, ])
# dge <- calcNormFactors(dge)
# vobjDream <- voomWithDreamWeights(dge, design, hautzInfo.f)
# fitmm <- dream(vobjDream, design, hautzInfo.f)
# fitmm <- eBayes(fitmm)
# # Viable results:
# viableResults <- topTable(fitmm, coef = "perfViable", n = nrow(dge)) %>%
#   filter(adj.P.Val < 0.05, !between(logFC, -0.25, 0.25)) %>%
#   mutate(gene = rownames(.))
# # Non-Viable results
# nonViableResults <- topTable(fitmm,coef = "perfNon-Viable", n = nrow(dge)) %>%
#   filter(adj.P.Val < 0.05, !between(logFC, -0.25, 0.25)) %>%
#   mutate(gene = rownames(.))
# # Read my own results:
# deg <- fread('/QRISdata/Q1144/Results/RNAseq/DEG/CIT_viable_significant_DEG.csv')
# sum(deg$gene %in% viableResults$gene)
# 
# # Circadian stuff:
# circGenes <- c( "CLOCK", "ARNTL", "DBP", "NPAS2", "PER1", "PER2", "PER3",
#                 "TEF", "HLF", "CRY1", "CRY2", "NR1D1")
# 
# CLOCK_GTEx <- read.table("/QRISdata/Q1144/Results/RNAseq/CLOCK/GTExMatrixRef.csv")
# matrix_levels <- c(
#   "ARNTL", "NPAS2", "CLOCK", "NFIL3", "CRY1",
#   "CRY2", "NR1D1", "NR1D2", "PER1", "PER2",
#   "PER3", "DBP", "TEF", "HLF"
# )
# matrix_levels <- str_to_title(matrix_levels)
# # Correlation vector:
# CLOCK_GTEx <- CLOCK_GTEx[order(rownames(CLOCK_GTEx)), ]
# rownames(CLOCK_GTEx) <- str_to_title(rownames(CLOCK_GTEx))
# GTEXcorMatrix <- cor(t(CLOCK_GTEx))
# # Get the correlation matrix for each of the experimental conditions:
# # Order the dataframe in the same way than the GTEx matrix
# # Generate the dataset of interest:
# cpmData <- cpm(hautz.f)
# rownames(cpmData) <- str_to_title(rownames(hautz.f))
# var <- apply(cpmData, 1, var)
# cpmData <- cpmData[which(var > 1), ]
# 
# cpmCIT <- cpmData[, hautzInfo.f$perf == "CIT"]
# cpmFailed <- cpmData[, hautzInfo.f$perf == "Non-Viable"]
# cpmSuccess <- cpmData[, hautzInfo.f$perf == "Viable"]
# # Order the dataframes:
# cpmCITordered <- cpmCIT[order(rownames(cpmCIT)), ]
# cpmFailedOrdered <- cpmFailed[order(rownames(cpmFailed)), ]
# cpmSuccessOrdered <- cpmSuccess[order(rownames(cpmSuccess)), ]
# # Create the correlation matrix:
# CITcorMatrix <- cor(t(cpmCITordered[which(rownames(cpmCITordered) %in% matrix_levels), ]))
# # Viable livers:
# FailedCorMatrix <- cor(t(cpmFailedOrdered[which(rownames(cpmFailedOrdered) %in% matrix_levels), ]))
# # Non-viable Livers:
# SuccessCorMatrix <- cor(t(cpmSuccessOrdered[which(rownames(cpmSuccessOrdered) %in% matrix_levels), ]))
# # Extract the correlation vector for each:
# GTEXcorVector <- GTEXcorMatrix[lower.tri(GTEXcorMatrix, diag = F)]
# CITcorVector <- CITcorMatrix[lower.tri(CITcorMatrix, diag = F)]
# FailedCorVector <- FailedCorMatrix[lower.tri(FailedCorMatrix, diag = F)]
# SuccessCorVector <- SuccessCorMatrix[lower.tri(SuccessCorMatrix, diag = F)]
# # Calculate the euclidean distance between the two:
# euclidean <- function(a, b) sqrt(sum((a - b)^2))
# CCDCIT <- euclidean(GTEXcorVector, CITcorVector)
# CCDFailed <- euclidean(GTEXcorVector, FailedCorVector)
# CCDSuccess <- euclidean(GTEXcorVector, SuccessCorVector)
# # Bootstrap distribution:
# bootstrapDistribution <- function(cpm, clockGenes, n = 1000, refCorVector) {
#   CCD <- c()
#   while (length(CCD) < n) {
#     meanVector <- rowMeans(cpm[which(rownames(cpm) %in% clockGenes), ])
#     # Randomly sample a set of the same size as the clock:
#     cpmSamplingSet <- cpm[-which(rownames(cpm) %in% clockGenes), ]
#     set <- sample(rownames(cpmSamplingSet), size = length(clockGenes))
#     # Calculate the mean expression of those genes:
#     setMeanVector <- rowMeans(cpmSamplingSet[set, ])
#     # If the genes have a similar mean to the CLOCK genes, accept the set:
#     if (t.test(meanVector, setMeanVector)$p.value > 0.05) {
#       # Accept the set and calculate the CCD:
#       SetCorMatrix <- cor(t(cpmSamplingSet[set, ]))
#       SetCorVector <- SetCorMatrix[lower.tri(SetCorMatrix, diag = F)]
#       CCD <- c(CCD, euclidean(refCorVector, SetCorVector))
#     }
#   }
#   return(CCD)
# }
# nperm <- 10000
# 
# CCDdistribCIT <- bootstrapDistribution(
#   cpm = cpmCITordered,
#   clockGenes = matrix_levels,
#   n = nperm, refCorVector = GTEXcorVector
# )
# CCDdisribFailed <- bootstrapDistribution(
#   cpm = cpmFailedOrdered,
#   clockGenes = matrix_levels,
#   n = nperm, refCorVector = GTEXcorVector
# )
# CCDdistribSuccess <- bootstrapDistribution(
#   cpm = cpmSuccessOrdered,
#   clockGenes = matrix_levels,
#   n = nperm, refCorVector = GTEXcorVector
# )
# # Calculate the p-value:
# pCIT <- pnorm(
#   q = CCDCIT, mean = mean(CCDdistribCIT, na.rm = T),
#   sd = sd(CCDdistribCIT, na.rm = T)
# )
# pFail <- pnorm(
#   q = CCDFailed, mean = mean(CCDdisribFailed, na.rm = T),
#   sd = sd(CCDdisribFailed, na.rm = T)
# )
# pSuccess <- pnorm(
#   q = CCDSuccess, mean = mean(CCDdistribSuccess, na.rm = T),
#   sd = sd(CCDdistribSuccess, na.rm = T)
# )
# pval <- p.adjust(c(pCIT, pSuccess, pFail), method = "bonfe")
# pval <- format(signif(pval, 3), scientific = T)
# 
# CCDdataFrame <- data.frame(
#   distance = c(CCDdistribCIT, CCDdisribFailed, CCDdistribSuccess),
#   condition = c(
#     rep("CIT", nperm),
#     rep("Non-Viable", nperm),
#     rep("Viable", nperm)
#   )
# )
# 
# title <- paste0(
#   "CCD CIT: ", round(CCDCIT, 2),
#   "; P-value:", pval[1], "\n",
#   "CCD Viable: ", round(CCDSuccess, 2),
#   "; P-value:", pval[2], "\n",
#   "CCD Non-Viable: ", round(CCDFailed, 2),
#   "; P-value:", pval[3], "\n"
# )
# 
# ggplot(CCDdataFrame, aes(x = distance, fill = condition)) +
#   geom_density(alpha = 0.5) +
#   geom_vline(xintercept = CCDCIT, col = "#BEDEE9", linewidth = 1) +
#   geom_vline(xintercept = CCDFailed, col = "#A6BCA2", linewidth = 1) +
#   geom_vline(xintercept = CCDSuccess, col = "#F0C56F", linewidth = 1) +
#   scale_fill_manual(
#     values = c("#BEDEE9", "#A6BCA2", "#F0C56F"),
#     name = "Empirical Distribution for:"
#   ) +
#   theme_minimal() +
#   xlab("Clock Correlation Distance") +
#   xlim(0, 7.5) +
#   ggtitle(title) -> ccdplot
# ggsave(ccdplot, filename = "/QRISdata/Q1144/Results/RNAseq/replication/CCDplot.jpeg", width = 10, height = 7, dpi = 300)
# 
# ### Try to run this locally:
# hautz <- fread('/Users/uqschauq/Documents/NMP/Data/hautzetalRNA/GSE263614_raw_counts.txt.gz')
# hautzInfo <- fread('/Users/uqschauq/Documents/NMP/Data/hautzetalRNA/hautzSampleInfo.csv')
# hautzInfo.f <- hautzInfo %>% filter(time %in% c('PRE','6H')) %>%
#   mutate(time = ifelse(time == 'PRE', 'CIT', '6H'),
#          perf = ifelse(perf == 'TP', 'Viable', 'Non-Viable'),
#          perf = ifelse(time == 'CIT', 'CIT', perf))  %>%
#   as.data.frame()
# 
# colID <- colnames(hautz) %in% hautzInfo.f$colID
# hautz.f <- hautz[,..colID]
# hautz.f <- as.data.frame(hautz.f)
# rownames(hautz.f) <- hautz$GENE
# # Heatmap:
# lvl <- c(
#   hautzInfo.f[which(hautzInfo.f$perf == "CIT"), ]$colID,
#   hautzInfo.f[ which(hautzInfo.f$perf == "Viable"),]$colID,
#   hautzInfo.f[which(hautzInfo.f$perf == "Non-Viable"),]$colID
# )
# # Annotation of the heatmap:
# rownames(hautzInfo.f) <- hautzInfo.f$colID
# sample_col <- data.frame(
#   sample = hautzInfo.f[lvl, ]$perf,
#   perfusion = as.vector(hautzInfo.f[lvl, ]$time)
# )
# row.names(sample_col) <- lvl
# 
# # Change the names for the figure:
# sample_col$sample[which(sample_col$sample == "Non-Viable")] <- "Non-Viable Liver"
# sample_col$sample[which(sample_col$sample == "Viable")] <- "Viable Liver"
# # Change the perfusion names:
# sample_col$perfusion[which(sample_col$perfusion == "6H")] <- "Perfusion: 6hours"
# # Min-max normalization for the heatmap:
# a <- -1
# b <- 1
# cal_norm <- function(x) {
#   a + ((x - min(x)) * (b - a) / (max(x) - min(x)))
# }
# hmp_count <- hautz.f[which(rownames(hautz.f) %in% unique(viableResults$gene)), ]
# hmp_count <- hmp_count[-which(rowSums(hmp_count) == 0),]
# hmp_norm <- t(apply(hmp_count, 1, cal_norm))
# hmp_norm <- hmp_norm[, lvl]
# # Create the heatmap:
# hmpColors <- met.brewer(name = "Ingres", n = 4)
# ann_colors <- list(
#   perfusion = c(
#     `CIT` = hmpColors[2],  `Perfusion: 6hours` = hmpColors[4]
#   ),
#   sample = c(
#     CIT = "light blue", `Viable Liver` = VolColorViable[2],
#     `Non-Viable Liver` = VolColorNonViable[2]
#   )
# )
# 
# # Create the heatmap:
# hmpAllClustering <- pheatmap(hmp_norm,
#                              annotation_col = sample_col, cluster_cols = T,
#                              show_colnames = F, annotation_colors = ann_colors,
#                              labels_row = F, show_rownames = F,
#                              treeheight_row = 0
# )
# ggsave(hmpAllClustering, filename = "/Users/uqschauq/Documents/NMP/Results/Replication/Heatmap.jpeg", width = 10, height = 7, dpi = 300)
# # Classifier:
# # X test:
# 
# sampleToRemoveHautz <- hautzInfo.f %>%
#   group_by(sampleID) %>%
#   mutate(n = n()) %>%
#   filter(n == 1) %>%
#   pull(colID)
# 
# Xtest <- cpm(hautz.f) %>%
#   as.data.frame() %>%
#   filter(rownames(.) %in% successfulSignature) %>%
#   as.matrix()
# Xtest <- Xtest[,!colnames(Xtest) %in% sampleToRemoveHautz]
# Xtest <- t(Xtest[-which(rowVars(Xtest) == 0), ])
# Xtest <- Xtest[, order(colnames(Xtest))]
# # Y test:
# hautzInfo.f
# 
# Ytest <- factor(hautzInfo.f %>% filter(!colID %in% sampleToRemoveHautz) %>%
#                   mutate(
#   perf = ifelse(perf == "Viable", "Perf_Successful", perf),
#   perf = ifelse(perf == "Non-Viable", "Perf_Unsuccessful", perf)) %>%
#     pull(perf),
#   levels = c("CIT", "Perf_Unsuccessful", "Perf_Successful"))
# 
# # Extract our data to create a classifier:
# Xtrain <- t(cpm(countData_sub)[
#   which(rownames(countData_sub) %in% colnames(Xtest)),
#   -which(colnames(countData_sub) %in% c("NMP4_3", "NMP4_4"))
# ])
# Xtrain <- Xtrain[order(rownames(Xtrain)), order(colnames(Xtrain))]
# # Remove the two misclassified samples:
# Ytrain <- sampleinfo_sub %>% filter(!colID %in% c("NMP4_3", "NMP4_4"))
# Ytrain <- Ytrain[order(Ytrain$colID), ]
# Ytrain$perf_success <- factor(Ytrain$perf_success, levels = c(
#   "CIT",
#   "Perf_Unsuccessful", "Perf_Successful"
# ))
# 
# sampleToRemove <- Ytrain %>%
#   group_by(sample) %>%
#   mutate(n = n()) %>%
#   filter(n == 1) %>%
#   pull(colID)
# Ytrain <- Ytrain %>% filter(!colID %in% sampleToRemove)
# Xtrain <- Xtrain[Ytrain$colID, ]
# # Build the classifier:
# # Create the factor for the sample origin:
# Ytrain$sample <- factor(Ytrain$sample)
# design <- data.frame(sample = Ytrain$sample)
# # Tune:
# list.keepX <- c(1:10, seq(20, 100, 5))
# # undergo the tuning process to determine the optimal number of variables
# tune.splsda <- tune.splsda(
#   X = Xtrain,
#   Y = Ytrain$perf_success,
#   ncomp = 3,
#   validation = "Mfold",
#   folds = 3, nrepeat = 100, # use repeated cross-validation
#   dist = "max.dist", # use max.dist measure
#   measure = "BER", # use balanced error rate of dist measure
#   test.keepX = list.keepX,
#   cpus = 4
# )
# # Extract parameters:
# optimal.ncomp <- tune.splsda$choice.ncomp$ncomp
# optimal.keepX <- tune.splsda$choice.keepX[1:optimal.ncomp]
# # Create the final model:
# final.splsda.multilevel <- splsda(Xtrain, Ytrain$perf_success,
#                                   multilevel = design,
#                                   ncomp = optimal.ncomp,
#                                   keepX = optimal.keepX
# )
# # Calculate the background:
# background <- background.predict(final.splsda.multilevel,
#                                  comp.predicted = 2,
#                                  dist = "mahalanobis.dist",
#                                  xlim = c(-12.5, 12.5), ylim = c(-10, 10)
# )
# plotIndiv(final.splsda.multilevel,
#           group = Ytrain$perf_success,
#           ind.names = design$sample,
#           pch = as.factor(design$sample), legend.title.pch = "Sample",
#           legend = TRUE, legend.title = "Time",
#           title = "Sample Plot of sPLS-DA", background = background
# )
# # Calculate the model stability:
# perf.final.multilevel <- perf(final.splsda.multilevel,
#                               folds = 3, nrepeat = 1000, validation = "loo",
#                               dist = "max.dist"
# )
# par(mfrow = c(1, 2))
# plot(perf.final.multilevel$features$stable[[1]],
#      type = "h", ylab = "Stability",
#      xlab = "Features", main = "(a) Comp 1", las = 2
# )
# plot(perf.final.multilevel$features$stable[[2]],
#      type = "h", ylab = "Stability",
#      xlab = "Features", main = "(b) Comp 2", las = 2
# )
# # Prediction:
# YtestSample <- data.frame(sample = factor((hautzInfo.f %>% filter(!colID %in% sampleToRemoveHautz))$sampleID))
# # YtestSample <- rbind(YtestSample, data.frame(sample = rep('nmp4', 2)))
# 
# predict.splsda <- predict(final.splsda.multilevel,
#                           newdata = Xtest,
#                           dist = "mahalanobis.dist",
#                           multilevel = YtestSample
# )
# # Prediction:
# predict.comp2 <- predict.splsda$class$mahalanobis.dist[, 2]
# table(factor(predict.comp2, levels = levels(Ytest)), Ytest)
# 
# ############### Plotting:
# PLStrain <- as.data.frame(final.splsda.multilevel$variates$X)
# PLStrain$set <- "train"
# PLStrain$class <- final.splsda.multilevel$Y
# PLStrain$sample <- rownames(PLStrain)
# # Test samples:
# PLStest <- as.data.frame(predict.splsda$variates)
# colnames(PLStest) <- c("comp1", "comp2")
# PLStest$sample <- rownames(PLStest)
# PLStest$class <- Ytest
# PLStest$set <- "test"
# 
# PLS <- rbind(PLStrain, PLStest)
# nameChange <- as.vector(PLS$class)
# nameChange[nameChange == "Perf_Successful"] <- "Viable Liver"
# nameChange[nameChange == "Perf_Unsuccessful"] <- "Non-Viable Liver"
# PLS$class <- factor(nameChange, levels = c(
#   "CIT", "Non-Viable Liver",
#   "Viable Liver"
# ))
# predicted <- data.frame(
#   sample = names(predict.comp2),
#   predicted = predict.comp2
# )
# predicted$predicted[
#   predicted$predicted == "Perf_Unsuccessful"
# ] <- "Non-Viable Liver"
# predicted$predicted[
#   predicted$predicted == "Perf_Successful"
# ] <- "Viable Liver"
# PLS <- left_join(PLS, predicted, by = "sample")
# PLS$predicted <- PLS$class == PLS$predicted
# PLS$set[PLS$predicted == F] <- "Misclassified"
# # PLS[str_detect(PLS$sample, 'NMP4'),]$set <- 'other'
# # scatter plot:
# col <- c("light blue", VolColorNonViable[2], VolColorViable[2])
# PLS.f <- PLS %>% filter(set %in% c("Misclassified", "test", "other"))
# ggplot(PLS.f) +
#   geom_point(aes(x = comp1, y = comp2, col = class, shape = set, group = class),
#              size = 4
#   ) +
#   theme_light() +
#   # stat_ellipse(aes(x = comp1, y = comp2, col = class), data = PLS %>%
#   #                filter(set == "test"), level = 0.95, lty = 2) +
#   scale_color_manual(values = col, name = "Real Label:") +
#   scale_shape_manual(values = c(3, 19, 18)) +
#   theme(legend.position = "bottom") +
#   xlab("Component 1 (55%)") +
#   ylab("Component 2 (4%)") +
#   geom_point(aes(x = Var1, y = Var2),
#              data = as.data.frame(background$CIT)
#   ) +
#   geom_point(aes(x = Var1, y = Var2),
#              data = as.data.frame(background$Perf_Successful)
#   ) +
#   geom_point(aes(x = Var1, y = Var2),
#              data = as.data.frame(background$Perf_Unsuccessful)
#   ) -> plsPlot
# ggExtra::ggMarginal(plsPlot, type = "boxplot", groupColour = T) -> plsPlot
# ggsave(plsPlot, filename = "/Users/uqschauq/Documents/NMP/Results/Replication/PLS.pdf", width = 10, height = 10, dpi = 300)
# # Get metrics:
# getMetrics <- function(predicted, true, label) {
#   # Predicted
#   predicted <- as.vector(predicted$class$mahalanobis.dist[, 2])
#   predicted[which(predicted != label)] <- "other"
#   # True
#   true <- as.vector(true)
#   true[which(true != label)] <- "other"
#   # Predictions:
#   TP <- sum(which(predicted == label) %in% which(true == label))
#   TN <- sum(which(predicted == "other") %in% which(true == "other"))
#   FP <- sum(which(predicted == "other") %in% which(true == label))
#   FN <- sum(which(predicted == label) %in% which(true == "other"))
#   # Metrics:
#   SPC <- TN / (FP + TN)
#   SENS <- TP / (TP + FN)
#   ACC <- (TP + TN) / (TP + TN + FP + FN)
#   PPV <- TP / (TP + FP)
#   return(c(SPC, SENS, PPV, ACC))
# }
# 
# auc_mannWhitney <- function(y, pred, class) {
#   y <- as.vector(Ytest)
#   pred <- as.vector(pred)
#   # Make it a multiclass classifier:
#   y[which(y != class)] <- "Other"
#   pred[which(pred != class)] <- "Other"
#   # Make y logical:
#   table <- table(y)
#   if (table[1] < table[2]) {
#     y <- y != class
#   } else {
#     y <- y == class
#   }
#   # Calculate the AUC using a Mann-Whitney test:
#   n1 <- sum(y)
#   n2 <- sum(!y)
#   R1 <- sum(rank(pred)[y])
#   U1 <- R1 - n1 * (n1 + 1) / 2
#   U1 / (n1 * n2)
# }
# 
# 
# aucCIT <- 1-auc_mannWhitney(Ytest, predict.comp2, "CIT")
# aucSuccess <- auc_mannWhitney(Ytest, predict.comp2, "Perf_Successful")
# aucFail <- auc_mannWhitney(Ytest, predict.comp2, "Perf_Unsuccessful")
# 
# 
# 
# # Get their classifier in ours:
# classifierGenes <- c('FAM43A', 'PKD2L1', 'CD274', 'PLD6', 'MLPH', 'DBP', 'WDR72')
# Ytest <- sampleinfo_sub %>% filter(perf_success %in% c('Perf_Successful', 
#                                                        'Perf_Unsuccessful')) %>%  dplyr::select(colID, perf_success)
# 
# logData <- log2(countData_sub+1)
# 
# classifierMatrix <- logData[rownames(logData) %in%  classifierGenes,colnames(logData) %in% Ytest$colID]
# 
# paperClassifier <- data.frame()
# for(i in colnames(classifierMatrix)) {
#   temp <- classifierMatrix[,i]
#   names(temp) <- rownames(classifierMatrix)
#   logP1P <- -0.385 + 
#     0.593*temp[names(temp)== 'FAM43A'] + 
#     0.917*temp[names(temp)== 'PKD2L1'] + 
#     0.081*temp[names(temp)== 'CD274'] - 
#     0.619*temp[names(temp)== 'PLD6'] - 
#     0.437*temp[names(temp)== 'MLPH'] + 
#     0.415*temp[names(temp)== 'DBP'] + 
#     0.592*temp[names(temp)== 'WDR72']
#   
#   Ytest  %>% filter(colID == i) %>% pull(perf_success) -> perf
#   pred = data.frame(colID = i, logP1P = logP1P, truth = perf) %>%
#     mutate(pred = ifelse(logP1P < 1, 'Perf_Successful', 'Perf_Unsuccessful'))
#   paperClassifier <- rbind(paperClassifier, pred)
# }
# table(paperClassifier$pred, paperClassifier$truth)
