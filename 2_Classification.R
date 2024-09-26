library(tidyverse)
library(MetBrewer)
library(pheatmap)
library(mixOmics)
library(gprofiler2)
library(edgeR)
library(matrixStats)


# Read the data:
VolColorViable <- c(met.brewer(name = "Renoir", n = 2))
VolColorNonViable <- c(met.brewer(name = "Hokusai1", n = 2))
oneDriveOutput <-
"/Users/uqschauq/Library/CloudStorage/OneDrive-TheUniversityofQueensland/NMP_liver/Results/"

# Read count data:
countData <- read.delim(
  "/Users/uqschauq/Documents/NMP/Data/RNAseq/featureCountsMatrixClean.csv",
  sep = ",", row.names = 1
)
# Read sample information:
sampleinfo <- readr::read_tsv(
  "/Users/uqschauq/Documents/NMP/Data/RNAseq/sampleinfo.tsv"
)
# Select only the sample of interest:
sampleinfo_sub <- sampleinfo[which(sampleinfo$perf_success %in%
  c(
    "CIT", "Perf_Unsuccessful",
    "Perf_Successful"
  )), ]
# Remove duplicated samples:
sampleinfo_sub <- sampleinfo_sub %>% filter(!colID %in% c(
  "NMP3_3", "NMP16_5",
  "NMP11_5", "NMP6_3",
  "NMP4_5", "NMP15_5"
))
sampleinfo_sub[sampleinfo_sub$colID == "NMP4_4", ]$perf <- "6h_perf"
sampleinfo_sub$RIN <- c(
  8.8, 7.2, 7.3, 8.7, 7.7, 7.7, 7.3, 8.2, 8.6, 9.1, 9.2, 9.2, 8.3, 8.8, 7.6,
  7.5, 8.7, 8.8, 8.6, 7.8, 8.6, 8.7, 7.9, 8.2, 8.8, 8.7, 8.5, 8.6, 8.8, 7.4,
  6.6, 9.2
)
# Order the samples:
sampleinfo_sub <- sampleinfo_sub[order(sampleinfo_sub$colID), ]
rownames(sampleinfo_sub) <- sampleinfo_sub$colID
# Order the count data and only keep the samples selected:
countData_sub <- countData[, which(colnames(countData) %in%
  sampleinfo_sub$colID)]
countData_sub <- countData_sub[, which(colnames(countData_sub) %in%
  (sampleinfo_sub$colID))]
countData_sub <- countData_sub[, order(colnames(countData_sub))]
colnames(countData_sub) == sampleinfo_sub$colID

# Read the DEG resutls:
viableResults <- read.csv(paste0(
  oneDriveOutput,
  "Supplementary/CIT_viable_significant_DEG.csv"
))
# Non viable:
nonViableResults <- read.csv(paste0(
  oneDriveOutput,
  "Supplementary/CIT_NonViable_significant_DEG.csv"
))
# Heatmap:
lvl <- c(
  sampleinfo_sub[which(sampleinfo_sub$perf_success == "CIT"), ]$colID,
  sampleinfo_sub[
    which(sampleinfo_sub$perf_success == "Perf_Unsuccessful"),
  ]$colID,
  sampleinfo_sub[
    which(sampleinfo_sub$perf_success == "Perf_Successful"),
  ]$colID
)
# Annotation of the heatmap:
rownames(sampleinfo_sub) <- sampleinfo_sub$colID
sample_col <- data.frame(
  sample = sampleinfo_sub[lvl, ]$perf_success,
  perfusion = as.vector(sampleinfo_sub[lvl, ]$perf)
)
row.names(sample_col) <- lvl

# Change the names for the figure:
sample_col$sample[which(sample_col$sample == "Perf_Unsuccessful")] <-
  "Non-Viable Liver"
sample_col$sample[which(sample_col$sample == "Perf_Successful")] <-
  "Viable Liver"
# Change the perfusion names:
sample_col$perfusion[
  which(sample_col$perfusion == "CIT_start")
] <- "CIT: start"
sample_col$perfusion[
  which(sample_col$perfusion == "CIT_end")
] <- "CIT: end"
sample_col$perfusion[
  which(sample_col$perfusion == "3h_perf")
] <- "Perfusion: 3hours"
sample_col$perfusion[
  which(sample_col$perfusion == "4h_perf")
] <- "Perfusion: 3hours"
sample_col$perfusion[
  which(sample_col$perfusion == "6h_perf")
] <- "Perfusion: 6hours"
# Min-max normalization for the heatmap:
a <- -1
b <- 1
cal_norm <- function(x) {
  a + ((x - min(x)) * (b - a) / (max(x) - min(x)))
}

hmp_count <- countData_sub[which(rownames(countData_sub) %in%
  unique(viableResults$gene)), ]
hmp_norm <- t(apply(hmp_count, 1, cal_norm))
hmp_norm <- hmp_norm[, lvl]
# Heatmap
hmpColors <- met.brewer(name = "Ingres", n = 4)
ann_colors <- list(
  perfusion = c(
    `CIT: start` = hmpColors[2], `CIT: end` = hmpColors[1],
    `Perfusion: 3hours` = hmpColors[3], `Perfusion: 6hours` = hmpColors[4]
  ),
  sample = c(
    CIT = "light blue", `Viable Liver` = VolColorViable[2],
    `Non-Viable Liver` = VolColorNonViable[2]
  )
)

# Create the heatmap:
hmpAllClustering <- pheatmap(hmp_norm,
  annotation_col = sample_col, cluster_cols = T,
  show_colnames = F, annotation_colors = ann_colors,
  labels_row = F, show_rownames = F,
  treeheight_row = 0
)
# Transfer the object as a ggplot object for saving purposes:
hmpAllClustering <- ggplotify::as.ggplot(hmpAllClustering)
# Save the heatmap:
ggsave(hmpAllClustering, filename = paste0(
  oneDriveOutput,
  "Figure2/Heatmap.pdf"
), width = 14, height = 14)
# INvestigate the pathways in NMP4
IRIgenes <- read.csv(paste0(oneDriveOutput,
  "Supplementary/IRIgenesCoefficient.csv")) %>% filter(IRI == 'Upregulated')

NMP4IRI <- cpmData_sub[IRIgenes %>% 
                         pull(gene) %>% 
                         unique,
                       str_detect(colnames(cpmData_sub),'NMP4')] %>%
  mutate(gene = rownames(.))
NMP4IRI <- NMP4IRI %>% pivot_longer(cols = c('NMP4_1', 'NMP4_3', 'NMP4_4'),
                                    values_to = 'counts', 
                                    names_to = 'sample') %>% 
  mutate(sample = factor(sample, levels = c('NMP4_1', 'NMP4_3', 'NMP4_4')))
  


ggplot(as.data.frame(NMP4IRI), 
       aes(x = as.integer(sample), counts, col = gene)) +
  geom_point() +
  geom_line() +
  stat_summary(
    fun.y = mean, fun.ymin = function(y) mean(y) - sd(y),
    fun.ymax = function(y) mean(y) + sd(y), color = "dark red",
    geom = "pointrange", show.legend = FALSE
  ) + theme_minimal() + scale_color_manual(values = rep('grey', 17)) +
  theme(legend.position = 'none') + xlab('Time of collection') +
  ylab('Gene expression Log2counts of IRI genes')

# Classification:
# Read the Raigani data:
files <- list.files("/Users/uqschauq/Documents/NMP/Data/RaiganiEtAl/",
  full.names = T
)[1:32]
countRaigan <- matrix(nrow = 60675)
for (i in files) {
  countRaigan <- cbind(countRaigan, read.table(i))
}
countRaigan <- countRaigan[, -1]
genes <- countRaigan[, seq(1, ncol(countRaigan), by = 2)]
countRaigan <- countRaigan[, seq(2, ncol(countRaigan), by = 2)]
rownames(countRaigan) <- genes[, 1]
colnames(countRaigan) <- str_remove(str_remove(
  files,
  "/Users/uqschauq/Documents/NMP/Data/RaiganiEtAl//"
), "_CountTable.txt.gz")
# Conditions:
conditions <- data.frame(
  colID = c(
    "GSM6124517_JB1", "GSM6124518_JB2", "GSM6124519_JB3", "GSM6124520_JB4",
    "GSM6124521_JB5", "GSM6124522_JB6", "GSM6124523_JB7", "GSM6124524_JB8",
    "GSM6124525_JB9", "GSM6124526_JB10", "GSM6124527_JB11",
    "GSM6124528_JB12", "GSM6124529_JB13", "GSM6124530_JB14",
    "GSM6124531_JB15", "GSM6124532_16", "GSM6124533_29", "GSM6124534_5",
    "GSM6124535_6", "GSM6124536_25", "GSM6124537_4", "GSM6124538_35",
    "GSM6124539_34", "GSM6124540_3", "GSM6124541_8", "GSM6124542_1",
    "GSM6124543_2", "GSM6124544_15", "GSM6124545_24", "GSM6124546_14",
    "GSM6124547_27", "GSM6124548_13"
  ),
  replicate = c(
    "replicate 1", "replicate 2", "replicate 3", "replicate 4",
    "replicate 5", "replicate 1", "replicate 3", "replicate 4",
    "replicate 1", "replicate 2", "replicate 3", "replicate 5",
    "replicate 2", "replicate 4", "replicate 5", "replicate 1",
    "replicate 1", "replicate 1", "replicate 2", "replicate 2",
    "replicate 2", "replicate 3", "replicate 3", "replicate 1",
    "replicate 1", "replicate 1", "replicate 2", "replicate 2",
    "replicate 2", "replicate 3", "replicate 3", "replicate 3"
  ),
  condition = c(
    "Emricasan", "Emricasan", "Emricasan", "Emricasan", "Emricasan",
    "Emricasan", "Emricasan", "Emricasan", "Emricasan", "Emricasan",
    "Emricasan", "Emricasan", "Emricasan", "Emricasan", "Emricasan",
    "Control", "Control", "Control", "Control", "Control", "Control",
    "Control", "Control", "Excluded Control", "Excluded Control",
    "Excluded Control", "Control", "Control", "Control", "Control",
    "Control", "Control"
  ),
  perf_success = c(
    "Adequate Hepatic Function Liver",
    "Adequate Hepatic Function Liver",
    "Inadequate Hepatic Function Liver",
    "Adequate Hepatic Function Liver",
    "Adequate Hepatic Function Liver",
    "Adequate Hepatic Function Liver",
    "Inadequate Hepatic Function Liver",
    "Adequate Hepatic Function Liver",
    "Adequate Hepatic Function Liver",
    "Adequate Hepatic Function Liver",
    "Inadequate Hepatic Function Liver",
    "Adequate Hepatic Function Liver",
    "Adequate Hepatic Function Liver",
    "Adequate Hepatic Function Liver",
    "Adequate Hepatic Function Liver",
    "Fatty Inadequate Hepatic Function Liver",
    "Fatty Inadequate Hepatic Function Liver",
    "Fatty Inadequate Hepatic Function Liver",
    "Fatty Inadequate Hepatic Function Liver",
    "Fatty Inadequate Hepatic Function Liver",
    "Fatty Inadequate Hepatic Function Liver",
    "Fatty Inadequate Hepatic Function Liver",
    "Fatty Inadequate Hepatic Function Liver",
    "Fatty Adequate Hepatic Function Liver",
    "Fatty Adequate Hepatic Function Liver",
    "Fatty Adequate Hepatic Function Liver",
    "Fatty Adequate Hepatic Function Liver",
    "Fatty Adequate Hepatic Function Liver",
    "Fatty Adequate Hepatic Function Liver",
    "Fatty Adequate Hepatic Function Liver",
    "Fatty Adequate Hepatic Function Liver",
    "Fatty Adequate Hepatic Function Liver"
  ),
  perf = c(
    "3 hours perfusion", "0 hours perfusion", "6 hours perfusion",
    "3 hours perfusion", "6 hours perfusion", "0 hours perfusion",
    "3 hours perfusion", "0 hours perfusion", "6 hours perfusion",
    "3 hours perfusion", "0 hours perfusion", "0 hours perfusion",
    "6 hours perfusion", "6 hours perfusion", "3 hours perfusion",
    "3 hours perfusion", "6 hours perfusion", "0 hours perfusion",
    "3 hours perfusion", "6 hours perfusion", "0 hours perfusion",
    "3 hours perfusion", "0 hours perfusion", "3 hours perfusion",
    "6 hours perfusion", "0 hours perfusion", "3 hours perfusion",
    "6 hours perfusion", "0 hours perfusion", "3 hours perfusion",
    "6 hours perfusion", "0 hours perfusion"
  )
)

conditions.f <- conditions %>% dplyr::filter(condition %in%
  c("Control", "Excluded Control"))
conditions.f$perf <- str_replace(
  conditions.f$perf,
  "0 hours perfusion", "CIT_end"
)
conditions.f$perf <- str_replace(
  conditions.f$perf,
  "3 hours perfusion", "3h_perf"
)
conditions.f$perf <- str_replace(
  conditions.f$perf,
  "6 hours perfusion", "6h_perf"
)
conditions.f$perf_success <- str_replace(
  conditions.f$perf_success,
  "Fatty Adequate Hepatic Function Liver", "Perf_Successful"
)
conditions.f$perf_success <- str_replace(
  conditions.f$perf_success,
  "Fatty Inadequate Hepatic Function Liver",
  "Perf_Unsuccessful"
)
conditions.f[which(conditions.f$colID %in% c(
  "GSM6124534_5",
  "GSM6124537_4", "GSM6124539_34",
  "GSM6124542_1", "GSM6124545_24",
  "GSM6124548_13"
)), ]$perf_success <- "CIT"
# Get the dataset:
countRaigan.f <- countRaigan[, conditions.f$colID]
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

# Discriminant analysis:
set.seed(5249)
# Get the signature:
successfulSignature <- unique(c(
  viableResults %>% pull(gene)
))


# Use the raigan data as a test set:
Xtest <- cpm(countRaigan.f[, 1:17]) %>%
  as.data.frame() %>%
  filter(rownames(.) %in% successfulSignature) %>%
  as.matrix()
Xtest <- t(Xtest[-which(rowVars(Xtest) == 0), ])
Xtest <- Xtest[, order(colnames(Xtest))]
Ytest <- factor((conditions.f %>% pull(perf_success)),
  levels = c("CIT", "Perf_Unsuccessful", "Perf_Successful")
)
# Extract our data to create a classifier:
Xtrain <- t(cpm(countData_sub)[
  which(rownames(countData_sub) %in% colnames(Xtest)),
  -which(colnames(countData_sub) %in% c("NMP4_3", "NMP4_4"))
])
Xtrain <- Xtrain[order(rownames(Xtrain)), order(colnames(Xtrain))]
# Remove the two misclassified samples:
Ytrain <- sampleinfo_sub %>% filter(!colID %in% c("NMP4_3", "NMP4_4"))
Ytrain <- Ytrain[order(Ytrain$colID), ]
Ytrain$perf_success <- factor(Ytrain$perf_success, levels = c(
  "CIT",
  "Perf_Unsuccessful", "Perf_Successful"
))

sampleToRemove <- Ytrain %>%
  group_by(sample) %>%
  mutate(n = n()) %>%
  filter(n == 1) %>%
  pull(colID)
Ytrain <- Ytrain %>% filter(!colID %in% sampleToRemove)
Xtrain <- Xtrain[Ytrain$colID, ]

# Add the misclassified sample to the test set:
# misclassifiedSample <- t(cpm(countData_sub)[
#   which(rownames(countData_sub) %in% colnames(Xtest)),
#   which(colnames(countData_sub) %in% c("NMP4_3", "NMP4_4"))
# ])
# Xtest <- rbind(Xtest, misclassifiedSample[
#   order(rownames(misclassifiedSample)), order(colnames(misclassifiedSample))
# ])
# Ytest[18:19] <- c("Perf_Unsuccessful", "Perf_Unsuccessful")

# Create the factor for the sample origin:
Ytrain$sample <- factor(Ytrain$sample)
design <- data.frame(sample = Ytrain$sample)
# Tune:
list.keepX <- c(1:10, seq(20, 100, 5))
# undergo the tuning process to determine the optimal number of variables
tune.splsda <- tune.splsda(
  X = Xtrain,
  Y = Ytrain$perf_success,
  ncomp = 3,
  validation = "Mfold",
  folds = 5, nrepeat = 50, # use repeated cross-validation
  dist = "max.dist", # use max.dist measure
  measure = "BER", # use balanced error rate of dist measure
  test.keepX = list.keepX,
  cpus = 4
)
# Extract parameters:
optimal.ncomp <- tune.splsda$choice.ncomp$ncomp
optimal.keepX <- tune.splsda$choice.keepX[1:optimal.ncomp]
# Create the final model:
final.splsda.multilevel <- splsda(Xtrain, Ytrain$perf_success,
  multilevel = design,
  ncomp = optimal.ncomp,
  keepX = optimal.keepX
)
# Calculate the background:
background <- background.predict(final.splsda.multilevel,
  comp.predicted = 2,
  dist = "mahalanobis.dist",
  xlim = c(-12.5, 12.5), ylim = c(-10, 10)
)
plotIndiv(final.splsda.multilevel,
  group = Ytrain$perf_success,
  ind.names = design$sample,
  pch = as.factor(design$sample), legend.title.pch = "Sample",
  legend = TRUE, legend.title = "Time",
  title = "Sample Plot of sPLS-DA", background = background
)
# Calculate the model stability:
perf.final.multilevel <- perf(final.splsda.multilevel,
  folds = 5, nrepeat = 1, validation = "loo",
  dist = "max.dist"
)
par(mfrow = c(1, 2))
plot(perf.final.multilevel$features$stable[[1]],
  type = "h", ylab = "Stability",
  xlab = "Features", main = "(a) Comp 1", las = 2
)
plot(perf.final.multilevel$features$stable[[2]],
  type = "h", ylab = "Stability",
  xlab = "Features", main = "(b) Comp 2", las = 2
)
# Predict:
YtestSample <- data.frame(sample = factor(conditions.f$replicate))
# YtestSample <- rbind(YtestSample, data.frame(sample = rep('nmp4', 2)))

predict.splsda <- predict(final.splsda.multilevel,
  newdata = Xtest,
  dist = "mahalanobis.dist",
  multilevel = YtestSample
)
# Prediction:
predict.comp2 <- predict.splsda$class$mahalanobis.dist[, 2]
table(factor(predict.comp2, levels = levels(Ytest)), Ytest)
# Try to plot the prediction:
# Plotting parameters:
# Get the PLS dataframe:
# Extract the PLS value for each of the sample
PLStrain <- as.data.frame(final.splsda.multilevel$variates$X)
PLStrain$set <- "train"
PLStrain$class <- final.splsda.multilevel$Y
PLStrain$sample <- rownames(PLStrain)
# Test samples:
PLStest <- as.data.frame(predict.splsda$variates)
colnames(PLStest) <- c("comp1", "comp2")
PLStest$sample <- rownames(PLStest)
PLStest$class <- Ytest
PLStest$set <- "test"

PLS <- rbind(PLStrain, PLStest)
nameChange <- as.vector(PLS$class)
nameChange[nameChange == "Perf_Successful"] <- "Viable Liver"
nameChange[nameChange == "Perf_Unsuccessful"] <- "Non-Viable Liver"
PLS$class <- factor(nameChange, levels = c(
  "CIT", "Non-Viable Liver",
  "Viable Liver"
))
predicted <- data.frame(
  sample = names(predict.comp2),
  predicted = predict.comp2
)
predicted$predicted[
  predicted$predicted == "Perf_Unsuccessful"
] <- "Non-Viable Liver"
predicted$predicted[
  predicted$predicted == "Perf_Successful"
] <- "Viable Liver"
PLS <- left_join(PLS, predicted, by = "sample")
PLS$predicted <- PLS$class == PLS$predicted
PLS$set[PLS$predicted == F] <- "Misclassified"
# PLS[str_detect(PLS$sample, 'NMP4'),]$set <- 'other'
# scatter plot:
col <- c("light blue", VolColorNonViable[2], VolColorViable[2])
PLS.f <- PLS %>% filter(set %in% c("Misclassified", "test", "other"))
ggplot(PLS.f) +
  geom_point(aes(x = comp1, y = comp2, col = class, shape = set, group = class),
    size = 4
  ) +
  theme_light() +
  # stat_ellipse(aes(x = comp1, y = comp2, col = class), data = PLS %>%
  #                filter(set == "test"), level = 0.95, lty = 2) +
  scale_color_manual(values = col, name = "Real Label:") +
  scale_shape_manual(values = c(3, 19, 18)) +
  theme(legend.position = "bottom") +
  xlab("Component 1 (55%)") +
  ylab("Component 2 (4%)") +
  geom_point(aes(x = Var1, y = Var2),
    data = as.data.frame(background$CIT)
  ) +
  geom_point(aes(x = Var1, y = Var2),
    data = as.data.frame(background$Perf_Successful)
  ) +
  geom_point(aes(x = Var1, y = Var2),
    data = as.data.frame(background$Perf_Unsuccessful)
  ) -> plsPlot
ggExtra::ggMarginal(plsPlot, type = "boxplot", groupColour = T) -> plsPlot

ggsave(plsPlot, filename = paste0(
  oneDriveOutput,
  "Figure2/PLS-DA_plot.pdf"
))

# Metrics extraction:
# Specificity
getMetrics <- function(predicted, true, label) {
  # Predicted
  predicted <- as.vector(predicted$class$mahalanobis.dist[, 2])
  predicted[which(predicted != label)] <- "other"
  # True
  true <- as.vector(true)
  true[which(true != label)] <- "other"
  # Predictions:
  TP <- sum(which(predicted == label) %in% which(true == label))
  TN <- sum(which(predicted == "other") %in% which(true == "other"))
  FP <- sum(which(predicted == "other") %in% which(true == label))
  FN <- sum(which(predicted == label) %in% which(true == "other"))
  # Metrics:
  SPC <- TN / (FP + TN)
  SENS <- TP / (TP + FN)
  ACC <- (TP + TN) / (TP + TN + FP + FN)
  PPV <- TP / (TP + FP)
  return(c(SPC, SENS, PPV, ACC))
}

auc_mannWhitney <- function(y, pred, class) {
  y <- as.vector(Ytest)
  pred <- as.vector(pred)
  # Make it a multiclass classifier:
  y[which(y != class)] <- "Other"
  pred[which(pred != class)] <- "Other"
  # Make y logical:
  table <- table(y)
  if (table[1] < table[2]) {
    y <- y != class
  } else {
    y <- y == class
  }
  # Calculate the AUC using a Mann-Whitney test:
  n1 <- sum(y)
  n2 <- sum(!y)
  R1 <- sum(rank(pred)[y])
  U1 <- R1 - n1 * (n1 + 1) / 2
  U1 / (n1 * n2)
}


aucCIT <- auc_mannWhitney(Ytest, predict.comp2, "CIT")
aucSuccess <- auc_mannWhitney(Ytest, predict.comp2, "Perf_Successful")
aucFail <- auc_mannWhitney(Ytest, predict.comp2, "Perf_Unsuccessful")


metrics <- data.frame(
  value = c(
    getMetrics(predict.splsda, Ytest, "CIT"),
    getMetrics(predict.splsda, Ytest, "Perf_Unsuccessful"),
    getMetrics(predict.splsda, Ytest, "Perf_Successful")
  ),
  metrics = rep(c("Specificity", "Sensitivity", "Precision", "Accuracy"), 3),
  condition = c(
    rep("CIT", 4), rep("Non-Viable Liver", 4),
    rep("Viable Liver", 4)
  )
)

# Test:

library(PRROC)
temp <- as.data.frame(predict.splsda$predict)
temp <- temp %>%
  mutate(
    CIT = CIT.dim1 + CIT.dim2,
    Perf_Unsuccessful = Perf_Unsuccessful.dim1 + Perf_Unsuccessful.dim2,
    Perf_Successful = Perf_Successful.dim1 + Perf_Successful.dim2
  ) %>%
  dplyr::select(CIT, Perf_Successful, Perf_Unsuccessful)

# Supplementary table 5 and 6:
# Basically, I wnat to include both the loading and how conserved the genes are:
# Final Loadings
loadings <- as.data.frame(final.splsda.multilevel$loadings$X)
loadings <- loadings %>%
  filter(comp1 != 0 | comp2 != 0) %>%
  mutate(Gene = rownames(.)) %>%
  pivot_longer(
    cols = c(comp1, comp2), values_to = "Loading",
    names_to = "Component"
  ) %>%
  filter(Loading != 0)
# Merge with the gene stability:
loadings <- left_join(loadings, data.frame(
  stability = c(
    perf.final.multilevel$features$stable[[1]],
    perf.final.multilevel$features$stable[[2]]
  ),
  Gene = c(
    names(perf.final.multilevel$features$stable[[1]]),
    names(perf.final.multilevel$features$stable[[2]])
  )
),
by = "Gene"
)

loadings$Gene <- factor(loadings$Gene,
  levels = loadings %>% group_by(Component) %>%
    arrange(Loading) %>% pull(Gene)
)


stability <- data.frame(
  Gene = c(
    names(perf.final.multilevel$features$stable[[1]]),
    names(perf.final.multilevel$features$stable[[2]])
  ),
  stability = c(
    perf.final.multilevel$features$stable[[1]],
    perf.final.multilevel$features$stable[[2]]
  ),
  component = c(
    rep(
      "Component1",
      length(perf.final.multilevel$features$stable[[1]])
    ),
    rep(
      "Component2",
      length(perf.final.multilevel$features$stable[[2]])
    )
  )
)

loadings <- left_join(loadings, stability, by = "Gene")
geneOrder <- loadings %>%
  group_by(component) %>%
  arrange(Loading) %>%
  pull(Gene)
loadings$Gene <- factor(loadings$Gene, levels = geneOrder)

# Boxplot:
ggplot(loadings) +
  geom_bar(aes(x = Loading, y = Gene, fill = Component), stat = "identity") +
  xlab("PLS-DA Gene Weight") +
  ylab("") +
  theme_minimal() +
  facet_wrap(~component, scales = "free_y") +
  scale_fill_manual(
    labels = c("Component 1", "Component 2"),
    values = met.brewer(n = 2, name = "Johnson"), name = ""
  ) +
  theme(legend.position = "top") -> LoadingPlot

ggsave(LoadingPlot,
  filename = paste0(
    oneDriveOutput,
    "Figure2/loadingPlot.pdf"
  ),
  width = 10, height = 10
)

write_csv(loadings,
  file = paste0(
    oneDriveOutput,
    "Supplementary/PLS_DA_loadings_stability.csv"
  ))

