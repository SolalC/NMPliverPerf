library(tidyverse)
library(MetBrewer)
library(data.table)
library(variancePartition)
library(mixOmics)
library(rlist)
library(plotROC)
library(cowplot)
# set the parralelization:
# param <- SnowParam(40, "SOCK", progressbar = TRUE)
# Read the data:
# nCores <- 40
mergedCount <- read.csv('/QRISdata/Q1144/Data/merged/mergedCountData.csv', row.names = 1)
colnames(mergedCount) <- str_remove(colnames(mergedCount), 'X')
colnames(mergedCount) <- str_replace(colnames(mergedCount), '\\.', '-')
# Read the sample information:
mergedInfo <- read.csv('/QRISdata/Q1144/Data/merged/mergedSampleInfo.csv', row.names = 1)
timelevels <- c('CIT', '1H', '3H', '6H', '12H', '24H', 'POST')
mergedInfo$time <- factor(mergedInfo$time, levels = timelevels)
# Remove a single gene that makes sPLSda fail:
mergedCount <- mergedCount[rownames(mergedCount) != 'MIR612',]
all(colnames(mergedCount) == mergedInfo$colID)
# Based on the common genes:
DEGallbatch <- readRDS('/scratch/project_mnt/S0007/solal/NMP/Results/DEG/viabilityDEG.rds')
viabilityDEG <- DEGallbatch[[1]]
# Get the DEG:
nGenes <- nrow(viabilityDEG)
viableDEG <- topTable(viabilityDEG, coef = 'viabilityViable', number = nGenes)
nonViableDEG <- topTable(viabilityDEG, coef = 'viabilityNon-Viable', number = nGenes)
# Get significant genes:
viableDEG.f <- viableDEG %>% filter(adj.P.Val < 0.05, !between(logFC, -1, 1)) %>% 
  mutate(gene = rownames(.))
nonViableDEG.f <- nonViableDEG %>% filter(adj.P.Val < 0.05, !between(logFC, -1, 1)) %>%
  mutate(gene = rownames(.))
# Define the sets of genes :
shared <-  intersect(viableDEG.f$gene, nonViableDEG.f$gene)
viableOnly <- setdiff(viableDEG.f$gene, nonViableDEG.f$gene)
nonViableOnly <- setdiff(nonViableDEG.f$gene, viableDEG.f$gene)
# Get the contrast genes:
contrastDEG <- DEGallbatch[[2]]
# Extract the value for all the genes:
constrastDEG <- topTable(contrastDEG, coef = 'viabilityContrast', number = nGenes)
# Get significant genes:
constrastDEG.f <- constrastDEG %>% filter(adj.P.Val < 0.05, !between(logFC, -1, 1)) %>% 
  mutate(gene = rownames(.))
# Genes:
# previous:
inHouseDEG <- fread('/scratch/project_mnt/S0007/solal/NMP/Data/PreviousDEG/CIT_viable_significant_DEG.csv') %>% pull(gene)
# allDEG:

DEGall <- unique(c(viableDEG.f$gene, nonViableDEG.f$gene, constrastDEG.f$gene))
# Filter based on hour:

# PCA function:
# Heatmap function:
PCAClassification <- function(count, colData, genes, filter =F, seed=3424, scale=F) {
  set.seed(seed)
  # Filter based on batch:
  if(is.character(filter)){
    colData <- colData %>% filter(batch %in% filter)
    count <- count[, colData$colID]
    print(all(colnames(count) == colData$colID))
  }
  if(all(colnames(count) == colData$colID) == F){
    print('Datasets in different order')
    break
  }
  # PCA:
  pcaCount <- count[c(rownames(count) %in% genes), ]
  pcaCount <- pcaCount[(rowSums(pcaCount) != 0),]
  pcaLog <- log2(count+1)
  pca <- prcomp(t(pcaLog), scale = scale)
  pcVar <- summary(pca)$importance[2,][1:5]
  df <- data.frame(PC1 = pca$x[,1], 
                   PC2 = pca$x[,2],
                   PC3 = pca$x[,3],
                   PC4 = pca$x[,4],
                   PC5 = pca$x[,5],
                   colID = names(pca$x[,1]))
  df <- df %>% left_join(., colData, by = 'colID')
  
  return(list(df, pcVar))
}
# PCA:
# # Correcat for batch effect:
# library(sva)
# # Set the Surrogate variable:
# mod1 <- model.matrix(~ viability, data = mergedInfoFiltered)
# # Null model matrix contains only the adjustment variables
# mod0 <- model.matrix(~ 1, data = mergedInfoFiltered)
# # filter low number genes:
# mergedCountFilterMatrix <- as.matrix(mergedCountFilter)
# svseq <- svaseq(mergedCountFilterMatrix, mod1, mod0)
# colnames(svseq$sv) <- paste0("SV", seq(1, ncol(svseq$sv)))
# # Add the SVA to the dataframe:
# mergedInfoFiltered <- cbind(mergedInfoFiltered, svseq$sv)
# # Create the experimental design (random effect, all known covariates
# # and surrogate variables:)
# covar_matrix <- model.matrix(~ SV1 + SV2 + SV3, data = mergedInfoFiltered)
# # Check that the columns and rows are equivalent:
# all(colnames(mergedCountFilter) == mergedInfoFiltered$colID)
# countBatchCorrected <- ComBat_seq(counts = as.matrix(mergedCountFilter), 
#                                   batch = (as.factor(mergedInfoFiltered$batch)),
#                                   group = (as.factor(mergedInfoFiltered$viability)),
#                                   covar_mod = covar_matrix)
# countBatchCorrected <- as.data.frame(countBatchCorrected)
# saveRDS(countBatchCorrected, 
#         '/scratch/project_mnt/S0007/solal/NMP/Results/Classification/countBatchCorrected.rds')
# countBatchCorrected <- readRDS('/scratch/project_mnt/S0007/solal/NMP/Results/Classification/countBatchCorrected.rds')
# # PCA on batch corrected counts:
# pcaAll <- PCAClassification(as.data.frame(countBatchCorrected), 
#                             mergedInfoFiltered, 
#                             DEGall, 
#                             scale = F,
#                             filter = 'raigani')
# # PC1/2
# (ggplot(pcaAll[[1]], aes(x = PC1, PC2, col = viability, shape = batch)) +
#   geom_point() +
#   xlab(paste0('PC1 (', round(pcaAll[[2]][1]*100, 2), '%)')) +
#   ylab(paste0('PC2 (', round(pcaAll[[2]][2]*100, 2), '%)')) -> pcPlot)
# # PC2/3
# ggplot(pcaAll[[1]], aes(x = PC3, PC4, col = viability, shape = batch)) +
#   geom_point() +
#   xlab(paste0('PC3 (', round(pcaAll[[2]][3]*100, 2), '%)')) +
#   ylab(paste0('PC4 (', round(pcaAll[[2]][4]*100, 2), '%)')) 
# ggsave(pcPlot, filename = '/scratch/project_mnt/S0007/solal/NMP/Results/Classification/PCA.png')
# 
# h4 <- heatmapClassification(countBatchCorrected, mergedInfoFiltered, DEGall)
# ggsave(h4, filename = '/scratch/project_mnt/S0007/solal/NMP/Results/Classification/DEGall/heatmapAllDF.png')
# sPLS-DA:

# TODO: apply their classifier to our samples
# TODO: Train a new classifier, apply it to our data.
# TODO: show the evolution of the score over time (even for samples that are POST the 3/6 hours)
# Leave-One-Dataset-Out Cross-Validation
# Normalise the data:
# Get the count corrected data:
# countBatchCorrected <- readRDS('/scratch/project_mnt/S0007/solal/NMP/Results/Classification/countBatchCorrected.rds')
# mergedInfo %>% filter(time %in% c('CIT', '3H', '6H')) -> mergedInfoFiltered
# countBatchCorrected[(rowSums(countBatchCorrected) > 10), 
#                     mergedInfoFiltered$colID] -> mergedCountFilter
mergedCountFilter <- mergedCountFilter[rownames(mergedCountFilter) != 'MIR612',]
all(colnames(mergedCountFilter) == mergedInfoFiltered$colID)
# Log2 transformation:
mergedLogfilter <- log2(mergedCountFilter+1)
# Split by data.frames to train the model:
# inHouse:
Y1 <- mergedInfoFiltered %>% filter(batch == 'inHouse') %>%
  group_by(sampleID) %>% filter(n() > 1) 
X1 <- mergedLogfilter[rownames(mergedLogfilter) %in% DEGall, Y1 %>% pull(colID)]

# Hautz:
Y2 <- mergedInfoFiltered %>% filter(batch == 'hautz') %>%
  group_by(sampleID) %>% filter(n() > 1) 
X2 <- mergedLogfilter[rownames(mergedLogfilter) %in% DEGall, Y2 %>% pull(colID)]
# Raigani:
Y3 <- mergedInfoFiltered %>% filter(batch == 'raigani')  %>%
  group_by(sampleID) %>% filter(n() > 1) 
X3 <- mergedLogfilter[rownames(mergedLogfilter) %in% DEGall, Y3 %>% pull(colID)]
# Use the logistic model regression that they used in the Hautz paper:
# Here are the weights:
hautzLogReg <- data.frame(
  gene = c('Intercept', 'FAM43A', 'PKD2L1', 'CD274', 'PLD6', 'MLPH', 'DBP', 'WDR72'),
  coef = c(-0.385, 0.593, 0.917, 0.081, -0.619, -0.437, 0.415, 0.592)
  )
# Get the data for the Hautz model:
hautzLogRegProba <- mergedCount[hautzLogReg$gene[-1],]
hautzLogRegProba <- log2(t(apply(hautzLogRegProba, 1, function(x)scale(x)))+1)
colnames(hautzLogRegProba) <- colnames(mergedCount)
rownames(hautzLogRegProba) <- hautzLogReg$gene[-1]
# Create a vector that we will use to multiply the genes.
coef_vector <- setNames(hautzLogReg$coef, hautzLogReg$gene)
# Calculate the result:
(hautzLogRegProba <- coef_vector['Intercept'] + 
    colSums(hautzLogRegProba * 
              coef_vector[rownames(hautzLogRegProba)]))

hautzLogRegProbaDF <- data.frame(
  colID = names(hautzLogRegProba), 
  logit = hautzLogRegProba) %>% mutate(pred = ifelse(logit > 0, 
                                              'Non-Viable',
                                              'Viable')) %>% 
  left_join(., (mergedInfoFiltered %>% dplyr::select(colID, batch, time, viability)), 
            by = 'colID')%>% filter(batch %in% c('inHouse', 'raigani'),
                                    time %in% c('6H'))
library(pROC)
# Generate ROC curve
roc_curve <- roc(hautzLogRegProbaDF$viability ~ hautzLogRegProbaDF$logit)
# Plot the curve
pdf("/scratch/project_mnt/S0007/solal/NMP/Results/Classification/rocCurveHautzClassifier.pdf",
    width = 5, height = 5)
plot(roc_curve, 
     main = "ROC Curve",
     col = "blue",
     lwd = 2)
auc_value <- auc(roc_curve)
legend("bottomright", 
       legend = paste("AUC =", round(auc_value, 3)),
       bty = "n")
dev.off()
# Let's ignore the regression and show the expression. That could also just work.
getlogFC <- function(logCount, info, genes, study, order){
  logCount <- logCount[genes,]
  logExp3hoursViable <- logCount[,info %>% 
                             filter(time == '3H', batch %in% study, 
                                    viability == 'Viable') %>%
                             pull(colID)]
  logExp3hoursViable <- rowMeans(logExp3hoursViable)
  
  logExp6hoursViable <- logCount[,info %>% 
                             filter(time == '6H', batch %in% study,
                                    viability == 'Viable') %>%
                             pull(colID)]
  logExp6hoursViable <- rowMeans(logExp6hoursViable)
  # Non-viable:
  logExp3hoursNonViable <- logCount[,info %>% 
                                   filter(time == '3H', batch %in% study,
                                          viability == 'Non-Viable') %>%
                                   pull(colID)]
  logExp3hoursNonViable <- rowMeans(logExp3hoursNonViable)
  
  logExp6hoursNonViable <- logCount[,info %>% 
                                   filter(time == '6H', batch %in% study,
                                          viability == 'Non-Viable') %>%
                                   pull(colID)]
  logExp6hoursNonViable <- rowMeans(logExp6hoursNonViable)
  
  
  logFCdf <- data.frame(
    NTP_TP_3H = logExp3hoursViable - logExp3hoursNonViable,
    NTP_TP_6H = logExp6hoursViable - logExp6hoursNonViable
  ) %>% mutate(gene = rownames(.)) %>% arrange(match(gene, order))
  
  return(logFCdf)
}
heatmapOrdering <- c('CD274', 'FAM43A', 'PKD2L1', 'DBP', 'PLD6', 'MLPH', 'WDR72')
logFCinHouse <- getlogFC(mergedLogfilter, mergedInfo, hautzLogReg$gene[-1], 
                         'inHouse', heatmapOrdering)
heatmapColor <- c("#0000FF", "#000099", "#000066", "#000033", 
                  "#000000", "#330000", "#660000", "#990000", 
                  "#FF0000")
htmp1 <- ggplotify::as.ggplot(pheatmap(logFCinHouse[,-3],
                                       breaks = seq(-2, 2, by = 0.5 ),
                                       legend_breaks = seq(-2, 2, by = 0.5 ),
                                       color = heatmapColor,
                                       cluster_rows = F, cluster_cols = F,
                                       main = 'inHouse dataset'))
logFCraigani <- getlogFC(mergedLogfilter, mergedInfo, hautzLogReg$gene[-1], 
                         'raigani', heatmapOrdering)
htmp2 <- ggplotify::as.ggplot(pheatmap(logFCraigani[,-3],
                                       breaks = seq(-2, 2, by = 0.5 ),
                                       legend_breaks = seq(-2, 2, by = 0.5 ),
                                       color = heatmapColor,
                                       cluster_rows = F, cluster_cols = F,
                                       main = 'raigani dataset'))

logFCboth <- getlogFC(mergedLogfilter, mergedInfo, hautzLogReg$gene[-1], 
                         c('inHouse', 'raigani'), heatmapOrdering)
htmp3 <- ggplotify::as.ggplot(pheatmap(logFCboth[,-3],
                                       breaks = seq(-2, 2, by = 0.5 ),
                                       legend_breaks = seq(-2, 2, by = 0.5 ),
                                       color = heatmapColor,
                                       cluster_rows = F, cluster_cols = F,
                                       main = 'both dataset'))
htmp <- plot_grid(htmp1, htmp2, htmp3, nrow = 1)
ggsave(htmp, filename = '/scratch/project_mnt/S0007/solal/NMP/Results/Classification/heatmapHautzClassifier.pdf', width = 8, height = 3)
# Classification:
X <- list(X1, X2, X3)
Y <- list(Y1, Y2, Y3)
list.keepX <- c(1:10, seq(20, 200, 5))
finalOutput <- list()
for(i in 1:3) {
  # Training:
  X_train <- t(do.call(cbind, X[-i]))
  Y_train <- do.call(rbind, Y[-i]) %>% 
    mutate(viability = factor(viability, c('CIT', 'Non-Viable', 'Viable')),
           sampleID = as.factor(sampleID))
  if(!all(rownames(X_train) == Y_train$colID)){
    message('Train dataset in different order')
    break
    }
  # Test:
  X_test <- t(X[[i]])
  Y_test <- Y[[i]] %>%
    mutate(viability = factor(viability, c('CIT', 'Non-Viable', 'Viable')))
  
  if(!all(rownames(X_test) == Y_test$colID)){
    message('Test dataset in different order')
    break
  }
  message("Start training of model: ", i)
  splsdaTuned <- tune.splsda(
    X = X_train,
    Y = Y_train$viability,
    ncomp = 3,
    validation = "Mfold",
    folds = 10, nrepeat = 100, # use repeated cross-validation
    dist = "max.dist", # use max.dist measure
    measure = "BER", # use balanced error rate of dist measure
    test.keepX = list.keepX
  )
  # Extract optimal parameters:
  optimal_ncomp <- 2
  optimal_keepX <- splsdaTuned$choice.keepX[1:optimal_ncomp]
  message('Optimal keepX:', optimal_keepX, '\n Optimal ncomp:', optimal_ncomp)
  #design:
  design <- data.frame(sample = factor(Y_train$sampleID))
  # Multilevel splsDA:
  final.splsda.multilevel <- splsda(X = X_train,
                                    Y = Y_train$viability,
                                    multilevel = design,
                                    ncomp = optimal_ncomp,
                                    keepX = optimal_keepX)
  # Calculate the model performance:
  perf.final.multilevel <- perf(final.splsda.multilevel,
                                folds = 5, nrepeat = 1, 
                                validation = "loo",
                                dist = "max.dist")
  
  predictions <- predict(final.splsda.multilevel,
                            newdata = X_test,
                            dist = "mahalanobis.dist",
                            multilevel = Y_test$sampleID)
  
  # Confusion matrix and accuracy
  confusion_matrix <- table(Predicted = predictions$class$mahalanobis.dist[, 2], 
                            Actual = Y_test$viability)
  iterationOutput <- list(final.splsda.multilevel, perf.final.multilevel,
                          predictions, confusion_matrix)
  names(iterationOutput) <- c('splsda.model', 'splsda.internal.perf',
                              'prediction.object', 'confusion_matrix')
  
  finalOutput <- list.append(finalOutput, iterationOutput)
}

saveRDS(finalOutput, '/scratch/project_mnt/S0007/solal/NMP/Results/Classification/splsdaOutput2Component.rds')
finalOutput <- readRDS('/scratch/project_mnt/S0007/solal/NMP/Results/Classification/splsdaOutput2Component.rds')
# Add iteration here:
perf <- data.frame()
metrics <- data.frame()
for(i in 1:3) {
  nComp <- finalOutput[[i]]$splsda.model$ncomp
  levels <- c("CIT", "Non-Viable", "Viable")
  predicted_value <- factor(finalOutput[[i]]$prediction.object$class$mahalanobis.dist[, nComp], levels)
  true_labels <- factor(Y[[i]]$viability, levels)
  
  for(positive_class in levels) {
    predDf <- data.frame(predicted_value, true_labels)
    predDf$binary <- as.integer(predDf$true_labels == positive_class)
    predDf$pred <- finalOutput[[i]]$prediction.object$predict[,,paste0('dim', nComp)][, positive_class]
    predDf$class <- positive_class
    predDf$batch <- unique(Y[[i]]$batch)
    perf <- rbind(perf, predDf)
    # Calculate confusion matrix
    confusion_matrix <- table(Predicted = predDf$predicted_value == positive_class, 
                              Actual = predDf$true_labels == positive_class)
    
    # Extract TP, FP, TN, FN
    TP <- confusion_matrix["TRUE", "TRUE"]
    FP <- confusion_matrix["TRUE", "FALSE"]
    TN <- confusion_matrix["FALSE", "FALSE"]
    FN <- confusion_matrix["FALSE", "TRUE"]
    
    # Metrics calculations
    sensitivity <- TP / (TP + FN)
    specificity <- TN / (TN + FP)
    precision <- TP / (TP + FP)
    NPV <- TN / (TN + FN)
    # Calculate AUC using the pROC package
    roc_obj <- roc(predDf$binary, predDf$pred)
    auc_value <- auc(roc_obj)
    # Append metrics to predDf
    roc_obj <- roc(predDf$binary, predDf$pred)
    auc_value <- auc(roc_obj)
    
    tempPerf <- data.frame(
      sensitivity = TP / (TP + FN),
      specificity = TN / (TN + FP),
      precision = TP / (TP + FP),
      NPV = TN / (TN + FN),
      batch = unique(Y[[i]]$batch),
      class = positive_class,
      AUC =  auc_value
    )
    metrics <- rbind(metrics, tempPerf)
  }
}
fwrite(metrics, '/scratch/project_mnt/S0007/solal/NMP/Results/Classification/LSOmetrics.csv')
# Plot the ROC curve
rocPlot <- ggplot(perf, aes(d = binary, m = pred, col = class)) +
  geom_roc(labels = FALSE, n.cuts = 0) +
  facet_wrap(~batch) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  style_roc(xlab = "1 - Specificity", ylab = 'Sensitivity') +
  scale_color_manual(values = c('#B0DAE7', '#17154F', '#F5BB50'))
ggsave(rocPlot, filename = '/scratch/project_mnt/S0007/solal/NMP/Results/Classification/LSOrocPlot.pdf', width = 8, height = 3)

# Extract the PLS values:
# pls <- data.frame()
# for(i in 1:3) {
#   as.data.frame(finalOutput[[i]]$prediction.object$variates) %>% 
#     mutate(viability = Y[[i]]$viability,
#            colID = rownames(.)) %>% 
#     pivot_longer(cols = c('dim1', 'dim2'), names_to = 'Comp', values_to = 'value') %>% 
#     mutate(batch = unique(Y[[i]]$batch)) -> plsTemp
#   pls <- rbind(pls, plsTemp)
# }
# 
# # inverse the value for the dim2 using the raigani as the test set:
# pls.f <- pls %>% mutate(value = ifelse((batch == 'raigani') & (Comp == 'dim2'), -value, value))
# # Plot:
# ggplot(pls.f, aes(x = batch, y = value, col = viability)) +
#   geom_boxplot(outliers = F) + 
#   geom_point(position = position_jitterdodge()) + 
#   facet_wrap(~Comp) +
#   scale_color_manual(values = c('#B0DAE7', '#17154F', '#F5BB50')) +
#   theme_minimal() -> dimPlot
# # Extract the stable feature and their effects:
# # Get the stability per model:
# extract_stable_features <- function(finalOutput, i) {
#   # Initialize an empty dataframe to store results
#   result_df <- data.frame(
#     list_name = character(),
#     feature_name = character(),
#     feature_value = numeric(),
#     stringsAsFactors = FALSE
#   )
#   # Iterate through each element of the list
#   for (list_index in seq_along(finalOutput[[i]]$splsda.internal.perf$features$stable)) {
#     # Get the current list element
#     current_list <- finalOutput[[i]]$splsda.internal.perf$features$stable[[list_index]]
#     # Convert the named vector to a dataframe
#     temp_df <- data.frame(
#       list_name = names(finalOutput[[i]]$splsda.internal.perf$features$stable)[list_index],
#       feature_value = current_list
#     )
#     # Bind to the result dataframe
#     result_df <- rbind(result_df, temp_df)
#   }
#   colnames(result_df) <- c('component', 'feature_name', 'feature_stability')
#   # Extract the loading:
#   loading <- as.data.frame(finalOutput[[i]]$splsda.model$loadings$X) %>% mutate(gene = rownames(.))
#   # Extract the different loadings:
#   LoadingComp1 <- loading %>% filter(gene %in% (result_df %>% 
#                                                filter(component == 'comp1') %>% 
#                                                pull(feature_name))) %>% 
#     mutate(component = 'comp1') %>%
#     dplyr::select(gene, comp1, component)
#   colnames(LoadingComp1) <- c('gene', 'loading', 'component')
#   # Comp2
#   LoadingComp2 <- loading %>% filter(gene %in% (result_df %>% 
#                                                filter(component == 'comp2') %>% 
#                                                pull(feature_name))) %>%
#     mutate(component = 'comp2') %>%
#     dplyr::select(gene, comp2, component)
#   colnames(LoadingComp2) <- c('gene', 'loading', 'component')
#   # Comp 3
#   # LoadingComp3 <- loading %>% filter(gene %in% (result_df %>% 
#   #                                              filter(component == 'comp3') %>% 
#   #                                              pull(feature_name))) %>%
#   #   mutate(component = 'comp3') %>%
#   #   dplyr::select(gene, comp3, component)
#   # colnames(LoadingComp3) <- c('gene', 'loading', 'component')
#   
#   loadingComp <- rbind(LoadingComp1, LoadingComp2)
#                        # LoadingComp3)
#   # Combine the results:
#   result_df <- left_join(result_df, loadingComp, by = c('feature_name' = 'gene', 'component'))
#   return(result_df)
# }
# stab1 <- extract_stable_features(finalOutput, 1) %>% mutate(batch = unique(Y[[1]]$batch))
# stab2 <- extract_stable_features(finalOutput, 2) %>% mutate(batch = unique(Y[[2]]$batch))
# stab3 <- extract_stable_features(finalOutput, 3) %>% mutate(batch = unique(Y[[3]]$batch))
# stab <- rbind(stab1, stab2, stab3)
# # something is wrong, I need to check the ITGAD gene:
# stab.f <- stab %>% filter(feature_stability > 0.8)
# # upset plot:
# create_upset_plot <- function(data, selected_component) {
#   # Filter data for the selected component
#   
#   presence_matrix <- data %>%
#     filter(component == selected_component) %>%
#     dplyr::select(feature_name, batch) %>%
#     mutate(present = 1) %>%
#     pivot_wider(names_from = batch,
#                 values_from = present,
#                 values_fill = 0) %>%
#     column_to_rownames("feature_name")
#   
#   # Create the upset plot
#   upset(as.data.frame(presence_matrix),
#         sets = colnames(presence_matrix),
#         order.by = "freq")
# }
# library(UpSetR)
# comp1FeatureStabPlot <- ggplotify::as.ggplot(create_upset_plot(stab.f, 'comp1'))
# comp2FeatureStabPlot <- ggplotify::as.ggplot(create_upset_plot(stab.f, 'comp2'))
# # comp3FeatureStabPlot <- ggplotify::as.ggplot(create_upset_plot(stab.f, 'comp3'))
# 
# featureStabPlot <- plot_grid(comp1FeatureStabPlot, comp2FeatureStabPlot, ncol = 1)
# # Feature loading:
# ggplot(stab.f, aes(x = loading, y = feature_name, col = batch)) +
#   geom_point() +
#   facet_wrap(~component, ncol = 1, scale = 'free_y') +  
#   theme_minimal() +
#   geom_vline(xintercept = 0, lty = 2) -> featureLoadingPlot
# # Correlation heatmap for all the genes idenfied:
# mergedLogfilter %>% filter(rownames(.) %in% stab.f$feature_name) -> mergedLogfilterStab
# # Calculate the correlation matrix:
# cor_matrix <- cor(t(mergedLogfilterStab), method = "pearson")
# 
# stab.ff <- stab.f %>% filter(!duplicated(feature_name))
# # Create the annotations:
# genes_in_cor <- rownames(cor_matrix)
# 
# annot_df <- stab.ff %>%
#   dplyr::select(feature_name, component) %>%
#   filter(feature_name %in% genes_in_cor) %>%
#   arrange(match(feature_name, genes_in_cor)) %>%
#   column_to_rownames("feature_name")
# 
# annot_col <- list(component = c(comp1 = "#E41A1C", comp2 = "#377EB8", comp3 = "#4DAF4A"))
# 
# pheatmap(cor_matrix,
#          annotation_row = annot_df,
#          annotation_col = annot_df,
#          annotation_colors = annot_col,
#          show_rownames = TRUE,
#          show_colnames = TRUE,
#          main = "Gene Expression Correlation Matrix with Component Annotation",
#          clustering_distance_rows = "correlation",
#          clustering_distance_cols = "correlation",
#          fontsize_row = 8,
#          fontsize_col = 8)
# Leave-One-Dataset-Out Cross-Validation
Xfull <- t(do.call(cbind, X))
Yfull <- do.call(rbind, Y) %>% mutate(viability = factor(viability, c('CIT', 'Non-Viable', 'Viable')),
                                     sampleID = as.factor(sampleID))
all(rownames(Xfull) == Yfull$colID)
list.keepX <- c(1:10, seq(20, 200, 5))

set.seed('65864')
design <- data.frame(sample = factor(Yfull$sampleID))
splsdaTuned <- tune.splsda(
  X = Xfull,
  Y = Yfull$viability,
  ncomp = 2,
  validation = "loo",
  dist = "max.dist", # use max.dist measure
  measure = "BER", # use balanced error rate of dist measure
  test.keepX = list.keepX,
  multilevel = design
)
# Extract optimal parameters:
optimal_ncomp <- 2
optimal_keepX <- splsdaTuned$choice.keepX[1:optimal_ncomp]
#design:
# Multilevel splsDA:
final.splsda.multilevel <- splsda(X = Xfull,
                                  Y = Yfull$viability,
                                  multilevel = design,
                                  ncomp = optimal_ncomp,
                                  keepX = optimal_keepX)
# Calculate the model performance:
perf.final.multilevel <- perf(final.splsda.multilevel,
                              nrepeat = 100, 
                              validation = "loo",
                              dist = "max.dist",
                              multilevel = design)

individualsPlot <- plotIndiv(final.splsda.multilevel, comp = c(1,2), 
          group = Yfull$viability, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE,
          title = 'sPLS-DA') 
individualsPlot <- individualsPlot$graph +
  scale_color_manual(values = c('#B0DAE7', '#17154F', '#F5BB50')) +
  theme(legend.position = 'bottom') + 
  xlab('Component1: 44%') + 
  ylab('Component 2: 4%')
# AUROC plot for the loo approach:
AUCplot <- auroc(final.splsda.multilevel, roc.comp = 2, print = F) 
AUCplot <- AUCplot$graph.Comp2 + 
    scale_color_manual(values = c('#B0DAE7', '#17154F', '#F5BB50')) +
  theme_minimal() +
  theme(legend.position = 'bottom')
AUCplotFull <- plot_grid(AUCplot, individualsPlot)
ggsave(AUCplotFull, filename = '/scratch/project_mnt/S0007/solal/NMP/Results/Classification/AUCplotLOO.pdf', width = 8, height = 4)
# Get metrics of interest for the final classifier:
nComp <- 2
levels <- c("CIT", "Non-Viable", "Viable")
predicted_value <- factor(as.data.frame(perf.final.multilevel$class$max.dist)[,2], levels)
true_labels <- factor(Yfull$viability, levels)
metrics <- data.frame()

for(positive_class in levels) {
  predDf <- data.frame(predicted_value, true_labels)
  predDf$binary <- as.integer(predDf$true_labels == positive_class)
  predDf$pred <- as.integer(predDf$predicted_value == positive_class)
  predDf$class <- positive_class
  # Calculate confusion matrix
  confusion_matrix <- table(Predicted = predDf$predicted_value == positive_class, 
                            Actual = predDf$true_labels == positive_class)
  
  # Extract TP, FP, TN, FN
  TP <- confusion_matrix["TRUE", "TRUE"]
  FP <- confusion_matrix["TRUE", "FALSE"]
  TN <- confusion_matrix["FALSE", "FALSE"]
  FN <- confusion_matrix["FALSE", "TRUE"]
  
  # Metrics calculations
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  precision <- TP / (TP + FP)
  NPV <- TN / (TN + FN)
  # Calculate AUC using the pROC package
  tempPerf <- data.frame(
    sensitivity = TP / (TP + FN),
    specificity = TN / (TN + FP),
    precision = TP / (TP + FP),
    NPV = TN / (TN + FN),
    class = positive_class
  )
  metrics <- rbind(metrics, tempPerf)
}
fwrite(metrics, '/scratch/project_mnt/S0007/solal/NMP/Results/Classification/LOOmetricsFull.csv')
# Extract the feature:
# Initialize an empty dataframe to store results
looDF <- data.frame(
  list_name = character(),
  feature_name = character(),
  feature_value = numeric(),
  stringsAsFactors = FALSE
)
# Iterate through each element of the list
for (list_index in seq_along(perf.final.multilevel$features$stable)) {
  # Get the current list element
  current_list <- perf.final.multilevel$features$stable[[list_index]]
  # Convert the named vector to a dataframe
  temp_df <- data.frame(
    list_name = names(perf.final.multilevel$features$stable)[list_index],
    feature_value = current_list
  )
  # Bind to the result dataframe
  looDF <- rbind(looDF, temp_df)
}
colnames(looDF) <- c('component', 'feature_name', 'feature_stability')
# Get the loadings:
loadingLoo <- as.data.frame(final.splsda.multilevel$loadings$X) %>% mutate(gene = rownames(.))
LoadingComp1 <- loadingLoo %>% filter(gene %in% (looDF %>% 
                                                filter(component == 'comp1') %>% 
                                                pull(feature_name))) %>% 
  mutate(component = 'comp1') %>% filter(comp1 != 0) %>% 
  dplyr::select(gene, comp1, component)
colnames(LoadingComp1) <- c('gene', 'loading', 'component')
# Second component:
LoadingComp2 <- loadingLoo %>% filter(gene %in% (looDF %>% 
                                                filter(component == 'comp2') %>% 
                                                pull(feature_name))) %>%
  mutate(component = 'comp2') %>% filter(comp2 != 0) %>% 
  dplyr::select(gene, comp2, component)
colnames(LoadingComp2) <- c('gene', 'loading', 'component')
# Merge the dataframe:
loadingLooComp <- rbind(LoadingComp1, LoadingComp2)
loadingLooComp <- left_join(loadingLooComp, looDF, by = c('gene'='feature_name', 'component'))
# Plot the feature:
order <- loadingLooComp %>% group_by(component) %>% arrange(loading) %>% pull(gene)
loadingLooComp$gene <- factor(loadingLooComp$gene, levels = order)
loadingLooComp <- loadingLooComp %>% 
  mutate(component = ifelse(component == 'comp1', 'Component 1', 'Component 2'))
# Plot the feature used:
ggplot(loadingLooComp) +
  geom_point(aes(x = loading, y = gene, col = feature_stability)) +
  facet_wrap(~component, scales = 'free', nrow = 1) +  
  theme_minimal() +
  geom_vline(xintercept = 0, lty = 2) +
  scale_color_gradient2() +
  xlab('Gene Loading') + ylab('') -> looFeaturePlot
ggsave(looFeaturePlot, filename = '/scratch/project_mnt/S0007/solal/NMP/Results/Classification/looFeaturePlot.pdf', width = 8, height = 3)
fwrite(loadingLooComp, '/scratch/project_mnt/S0007/solal/NMP/Results/Classification/looFeature.csv')
# Heatmap:
# Heatmap function:
heatmapClassification <- function(count, colData, genes, a=-1, b=1, filter =F) {
  require(pheatmap)
  # Filter based on batch:
  if(is.character(filter)){
    colData <- colData %>% filter(batch %in% filter)
    count <- count[, colData$colID]
    print(all(colnames(count) == colData$colID))
  }
  if(all(colnames(count) == colData$colID) == F){
    print('Datasets in different order')
    break
  }
  # Harmonise this shit because kill me.
  genes <- genes[genes %in% rownames(count)]
  levels <- c(
    colData[which(colData$viability == "CIT"), ]$colID,
    colData[which(colData$viability == "Non-Viable"),]$colID,
    colData[which(colData$viability == "Viable"),]$colID
  )
  rownames(colData) <- colData$colID
  sample_col <- data.frame(
    Viability = colData[levels, ]$viability,
    Time = colData[levels, ]$time,
    Batch = colData[levels, ]$batch
  )
  rownames(sample_col) <- levels
  cal_norm <- function(x) {a + ((x - min(x)) * (b - a) / (max(x) - min(x)))}
  # Remove some genes that might be not expressed in one specific batch:
  hmp_count <- hmp_count[(rowSums(hmp_count) != 0),]
  # Min-max normalization for the heatmap:
  # Per study normalisation:
  studies <- unique(colData$batch)
  hmp_norm <- as.data.frame(matrix(nrow = length(genes)))
  for(batch in studies){
    hmp_study_specific <- count[which(rownames(count) %in% unique(genes)), colData[which(colData$batch == batch),]$colID]
    hmp_study_specific <- t(apply(hmp_study_specific, 1, cal_norm))
    hmp_norm <- cbind(hmp_norm, hmp_study_specific)
  }
  hmp_norm <- hmp_norm[, -1]
  hmp_norm <- hmp_norm[,levels]
  # Annotation colors:
  # Color:
  timeColors <- met.brewer(name = "Ingres", n = length(table(colData$time)))
  batchColors <- met.brewer(name = "VanGogh2", n = 3)
  time <- levels(colData$time)
  ann_colors <- list(
    Time = c(timeColors),
    Viability = c(CIT = "light blue",
                  `Viable` = '#f5bb50',
                  `Non-Viable` = '#94b594'),
    Batch = c(inHouse = batchColors[1],
              hautz = batchColors[2],
              raigani = batchColors[3]))
  names(ann_colors[[1]]) <- time
  
  hmpClustering <- pheatmap(hmp_norm,
                            annotation_col = sample_col,
                            cluster_cols = T,
                            show_colnames = T,
                            annotation_colors = ann_colors,
                            labels_row = F,
                            show_rownames = F,
                            treeheight_row = 0)
  hmpClustering <- ggplotify::as.ggplot(hmpClustering)
  return(hmpClustering)
}

mergedInfo %>% filter(time %in% c('CIT', '3H', '6H')) -> mergedInfoFiltered
mergedCount[(rowSums(mergedCount) > 10), mergedInfoFiltered$colID] -> mergedCountFilter
all(colnames(mergedCountFilter) == mergedInfoFiltered$colID)
# Classification:
# Previous genes:
h1 <- ggplotify::as.ggplot(heatmapClassification(t(Xfull), Yfull, DEGall))
ggsave(h1, filename = '/scratch/project_mnt/S0007/solal/NMP/Results/Classification/DEGinHouse/heatmapDEGall.pdf',
       width = 8, height = 6)