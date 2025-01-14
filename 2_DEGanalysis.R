library(tidyverse)
library(sva)
library(variancePartition)
library(BiocParallel)
library(gprofiler2)
library(eulerr)
library(MetBrewer)
library(lme4)
library(performance)
library(ggrepel)
library(cowplot)
# set the parralelization:
param <- SnowParam(40, "SOCK", progressbar = TRUE)
# Read the data:
mergedCount <- read.csv('/QRISdata/Q1144/Data/merged/mergedCountData.csv', row.names = 1)
colnames(mergedCount) <- str_remove(colnames(mergedCount), 'X')
colnames(mergedCount) <- str_replace(colnames(mergedCount), '\\.', '-')
# Read the sample information:
mergedInfo <- read.csv('/QRISdata/Q1144/Data/merged/mergedSampleInfo.csv', row.names = 1)
all(colnames(mergedCount) == mergedInfo$colID)
# DEG analysis:
performDEG <- function(count, info, filtCount = 5, param,
                       multiBatch){
  # Double check everything is ordered correctly
  if(!all(colnames(count) == info$colID)){
    print('Datasets in different order')
    break
  }
  if (multiBatch) {
    mod1 <- model.matrix(~ viability, data = info)
    # Null model matrix contains only the adjustment variables
    mod0 <- model.matrix(~ 1, data = info)
    # filter low number genes:
    filt <- rowSums(count) > filtCount
    filteredCount <- as.matrix(count[filt, ])
    # Run SVA for sequencing data - restrict to only 2 SVs
    print('Processing sva:')
    svseq <- svaseq(filteredCount, mod1, mod0)
    colnames(svseq$sv) <- paste0("SV", seq(1, ncol(svseq$sv)))
    # Add the SVA to the dataframe:
    info <- cbind(info, svseq$sv)
    # Create the experimental design (random effect, all known covariates
    # and surrogate variables:)
    nSV <- sum(str_detect(colnames(info), 'SV'))
    design <- paste0('~ (1 | batch/sampleID) + ', paste0('SV', 1:nSV, collapse = ' + '), ' + viability')
    # Dream:
    print('Processing Voom:')
    filteredCount <- as.data.frame(filteredCount)
    vobjDream <- voomWithDreamWeights(filteredCount, design, info)
    print('Processing Dream:')
    fitmm <- dream(vobjDream, design, info, BPPARAM = param)
    fitmm <- eBayes(fitmm)
    # Run dream with a contrast for viable vs non-viable:
    L <- makeContrastsDream(design, info,
                            contrasts = c(viabilityContrast = "`viabilityNon-Viable`- viabilityViable"))
    fitContrast <- dream(vobjDream, design, info, L, BPPARAM = param)
    fitContrast <- eBayes(fitContrast)
    
    return(list(fitmm, fitContrast))
  }
  if (!multiBatch) {
    mod1 <- model.matrix(~ viability, data = info)
    # Null model matrix contains only the adjustment variables
    mod0 <- model.matrix(~ 1, data = info)
    # filter low number genes:
    filt <- rowSums(count) > filtCount
    filteredCount <- as.matrix(count[filt, ])
    # Run SVA for sequencing data - restrict to only 2 SVs
    print('Processing sva:')
    svseq <- svaseq(filteredCount, mod1, mod0)
    colnames(svseq$sv) <- paste0("SV", seq(1, ncol(svseq$sv)))
    # Add the SVA to the dataframe:
    info <- cbind(info, svseq$sv)
    # Create the experimental design (random effect, all known covariates
    # and surrogate variables:)
    nSV <- sum(str_detect(colnames(info), 'SV'))
    design <- paste0('~ (1 | sampleID) + ', paste0('SV', 1:nSV, collapse = ' + '), ' + viability')
    # Dream:
    print('Processing Voom:')
    filteredCount <- as.data.frame(filteredCount)
    vobjDream <- voomWithDreamWeights(filteredCount, design, info)
    print('Processing Dream:')
    fitmm <- dream(vobjDream, design, info, BPPARAM = param)
    fitmm <- eBayes(fitmm)
    
    L <- makeContrastsDream(design, info,
                            contrasts = c(viabilityContrast = "`viabilityNon-Viable`- viabilityViable"))
    fitContrast <- dream(vobjDream, design, info, L, BPPARAM = param)
    fitContrast <- eBayes(fitContrast)
    
    return(list(fitmm, fitContrast))
  }
}
# Merges together analysis:
viabilityInfo <- mergedInfo %>% filter(time %in% c('CIT', '3H', '6H'))
viabilityCount <- mergedCount[,viabilityInfo$colID]
# Merged analysis:
viabilityDEG <- performDEG(viabilityCount, viabilityInfo, 5, param,
                           multiBatch = T)
saveRDS(viabilityDEG, "/scratch/project_mnt/S0007/solal/NMP/Results/DEG/viabilityDEG.rds")
# Inhouse:
# Merges together analysis:
inHouseInfo <- mergedInfo %>% filter(batch == 'inHouse',
                                     time %in% c('CIT', '3H', '6H'))
inHouseCount <- mergedCount[,inHouseInfo$colID]
inHouseDEG <- performDEG(inHouseCount, inHouseInfo, 5, param,
                         multiBatch = F)
saveRDS(inHouseDEG, "/scratch/project_mnt/S0007/solal/NMP/Results/DEG/inHouseDEG.rds")
# Hautz:
# Merges together analysis:
hautzInfo <- mergedInfo %>% filter(batch == 'hautz',
                                   time %in% c('CIT', '3H', '6H'))
hautzCount <- mergedCount[,hautzInfo$colID]
hautzDEG <- performDEG(hautzCount, hautzInfo, 5, param,
                       multiBatch = F)
saveRDS(hautzDEG, "/scratch/project_mnt/S0007/solal/NMP/Results/DEG/hautzDEG.rds")
# Raigani:
# Merges together analysis:
raiganiInfo <- mergedInfo %>% filter(batch == 'raigani',
                                     time %in% c('CIT', '3H', '6H'))
raiganiCount <- mergedCount[,raiganiInfo$colID]
raiganiDEG <- performDEG(raiganiCount, raiganiInfo, 5, param,
                         multiBatch = F)
saveRDS(raiganiDEG, "/scratch/project_mnt/S0007/solal/NMP/Results/DEG/raiganiDEG.rds")

### Plotting purposes:
# Color variable for plot:
col <- met.brewer('Renoir', n = 3)
# Read the DEG from all the combined analysis:
DEGallbatch <- readRDS("/scratch/project_mnt/S0007/solal/NMP/Results/DEG/viabilityDEG.rds")
viabilityDEG <- DEGallbatch[[1]]
# Extract the value for all the genes:
nGenes <- nrow(viabilityDEG)
geneBackground <- rownames(viabilityDEG)
viableDEG <- topTable(viabilityDEG, coef = 'viabilityViable', number = nGenes)
nonViableDEG <- topTable(viabilityDEG, coef = 'viabilityNon-Viable', number = nGenes)
# Get significant genes:
viableDEG.f <- viableDEG %>% filter(adj.P.Val < 0.05, !between(logFC, -1, 1)) %>% 
  mutate(gene = rownames(.))
write.csv(viableDEG.f, '/scratch/project_mnt/S0007/solal/NMP/Results/DEG/viableDEG.csv')
nonViableDEG.f <- nonViableDEG %>% filter(adj.P.Val < 0.05, !between(logFC, -1, 1)) %>%
  mutate(gene = rownames(.))
write.csv(nonViableDEG.f, '/scratch/project_mnt/S0007/solal/NMP/Results/DEG/nonViableDEG.csv')
# Define the sets of genes :
shared <-  intersect(viableDEG.f$gene, nonViableDEG.f$gene)
viableOnly <- setdiff(viableDEG.f$gene, nonViableDEG.f$gene)
nonViableOnly <- setdiff(nonViableDEG.f$gene, viableDEG.f$gene)
# Compare the sign for all genes:
allGenes <- left_join(viableDEG %>% mutate(gene = rownames(.)) %>% select(gene, logFC), 
                      nonViableDEG %>% mutate(gene = rownames(.)) %>% select(gene, logFC),
                      by = 'gene', suffix = c('_Viable', '_NonViable')) %>% 
  filter(gene %in% unique(c(shared, viableOnly, nonViableOnly)))

# Compare the sign for the shared genes:
sharedGenes <- allGenes %>% filter(gene %in% shared)
sharedGenedLm <- lm(logFC_NonViable ~ logFC_Viable, data = sharedGenes)
sharedGenedLm <- round(coef(sharedGenedLm), 3)  # rounds to 3 decimal places
# Add to plot
ggplot(sharedGenes, aes(x = logFC_Viable, y = logFC_NonViable)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  annotate('text', x = -2,y = 4, label = paste("y=",sharedGenedLm[1], '+', sharedGenedLm[2], 'x')) +
  theme_minimal() +
  xlim(-3,6.5) + ylim(-3,6.5) +
  ggtitle('Share Genes') +
  geom_abline(slope = 1, intercept = 0, lty = 1, col = 'darkred') -> sharedGenesComparisonPlot
# Non-viable:
nonViableGenes <- allGenes %>% filter(gene %in% nonViableOnly)
nonViableGenes <- nonViableGenes %>% mutate(sign = ifelse(sign(logFC_Viable) == sign(logFC_NonViable), T, F))
nonViableGenesLm <- lm(logFC_NonViable ~ logFC_Viable, data = nonViableGenes)
nonViableGenesLm <- round(coef(nonViableGenesLm), 3)  # rounds to 3 decimal places
# Add to plot
ggplot(nonViableGenes, aes(x = logFC_Viable, y = logFC_NonViable)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  annotate('text', x = -2,y = 4, label = paste("y=",nonViableGenesLm[1], '+', nonViableGenesLm[2], 'x')) +
  theme_minimal() +
  ggtitle('Non-Viable only') +
  xlim(-3,6.5) + ylim(-3,6.5) +
  geom_abline(slope = 1, intercept = 0, lty = 1, col = 'darkred') +
  geom_text_repel(aes(label = gene), data = nonViableGenes %>% filter(!sign)) -> nonViableComparisonPlot
# Viable:
viableGenes <- allGenes %>% filter(gene %in% viableOnly)
viableGenes <- viableGenes %>% mutate(sign = ifelse(sign(logFC_Viable) == sign(logFC_NonViable), T, F))
viableGenesLm <- lm(logFC_NonViable ~ logFC_Viable, data = viableGenes)
viableGenesLm <- round(coef(viableGenesLm), 3)  # rounds to 3 decimal places
# Add to plot
ggplot(viableGenes, aes(x = logFC_Viable, y = logFC_NonViable)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  annotate('text', x = -2,y = 4, label = paste("y=",viableGenesLm[1], '+', viableGenesLm[2], 'x')) +
  theme_minimal() +
  ggtitle('Viable only') +
  xlim(-3,6.5) + ylim(-3,6.5) +
  geom_abline(slope = 1, intercept = 0, lty = 1, col = 'darkred') +
  geom_text_repel(aes(label = gene), data = viableGenes %>% filter(!sign)) -> viableComparisonPlot
DEGcomparison <- plot_grid(sharedGenesComparisonPlot, nonViableComparisonPlot, viableComparisonPlot, nrow = 1)
ggsave(DEGcomparison, filename = '/scratch/project_mnt/S0007/solal/NMP/Results/DEG/Supplementary/SupplementaryFigure_SharedDEGcomparison.pdf')
# venn diagram:
# List of items
vennList <- c(
  Viable = length(viableOnly),
  `Non-viable` = length(nonViableOnly),
  `Viable&Non-viable` = length(shared)
)
pdf("/scratch/project_mnt/S0007/solal/NMP/Results/DEG/vennDiagram.pdf")
plot(euler(vennList),
     fill = c('#F5BB50', '#16144F'),
     quantities = T
)
dev.off()

# Make a constrast between viable and non-viable:
# Pathways Analysis:
# Viable genes:
set.seed(9763543)
pathwayDB <- 'KEGG'
viableAllPathways <- gost(query = viableDEG.f$gene, 
                          organism = "hsapiens",significant = TRUE, 
                          user_threshold = 0.05, correction_method = "g_SCS", 
                          domain_scope = "custom", custom_bg = geneBackground, 
                          sources = pathwayDB)$result
viableAllPathways <- viableAllPathways %>% mutate(condition = 'Viable - All Genes')
# Non-viable genes:
nonViableAllPathways <- gost(query = nonViableDEG.f$gene, 
                             organism = "hsapiens",significant = TRUE, 
                             user_threshold = 0.05, correction_method = "g_SCS", 
                             domain_scope = "custom", custom_bg = geneBackground, 
                             sources = pathwayDB)$result
nonViableAllPathways <- nonViableAllPathways %>% mutate(condition = 'Non-Viable - All Genes')
# Intersect:
sharedPathways <-  gost(query = shared, 
                        organism = "hsapiens",significant = TRUE, 
                        user_threshold = 0.05, correction_method = "g_SCS", 
                        domain_scope = "custom", custom_bg = geneBackground, 
                        sources = pathwayDB)$result
sharedPathways <- sharedPathways %>% mutate(condition = 'Shared Genes')

# viable only:
viableOnlyPathways <-  gost(query = viableOnly, 
                            organism = "hsapiens",significant = TRUE, 
                            user_threshold = 0.05, correction_method = "g_SCS", 
                            domain_scope = "custom", custom_bg = geneBackground, 
                            sources = pathwayDB)$result
viableOnlyPathways <- viableOnlyPathways %>% mutate(condition = 'Viable - Only genes')
# Non-viable only:
nonViableOnlyPathways <-  gost(query = nonViableOnly, 
                               organism = "hsapiens",significant = TRUE, 
                               user_threshold = 0.05, correction_method = "g_SCS", 
                               domain_scope = "custom", custom_bg = geneBackground, 
                               sources = pathwayDB)$result
nonViableOnlyPathways <- nonViableOnlyPathways %>% mutate(condition = 'Non-Viable - Only genes')
# Merge all the results:
allPathways <- bind_rows(viableAllPathways, nonViableAllPathways, sharedPathways, viableOnlyPathways, nonViableOnlyPathways)
data.table::fwrite(as.data.frame(allPathways),  file = '/scratch/project_mnt/S0007/solal/NMP/Results/DEG/allPathways.csv')
# Plot the results:
allPathways.f <- allPathways %>% 
  filter(!condition %in% c('Viable - All Genes', 'Non-Viable - All Genes')) %>%
  mutate(term_name = reorder(term_name, -log10(p_value))) %>% 
  group_by(term_name) %>%  mutate(n = as.character(n()))
# Plot:
ggplot(allPathways.f, aes(x = condition, y = term_name, 
                          size = -log10(p_value), col = condition)) +
  geom_point() +
  scale_color_manual(values = col) +
  theme_light() +
  facet_wrap(~condition, scale = 'free_y') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> KEGGplot
ggsave(KEGGplot, filename = '/scratch/project_mnt/S0007/solal/NMP/Results/DEG/KEGGplot.pdf')
# Viable versus non-viable
contrastDEG <- DEGallbatch[[2]]
# Extract the value for all the genes:
constrastDEG <- topTable(contrastDEG, coef = 'viabilityContrast', number = nGenes)
# Get significant genes:
constrastDEG.f <- constrastDEG %>% filter(adj.P.Val < 0.05, !between(logFC, -1, 1)) %>% 
  mutate(gene = rownames(.))
data.table::fwrite(constrastDEG.f, '/scratch/project_mnt/S0007/solal/NMP/Results/DEG/constrastDEG.csv')
ggplot(constrastDEG.f %>% mutate(gene = fct_reorder(gene, logFC))) + 
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(aes(x = logFC, y = gene)) +
  theme_minimal()

t(log2(viabilityCount+1)[c('CIDEC', 'PPM1E', 'EGR1', 'PRSS22'),]) %>% as.data.frame() %>%  mutate(colID = rownames(.)) %>%
  left_join(., viabilityInfo, by = 'colID') %>% pivot_longer(cols = c('CIDEC', 'PPM1E', 'EGR1', 'PRSS22'), 
                                                             names_to = 'gene', 
                                                             values_to = 'count') %>%
  ggplot(aes(x = gene, y = count, fill = viability)) + 
  geom_boxplot()

# Pathway analysis:
contrastPathways <-  gost(query = constrastDEG.f$gene, 
                               organism = "hsapiens",significant = TRUE, 
                               user_threshold = 0.05, correction_method = "g_SCS", 
                               domain_scope = "custom", custom_bg = geneBackground, 
                               sources = pathwayDB)$result

# IRI analysis:
IRIgenes <- c(
  "SGPL1", "DLG1", "RLN1", "MLH3", "TFAP2A", "TTPA", "LPA", "KLHDC10",
  "ZNF710", "LRRN3", "EMP1", "INSIG1", "PPP1R15A", "TP53BP2", "MCL1",
  "WEE1", "GADD45A", "ETS1", "TRAM1", "LIFR", "TRAF4", "KLF5", "AKAP12",
  "MMP19", "PTPN1", "TACSTD2", "PPP2R2A", "RORA", "ADM", "MAFF",
  "TMEM184B", "DUSP5", "BAG3", "IL18RAP", "NFKB1", "PLAUR", "SERPINE1",
  "C5AR1", "B4GALT6", "ALAS1", "CORO1C", "HBB", "SLC20A1", "SLC7A5",
  "NABP1", "SMIM13", "APOLD1", "TIPARP", "ISG20L2"
)
# Check how many genes are enriched within DEGs sets:
print(sum(viableDEG.f$gene %in% IRIgenes) / length(IRIgenes))
print(sum(nonViableDEG.f$gene %in% IRIgenes) / length(IRIgenes))
IRIgenes <- intersect(IRIgenes, rownames(viabilityCount))
IRItotest <- intersect(intersect(IRIgenes, viableDEG.f$gene), nonViableDEG.f$gene)
# Create the dataframe needed to test gene progression:
viabilityCPM <- as.data.frame(edgeR::cpm(viabilityCount))[IRIgenes,]
CITsamples <- viabilityInfo %>% filter(time == 'CIT') %>% pull(colID)
viable3Hsamples <- viabilityInfo %>% filter(viability == 'Viable', time == '3H') %>% pull(colID)
viable6Hsamples <- viabilityInfo %>% filter(viability == 'Viable', time == '6H') %>% pull(colID)

nonViable3Hsamples <- viabilityInfo %>% filter(viability == 'Non-Viable', time == '3H') %>% pull(colID)
nonViable6Hsamples <- viabilityInfo %>% filter(viability == 'Non-Viable', time == '6H') %>% pull(colID)

viabilityCPMCIT <- as.data.frame(t(viabilityCPM[,CITsamples] )) %>% mutate(condition = 'CIT')
viabilityCPMviable3H <- as.data.frame(t(viabilityCPM[,viable3Hsamples] )) %>% mutate(condition = 'Viable: 3hours')
viabilityCPMviable6H <- as.data.frame(t(viabilityCPM[,viable6Hsamples] )) %>% mutate(condition = 'Viable: 6hours')
viabilityCPMnonViable3H <- as.data.frame(t(viabilityCPM[,nonViable3Hsamples] )) %>% mutate(condition = 'Non-Viable: 3hours')
viabilityCPMnonViable6H <- as.data.frame(t(viabilityCPM[,nonViable6Hsamples] )) %>% mutate(condition = 'Non-Viable: 6hours')
# Bind the dataframe:
IRIdf <- rbind(viabilityCPMCIT, viabilityCPMviable3H, viabilityCPMviable6H, 
               viabilityCPMnonViable3H,viabilityCPMnonViable6H) %>% 
  mutate(condition = factor(condition,
                            levels = c('CIT', 'Viable: 3hours', 'Viable: 6hours',
                                       'Non-Viable: 3hours', 'Non-Viable: 6hours')),
         colID = rownames(.)) %>%
  left_join(., viabilityInfo, by = 'colID')

# Mixed model regression for all the IRI genes:
IRImixedModel <- data.frame()
for (i in IRIgenes) {
  designLMM <- paste0(i, "~ (1|batch/sampleID) +  condition")
  lmm <- lmer(designLMM, data = IRIdf, lmerControl(optimizer = "nloptwrap",calc.derivs = F),
              REML = T)
  # sum <- summary(lmm, ddf = "Kenward-Roger")
  sum <- summary(lmm)
  perf <- model_performance(lmm)
  report <- as.data.frame(sum$coefficients)
  report$gene <- i
  report$R2_conditional <- perf$R2_conditional
  report$R2_marginal <- perf$R2_marginal
  IRImixedModel <- rbind(IRImixedModel, report)
}
# Wrangle the regression results:
downregulatedIRIgenes <- c("SGPL1", "DLG1", "RLN1", "MLH3",
                           "TFAP2A", "TTPA", "LPA", "KLHDC10", "ZNF710", "LRRN3",
                           "BLTP2")
upregulatedIRIgenes <- c(  "EMP1", "INSIG1", "PPP1R15A",
                           "TP53BP2", "MCL1", "WEE1", "GADD45A", "EMP1", "TRAM1",
                           "LIFR", "TRAF4", "KLF5", "AKAP12", "MMP19", "PTPN1",
                           "TACSTD2", "PPP2R2A", "RORA", "ADM", "MAFF", "TMEM184B",
                           "ETS1", "DUSP5", "BAG3", "IL18RAP", "NFKB1", "PLAUR",
                           "SERPINE1", "C5AR1", "B4GALT6", "ALAS1", "NABP1",
                           "SMIM13", "APOLD1", "TIPARP", "ISG20L2", "CORO1C", "HBB",
                           "SLC20A1", "SLC7A5")

IRImixedModel.f <- IRImixedModel %>% mutate(comparison = str_remove(rownames(.), "condition"),
                                            comparison = gsub("[0-9]+$", "", comparison),
                                            comparison = factor(comparison, 
                                                                levels(IRIdf$condition)),
                                            viability = as.integer(!str_detect(comparison, 
                                                                               'Non-Viable')),
                                            sign = sign(Estimate),
                                            IRI = ifelse(gene %in% downregulatedIRIgenes,
                                                         "Downregulated", NA),
                                            IRI = ifelse(gene %in% upregulatedIRIgenes,
                                                         "Upregulated", IRI)) %>% 
  filter(!is.na(comparison)) 
# Plot:
se <- function(y) {
  sd(y) / sqrt(length(y))
}

IRImixedModel.ff <- IRImixedModel.f %>% select(comparison, Estimate, gene, `Std. Error`, IRI)
colnames(IRImixedModel.ff)[4] <- 'se'
IRImixedModel.ff$viability <- unlist(map(str_split(IRImixedModel.ff$comparison, ':'), 1))
data.table::fwrite(IRImixedModel.ff %>% filter(comparison != 'CIT'), 
                   file = '/scratch/project_mnt/S0007/solal/NMP/Results/DEG/IRImixedModel.csv')
CITforPlotting <- data.frame(comparison=rep('CIT', 4),
                  Estimate = rep(0, 4),
                  gene = rep(c('EMP1','SGPL1'), 2),
                  se = rep(0, 4),
                  IRI = rep(c('Upregulated', 'Downregulated'), 2),
                  viability = c('Viable', 'Viable', 'Non-Viable', 'Non-Viable'))

IRImixedModel.ff <- rbind(IRImixedModel.ff, CITforPlotting)
IRImixedModel.ff <- IRImixedModel.ff %>% 
  mutate(viability = factor(viability, levels = c('Viable', 'Non-Viable')))

# Extract the coefficients:
coefs <- IRImixedModel.ff %>% filter(comparison != 'CIT') %>% 
  group_by(IRI, comparison) %>% 
  summarise(mean = round(mean(Estimate), 2), 
            se = round(se(Estimate), 2))


ggplot(IRImixedModel.ff %>% filter(IRI == 'Upregulated'),aes(x = comparison, y = Estimate)) +
  facet_wrap(~ viability, scale = 'free_x') +
  geom_point(col = 'lightgrey', alpha = 0.8, position = 
               position_jitter(width = 0.1, height = 0.1)) +
  # geom_text_repel(data = IRImixedModel.ff %>% filter(abs(Estimate) > 1),
  #                 aes(label = gene), size = 3) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_minimal() +
  # scale_color_manual(values = rep("grey", 40)) +
  theme(legend.position = "none") +
  stat_summary(
    fun = mean,
    fun.min = function(y) mean(y) - se(y),
    fun.max = function(y) mean(y) + se(y),
    color = "dark red",
    geom = "pointrange",
    show.legend = FALSE
  ) +
  xlab("Perfusion Time") +
  ylab("Gene expression linear coefficient") +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank()) -> IRIup

ggplot(IRImixedModel.ff %>% filter(IRI == 'Downregulated'),aes(x = comparison, y = Estimate)) +
  facet_wrap(~ viability, scale = 'free_x') +
  geom_point(col = 'lightgrey', alpha = 0.8, position = 
               position_jitter(width = 0.1, height = 0.1)) +
  # geom_text_repel(data = IRImixedModel.ff %>% filter(abs(Estimate) > 1),
  #                 aes(label = gene), size = 3) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_minimal() +
  # scale_color_manual(values = rep("grey", 40)) +
  theme(legend.position = "none") +
  stat_summary(
    fun = mean,
    fun.min = function(y) mean(y) - se(y),
    fun.max = function(y) mean(y) + se(y),
    color = "dark red",
    geom = "pointrange",
    show.legend = FALSE
  ) +
  xlab("Perfusion Time") +
  ylab("Gene expression linear coefficient") -> IRIdown

IRIplot <- cowplot::plot_grid(IRIup, IRIdown, ncol = 1)
ggsave(IRIplot, filename = '/scratch/project_mnt/S0007/solal/NMP/Results/DEG/IRIplot.pdf')

# clock genes DEG:
clockGenes <- c('PER1', 'PER2', 'PER3', 'CRY1', 'CRY2','CLOCK', 
                'NPAS2',  'BHLHE40','BHLHE41', 'NR1D1', 
                'NR1D2', 'RORA', 'RORC', 'NFIL3', 'DBP',
                'TEF', 'HLF')

viableDEG.f %>% filter(gene %in% clockGenes)
nonViableDEG.f %>% filter(gene %in% clockGenes)
constrastDEG.f %>% mutate(gene = rownames(.)) %>% filter(gene %in% clockGenes)


###  Interaction term:
viabilityInfo <- mergedInfo %>% filter(time %in% c('CIT', '3H', '6H')) %>% 
  mutate(time = ifelse(time == 'CIT', 0, time),
         time = str_remove(time, 'H') %>% 
           as.numeric()) %>%
  group_by(sampleID) %>%
  mutate(viability = ifelse(viability == "CIT",
                            first(viability[viability != "CIT"], 
                                  default = "Non-Viable"), viability)) %>%
  ungroup()

viabilityCount <- mergedCount[,viabilityInfo$colID]

mod1 <- model.matrix(~ time + viability, data = viabilityInfo)
# Null model matrix contains only the adjustment variables
mod0 <- model.matrix(~ 1, data = viabilityInfo)
# filter low number genes:
filt <- rowSums(viabilityCount) > 5
viabilityCountFiltered <- as.matrix(viabilityCount[filt, ])
# Run SVA for sequencing data - restrict to only 2 SVs
svseq <- svaseq(viabilityCountFiltered, mod1, mod0)
colnames(svseq$sv) <- paste0("SV", seq(1, ncol(svseq$sv)))
# Add the SVA to the dataframe:
viabilityInfo <- cbind(viabilityInfo, svseq$sv)
# Create the experimental design (random effect, all known covariates
# and surrogate variables:)
nSV <- sum(str_detect(colnames(viabilityInfo), 'SV'))
design <- paste0('~ (1 | batch/sampleID) + ', 
                 paste0('SV', 1:nSV, collapse = ' + '), 
                 ' + time+ viability + time:viability')
# Dream:
print('Processing Voom:')
viabilityCountFiltered <- as.data.frame(viabilityCountFiltered)
vobjDream <- voomWithDreamWeights(viabilityCountFiltered, design, viabilityInfo)
print('Processing Dream:')
fitmm <- dream(vobjDream, design, viabilityInfo, BPPARAM = param)
fitmm <- eBayes(fitmm)

saveRDS(fitmm, "/scratch/project_mnt/S0007/solal/NMP/Results/DEG/interactionTermDEG.rds")

nGenes <- nrow(fitmm)
geneBackground <- rownames(fitmm)
viableDEG <- topTable(fitmm, coef = 'viabilityViable', number = nGenes)
timeViableDEG <- topTable(fitmm, coef = 'time:viabilityViable', number = nGenes)

plot_gene_expression <- function(gene_name, viabilityCountFiltered, viabilityInfo) {
  gene_expression <- viabilityCountFiltered[gene_name, ]
  plot_data <- data.frame(
    time = viabilityInfo$time,
    viability = viabilityInfo$viability,
    expression = log(unlist(c(gene_expression))+1)
  )
  
  summary_data <- plot_data %>%
    group_by(time, viability) %>%
    summarize(
      mean_expression = mean(expression),
      se_expression = sd(expression) / sqrt(n()),
      .groups = "drop"
    )
  # Generate the plot
  p <- ggplot(plot_data, aes(x = time, y = expression, color = viability)) +
    geom_point(size = 1) +
    geom_line(aes(x = time, y = mean_expression, group = viability), data = summary_data) +
    geom_errorbar(
      data = summary_data,
      aes(y = mean_expression, 
          ymin = mean_expression - se_expression, 
          ymax = mean_expression + se_expression),
      width = 0.2) +
    geom_point(data = summary_data,
               aes(y = mean_expression, shape = viability),
               size = 1) +
    labs(
      title = paste("Gene Expression of", gene_name),
      x = "Time",
      y = "Expression Level",
      color = "Viability",
      shape = "Viability"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  return(p)
}

p1 <- plot_gene_expression('PPM1E', viabilityCountFiltered, viabilityInfo)
p2 <- plot_gene_expression('CALCRL', viabilityCountFiltered, viabilityInfo)
p3 <- plot_gene_expression('ZNF34', viabilityCountFiltered, viabilityInfo)
p4 <- plot_gene_expression('KSR1', viabilityCountFiltered, viabilityInfo)
p5 <- plot_gene_expression('GJC2', viabilityCountFiltered, viabilityInfo)
p6 <- plot_gene_expression('ITGA9', viabilityCountFiltered, viabilityInfo)
p <- cowplot::plot_grid(p1, p2, p3, p4, p5, p6)
ggsave(p, filename = '/scratch/project_mnt/S0007/solal/NMP/Results/DEG/interactionPlot.pdf',
       width = 14, height = 14, dpi = 300)







plink --ped /QRISdata/Q7912/Data/Genotype/PLINK_071221_0922/CLOCKOnlyReports.ped --extract /QRISdata/Q7912/Data/Genotype/SNPs/snp-names.txt --recodeA --out /QRISdata/Q7912/Data/Genotype/SNPs/genotypesExtracted --map /QRISdata/Q7912/Data/Genotype/PLINK_071221_0922/CLOCKOnlyReports.map


genotypes <- meffil.extract.genotypes('/QRISdata/Q7912/Data/Genotype/methylationSNPs/CLOCKmethylationSNPs.raw')
