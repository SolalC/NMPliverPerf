# Required packages
library(DESeq2)
library(tidyverse)
library(MASS)
library(lme4)
library(parallel)
library(MetBrewer)
library(cowplot)
source('/scratch/project_mnt/S0007/solal/NMP/Code/script/CreateMatrixDesign.R')
source('/scratch/project_mnt/S0007/solal/NMP/Code/script/FitNegativeBinomia.R')
source('/scratch/project_mnt/S0007/solal/NMP/Code/script/ModelWeight.R')
source('/scratch/project_mnt/S0007/solal/NMP/Code/script/ResultFunctions.R')
cores <- 40
mergedCount <- read.csv('/QRISdata/Q1144/Data/merged/mergedCountData.csv', 
                        row.names = 1)
colnames(mergedCount) <- str_remove(colnames(mergedCount), 'X')
colnames(mergedCount) <- str_replace(colnames(mergedCount), '\\.', '-')
# Read the sample information:
mergedInfo <- read.csv('/QRISdata/Q1144/Data/merged/mergedSampleInfo.csv', 
                       row.names = 1)
all(colnames(mergedCount) == mergedInfo$colID)
# Read the sample information:
singleton <- mergedInfo %>% group_by(sampleID) %>% summarise(n = n()) %>% 
  filter(n <= 2) %>% pull(sampleID)
# Filter the data:
mergedInfoHarmonic <- mergedInfo %>% 
  filter(time %in% c('CIT','1H', '3H', '6H', '12H', '24H'),
         !sampleID %in% singleton,
         !colID %in% c('NMP5_1', 'NMP8_2'))
# Change the CIT label to their viability status:
update_cit_viability <- function(df) {
  time_order <- c("CIT", "1H", "3H", "6H", "12H", "24H")
  df$time <- factor(df$time, levels = time_order, ordered = TRUE)
  df <- df %>%
    arrange(sampleID, time)
  df <- df %>%
    group_by(sampleID) %>%
    mutate(
      viability = case_when(
        viability == "CIT" ~ lead(viability), TRUE ~ viability)) %>%
    ungroup()
  return(df)
}
mergedInfoHarmonic <- update_cit_viability(mergedInfoHarmonic)
# Express time in radian:
time_to_radians <- function(time) {
  2 * pi * time / 24
}

mergedInfoHarmonic <- mergedInfoHarmonic %>% mutate(
  timeInt = as.numeric(str_remove(time, "H")),
  timeInt = ifelse(time == "CIT", 0, timeInt),
  phi = time_to_radians(timeInt))

mergedInfoInput <- mergedInfoHarmonic %>% as.data.frame()
rownames(mergedInfoInput) <- mergedInfoInput$colID

mergedInfoInput <- mergedInfoInput %>% 
  dplyr::select(viability, batch, timeInt, sampleID)

colnames(mergedInfoInput) <- c('condition', 'batch', 'time', 'sample')
countHarmonic <- mergedCount[ ,rownames(mergedInfoInput)]
# Add an expression filter: 
minCount <- 50
sumNonViable <- rowSums(countHarmonic[, mergedInfoInput %>%
                                        filter(condition == 'Non-Viable') %>%
                                        rownames(.)])
sumViable <- rowSums(countHarmonic[, mergedInfoInput %>%
                                     filter(condition == 'Viable') %>%
                                     rownames(.)])
keepGenes <- intersect(names(sumViable[sumViable > minCount]),
                       names(sumNonViable[sumNonViable > minCount]))

countHarmonicFiltered <- countHarmonic[which(rownames(countHarmonic) %in% keepGenes),]
all(colnames(countHarmonicFiltered) == rownames(mergedInfoInput))
# Plot the core clock genes:
coreClock <- c('PER2','CRY1')
countCoreClock <- t(mergedCount[coreClock,rownames(mergedInfoInput)]) %>% 
  as.data.frame() %>% mutate(colID = rownames(.)) %>% 
  left_join(., mergedInfoInput %>% mutate(colID = rownames(.)), by = 'colID') %>%
  mutate(timeR = (time/24)*2*pi) %>%
  pivot_longer(cols = coreClock, names_to = 'gene', values_to = 'count') %>% 
  mutate(logCount = log2(count+1),
         cpm = edgeR::cpm(count),
         condition = factor(condition, levels = c('Viable', 'Non-Viable')))

radianBreaks <- (c(0, 6, 12, 18, 24)/24*2*pi)
ggplot(countCoreClock, aes(x = timeR, y = logCount, col = condition)) +
  geom_point() +
  facet_wrap(~gene+condition) +
  geom_smooth(method = 'lm', formula = y~cos(x)+sin(x), col = 'darkred') +
  geom_smooth(method = 'lm', formula = y~x, col = 'black') +
  scale_color_manual(values = c('#F5BB50' , '#16144F'), name ='') +
  theme_minimal() +
  xlab('Perfusion Time (h)') +
  ylab('log2(Count + 1)') +
  scale_x_continuous(breaks = radianBreaks,
                     labels = c(0,6,12,18,24)) +
  theme(legend.position = 'bottom') -> coreClockPlot
ggsave(coreClockPlot, filename = '/scratch/project_mnt/S0007/solal/NMP/Results/Circadian/coreClockPlot.pdf', width = 7, height = 7, dpi = 300)
# Log likelihood selection for the two genes in viable and non-viable:
library(lmtest)
# Viable:
linearCRY1Viable <- lm(logCount ~ timeR, data = countCoreClock%>% filter(gene == 'CRY1', condition == 'Viable'))
harmonicCRY1Viable <- lm(logCount ~ cos(timeR)+sin(timeR), data = countCoreClock %>% filter(gene == 'CRY1', condition == 'Viable'))
lrtest(linearCRY1Viable, harmonicCRY1Viable)
linearPER2Viable <- lm(logCount ~ timeR, data = countCoreClock%>% filter(gene == 'PER2', condition == 'Viable'))
harmonicPER2Viable <- lm(logCount ~ cos(timeR)+sin(timeR), data = countCoreClock %>% filter(gene == 'PER2', condition == 'Viable'))
lrtest(linearPER2Viable, harmonicPER2Viable)
# Non-Viable:
linearCRY1NonViable <- lm(logCount ~ timeR, data = countCoreClock%>% filter(gene == 'CRY1', condition == 'Non-Viable'))
harmonicCRY1NonViable <- lm(logCount ~ cos(timeR)+sin(timeR), data = countCoreClock %>% filter(gene == 'CRY1', condition == 'Non-Viable'))
lrtest(linearCRY1NonViable, harmonicCRY1NonViable)
linearPER2NonViable <- lm(logCount ~ timeR, data = countCoreClock%>% filter(gene == 'PER2', condition == 'Non-Viable'))
harmonicPER2NonViable <- lm(logCount ~ cos(timeR)+sin(timeR), data = countCoreClock %>% filter(gene == 'PER2', condition == 'Non-Viable'))
lrtest(linearPER2NonViable, harmonicPER2NonViable)
# Show those genes based on real time:
inHouseSampleInfo <- data.table::fread('/QRISdata/Q1144/Data/inHouse/NMPsampleCovariates.csv')
realTimeCoreClock <- countCoreClock %>% filter(batch == 'inHouse') %>% 
  left_join(., inHouseSampleInfo %>% mutate(realTime = time) %>% 
              dplyr::select(colID, realTime),
            by = 'colID') %>% mutate(time_parts = as.POSIXlt(strptime(
              realTime, format = "%H:%M:%S")),
              decimal_hours = time_parts$hour + time_parts$min/60 + time_parts$sec/3600,
              realTimeRadian = decimal_hours * 2 * pi / 24)
# Plot based on time of the day.
ggplot(realTimeCoreClock, aes(x = realTimeRadian, y = logCount, col = condition)) +
  geom_point() +
  facet_wrap(~gene+condition) +
  scale_color_manual(values = c('#F5BB50' , '#16144F'), name ='') +
  theme_minimal() +
  xlab('Sample time (h)') +
  ylab('log2(Count + 1)') +
  scale_x_continuous(breaks = radianBreaks,
                     labels = c(0,6,12,18,24)) +
  theme(legend.position = 'bottom') -> inHouseRealTimeClock

ggplot(realTimeCoreClock, aes(x = timeR, y = logCount, col = condition)) +
  geom_point() +
  facet_wrap(~gene+condition) +
  scale_color_manual(values = c('#F5BB50' , '#16144F'), name ='') +
  theme_minimal() +
  xlab('Perfusion Time (h)') +
  ylab('log2(Count + 1)') +
  scale_x_continuous(breaks = radianBreaks,
                     labels = c(0,6,12,18,24)) +
  theme(legend.position = 'bottom') -> inHousePerfusionTimeClock
ggsave(plot_grid(inHouseRealTimeClock, inHousePerfusionTimeClock, ncol = 2), 
       filename = '/scratch/project_mnt/S0007/solal/NMP/Results/Circadian/coreClockPlotRealTime.pdf', 
       width = 8, height = 5.85, dpi = 300)


# Fit the mixed effect:
fitAll_mixed <- fit_temporal_expression_mixed(countHarmonicFiltered, 
                                              mergedInfoInput, 
                                              period = 24, 
                                              n_cores=cores)
# Save the results to the disk:
saveRDS(fitAll_mixed, '/scratch/project_mnt/S0007/solal/NMP/Results/Circadian/AllGenesRhythmicity_mixed.rds')


step1_fits <- readRDS('/scratch/project_mnt/S0007/solal/NMP/Results/Circadian/intermediary/step1_fits.rds')
# Add the names to the results:
for (i in names(step1_fits)) {
  names(step1_fits[i][[1]]) <-   (c(i, 'best_model', 'weights'))
}
# Read the results:
results <- extract_parameters(step1_fits)
results.f <- lapply(results, function(x) x %>% filter(BICW >= 0.95))
results.f[[1]] <- results.f[[1]] %>% filter(amplitude > 0.25)
results.f[[3]] <- results.f[[3]] %>% filter(amplitude > 0.25)
results.f[[4]] <- results.f[[4]] %>% filter(amplitude > 0.25)
# Put the phase on a 24 hours window:
results.f[[1]] <- results.f[[1]] %>% mutate(phase = phase + 12)
results.f[[3]] <- results.f[[3]] %>% mutate(phase = phase + 12)
results.f[[4]] <- results.f[[4]] %>% mutate(phase = phase + 12)
results.f[[5]] <- results.f[[5]] %>% mutate(phase1 = phase1 + 12, 
                                            phase2 = phase2 + 12)
# Write all of the circadian genes to the disk:
results.f[[1]]  %>% write.csv('/scratch/project_mnt/S0007/solal/NMP/Results/Circadian/GLMM/csv/sharedRhythmicity.csv')
results.f[[3]]  %>% write.csv('/scratch/project_mnt/S0007/solal/NMP/Results/Circadian/GLMM/csv/nonViableRhythmicity.csv')
results.f[[4]]  %>% write.csv('/scratch/project_mnt/S0007/solal/NMP/Results/Circadian/GLMM/csv/viableRhythmicity.csv')
results.f[[5]]  %>% write.csv('/scratch/project_mnt/S0007/solal/NMP/Results/Circadian/GLMM/csv/independentRhythmicity.csv')

col <- MetBrewer::met.brewer('Hokusai3', n = 5)
rhythmicity <- data.frame(rhythmicity = unlist(lapply(results.f, nrow)),
           model = names(results.f)) %>% 
  mutate(model = ifelse(model == 'condition1', 'Non-Viable', model),
         model = ifelse(model == 'condition2', 'Viable', model),
         model = str_to_title(model),
         model = factor(model, 
                        rev(c('None', 'Shared', 'Viable', 'Non-Viable', 'Independent'))),
         percent = rhythmicity/nrow(countHarmonicFiltered)) %>%
  filter(model != 'None')

ggplot(rhythmicity, aes(x = rhythmicity,y = model, fill = model)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = paste0(rhythmicity, '\n',
                               scales::percent(percent, accuracy = 0.01))), hjust = 0.2) +
  theme_light() +
  theme(legend.position = 'none') +  
  scale_fill_manual(values = c('#ADA43B','#16144F','#F5BB50','#B0799A')) -> barPlotRhythmicity

ggsave(barPlotRhythmicity, 
       filename = '/scratch/project_mnt/S0007/solal/NMP/Results/Circadian/GLMM/barPlotRhythmicity.pdf')

# Viable rhythmicity:
rhythmicityTest <- data.frame(nGene = unlist(lapply(results.f, nrow)),
                          condition = names(unlist(lapply(results.f, nrow))))
# Percentage of rhythmicity in the non-viable samples:
nonViableNrhythm <- rhythmicityTest %>% filter(condition %in% c('shared', 'condition1',
                                                             'independent')) %>%
  pull(nGene) %>% sum()/nrow(mergedCount)
print(round(nonViableNrhythm, 2))
# Percentage of rhythmicity in the viable samples:
viableNrhythm <- rhythmicityTest %>% filter(condition %in% c('shared', 'condition2',
                                                             'independent')) %>%
  pull(nGene) %>% sum()/nrow(mergedCount)
print(round(viableNrhythm, 2))
# Amplitude plot:
ampPlot <- amplitudePlot(results.f, 4, c('#F5BB50', '#16144F'))
ggsave(ampPlot, filename = '/scratch/project_mnt/S0007/solal/NMP/Results/Circadian/GLMM/amplitudePlot.pdf', width = 7, height = 7)

# Phase Plot:
phsPlot <- phasePlot(results.f, c('#F5BB50', '#16144F'))
ggsave(phsPlot, filename = '/scratch/project_mnt/S0007/solal/NMP/Results/Circadian/GLMM/phasePlot.pdf', width = 7, height = 7)

# Plot all circadian genes:
clockGenes <- c('PER1', 'PER2', 'PER3', 'CRY1', 'CRY2','CLOCK', 
                'NPAS2',  'BHLHE40','BHLHE41', 'NR1D1', 
                'NR1D2', 'RORA', 'RORC', 'NFIL3', 'DBP',
                'TEF', 'HLF')
rhythmClockGenes <- do.call(rbind, lapply(results, function(x) x %>% filter(gene %in% clockGenes) %>%
                                            dplyr::select(gene, model, BICW))) 

rhythmClockGenes$gene <- factor(rhythmClockGenes$gene, levels = clockGenes)
rhythmClockGenes <- rhythmClockGenes %>% mutate(model = ifelse(model == 'condition1', 'Non-Viable', model),
                                                model = ifelse(model == 'condition2', 'Viable', model),
                                                model = str_to_title(model))

ggplot(rhythmClockGenes) +
  geom_tile(aes(x = gene, 
                y = 1, 
                fill = model),
            col = 'black') + 
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = 'top') +
  coord_equal() +
  geom_segment(aes(x = 0.6, xend = 5.4, y = 1.6, yend = 1.6), col = 'black') +
  geom_segment(aes(x = 5.6, xend = 11.4, y = 1.6, yend = 1.6), col = 'black') +
  geom_segment(aes(x = 11.6, xend = 15.4, y = 1.6, yend = 1.6), col = 'black') +
  geom_text(aes(x = 3, y = 2, label = 'Core loop'), col = 'black') +
  geom_text(aes(x = 8.5, y = 2, label = 'Interlocking loops'), col = 'black') +
  geom_text(aes(x = 13.5, y = 2, label = 'Output'), col = 'black') +
  geom_segment(aes(x = 0.6, xend = 15.4, y = 2.5, yend = 2.5), col = 'black') +
  scale_fill_manual(values = c('#ADA43B','#16144F','black','#B0799A', '#F5BB50')) -> rhythmicityPlotModel

ggplot(rhythmClockGenes) +
  geom_tile(aes(x = gene, 
                y = 1, 
                fill = BICW),
            col = 'black') + 
  geom_text(aes(x = gene, y = 1, label = ifelse(BICW > 0.90, "*", "")),
            color = "white",
            size = 5, 
            vjust = 0.5) +
  theme_void() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.51),
        legend.position = 'bottom') +
  coord_equal() +
  scale_fill_gradient(low = "#373595", high = "#98423C") -> rhythmicityPlotBICW
rhythmicityPlot <- plot_grid(rhythmicityPlotModel, rhythmicityPlotBICW, ncol = 1)

ggsave(rhythmicityPlot, filename = '/scratch/project_mnt/S0007/solal/NMP/Results/Circadian/GLMM/GeneRhythmicity.pdf')

###########################
# temporal enrichment:
temporalEnrichmentShared <- tempEnrichment(results.f, names(step1_fits), condition = 'Viable')
# temporalEnrichmentNonViable <- tempEnrichment(results.f, names(step1_fits), condition = 'Non-Viable')

# temporalEnrichmentViable$condition <- 'Viable'
# temporalEnrichmentNonViable$condition <- 'Non-Viable'
temporalEnrichment <- temporalEnrichmentShared
pathwaysViable <- temporalEnrichment %>% 
  filter(significant == T) %>% 
  pull(term_name) %>% unique()
temporalEnrichment.f <- temporalEnrichment %>% filter(term_name %in% pathwaysViable)
fwrite(temporalEnrichment.f, '/scratch/project_mnt/S0007/solal/NMP/Results/Circadian/GLMM/csv/temporalEnrichmentViable.csv')
temporalEnrichment.f
# Categories:
# Carbohydrate Metabolism
carboMet <- c('Carbon metabolism', 'Glycolysis / Gluconeogenesis', 
              'Fructose and mannose metabolism',
              'Pentose and glucuronate interconversions',
              'Pyruvate metabolism', 'Citrate cycle (TCA cycle)')
# Amino Acid metabolism
AAmet <- c('Biosynthesis of amino acids', 
           # 'Valine, leucine and isoleucine degradation', 
           # 'Tryptophan metabolism', 'Tyrosine metabolism',
           'Cysteine and methionine metabolism', 
           'Histidine metabolism', 'beta-Alanine metabolism')
# Fatty Acid and Lipid Metabolism   
lipidMet <- c('Fatty acid degradation',
              # 'Butanoate metabolism',
              'PPAR signaling pathway','Retinol metabolism')
# Cytochrome P450
cytP450 <- c('Drug metabolism - cytochrome P450',
             'Metabolism of xenobiotics by cytochrome P450' )
# Apoptosis and Immune Responses     
apopImmune <- c('Apoptosis','NF-kappa B signaling pathway',
                'Complement and coagulation cascades',
                'Neutrophil extracellular trap formation')
# Cellular Signaling:
cellSig <- c('Oxytocin signaling pathway',
             'Parathyroid hormone synthesis, secretion and action',
             'GnRH secretion','Serotonergic synapse','Long-term potentiation')


temporalEnrichment.f <- temporalEnrichment.f %>% 
  mutate(term_category = ifelse(term_name %in%
                                  carboMet, 'Carbohydrate Metabolism', NA),
         term_category = ifelse(term_name %in% 
                                  AAmet, 'Amino Acid metabolism', term_category),
         term_category = ifelse(term_name %in% 
                                  lipidMet, 
                                'Fatty Acid and Lipid Metabolism', term_category),
         term_category = ifelse(term_name %in% 
                                  cytP450, 'Cytochrome P450', term_category),
         term_category = ifelse(term_name %in% 
                                  apopImmune, 'Apoptosis and Immune Responses',
                                term_category),
         term_category = ifelse(term_name %in% cellSig, 
                                'Cellular Signaling', term_category))

temporalEnrichment.f <- temporalEnrichment.f %>% filter(!is.na(term_category))
# lvls <- temporalEnrichment.f %>% filter(significant == T) %>% arrange(phase) %>% pull(term_name) %>% unique()
lvlsCategory <- temporalEnrichment.f %>% filter(significant ==T) %>% arrange(phase) %>% pull(term_category) %>% unique()

lvlsTermName <- temporalEnrichment.f %>% filter(significant == T) %>%
  arrange(term_category) %>% filter(!duplicated(term_name)) %>% pull(term_name)

temporalEnrichment.f <- temporalEnrichment.f %>% 
  mutate(term_category = factor(term_category, levels = lvlsCategory),
         term_name = factor(term_name, levels = lvlsTermName))

# Add the p-value of zero for plotting purposes:
temporalEnrichment.f <- rbind(temporalEnrichment.f,
                              temporalEnrichment.f %>% 
                                filter(term_name == 'Fatty acid degradation',
                                       !duplicated(term_name)) %>% 
                                mutate(phase = 15,p_value = 1))
# Get a similar color by category:
colorNumber <- temporalEnrichment.f %>% filter(!duplicated(term_name)) %>% group_by(term_category) %>% summarise(n=n()) %>% pull(n)

colors <-  MetBrewer::colorblind_palettes

col1 <- met.brewer(colors[1], n = colorNumber[1])
col2 <- met.brewer(colors[2], n = colorNumber[2])
col3 <- met.brewer(colors[3], n = colorNumber[3])
col4 <- met.brewer(colors[4], n = colorNumber[4])
col5 <- met.brewer(colors[5], n = colorNumber[5])
col6 <- met.brewer(colors[6], n = colorNumber[6])

ggplot(temporalEnrichment.f,
       aes(x = phase, y = -log10(p_value), col = term_name)) +
  geom_polygon(aes(group = term_name, fill = term_name), col = 'black', 
               alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  theme_minimal() +
  facet_wrap(~term_category) + 
  coord_radial(expand = F, r.axis.inside = T, start = 0) + 
  scale_fill_manual(values = c(col1, col2, col3, col4, col5, col6)) + 
  scale_x_continuous(
    limits = c(0, 24),
    breaks = seq(0, 24, 2),
    labels = seq(0, 24, 2)
  ) +
  labs(x = 'Perfusion Time (h)', y = '-log10(p-value)') +
  theme(legend.position = 'bottom') -> RadialPlots 
ggsave(RadialPlots, filename = '/scratch/project_mnt/S0007/solal/NMP/Results/Circadian/GLMM/EnrichmentRadialPlot.pdf', width = 20, height = 20)






# 
# ### Comparison between the two analysis:
# # 
# geneA <- results.f[['shared']] %>% pull(gene)
# geneB <- dResults.f %>% filter(chosen_model_stringent == '4') %>% mutate(gene = rownames(.)) %>% pull(gene)
# rhythmicGenes <- intersect(geneA, geneB)
# length(rhythmicGenes)/length(geneA)
# length(rhythmicGenes)/length(geneB)
# # Compare all genes:
# geneSet <- intersect(unlist(lapply(results.f, function(x) x %>% pull(gene))), rownames(dResults.f))
# 
# do.call(rbind, lapply(results.f, function(x) x %>% dplyr::select(gene, model))) %>% mutate(
#   chosen_model_inHouse = ifelse(model == 'none', 1, model),
#   chosen_model_inHouse = ifelse(chosen_model_inHouse == 'condition1', 2, chosen_model_inHouse),
#   chosen_model_inHouse = ifelse(chosen_model_inHouse == 'condition2', 3, chosen_model_inHouse),
#   chosen_model_inHouse = ifelse(chosen_model_inHouse == 'shared', 4, chosen_model_inHouse),
#   chosen_model_inHouse = ifelse(chosen_model_inHouse == 'independent', 5, chosen_model_inHouse)) %>%
#   filter(gene %in% geneSet) -> inHouseComparison
# 
# # Merge the two dataframes:
# comp <- dResults.f %>% mutate(gene = rownames(.),
#                               chosen_model_dryR = chosen_model_stringent,
#                               chosen_model_dryR =ifelse(chosen_model_dryR == 'Ambiguous',NA, chosen_model_dryR)) %>%
#   dplyr::select(gene, chosen_model_dryR) %>%
#   filter(gene %in% geneSet) %>%
#   left_join(., inHouseComparison, by = 'gene')
# # Compare the assignment:
# table(comp$chosen_model_inHouse, comp$chosen_model_dryR, useNA = 'always')
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # # Now plot the rhythmicity for genes:
# # countClock <- countHarmonicFiltered[rownames(countHarmonicFiltered) %in% clockGenes,]
# # logClock <- t(log2(countClock+1)) %>% as.data.frame() %>% mutate(colID = rownames(.)) %>%
# #   left_join(., mergedInfoHarmonic, by = 'colID') %>% 
# #   pivot_longer(cols = clockGenes, names_to = 'gene', values_to = 'log')
# # 
# # lmHarmo <- y ~ sin(x) + cos(x)
# # ggplot(logClock, aes(x = ((timeInt/24)*2*pi), log, col = viability)) + 
# #   geom_point() +
# #   facet_wrap(~gene) +
# #   geom_smooth(formula = y ~ sin(x) + cos(x), method = "lm")


# # Get the genes which are DEG as well as rhythmic:
# DEGallbatch <- readRDS("/scratch/project_mnt/S0007/solal/NMP/Results/DEG/viabilityDEG.rds")
# # Overlap with the DEG? Not sure how to put that forward yet. 
# viabilityDEG <- DEGallbatch[[1]]
# # Extract the value for all the genes:
# nGenes <- nrow(viabilityDEG)
# geneBackground <- rownames(viabilityDEG)
# viableDEG <- topTable(viabilityDEG, coef = 'viabilityViable', number = nGenes)
# nonViableDEG <- topTable(viabilityDEG, coef = 'viabilityNon-Viable', number = nGenes)
# # Get significant genes:
# viableDEG.f <- viableDEG %>% filter(adj.P.Val < 0.05, !between(logFC, -1, 1)) %>% 
#   mutate(gene = rownames(.))
# nonViableDEG.f <- nonViableDEG %>% filter(adj.P.Val < 0.05, !between(logFC, -1, 1)) %>%
#   mutate(gene = rownames(.))
# # Define the sets of genes :
# shared <-  intersect(viableDEG.f$gene, nonViableDEG.f$gene)
# viableOnly <- setdiff(viableDEG.f$gene, nonViableDEG.f$gene)
# nonViableOnly <- setdiff(nonViableDEG.f$gene, viableDEG.f$gene)
# 
# unlist(lapply(results.f, function(x) x %>% filter(gene %in% shared) %>% nrow()) )
# unlist(lapply(results.f, function(x) x %>% filter(gene %in% viableOnly) %>% nrow()))
# unlist(lapply(results.f, function(x) x %>% filter(gene %in% nonViableOnly) %>% nrow()))

