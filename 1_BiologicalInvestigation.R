library(gprofiler2)
library(variancePartition)
library(tidyverse)
library(MetBrewer)
library(sva)
library(edgeR)
library(EnhancedVolcano)
library(qqman)
library(performance)
library(lmerTest)
library(lubridate)
library(rlist)

# Color variables:
VolColorViable <- c(met.brewer(name = "Renoir", n = 2))
VolColorNonViable <- c(met.brewer(name = "Hokusai1", n = 2))
oneDriveOutput <- "/Users/uqschauq/Library/CloudStorage/OneDrive-TheUniversityofQueensland/NMP_liver/Results/"
# Set the directory:
setwd("/Users/uqschauq/Documents/NMP/Code")
# Read count data:
countData <- read.delim("../Data/RNAseq/featureCountsMatrixClean.csv",
  sep = ",", row.names = 1
)

# Read sample information:
sampleinfo <- readr::read_tsv("../Data/RNAseq/sampleinfo.tsv")
# Select only the sample of interest:
sampleinfo_sub <- sampleinfo[which(sampleinfo$perf_success %in%
  c("CIT", "Perf_Unsuccessful", "Perf_Successful")), ]
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
# Create the surrogate variables:
# Create the DGElist object
dge <- DGEList(
  counts = countData_sub, samples = sampleinfo_sub,
  genes = rownames(countData)
)
# Double check everything is ordered correctly
all(colnames(dge$counts) == rownames(dge$samples))
# Set null and alternative models (ignore batch)
# full model matrix - including both the adjustment variables and the
# variable of interest
mod1 <- model.matrix(~ sex + bmi + age + death + RIN + perf_success,
  data = dge$samples
)
# Null model matrix contains only the adjustment variables
mod0 <- model.matrix(~ sex + bmi + age + death + RIN, data = dge$samples)
# filter low number genes:
filt <- apply(dge$counts, 1, function(x) length(x[x > 5]) >= 2)
filtered <- as.matrix(dge$counts[filt, ])
# Run SVA for sequencing data - restrict to only 2 SVs
svseq <- svaseq(filtered, mod1, mod0)
colnames(svseq$sv) <- paste0("SV", seq(1, ncol(svseq$sv)))
# Add the SVA to the dataframe:
sampleinfo_sub <- cbind(sampleinfo_sub, svseq$sv)
# Create the experimental design (random effect, all known covariates
# and surrogate variables:)
design <- ~ (1 | sample) + sex + bmi + age + death + RIN + SV1 +
  SV2 + SV3 + perf_success
# filter genes by number of counts
isexpr <- rowSums(cpm(countData_sub) > 0.1) >= 5
# Standard usage of limma/voom
dge <- DGEList(countData_sub[isexpr, ])
dge <- calcNormFactors(dge)
vobjDream <- voomWithDreamWeights(dge, design, sampleinfo_sub)
fitmm <- dream(vobjDream, design, sampleinfo_sub)
fitmm <- eBayes(fitmm)
# Write the genes used for DEG testing to the disk:
data.frame(genes = rownames(dge)) %>% write.csv(
  paste0(
    oneDriveOutput,
    "Supplementary/DEGbackgroundGenesList.csv"
  ),
  row.names = F
)

# Extract the results:
viableResults <- topTable(fitmm,
  coef = "perf_successPerf_Successful",
  number = 300000
)
nonViableResults <- topTable(fitmm,
  coef = "perf_successPerf_Unsuccessful",
  number = 300000
)
# Check the QQ plots for possible statistical inflation:
qqman::qq(-2 * log(viableResults$P.Value))
qqman::qq(-2 * log(nonViableResults$P.Value))

# extract significant results:
# For viable livers:
fitperfViable <- viableResults %>%
  filter(adj.P.Val < 0.05, !between(logFC, -0.25, 0.25)) %>%
  mutate(gene = rownames(.)) %>%
  pull(gene)
# And non-viables:
fitperfNonViable <- nonViableResults %>%
  filter(adj.P.Val < 0.05, !between(logFC, -0.25, 0.25)) %>%
  mutate(gene = rownames(.)) %>%
  pull(gene)
# Write the two genesets to the disk:
viableResults %>%
  filter(adj.P.Val < 0.05, !between(logFC, -0.25, 0.25)) %>%
  mutate(gene = rownames(.)) %>%
  write.csv(
    paste0(
      oneDriveOutput,
      "Supplementary/CIT_viable_significant_DEG.csv"
    ),
    row.names = F
  )
# Non viable:
nonViableResults %>%
  filter(adj.P.Val < 0.05, !between(logFC, -0.25, 0.25)) %>%
  mutate(gene = rownames(.)) %>%
  write.csv(paste0(oneDriveOutput, "Supplementary/CIT_NonViable_significant_DEG.csv"),
    row.names = F
  )

# Write the venn diagram to the disk:
library(eulerr)
vennList <- c(
  Viable = length(fitperfViable),
  `Non-viable` = length(fitperfNonViable),
  `Viable&Non-viable` = length(intersect(
    fitperfViable,
    fitperfNonViable
  ))
)
pdf(paste0(oneDriveOutput, "Figure1/VennDiagram.pdf"))
plot(euler(vennList),
  fill = c(VolColorViable[2], VolColorNonViable[2]),
  quantities = T
)
dev.off()

### Correlation between shared genes:
# Write the two gene sets to the disk:
intersectGenes <- intersect(fitperfViable, fitperfNonViable)

set1 <- viableResults %>%
  filter(adj.P.Val < 0.05, !between(logFC, -0.25, 0.25)) %>%
  mutate(gene = rownames(.)) %>%
  filter(gene %in% intersectGenes) %>%
  arrange(gene)
# Non viable:
set2 <- nonViableResults %>%
  filter(adj.P.Val < 0.05, !between(logFC, -0.25, 0.25)) %>%
  mutate(gene = rownames(.)) %>%
  filter(gene %in% intersectGenes) %>%
  arrange(gene)

summary(lm(set1$logFC ~ set2$logFC))

# Volcano plots:
# Volcano plots using the enhanced volcano library:
# Viable:
ViableVolcano <- EnhancedVolcano(viableResults,
  lab = rownames(viableResults), x = "logFC",
  y = "P.Value", title = "Successful Perfusion",
  pCutoffCol = "adj.P.Val", pCutoff = 0.05,
  col = c(
    "grey30", "#769FB6", "#0094C6",
    VolColorViable[2]
  ), FCcutoff = 0.5
)

ggsave(ViableVolcano, filename = paste0(
  oneDriveOutput,
  "Figure1/volcanoPlotViable.pdf"
))
# dev.off()
# Non viable:
NonViableVolcano <- EnhancedVolcano(nonViableResults,
  lab = rownames(nonViableResults),
  x = "logFC", y = "P.Value",
  title = "Successful Perfusion",
  pCutoffCol = "adj.P.Val", pCutoff = 0.05,
  col = c(
    "grey30", "#F5E5FC", "#B28B84",
    VolColorNonViable[2]
  ), FCcutoff = 0.5
)
ggsave(NonViableVolcano, filename = paste0(
  oneDriveOutput,
  "Figure1/volcanoPlotNonViable.pdf"
))

# Pathway enrichment:
set.seed(246787)
# Enrichment was performed using the enrichR website.

# The genes investigated for DEG were used as the background.
enrichViable <- read.csv(paste0(
  oneDriveOutput,
  "Supplementary/PathwayViableBackground.csv"
))
# multiple testing correction:
# enrichViable$Adjusted.P.value <- p.adjust(enrichViable$P.value, 'BH')
enrichViable_f <- enrichViable %>% filter(Adjusted.P.value < 0.05)

enrichNonViable <- read.csv(paste0(
  oneDriveOutput,
  "Supplementary/PathwayNonViableBackground.csv"
))
# enrichNonViable$Adjusted.P.value <- p.adjust(enrichNonViable$P.value, 'BH')
enrichNonViable_f <- enrichNonViable %>% filter(Adjusted.P.value < 0.05)
# Bind the results:
enrichRresults <- rbind(enrichViable_f, enrichNonViable_f)
enrichRresults$Term <- str_remove(enrichRresults$Term, " WP.*")
# Plot:
# EnrichR:
enrichRresults$Term <- factor(enrichRresults$Term,
  levels = unique(enrichRresults$Term[order(enrichRresults$P.value)])
)
# Plot:
ggplot(enrichRresults %>% filter(Set == "WP"), aes(
  x = Term, y = Test,
  size = -log10(Adjusted.P.value), col = Test
)) +
  geom_point() +
  scale_color_manual(values = c(
    VolColorNonViable[2],
    VolColorViable[2]
  )) +
  theme_minimal() +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5, hjust = 1
  )) ->
enrichmentPlot
# Save the plot:
ggsave(enrichmentPlot, filename = paste0(
  oneDriveOutput,
  "Figure1/PathwaysPlot.pdf"
), width = 14, height = 7)

# IRI genes:
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
geneset <- c(fitperfViable)
print(sum(fitperfViable %in% IRIgenes) / length(IRIgenes))
print(sum(fitperfNonViable %in% IRIgenes) / length(IRIgenes))

IRItotest <- intersect(intersect(IRIgenes, fitperfViable), fitperfNonViable)
# Create the dataframe needed to test gene progression:
cpmData_sub <- as.data.frame(cpm(countData_sub))

CITstartIRI <- as.data.frame(t(cpmData_sub
[IRIgenes, sampleinfo_sub
  [which(sampleinfo_sub$perf == "CIT_start"), ]
  $colID]))
CITstartIRI$condition <- "CIT:start"
# End
CITendIRI <- as.data.frame(t(cpmData_sub
[IRIgenes, sampleinfo_sub
  [which(sampleinfo_sub$perf == "CIT_end"), ]
  $colID]))
CITendIRI$condition <- "CIT:end"
# NonViable:
# 3 hours:
NonViable3hIRI <- as.data.frame(t(cpmData_sub[
  IRIgenes,
  sampleinfo_sub
  [which(sampleinfo_sub$perf_success
    == "Perf_Unsuccessful"), ]
  [sampleinfo_sub[which(
      sampleinfo_sub$perf_success
      == "Perf_Unsuccessful"
    ), ]$perf
      %in% c("3h_perf", "4h_perf"), ]
  $colID
]))
NonViable3hIRI$condition <- "Non-Viable: 3hours"
# 6 hours:
NonViable6hIRI <- as.data.frame(t(cpmData_sub[
  IRIgenes,
  sampleinfo_sub
  [which(sampleinfo_sub$perf_success
    == "Perf_Unsuccessful"), ]
  [sampleinfo_sub[which(
      sampleinfo_sub$perf_success ==
        "Perf_Unsuccessful"
    ), ]$perf
      %in% c("6h_perf"), ]$colID
]))
NonViable6hIRI$condition <- "Non-Viable: 6hours"
# Viable
# 3 hours:
Viable3hIRI <- as.data.frame(t(cpmData_sub[
  IRIgenes,
  sampleinfo_sub
  [which(sampleinfo_sub$perf_success
    == "Perf_Successful"), ]
  [sampleinfo_sub[
      which(sampleinfo_sub$perf_success
      == "Perf_Successful"),
    ]
    $perf %in% c("3h_perf"), ]$colID
]))
Viable3hIRI$condition <- "Viable: 3hours"
# 6 hours:
Viable6hIRI <- as.data.frame(t(cpmData_sub[
  IRIgenes,
  sampleinfo_sub[
    which(sampleinfo_sub$perf_success
    == "Perf_Successful"),
  ]
  [sampleinfo_sub[
      which(sampleinfo_sub$perf_success
      == "Perf_Successful"),
    ]
    $perf %in% c("6h_perf"), ]$colID
]))
Viable6hIRI$condition <- "Viable: 6hours"
# Bind the dataframe:
IRIdf <- rbind(
  CITstartIRI, CITendIRI, NonViable3hIRI, NonViable6hIRI,
  Viable3hIRI, Viable6hIRI
)
# IRIgenesLonger <- IRIgenes %>% pivot_longer(cols = !condition)
IRIdf$condition <- factor(IRIdf$condition,
  levels = c(
    "CIT:start", "CIT:end",
    "Non-Viable: 3hours", "Non-Viable: 6hours",
    "Viable: 3hours", "Viable: 6hours"
  )
)
IRIdf$colID <- rownames(IRIdf)
IRIdf <- left_join(IRIdf, sampleinfo_sub, by = "colID")
# Mixed model regression for all the IRI genes:
IRImixedModel <- data.frame()
for (i in IRIgenes) {
  designLMM <- paste0(i, "~ (1|sample) +  condition")
  lmm <- lmer(designLMM,
    data = IRIdf, lmerControl(
      optimizer = "nloptwrap",
      calc.derivs = F
    ),
    REML = T
  )
  # sum <- summary(lmm, ddf = "Kenward-Roger")
  sum <- summary(lmm)
  perf <- model_performance(lmm)
  report <- as.data.frame(sum$coefficients)
  report$gene <- i
  report$R2_conditional <- perf$R2_conditional
  report$R2_marginal <- perf$R2_marginal
  IRImixedModel <- rbind(IRImixedModel, report)
}

length(unique(IRImixedModel[is.na(IRImixedModel$R2_conditional), ]$gene))
# Get only the ones that are identified in viable livers:
IRImixedModel.f <- IRImixedModel %>%
  filter(gene %in% fitperfViable) %>%
  mutate(comparison = str_remove(rownames(.), "condition"))



IRImixedModel.f <- IRImixedModel %>%
  filter(!is.na(R2_conditional)) %>%
  mutate(comparison = str_remove(rownames(.), "condition"))


IRImixedModel_f <- IRImixedModel.f[!str_detect(
  IRImixedModel.f$comparison,
  "Intercept"
), ]
IRImixedModel_f$comparison <- gsub("[0-9]+$", "", IRImixedModel_f$comparison)
IRImixedModel_f$comparison <- factor(IRImixedModel_f$comparison,
  levels = c(
    "CIT:end", "Non-Viable: 3hours",
    "Non-Viable: 6hours",
    "Viable: 3hours",
    "Viable: 6hours"
  )
)
IRImixedModel_f$gene <- factor(IRImixedModel_f$gene, levels = IRIgenes)


IRImixedModel_f$success <- unlist(map(str_split(
  IRImixedModel_f$comparison,
  ":"
), 1))
IRImixedModel_f$time <- unlist(map(str_split(
  IRImixedModel_f$comparison,
  ":"
), 2))
IRImixedModel_f$time <- as.integer(factor(IRImixedModel_f$time,
  levels = c(
    "end",
    " 3hours",
    " 6hours"
  )
))
IRImixedModel_f <- rbind(IRImixedModel_f, (IRImixedModel_f %>%
  filter(success == "CIT") %>%
  mutate(success = "Viable")))

IRImixedModel_f <- rbind(IRImixedModel_f, (IRImixedModel_f %>%
  filter(success == "CIT") %>%
  mutate(success = "Non-Viable")))
IRImixedModel_f <- IRImixedModel_f %>% mutate(sign = sign(Estimate))

# Ok, so now, I need to split by up and downregulated:
IRImixedModel_f$IRI <- ""
IRImixedModel_f[which(IRImixedModel_f$gene %in% c(
  "SGPL1", "DLG1", "RLN1", "MLH3",
  "TFAP2A", "TTPA", "LPA", "KLHDC10", "ZNF710", "LRRN3",
  "BLTP2"
)), ]$IRI <- "Downregulated"

IRImixedModel_f[which(IRImixedModel_f$gene %in% c(c(
  "EMP1", "INSIG1", "PPP1R15A",
  "TP53BP2", "MCL1", "WEE1", "GADD45A", "EMP1", "TRAM1",
  "LIFR", "TRAF4", "KLF5", "AKAP12", "MMP19", "PTPN1",
  "TACSTD2", "PPP2R2A", "RORA", "ADM", "MAFF", "TMEM184B",
  "ETS1", "DUSP5", "BAG3", "IL18RAP", "NFKB1", "PLAUR",
  "SERPINE1", "C5AR1", "B4GALT6", "ALAS1", "NABP1",
  "SMIM13", "APOLD1", "TIPARP", "ISG20L2", "CORO1C", "HBB",
  "SLC20A1", "SLC7A5"
))), ]$IRI <- "Upregulated"

IRIplot <- IRImixedModel_f %>%
  dplyr::select(
    success, gene, comparison, time, sign,
    IRI, Estimate
  ) %>%
  filter(success != "CIT")

IRIgenesNotToPlot <- IRIplot %>%
  filter(!between(Estimate, -300, 300)) %>%
  pull(gene) %>%
  unique()

IRIplot <- IRIplot %>% filter(!gene %in% IRIgenesNotToPlot)
IRIplot <- IRIplot %>% filter(gene %in% fitperfViable)

se <- function(y) {
  sd(y) / sqrt(length(y))
}
IRIplot$success <- factor(IRIplot$success, level = c('Viable', 'Non-Viable'))


ggplot(
  IRIplot %>% filter(gene %in% IRItotest),
  aes(x = time, y = Estimate, col = gene)
) +
  facet_wrap(~ IRI + success) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0, lty = 2) +
  theme_minimal() +
  scale_color_manual(values = rep("grey", 40)) +
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
  ylab("Gene expression linear coefficient") -> estimatePlot

ggsave(estimatePlot, filename = paste0(
  oneDriveOutput,
  "Figure1/estimatePlot.pdf"
), width = 7, height = 4)

IRIplot %>%
  filter(gene %in% IRItotest) %>%
  write_csv(, file = paste0(
    oneDriveOutput,
    "Supplementary/IRIgenesCoefficient.csv"
  ))


# Viable values:
IRIplot %>%
  filter(
    gene %in% IRItotest,
    comparison == "Viable: 3hours"
  ) %>%
  filter(IRI == "Upregulated") %>%
  pull(Estimate) -> ViableUp3hours

IRIplot %>%
  filter(
    gene %in% IRItotest,
    comparison == "Viable: 6hours"
  ) %>%
  filter(IRI == "Upregulated") %>%
  pull(Estimate) -> ViableUp6hours
# Non-Viable values:
IRIplot %>%
  filter(
    gene %in% IRItotest,
    comparison == "Non-Viable: 3hours"
  ) %>%
  filter(IRI == "Upregulated") %>%
  pull(Estimate) -> NonViableUp3hours

IRIplot %>%
  filter(
    gene %in% IRItotest,
    comparison == "Non-Viable: 6hours"
  ) %>%
  filter(IRI == "Upregulated") %>%
  pull(Estimate) -> NonViableUp6hours

# Viable:
round(mean(ViableUp3hours), 2)
round(sd(ViableUp3hours), 2)/sqrt(length(ViableUp3hours))
round(mean(ViableUp6hours), 2)
round(sd(ViableUp6hours), 2)/sqrt(length(ViableUp6hours))
# Non-Viable:
round(mean(NonViableUp3hours), 2)
round(sd(NonViableUp3hours), 2)/sqrt(length(NonViableUp3hours))
round(mean(NonViableUp6hours), 2)
round(sd(NonViableUp6hours), 2)/sqrt(length(NonViableUp6hours))


# Sample characteristics:
library(coin)

sampleinfoDist <- sampleinfo_sub %>%
  filter(!duplicated(sample), success %in% c("Successful", "Unsuccessful")) %>%
  mutate(Is.F = ifelse(sex == "F", "TRUE", "FALSE")) %>%
  mutate(Is.DCD = ifelse(death == "DCD", "TRUE", "FALSE")) %>%
  group_by(success) %>%
  dplyr::select(age, Is.F, bmi, Is.DCD, RIN) %>%
  summarise(
    across(where(is.numeric), .fns = list(
      mean = ~ mean(.x, na.rm = T),
      sd = ~ sd(.x, na.rm = T),
      median = ~ median(.x, na.rm = T)
    )),
    across(!where(is.numeric), .fns = list(
      N = ~ length(which(.x == "TRUE")),
      perc = ~ (length(which(.x == "TRUE")) / n()) * 100
    )),
    Total_N = n()
  ) %>%
  mutate(Total_perc = (Total_N / sum(Total_N)) * 100) %>%
  # Round all numbers
  mutate(across(where(is.numeric), ~ round(.x, digits = 1))) %>%
  # Pivot it longer
  pivot_longer(
    cols = !success, names_to = c("Variable", ".value"),
    names_pattern = "^(.+)_([[:alnum:]]+)$"
  ) %>%
  # For the mean and SD, combine
  mutate(Value = case_when(
    is.na(N) ~ paste0(mean, " (", sd, ")"),
    Variable == "Total" ~ as.character(N),
    TRUE ~ as.character(perc)
  )) %>%
  dplyr::select(success, Variable, Value) %>%
  pivot_wider(id_cols = Variable, names_from = success, values_from = Value) %>%
  arrange(Variable)

# Significance df - t-test
sampleinfoStat <- sampleinfo_sub %>%
  filter(!duplicated(sample), success %in% c("Successful", "Unsuccessful")) %>%
  mutate(Is.F = ifelse(sex == "F", "TRUE", "FALSE")) %>%
  mutate(Is.DCD = ifelse(death == "DCD", "TRUE", "FALSE")) %>%
  dplyr::select(success, age, Is.F, bmi, Is.DCD, RIN) %>%
  # Make all the characters into factors
  mutate(across(where(is.character), as.factor)) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        pvalue = ~ pvalue(wilcox_test(.x ~ success)),
        rvalue = (~ statistic(
          wilcox_test(.x ~ success),
          "test"
        ) / sqrt(n())),
        nmissing = ~ sum(is.na(.x))
      )
    ),
    across(
      (!where(is.numeric) & !matches("success")),
      list(
        phivalue = ~ sqrt(statistic(chisq_test(.x ~ success,
          distribution = "exact"
        )) / n()),
        pvalue = ~ pvalue(chisq_test(.x ~ success,
          distribution = "exact"
        )),
        nmissing = ~ sum(is.na(.x))
      )
    )
  ) %>%
  mutate(across(everything(), as.numeric)) %>%
  # Pivot longer
  pivot_longer(everything(),
    names_to = c("Variable", ".value"),
    names_pattern = "^(.+)_([[:alnum:]]+)$"
  ) %>%
  # Add the test type
  mutate(Test = case_when(
    is.na(rvalue) ~ "Chi-squared",
    is.na(phivalue) ~ "Wilcoxon rank-sum (r value)"
  )) %>%
  # Add the BH p-value?
  mutate(BH_pvalue = p.adjust(pvalue, method = "BH")) %>%
  # Format the p-value part
  mutate("P_value" = case_when(
    pvalue < 0.001 ~
      as.character(scales::scientific(pvalue)),
    pvalue >= 0.001 ~
      as.character(signif(pvalue, digits = 1))
  )) %>%
  mutate(BH_pvalue = signif(BH_pvalue, digits = 2)) %>%
  # Format the BH p-vals
  mutate("P_value_BH" = case_when(
    BH_pvalue < 0.01 ~
      as.character(scales::scientific(BH_pvalue)),
    BH_pvalue >= 0.01 ~
      as.character(signif(BH_pvalue, digits = 2))
  )) %>%
  # Format the test statistics
  mutate(Statistic = case_when(
    !is.na(rvalue) ~ round(rvalue, digits = 2),
    !is.na(phivalue) ~ round(phivalue, digits = 2)
  )) %>%
  mutate(N_Missing = nmissing) %>%
  dplyr::select(Variable, P_value, P_value_BH, Statistic, nmissing, Test) %>%
  arrange(Variable)

sampleinfoStat <- left_join(sampleinfoDist, sampleinfoStat, by = "Variable")

# CLock analysis:
circGenes <- c(
  "CLOCK", "ARNTL", "DBP", "NPAS2", "PER1", "PER2", "PER3",
  "TEF", "HLF", "CRY1", "CRY2", "NR1D1"
)
# Mixed linear regression:
cpmData_sub <- as.data.frame(cpm(countData_sub))

CITstartCLOCK <- as.data.frame(t(cpmData_sub
[circGenes, sampleinfo_sub
  [which(sampleinfo_sub$perf == "CIT_start"), ]
  $colID]))
CITstartCLOCK$condition <- "CIT:start"
# End
CITendCLOCK <- as.data.frame(t(cpmData_sub
[circGenes, sampleinfo_sub
  [which(sampleinfo_sub$perf == "CIT_end"), ]
  $colID]))
CITendCLOCK$condition <- "CIT:end"
# NonViable:
# 3 hours:
NonViable3hCLOCK <- as.data.frame(t(cpmData_sub[
  circGenes,
  sampleinfo_sub
  [which(sampleinfo_sub$perf_success
    == "Perf_Unsuccessful"), ]
  [sampleinfo_sub[which(
      sampleinfo_sub$perf_success
      == "Perf_Unsuccessful"
    ), ]$perf
      %in% c("3h_perf", "4h_perf"), ]
  $colID
]))
NonViable3hCLOCK$condition <- "Non-Viable: 3hours"
# 6 hours:
NonViable6hCLOCK <- as.data.frame(t(cpmData_sub[
  circGenes,
  sampleinfo_sub
  [which(sampleinfo_sub$perf_success
    == "Perf_Unsuccessful"), ]
  [sampleinfo_sub[which(
      sampleinfo_sub$perf_success ==
        "Perf_Unsuccessful"
    ), ]$perf
      %in% c("6h_perf"), ]$colID
]))
NonViable6hCLOCK$condition <- "Non-Viable: 6hours"
# Viable
# 3 hours:
Viable3hCLOCK <- as.data.frame(t(cpmData_sub[
  circGenes,
  sampleinfo_sub
  [which(sampleinfo_sub$perf_success
    == "Perf_Successful"), ]
  [sampleinfo_sub[
      which(sampleinfo_sub$perf_success
      == "Perf_Successful"),
    ]
    $perf %in% c("3h_perf"), ]$colID
]))
Viable3hCLOCK$condition <- "Viable: 3hours"
# 6 hours:
Viable6hCLOCK <- as.data.frame(t(cpmData_sub[
  circGenes,
  sampleinfo_sub[
    which(sampleinfo_sub$perf_success
    == "Perf_Successful"),
  ]
  [sampleinfo_sub[
      which(sampleinfo_sub$perf_success
      == "Perf_Successful"),
    ]
    $perf %in% c("6h_perf"), ]$colID
]))
Viable6hCLOCK$condition <- "Viable: 6hours"
# Bind the dataframe:
CLOCKdf <- rbind(
  CITstartCLOCK, CITendCLOCK, NonViable3hCLOCK, NonViable6hCLOCK,
  Viable3hCLOCK, Viable6hCLOCK
)
# IRIgenesLonger <- IRIgenes %>% pivot_longer(cols = !condition)
CLOCKdf$condition <- factor(CLOCKdf$condition,
  levels = c(
    "CIT:start", "CIT:end",
    "Non-Viable: 3hours", "Non-Viable: 6hours",
    "Viable: 3hours", "Viable: 6hours"
  )
)
CLOCKdf$colID <- rownames(CLOCKdf)
CLOCKdf <- left_join(CLOCKdf, sampleinfo_sub, by = "colID")
# Mixed model regression for all the IRI genes:
CLOCKmixedModel <- data.frame()

for (i in circGenes) {
  designLMM <- paste0(i, "~ (1|sample) +  condition")
  lmm <- lmer(designLMM, data = CLOCKdf)
  sum <- summary(lmm, ddf = "Kenward-Roger")
  perf <- model_performance(lmm)
  report <- as.data.frame(sum$coefficients)
  report$gene <- i
  report$R2_conditional <- perf$R2_conditional
  report$R2_marginal <- perf$R2_marginal
  CLOCKmixedModel <- rbind(CLOCKmixedModel, report)
}

CLOCKmixedModel_f <- CLOCKmixedModel %>%
  mutate(comparison = str_remove(rownames(.), "condition"))
CLOCKmixedModel_f <- CLOCKmixedModel_f[!str_detect(
  CLOCKmixedModel_f$comparison,
  "Intercept"
), ]
CLOCKmixedModel_f$comparison <- gsub(
  "[0-9]+$", "",
  CLOCKmixedModel_f$comparison
)
CLOCKmixedModel_f$comparison <- factor(CLOCKmixedModel_f$comparison,
  levels = c(
    "CIT:end", "Non-Viable: 3hours",
    "Non-Viable: 6hours",
    "Viable: 3hours",
    "Viable: 6hours"
  )
)
CLOCKmixedModel_f$gene <- factor(CLOCKmixedModel_f$gene, levels = circGenes)


CLOCKmixedModel_f$success <- unlist(map(str_split(
  CLOCKmixedModel_f$comparison,
  ":"
), 1))
CLOCKmixedModel_f$time <- unlist(map(str_split(
  CLOCKmixedModel_f$comparison,
  ":"
), 2))
CLOCKmixedModel_f$time <- as.integer(factor(CLOCKmixedModel_f$time,
  levels = c(
    "end",
    " 3hours",
    " 6hours"
  )
))
CLOCKmixedModel_f <- rbind(CLOCKmixedModel_f, (CLOCKmixedModel_f %>%
  filter(success == "CIT") %>%
  mutate(success = "Viable")))

CLOCKmixedModel_f <- rbind(CLOCKmixedModel_f, (CLOCKmixedModel_f %>%
  filter(success == "CIT") %>%
  mutate(success = "Non-Viable")))
CLOCKmixedModel_f <- CLOCKmixedModel_f %>% mutate(sign = sign(Estimate))

# Split genes based on day and night:
CLOCKmixedModel_f$CLOCK <- "Night"

CLOCKmixedModel_f[CLOCKmixedModel_f$gene %in% c(
  "PER1", "PER2", "PER3",
  "TEF", "HLF", "CRY1",
  "CRY2", "NR1D1"
), ]$CLOCK <- "Day"

# Night genes:
c("CLOCK", "ARNTL", "DBP", "NPAS2")
# Day genes
c("PER1", "PER2", "PER3", "TEF", "HLF", "CRY1", "CRY2", "NR1D1")

CLOCKplot <- CLOCKmixedModel_f %>%
  dplyr::select(
    success, gene, comparison, time, sign, Estimate, CLOCK
  ) %>%
  filter(success != "CIT")

CLOCKmixedModel_f %>% filter(`Pr(>|t|)` < 0.05)
ggplot(
  CLOCKplot %>% filter(gene %in% c("PER1", "PER2", "NR1D1")),
  aes(x = time, y = Estimate, col = gene)
) +
  facet_wrap(~ gene + success, nrow = 3) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0, lty = 2) +
  theme_minimal() +
  scale_color_manual(values = rep("grey", 40)) +
  theme(legend.position = "none") +
  stat_summary(
    fun.y = mean, fun.ymin = function(y) mean(y) - sd(y),
    fun.ymax = function(y) mean(y) + sd(y), color = "dark red",
    geom = "pointrange", show.legend = FALSE
  ) +
  xlab("Perfusion Time") +
  ylab("Gene expression linear coefficient") -> clockPlot

ggsave(clockPlot, filename = paste0(
  oneDriveOutput,
  "Figure1/CLOCKestimatePlot.pdf"
), width = 7, height = 7)

ggplot(
  CLOCKplot,
  aes(x = time, y = Estimate, col = gene)
) +
  facet_wrap(~ gene + success, nrow = 3) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0, lty = 2) +
  theme_minimal() +
  scale_color_manual(values = rep("grey", 40)) +
  theme(legend.position = "none") +
  stat_summary(
    fun.y = mean, fun.ymin = function(y) mean(y) - sd(y),
    fun.ymax = function(y) mean(y) + sd(y), color = "dark red",
    geom = "pointrange", show.legend = FALSE
  ) +
  xlab("Perfusion Time") +
  ylab("Gene expression linear coefficient") -> allClockPlot

ggsave(allClockPlot, filename = paste0(
  oneDriveOutput,
  "Figure1/CLOCKallGenesestimatePlot.pdf"
), width = 14, height = 14)


# CCD analysis:
library(nCV)
sampleinfo_sub[str_detect(sampleinfo_sub$sample, "nmp4"), ]$perf_success[2:3] <- "Perf_Unsuccessful"

# Normalised connection of variation:
homemadenCVnet <- function(inputD, benchD, nperm = 1000, seedN = 10) {
  inputD <- as.data.frame(inputD)
  benchD <- as.data.frame(benchD)
  colnames(inputD)[1] <- colnames(benchD)[1] <- "geneSym"
  if (nrow(inputD) == length(unique(inputD$geneSym))) {
    rownames(inputD) <- inputD$geneSym
    rownames(benchD) <- benchD$geneSym
    bothID <- inner_join(data.frame(id = benchD$geneSym),
      data.frame(id = inputD$geneSym),
      by = "id"
    )

    corR <- as.matrix(benchD[bothID$id, bothID$id])
    corD <- cor(t(inputD[bothID$id, -1]), method = "spearman")
    rownames(corD) <- colnames(corD) <- bothID$id
    simD <- ape::mantel.test(corD, corR, nperm = 1)
    set.seed(seedN)
    indexM <- NULL
    for (np in 1:nperm) {
      indexM <- rbind(indexM, sample(
        1:nrow(inputD),
        nrow(bothID)
      ))
    }
    pstat <- apply(indexM, 1, function(indexv) {
      corP <- cor(t(inputD[indexv, -1]), method = "spearman")
      tepD <- ape::mantel.test(corP, corR, nperm = 1)
      return(tepD$z.stat)
    })
    pva <- sum(pstat > simD$z.stat) / nperm
    return(list(
      zstat = data.frame(
        tag = "nCVnet", geneNumber = nrow(bothID),
        zstat = simD$z.stat, pvalue = pva
      ), npermV = pstat,
      cmatrix = corD
    ))
  }
}

cpmData <- cpm(countData_sub)
var <- apply(cpmData, 1, var)
cpmData <- cpmData[which(var > 1), ]
rownames(cpmData) <- str_to_title(rownames(cpmData))
cpmData <- as.data.frame(cpmData)
cpmData$geneSym <- rownames(cpmData)
cpmData <- cpmData[, c("geneSym", colnames(cpmData)[-45])]
# CIT:
# Successful:
# Failed:
cpmCIT <- cpmData[, -1][, c(which(
  sampleinfo_sub$perf_success == "CIT"
))]
cpmSuccess <- cpmData[, -1][, c(which(
  sampleinfo_sub$perf_success == "Perf_Successful"
))]
cpmFailed <- cpmData[, -1][, c(which(
  sampleinfo_sub$perf_success == "Perf_Unsuccessful"
))]
# Switch NMP4 from part:
cpmFailed <- cbind(cpmFailed, cpmSuccess[, str_detect(
  colnames(cpmSuccess), "NMP4"
)])
cpmSuccess <- cpmSuccess[, !str_detect(
  colnames(cpmSuccess), "NMP4"
)]

cpmCIT <- cbind(cpmData[, 1], cpmCIT)
cpmSuccess <- cbind(cpmData[, 1], cpmSuccess)
cpmFailed <- cbind(cpmData[, 1], cpmFailed)

mClockD <- nCV::mClockD
nperm <- 10000
nCVcit <- homemadenCVnet(cpmCIT, mClockD, nperm = nperm, seedN = 154)
nCVfail <- homemadenCVnet(cpmFailed, mClockD, nperm = nperm, seedN = 154)
nCVsuccess <- homemadenCVnet(cpmSuccess, mClockD, nperm = nperm, seedN = 154)
# nCVGenes:

cgenes <- str_to_title(c(
  "ARNTL", "CLOCK", "NPAS2", "CRY1",
  "NR1D1", "CIART", "DBP", "PER1", "CRY2", "PER2"
))
nCIT <- nCVgene(inputD = cpmCIT, cgenes)
nFail <- nCVgene(inputD = cpmFailed, cgenes)
nSuccess <- nCVgene(inputD = cpmSuccess, cgenes)

mean(nCIT$nCV)
mean(nSuccess$nCV)
mean(nFail$nCV)

# CCD:
# Grab the GTEx vector:
CLOCK_GTEx <- read.table("../../CLOCK/Results/CLOCK/GTEx/GTExMatrixRef.csv")
matrix_levels <- c(
  "ARNTL", "NPAS2", "CLOCK", "NFIL3", "CRY1",
  "CRY2", "NR1D1", "NR1D2", "PER1", "PER2",
  "PER3", "DBP", "TEF", "HLF"
)
matrix_levels <- str_to_title(matrix_levels)
# Correlation vector:
CLOCK_GTEx <- CLOCK_GTEx[order(rownames(CLOCK_GTEx)), ]
rownames(CLOCK_GTEx) <- str_to_title(rownames(CLOCK_GTEx))
GTEXcorMatrix <- cor(t(CLOCK_GTEx))
# Get the correlation matrix for each of the experimental conditions:
# Order the dataframe in the same way than the GTEx matrix
cpmCITordered <- cpmCIT[order(rownames(cpmCIT)), -1]
cpmFailedOrdered <- cpmFailed[order(rownames(cpmFailed)), -1]
cpmSuccessOrdered <- cpmSuccess[order(rownames(cpmSuccess)), -1]
# Create the correlation matrix:
CITcorMatrix <- cor(t(cpmCITordered[which(rownames(cpmCITordered) %in%
  matrix_levels), ]))
# Viable livers:
FailedCorMatrix <- cor(t(cpmFailedOrdered[which(rownames(cpmFailedOrdered) %in%
  matrix_levels), ]))
# Non-viable Livers:
SuccessCorMatrix <- cor(t(cpmSuccessOrdered[which(
  rownames(cpmSuccessOrdered) %in% matrix_levels
), ]))
# Extract the correlation vector for each:
GTEXcorVector <- GTEXcorMatrix[lower.tri(GTEXcorMatrix, diag = F)]
CITcorVector <- CITcorMatrix[lower.tri(CITcorMatrix, diag = F)]
FailedCorVector <- FailedCorMatrix[lower.tri(FailedCorMatrix, diag = F)]
SuccessCorVector <- SuccessCorMatrix[lower.tri(SuccessCorMatrix, diag = F)]
# Calculate the euclidean distance between the two:
euclidean <- function(a, b) sqrt(sum((a - b)^2))
CCDCIT <- euclidean(GTEXcorVector, CITcorVector)
CCDFailed <- euclidean(GTEXcorVector, FailedCorVector)
CCDSuccess <- euclidean(GTEXcorVector, SuccessCorVector)
# Get the standard error for each value:
se <- function(y) {sd(y) / sqrt(length(y))}
se(CITcorVector)
se(FailedCorVector)
se(SuccessCorVector)

# Calculate the distribution of the distance within each condition.
bootstrapDistribution <- function(cpm, clockGenes, n = 1000, refCorVector) {
  CCD <- c()
  while (length(CCD) < n) {
    meanVector <- rowMeans(cpm[which(rownames(cpm) %in% clockGenes), ])
    # Randomly sample a set of the same size as the clock:
    cpmSamplingSet <- cpm[-which(rownames(cpm) %in% clockGenes), ]
    set <- sample(rownames(cpmSamplingSet), size = length(clockGenes))
    # Calculate the mean expression of those genes:
    setMeanVector <- rowMeans(cpmSamplingSet[set, ])
    # If the genes have a similar mean to the CLOCK genes, accept the set:
    if (t.test(meanVector, setMeanVector)$p.value > 0.05) {
      # Accept the set and calculate the CCD:
      SetCorMatrix <- cor(t(cpmSamplingSet[set, ]))
      SetCorVector <- SetCorMatrix[lower.tri(SetCorMatrix, diag = F)]
      CCD <- c(CCD, euclidean(refCorVector, SetCorVector))
    }
  }
  return(CCD)
}


CCDdistribCIT <- bootstrapDistribution(
  cpm = cpmCITordered,
  clockGenes = matrix_levels,
  n = nperm, refCorVector = GTEXcorVector
)
CCDdisribFailed <- bootstrapDistribution(
  cpm = cpmFailedOrdered,
  clockGenes = matrix_levels,
  n = nperm, refCorVector = GTEXcorVector
)
CCDdistribSuccess <- bootstrapDistribution(
  cpm = cpmSuccessOrdered,
  clockGenes = matrix_levels,
  n = nperm, refCorVector = GTEXcorVector
)
# Calculate the p-value:
pCIT <- pnorm(
  q = CCDCIT, mean = mean(CCDdistribCIT),
  sd = sd(CCDdistribCIT)
)
pFail <- pnorm(
  q = CCDFailed, mean = mean(CCDdisribFailed),
  sd = sd(CCDdisribFailed)
)
pSuccess <- pnorm(
  q = CCDSuccess, mean = mean(CCDdistribSuccess),
  sd = sd(CCDdistribSuccess)
)
pval <- p.adjust(c(pCIT, pSuccess, pFail), method = "bonfe")
pval <- format(signif(pval, 3), scientific = T)


CCDdataFrame <- data.frame(
  distance = c(CCDdistribCIT, CCDdisribFailed, CCDdistribSuccess),
  condition = c(
    rep("CIT", nperm),
    rep("Non-Viable", nperm),
    rep("Viable", nperm)
  )
)

title <- paste0(
  "CCD CIT: ", round(CCDCIT, 2),
  "; P-value:", pval[1], "\n",
  "CCD Viable: ", round(CCDSuccess, 2),
  "; P-value:", pval[2], "\n",
  "CCD Non-Viable: ", round(CCDFailed, 2),
  "; P-value:", pval[3], "\n"
)


library(ggpattern)

ggplot(CCDdataFrame, aes(x = distance)) +
  geom_density(aes(fill = condition), alpha = 0.5, col = "NA") +
  geom_vline(xintercept = CCDCIT, col = "#BEDEE9", linewidth = 1) +
  geom_vline(xintercept = CCDFailed, col = "#A6BCA2", linewidth = 1) +
  geom_vline(xintercept = CCDSuccess, col = "#F0C56F", linewidth = 1) +
  scale_fill_manual(
    values = c("#BEDEE9", "#A6BCA2", "#F0C56F"),
    name = "Empirical Distribution for:"
  ) +
  theme_minimal() +
  xlab("Clock Correlation Distance") +
  xlim(0, 7.5) +
  ggtitle(title) -> ccdPlot

ggsave(ccdPlot, filename = paste0(
  oneDriveOutput,
  "Figure1/CCDplot.pdf"
), width = 7, height = 4)

# Metabolite plots:
create_df <- function(nmp_df, sample) {
  temp <- as.data.frame(t(nmp_df))
  colnames(temp) <- temp[1, ]
  temp <- temp[-1, ]
  temp$sample <- sample
  return(temp)
}

breaks <- seq(-100, 350, by = 50)
midpoints <- function(x, dp = 2) {
  lower <- as.numeric(gsub(",.*", "", gsub("\\(|\\[|\\)|\\]", "", x)))
  upper <- as.numeric(gsub(".*,", "", gsub("\\(|\\[|\\)|\\]", "", x)))
  return(round(lower + (upper - lower) / 2, dp))
}



# Successful NMP:
metPath <- "/Users/uqschauq/Documents/NMP/Data/Metabolites/"
nmp2 <- create_df(read_tsv(
  file = paste0(metPath, "nmp2.tsv"),
  col_names = F
), "nmp2")
nmp5 <- create_df(read_tsv(
  file = paste0(metPath, "nmp5.tsv"),
  col_names = F
), "nmp5")
nmp11 <- create_df(read_tsv(
  file = paste0(metPath, "nmp11.tsv"),
  col_names = F
), "nmp11")
nmp12 <- create_df(read_tsv(
  file = paste0(metPath, "nmp12.tsv"),
  col_names = F
), "nmp12")
nmp15 <- create_df(read_tsv(
  file = paste0(metPath, "nmp15.tsv"),
  col_names = F
), "nmp15")
nmp16 <- create_df(read_tsv(
  file = paste0(metPath, "nmp16.tsv"),
  col_names = F
), "nmp16")
# Unsuccessful NMP:
nmp3 <- create_df(read_tsv(
  file = paste0(metPath, "nmp3.tsv"),
  col_names = F
), "nmp3")
nmp4 <- create_df(read_tsv(
  file = paste0(metPath, "nmp4.tsv"),
  col_names = F
), "nmp4")
nmp8 <- create_df(read_tsv(
  file = paste0(metPath, "nmp8.tsv"),
  col_names = F
), "nmp8")
nmp9 <- create_df(read_tsv(
  file = paste0(metPath, "nmp9.tsv"),
  col_names = F
), "nmp9")
# Columns names:
nmp2$perf <- difftime(
  time1 = ymd_hms(paste0(
    "2022-02-10 ",
    nmp2$`Time:`, ":00"
  )),
  time2 = ymd_hms("2022-02-10 10:00:00")
)

nmp3$perf <- difftime(
  time1 = ymd_hms(
    paste0(
      "2022-02-10 ",
      nmp3$`Time:`[1:10], ":00"
    ),
    paste0(
      "2022-02-11 ",
      nmp3$`Time:`[11:13], ":00"
    )
  ),
  time2 = ymd_hms("2022-02-10 21:25:00")
)

nmp4$perf <- difftime(
  time1 = ymd_hms(
    paste0(
      "2022-02-10 ",
      nmp4$`Time:`[1], ":00"
    ),
    paste0(
      "2022-02-11 ",
      nmp4$`Time:`[2:16], ":00"
    )
  ),
  time2 = ymd_hms("2022-02-11 01:10:00")
)

nmp5$perf <- difftime(
  time1 = ymd_hms(c(
    paste0(
      "2022-02-10 ",
      nmp5$`Time:`[1:5], ":00"
    ),
    paste0(
      "2022-02-11 ",
      nmp5$`Time:`[6:72], ":00"
    )
  )),
  time2 = ymd_hms("2022-02-10 09:00:00")
)

nmp8$perf <- difftime(
  time1 = ymd_hms(
    paste0(
      "2022-02-10 ",
      nmp8$`Time:`[1:19], ":00"
    ),
    paste0(
      "2022-02-11 ",
      nmp8$`Time:`[20:22], ":00"
    )
  ),
  time2 = ymd_hms("2022-02-10 18:55:00")
)

nmp9$perf <- difftime(
  time1 = ymd_hms(paste0(
    "2022-02-10 ",
    nmp9$`Time:`, ":00"
  )),
  time2 = ymd_hms("2022-02-10 17:30:00")
)

nmp11$perf <- difftime(
  time1 = ymd_hms(
    paste0(
      "2022-02-10 ",
      nmp11$`Time:`[1:2], ":00"
    ),
    paste0(
      "2022-02-11 ",
      nmp11$`Time:`[3:17], ":00"
    )
  ),
  time2 = ymd_hms("2022-02-11 00:10:00")
)

nmp12$perf <- difftime(
  time1 = ymd_hms(paste0(
    "2022-02-10 ",
    nmp12$`Time:`, ":00"
  )),
  time2 = ymd_hms("2022-02-10 08:05:00")
)

nmp15$perf <- difftime(
  time1 = ymd_hms(
    paste0(
      "2022-02-10 ",
      nmp15$`Time:`[1:4], ":00"
    ),
    paste0(
      "2022-02-11 ",
      nmp15$`Time:`[5:26], ":00"
    )
  ),
  time2 = ymd_hms("2022-02-11 00:45:00")
)

nmp16$perf <- difftime(
  time1 = ymd_hms(
    paste0(
      "2022-02-10 ",
      nmp16$`Time:`[1:18], ":00"
    ),
    paste0(
      "2022-02-11 ",
      nmp16$`Time:`[19:20], ":00"
    )
  ),
  time2 = ymd_hms("2022-02-10 18:45:00")
)

col <- c(
  "Serie", "Time", "pH", "pCO2", "pO2", "HCO3", "BE", "Hct",
  "tHb", "sO2", "FO2Hb", "FCOHb", "FmetHb", "FHHb", "Na", "K",
  "Ca", "Cl", "AnGap", "Glu", "Lac", "sample", "perf"
)

breaks <- seq(-100, 350, by = 50)
# Successful NMP:
metaSuccessful <- rbind(nmp2, nmp5, nmp11, nmp12, nmp15, nmp16, nmp4)
colnames(metaSuccessful) <- col
metaSuccessful$Glu <- as.double(metaSuccessful$Glu)
metaSuccessful$Lac <- as.double(metaSuccessful$Lac)
metaSuccessful$bin <- cut(as.double(metaSuccessful$perf), breaks)
metaSuccessful$midpoint <- midpoints(metaSuccessful$bin)
metaSuccessful$NMP <- "Successful"
metaSuccessful$pH <- as.double(metaSuccessful$pH)
# Unsuccessful NMP:
metaUnsuccessful <- rbind(nmp3, nmp8, nmp9)
colnames(metaUnsuccessful) <- col
metaUnsuccessful$Glu <- as.double(metaUnsuccessful$Glu)
metaUnsuccessful$Lac <- as.double(metaUnsuccessful$Lac)
metaUnsuccessful$bin <- cut(as.double(metaUnsuccessful$perf), breaks)
metaUnsuccessful$midpoint <- midpoints(metaUnsuccessful$bin)
metaUnsuccessful$NMP <- "Not Successful"
metaUnsuccessful$pH <- as.double(metaUnsuccessful$pH)
metaFull <- rbind(metaSuccessful, metaUnsuccessful)
metaFull$NMP <- factor(metaFull$NMP, levels = c("Successful", "Not Successful"))
# Glucose Plot:
ggplot(metaFull %>% group_by(bin, NMP) %>% mutate(
  mglu = mean(Glu, na.rm = T),
  sdglu = sd(Glu, na.rm = T)
)) +
  geom_point(aes(x = perf, y = Glu, color = sample), alpha = .4) +
  geom_line(aes(x = perf, y = Glu, color = sample), alpha = .4) +
  geom_pointrange(aes(
    x = midpoint, y = mglu, ymin = mglu - sdglu,
    ymax = mglu + sdglu
  )) +
  geom_line(aes(x = midpoint, y = mglu)) +
  geom_vline(aes(xintercept = 0), lty = "dashed", col = "dark red") +
  theme_linedraw() +
  xlim(-100, 400) +
  ylim(0, 50) +
  xlab("Perfusion time (minutes)") +
  ylab("Glucose (mmol/L)") +
  facet_wrap(~NMP) +
  scale_color_manual(values = rep("grey", 10)) +
  theme(legend.position = "none") -> GluPlot
# Lactate Plot:
ggplot(metaFull %>% group_by(bin, NMP) %>% mutate(
  mlac = mean(Lac),
  sdlac = sd(Lac)
)) +
  geom_point(aes(x = perf, y = Lac, color = sample), alpha = .4) +
  geom_line(aes(x = perf, y = Lac, color = sample), alpha = .4) +
  geom_pointrange(aes(
    x = midpoint, y = mlac, ymin = mlac - sdlac,
    ymax = mlac + sdlac
  )) +
  geom_line(aes(x = midpoint, y = mlac)) +
  geom_vline(aes(xintercept = 0), lty = "dashed", col = "dark red") +
  theme_linedraw() +
  xlim(-100, 400) +
  xlab("Perfusion time (minutes)") +
  ylab("Lactose (mmol/L)") +
  scale_color_manual(values = rep("grey", 10)) +
  facet_wrap(~NMP) +
  geom_hline(yintercept = 2.8, col = "lightgreen") +
  geom_hline(yintercept = 2.5, col = "darkred") +
  geom_hline(yintercept = 1.7, col = "darkblue") +
  theme(legend.position = "none") -> LacPlot
# ph plot:
ggplot(metaFull %>% group_by(bin, NMP) %>% mutate(
  mglu = mean(pH, na.rm = T),
  sdglu = sd(pH, na.rm = T)
)) +
  geom_point(aes(x = perf, y = pH, color = sample), alpha = .4) +
  geom_line(aes(x = perf, y = pH, color = sample), alpha = .4) +
  geom_pointrange(aes(
    x = midpoint, y = mglu, ymin = mglu - sdglu,
    ymax = mglu + sdglu
  )) +
  geom_line(aes(x = midpoint, y = mglu)) +
  geom_vline(aes(xintercept = 0), lty = "dashed", col = "dark red") +
  theme_linedraw() +
  xlim(-100, 400) +
  ylim(6.5, 8) +
  xlab("Perfusion time (minutes)") +
  ylab("pH") +
  facet_wrap(~NMP) +
  scale_color_manual(values = rep("grey", 10)) +
  geom_hline(yintercept = 7.2, col = "lightgreen") +
  geom_hline(yintercept = 7.3, col = "darkred") +
  geom_hline(yintercept = 7.35, col = "darkblue") +
  theme(legend.position = "none") -> phPlot


plot <- cowplot::plot_grid(phPlot, LacPlot, GluPlot, ncol = 1)
ggsave(plot,
  filename = "/Users/uqschauq/Library/CloudStorage/OneDrive-TheUniversityofQueensland/NMP_liver/Results/MetabolitePlot.pdf", width = 8, height = 11
)


ggplot(metaFull %>% filter(sample == "nmp4") %>% group_by(bin) %>%
  mutate(mglu = mean(pH, na.rm = T), sdglu = sd(pH, na.rm = T))) +
  geom_point(aes(x = perf, y = pH, color = sample), alpha = .4) +
  geom_line(aes(x = perf, y = pH, color = sample), alpha = .4) +
  geom_pointrange(aes(
    x = midpoint, y = mglu, ymin = mglu - sdglu,
    ymax = mglu + sdglu
  )) +
  geom_line(aes(x = midpoint, y = mglu)) +
  geom_vline(aes(xintercept = 0), lty = "dashed", col = "dark red") +
  theme_linedraw() +
  xlim(-100, 400) +
  ylim(6.5, 8) +
  xlab("Perfusion time (minutes)") +
  ylab("pH") +
  facet_wrap(~NMP) +
  scale_color_manual(values = rep("grey", 10)) +
  geom_hline(yintercept = 7.2, col = "lightgreen") +
  geom_hline(yintercept = 7.3, col = "darkred") +
  geom_hline(yintercept = 7.35, col = "darkblue") +
  theme(legend.position = "none") -> nmp4phPlot

ggplot(metaFull %>% filter(sample == "nmp4") %>% group_by(bin) %>%
  mutate(mlac = mean(Lac), sdlac = sd(Lac))) +
  geom_point(aes(x = perf, y = Lac, color = sample), alpha = .4) +
  geom_line(aes(x = perf, y = Lac, color = sample), alpha = .4) +
  geom_pointrange(aes(
    x = midpoint, y = mlac, ymin = mlac - sdlac,
    ymax = mlac + sdlac
  )) +
  geom_line(aes(x = midpoint, y = mlac)) +
  geom_vline(aes(xintercept = 0), lty = "dashed", col = "dark red") +
  theme_linedraw() +
  xlim(-100, 400) +
  xlab("Perfusion time (minutes)") +
  ylab("Lactose (mmol/L)") +
  scale_color_manual(values = rep("grey", 10)) +
  facet_wrap(~NMP) +
  geom_hline(yintercept = 2.8, col = "lightgreen") +
  geom_hline(yintercept = 2.5, col = "darkred") +
  geom_hline(yintercept = 1.7, col = "darkblue") +
  theme(legend.position = "none") -> npm4LacPlot

ggplot(metaFull %>% filter(sample == "nmp4") %>% group_by(bin) %>%
  mutate(mglu = mean(Glu, na.rm = T), sdglu = sd(Glu, na.rm = T))) +
  geom_point(aes(x = perf, y = Glu, color = sample), alpha = .4) +
  geom_line(aes(x = perf, y = Glu, color = sample), alpha = .4) +
  geom_pointrange(aes(
    x = midpoint, y = mglu, ymin = mglu - sdglu,
    ymax = mglu + sdglu
  )) +
  geom_line(aes(x = midpoint, y = mglu)) +
  geom_vline(aes(xintercept = 0), lty = "dashed", col = "dark red") +
  theme_linedraw() +
  xlim(-100, 400) +
  ylim(0, 50) +
  xlab("Perfusion time (minutes)") +
  ylab("Glucose (mmol/L)") +
  facet_wrap(~NMP) +
  scale_color_manual(values = rep("grey", 10)) +
  theme(legend.position = "none") -> nmp4GluPlot

### Figure:

viableDownsampleResults <- read_rds("/Users/uqschauq/Documents/NMP/Results/Figure 1/viableDownsampleResults.rds")

fitperfViable <- read.csv(paste0(
  oneDriveOutput,
  "Supplementary/CIT_viable_significant_DEG.csv"
)) %>% pull(gene)

fitperfNonViable <- read.csv(paste0(
  oneDriveOutput,
  "Supplementary/CIT_NonViable_significant_DEG.csv"
)) %>% pull(gene)

# Get the number of genes:
downSamplePlotDF <- data.frame(downSample = unlist(
  lapply(viableDownsampleResults, nrow)
)) %>%
  summarise(
    sd = sd(downSample),
    mean = mean(downSample),
    group = "Downsampled"
  )

# Number of genes:
ggplot(downSamplePlotDF) +
  geom_bar(aes(x = group, y = mean),
    stat = "identity",
    fill = VolColorViable[1], alpha = 0.9
  ) +
  geom_errorbar(aes(x = group, ymin = mean - sd, ymax = mean + sd),
    colour = "black"
  ) +
  geom_bar(aes(x = group, y = mean),
    stat = "identity",
    fill = VolColorViable[2],
    data = data.frame(
      group = "Viable",
      mean = 2413
    )
  ) +
  geom_bar(aes(x = group, y = mean),
    stat = "identity",
    fill = VolColorNonViable[2],
    data = data.frame(
      group = "Non-Viable",
      mean = 572
    )
  ) +
  theme_minimal() +
  scale_x_discrete(limits = c("Viable", "Downsampled", "Non-Viable")) +
  xlab("") +
  ylab("Number of DEGs")



# Calculate correlation between the downsampled genes:
unlist(lapply(viableDownsampleResults, function(x) {
  temp <- left_join(x, fitperfViable, by = "gene") %>%
    select(gene, logFC.x, logFC.y) %>%
    filter(!is.na(logFC.y))
  cor(temp$logFC.x, temp$logFC.y)
})) -> downSampleCor
