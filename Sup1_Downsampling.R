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
library(data.table)

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
                                     c("CIT", "Perf_Unsuccessful", 
                                       "Perf_Successful")), ]
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
##################################
######### Downsample analysis DEG:
##################################
### Downsampling:
# sample selection:
# Get all the different
sampleinfo_downsampled <- rbind(
  sampleinfo_sub %>% filter(
    perf_success == "Perf_Unsuccessful",
    sample != "nmp4"
  ),
  sampleinfo_sub %>% filter(
    perf_success == "CIT",
    sample != "nmp4"
  ),
  sampleinfo_sub %>% filter(perf_success == "Perf_Successful") %>%
    sample_n(5)
)

successSampleToUse <- combn(sampleinfo_sub %>%
                              filter(perf_success == "Perf_Successful") %>%
                              pull(sample) %>% unique(), 3, simplify = F)

# Get the genes removed from the first analysis:
isexpr <- rowSums(cpm(countData_sub) > 0.1) >= 5
viableDownsampleResults <- list()
for (i in successSampleToUse) {
  sampleinfo_downsampled <- rbind(
    sampleinfo_sub %>% filter(
      perf_success == "Perf_Unsuccessful",
      sample != "nmp4"
    ),
    sampleinfo_sub %>% filter(
      perf_success == "CIT",
      sample != "nmp4"
    ),
    sampleinfo_sub %>% filter(
      perf_success == "Perf_Successful",
      sample %in% i
    )
  )
  countData_downsampled <- countData_sub[sampleinfo_downsampled$colID]
  
  dge <- DGEList(
    counts = countData_downsampled, samples = sampleinfo_downsampled,
    genes = rownames(countData)
  )
  # Double check everything is ordered correctly
  all(colnames(dge$counts) == rownames(dge$samples))
  # Create the experimental design (random effect, all known covariates
  # and surrogate variables:)
  design <- ~ (1 | sample) + sex + bmi + age + death + RIN + SV1 +
    SV2 + SV3 + perf_success
  # Standard usage of limma/voom
  dge <- DGEList(countData_downsampled[isexpr, ])
  dge <- calcNormFactors(dge)
  vobjDream <- voomWithDreamWeights(dge, design, sampleinfo_downsampled)
  fitmm <- dream(vobjDream, design, sampleinfo_downsampled)
  fitmm <- eBayes(fitmm)
  
  viableResultsDownsampled <- topTable(fitmm,
                                       coef = "perf_successPerf_Successful",
                                       number = 300000
  )
  
  fitperfViableDownsampled <- viableResultsDownsampled %>%
    filter(adj.P.Val < 0.05, !between(logFC, -0.25, 0.25)) %>%
    mutate(gene = rownames(.))
  
  viableDownsampleResults <- list.append(
    viableDownsampleResults,
    fitperfViableDownsampled
  )
}

write_rds(viableDownsampleResults,
          file = "/Users/uqschauq/Documents/NMP/Results/Figure 1/viableDownsampleResults.rds"
)

##################################
######### Downsample analysis IRI:
##################################
IRIinitialAnalysis <- fread(paste0(oneDriveOutput,
             "Supplementary/IRIgenesCoefficient.csv"))
# Check how many genes are enriched within DEGs sets:
IRIgenes <- c("EMP1", "TRAF4", "MMP19", "PTPN1", "PLAUR", "C5AR1")
IRItotest <- unique(IRIinitialAnalysis$Genes)

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
downSampledIRIanalysis <- list()
fold <- 1
for (set in successSampleToUse){
  # Downsample 3 hours:
  downSampleViable3h <- sampleinfo_sub %>% filter(
    perf_success == "Perf_Successful",
    sample %in% set,
    perf %in% c("3h_perf", "4h_perf")) %>% 
    pull(colID)
  Viable3hIRI <- as.data.frame(t(cpmData_sub[IRIgenes,downSampleViable3h]))
  Viable3hIRI$condition <- "Viable: 3hours"
  # downsample 6 hours:
  downSampleViable6h <- sampleinfo_sub %>% filter(
    perf_success == "Perf_Successful",
    sample %in% set,
    perf %in% c("6h_perf")) %>% 
    pull(colID)
  Viable6hIRI <- as.data.frame(t(cpmData_sub[IRIgenes,downSampleViable6h]))
  Viable6hIRI$condition <- "Viable: 6hours"
  # Bind the dataframe together:
  IRIdf <- rbind(
    CITstartIRI, CITendIRI, Viable3hIRI, Viable6hIRI
  )
  # IRIgenesLonger <- IRIgenes %>% pivot_longer(cols = !condition)
  IRIdf$condition <- factor(IRIdf$condition,
                            levels = c(
                              "CIT:start", "CIT:end",
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
    report$downsample <- fold
    IRImixedModel <- rbind(IRImixedModel, report)
  }

  downSampledIRIanalysis <- list.append(downSampledIRIanalysis, IRImixedModel)
  fold <- fold + 1 
  }
downSampledIRIanalysis <- do.call(rbind, downSampledIRIanalysis)

se <- function(y) {
  sd(y) / sqrt(length(y))
}

downSampledIRIanalysis %>% 
  mutate(condition = rownames(.), 
         condition = ifelse(str_detect(condition, 'conditionCIT'),
                            'conditionCIT: end',
                            condition),
         condition = ifelse(str_detect(condition, 'conditionViable: 3hours'),
                            'conditionViable: 3hours',
                            condition),
         condition = ifelse(str_detect(condition, 'conditionViable: 6hours'),
                            'conditionViable: 6hours',
                            condition)) %>%
  filter(!str_detect(condition, 'Intercept')) %>%
  group_by(gene, condition) %>% 
  summarise(meanEstimate = mean(Estimate),
            seEstimate = se(Estimate)) -> downSampledIRIanalysisSummarised
# Downsampled figure:


ggplot(
  downSampledIRIanalysisSummarised %>% filter(gene %in% IRItotest),
  aes(x = condition, y = meanEstimate, col = gene)
) +
  geom_point() +
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
  ylab("Gene expression linear coefficient") -> estimatePlotDownsampled


ggsave(estimatePlotDownsampled, filename = paste0(
  oneDriveOutput,
  "Figure1/estimatePlotDownsampled.pdf"
), width = 7, height = 4)

downSampledIRIanalysis %>% 
  write_csv(, file = paste0(
    oneDriveOutput,
    "Supplementary/IRIgenesCoefficientDownsampled.csv"
  ))

downSampledIRIanalysisSummarised %>% group_by(condition) %>% 
  summarise(mean = mean(meanEstimate),
         se = se(meanEstimate))

# Get the coefficients:
downSampledIRIanalysis %>%
##################################
######### Downsample analysis CLOCK:
##################################
library(nCV)
sampleinfo_sub[str_detect(sampleinfo_sub$sample, "nmp4"), ]$perf_success[2:3] <- "Perf_Unsuccessful"
# Subsampling variable for CCD:
# One sample is reclassed as non-viable, we tehrefore need a  subsampling of 4:
successSampleToUseCLOCK <- combn(sampleinfo_sub %>%
                              filter(perf_success == "Perf_Successful") %>%
                              pull(sample) %>% unique(), 4, simplify = F)

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
# Get the GTEx reference:
CLOCK_GTEx <- read.table("../../CLOCK/Results/CLOCK/GTEx/GTExMatrixRef.csv")
matrix_levels <- c(
  "ARNTL", "NPAS2", "CLOCK", "NFIL3", "CRY1",
  "CRY2", "NR1D1", "NR1D2", "PER1", "PER2",
  "PER3", "DBP", "TEF", "HLF"
)
matrix_levels <- str_to_title(matrix_levels)
CLOCK_GTEx <- CLOCK_GTEx[order(rownames(CLOCK_GTEx)), ]
rownames(CLOCK_GTEx) <- str_to_title(rownames(CLOCK_GTEx))
GTEXcorMatrix <- cor(t(CLOCK_GTEx))
GTEXcorVector <- GTEXcorMatrix[lower.tri(GTEXcorMatrix, diag = F)]

# Downsample: 

CCDdownsampled <- c()
for (set in successSampleToUseCLOCK) {
  downSampleViable <- sampleinfo_sub %>% filter(
    perf_success == "Perf_Successful",
    sample %in% set) %>% 
    pull(colID)
  
  cpmSuccess <- cpmData[, -1][,downSampleViable]
  cpmSuccess <- cbind(cpmData[, 1], cpmSuccess)
  
  # mClockD <- nCV::mClockD
  # nperm <- 10000
  # nCVsuccess <- homemadenCVnet(cpmSuccess, mClockD, nperm = nperm, seedN = 154)
  # Correlation vector:
  # Order the dataframe in the same way than the GTEx matrix
  cpmSuccessOrdered <- cpmSuccess[order(rownames(cpmSuccess)), -1]
  # Viable Livers:
  SuccessCorMatrix <- cor(t(cpmSuccessOrdered[which(
    rownames(cpmSuccessOrdered) %in% matrix_levels
  ), ]))
  # Extract the correlation vector for each:
  SuccessCorVector <- SuccessCorMatrix[lower.tri(SuccessCorMatrix, diag = F)]
  # Calculate the euclidean distance between the two:
  euclidean <- function(a, b) sqrt(sum((a - b)^2))
  CCDSuccess <- euclidean(GTEXcorVector, SuccessCorVector)
  CCDdownsampled <- c(CCDdownsampled, CCDSuccess)
  }
  
  
  

cpmSuccess <- cbind(cpmData[, 1], cpmSuccess)

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
cpmSuccessOrdered <- cpmSuccess[order(rownames(cpmSuccess)), -c(1,2)]
# Create the correlation matrix:
CITcorMatrix <- cor(t(cpmCITordered[which(rownames(cpmCITordered) %in%
  matrix_levels), ]))
# Viable livers:
FailedCorMatrix <- cor(t(cpmFailedOrdered[which(rownames(cpmFailedOrdered) %in%
  matrix_levels), ]))
# Non-viable Livers:

SuccessCorMatrix <- cor(t(cpmSuccessOrdered[which(rownames(cpmSuccessOrdered) %in% matrix_levels), ]))
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
# SE:
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