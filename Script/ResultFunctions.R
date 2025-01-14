extract_parameters <- function(listResult, period = 24){
  require(rlist)
  amplitude <- function(a, b) {2*sqrt(a^2+b^2)}
  phase <- function(a, b) {period/(2*pi)*atan2(b,a)}
  geneTested <- data.frame(best_model = unlist(lapply(listResult, function(x) x$best_model)),
                           gene = names(unlist(lapply(listResult, function(x) x$best_model))))
  listOutput <- list()
  rhythm <- unique(geneTested$best_model)
  for(model in rhythm) {
    modelGene <- geneTested %>% filter(best_model == model) %>% pull(gene) %>% unique()
    print(paste0(model, ': ' ,length(modelGene), ' Genes'))
    df <- data.frame()
    pb <- txtProgressBar(min = 0, max = length(modelGene), style = 3)
    cnt <- 1
    for(gene in modelGene) {
      
      coef = fixef(listResult[[gene]][[1]][[model]]$fit)
      BICW <- max(step1_fits[[gene]]$weights$BICW)
      temp <- data.frame(term = names(coef), coef = coef, 
                         gene = gene, model = model, BICW=BICW)
      temp <- temp %>% pivot_wider(names_from=term, values_from=coef)
      df <- rbind(df, temp)
      cnt <- cnt+1
      setTxtProgressBar(pb, cnt)
    }
    colnames(df) <- str_replace(colnames(df), 'cos', 'a')
    colnames(df) <- str_replace(colnames(df), 'sin', 'b')
    if(model == 'none'){
      df$amplitude <- NA
      df$phase <- NA
    }
    if(model %in% c('shared', 'condition1', 'condition2')) {
      df$amplitude <- amplitude(df$a, df$b)
      df$phase <- phase(df$a, df$b)
    } 
    if(model == 'independent') {
      df$amplitude1 <- amplitude(df$a1, df$b1)
      df$phase1 <- phase(df$a1, df$b1)
      df$amplitude2 <- amplitude(df$a2, df$b2)
      df$phase2 <- phase(df$a2, df$b2)
    }
    
    listOutput <- list.append(listOutput, df)
  } 
  names(listOutput) <- rhythm
  return(listOutput)
}

# Add the code for the amplitude plot

amplitudePlot <- function(result, thresh, col = met.brewer(n = 2, 'Juarez')) {
  ampThreshold <- seq(0.25, thresh, by = 0.1)
  # Viable:
  rbind(result[['shared']] %>% dplyr::select(gene, amplitude),
        result[['condition2']] %>% dplyr::select(gene, amplitude),
        result[['independent']] %>% mutate(amplitude = amplitude2) %>% 
          dplyr::select(gene, amplitude)) -> ampViable
  # Viable count:
  viableCount <- sapply(ampThreshold, 
                        function(x) sum(ampViable$amplitude > x, na.rm=TRUE))
  # Non-viable:
  rbind(result[['shared']] %>% dplyr::select(gene, amplitude),
        result[['condition1']] %>% dplyr::select(gene, amplitude),
        result[['independent']] %>% mutate(amplitude = amplitude1) %>% 
          dplyr::select(gene, amplitude)) -> ampNonViable
  nonViableCount <- sapply(ampThreshold, 
                           function(x) sum(ampNonViable$amplitude > x, na.rm=TRUE))
  # Percentage loss:
  ampDf <- data.frame(ampThreshold, viableCount, nonViableCount) %>%
    mutate(percentLoss = (viableCount - nonViableCount) / viableCount * 100)
  # Longer plot:
  ampDf.longer <- ampDf %>%
    pivot_longer(cols = c('viableCount', 'nonViableCount'),
                 names_to = 'condition',
                 values_to = 'count') %>%
    mutate(condition = factor(condition, levels = c('viableCount', 'nonViableCount')))
  # Plot:
  ggplot(ampDf.longer, aes(x = ampThreshold, y = log10(count), col= condition)) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = 'top') + 
    ylab('Log10 \nNumber of \ngenes') +
    xlab('Amplitude (log2FC)') +
    xlim(0,4) +
    scale_color_manual(values = col,
                       label = c('Viable', 'Non-Viable'),
                       '') -> thresholdPlot
  # Loss percentage:
  ggplot(ampDf, aes(x = ampThreshold, y = percentLoss)) +
    geom_area() +
    ylab('Loss of rhythm \n(Viable Vs. Non-Viable)\n%of genes') +
    xlab('Amplitude (log2FC)') +
    xlim(0,4) +
    theme_minimal() -> lossPlot
  ampPlot <- cowplot::plot_grid(thresholdPlot, lossPlot, ncol = 1)
  return(ampPlot)
}

# Add the code for the phase plot

phasePlot <- function(result, col = met.brewer(n = 2, 'Juarez')){
  # Viable:
  rbind(result[['shared']] %>% dplyr::select(gene, phase, model),
        result[['condition2']] %>% dplyr::select(gene, phase, model),
        result[['independent']] %>% mutate(phase = phase2, model) %>% 
          dplyr::select(gene, phase, model)) -> phaseViable
  # Non-viable:
  rbind(result[['shared']] %>% dplyr::select(gene, phase, model),
        result[['condition1']] %>% dplyr::select(gene, phase, model),
        result[['independent']] %>% mutate(phase = phase1, model) %>% 
          dplyr::select(gene, phase, model)) -> phaseNonViable
  # Viable:
  ggplot(phaseViable, aes(x = phase)) +
    geom_histogram(fill = col[1], binwidth = 0.05) +
    theme_minimal() +
    ylab('# of gene') +
    scale_x_continuous(breaks = c(0,6,12,18,24),
                       labels = c(0,6,12,18,24)) +
    xlab('Perfusion time (h)') -> viablePhase
  # Viable phase plot:
  ggplot(phaseNonViable, aes(x = phase)) +
    geom_histogram(fill = col[2], binwidth = 0.05) +
    theme_minimal() + 
    ylab('# of gene') +
    scale_x_continuous(breaks = c(0,6,12,18,24),
                       labels = c(0,6,12,18,24)) +
    xlab('Perfusion time (h)') -> nonViablePhase
  
  cowplot::plot_grid(viablePhase, nonViablePhase, ncol = 1) -> phasePlot
  return(phasePlot)
}


tempEnrichment <- function(result, geneBackground, condition, pathwayDB='KEGG', seed = 13845) {
  require(gprofiler2)
  set.seed(seed)
  # Get the gene and their phase:
  if(condition == 'Viable'){
    # phase <- rbind(result[['shared']] %>% dplyr::select(gene, phase, model),
    #                result[['condition2']] %>% dplyr::select(gene, phase, model),
    #                result[['independent']] %>% mutate(phase = phase2, model) %>% 
    #                  dplyr::select(gene, phase, model))
    phase <- result[['shared']] %>% dplyr::select(gene, phase, model)
  }
  if(condition == 'Non-Viable'){
    phase <- rbind(result[['shared']] %>% dplyr::select(gene, phase, model),
                   result[['condition1']] %>% dplyr::select(gene, phase, model),
                   result[['independent']] %>% mutate(phase = phase1, model) %>% 
                     dplyr::select(gene, phase, model))
  }
  # Get the phase breaks for enrichment:
  phaseBreaks <- seq(0, 24, by = 2)
  tempEnrich <- data.frame()
  # Perform the enrichment by breaks:
  pb <- txtProgressBar(min = 0, max = length(phaseBreaks), style = 3)
  ct <- 1
  for(i in phaseBreaks){
    geneSet <- phase %>% filter(
      dplyr::between(phase, left = i-1, right = i+1)) %>% 
      pull(gene)
    # Perform enrichment:
    pathways <-  gost(query = geneSet, 
                      organism = "hsapiens",significant = F, 
                      user_threshold = 0.05, correction_method = "g_SCS", 
                      domain_scope = "custom", custom_bg = geneBackground, 
                      sources = pathwayDB)$result
    pathways <- pathways %>% mutate(phase = i)
    tempEnrich <- rbind(tempEnrich, pathways)
    setTxtProgressBar(pb, ct)
    ct <- ct +1 
  }
  return(tempEnrich)
}

