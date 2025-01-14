# GLM negative binomial model:
fit_nb_model <- function(counts, design, size_factors) {
  counts <- unlist(counts)
  fit <- glm.nb(counts ~ offset(log(size_factors)) + design - 1, link = log)
  return(list(
    coefficients = coef(fit),
    loglik = logLik(fit),
    df = length(coef(fit)),
    fit = fit
  ))
}
# GLMM negative binomial model:
fit_nb_mixed <- function(counts, design, size_factors) {
  # Extract the fixed effect
  fitDF <- design
  effects <- colnames(fitDF)
  fixedEffects <- effects[-length(effects)]
  # Extract the random effect
  randomEffects <- effects[length(effects)] 
  # Create the formula:
  formula <- paste("counts ~", paste(fixedEffects, collapse = " + "))
  # Add random effect as well as the size factor:
  formula <- paste0(formula, ' + ', paste0("(1|", randomEffects, ") + offset(log(size_factors))"))
  # Add the count and size factor columns:
  counts <- unlist(counts)
  fitDF <- cbind(fitDF, counts)
  fitDF <- cbind(fitDF, size_factors)
  fitDF <- as.data.frame(fitDF)
  # Fit the glmm negative binomial model:
  fit <- glmer.nb(as.formula(formula), data = (fitDF))
  return(
    list(
      coefficients = fixef(fit),
      loglik = logLik(fit),
      df = length(fixef(fit)),
      fit = fit)
  )
}
# Main function to call:
fit_temporal_expression <- function(counts, colData, period = 24) {
  
  print(paste0('Condition1: ', unique(colData$condition)[1]))
  print(paste0('Condition2: ', unique(colData$condition)[2]))
  # Validate input
  if(length(unique(colData$condition)) != 2) {
    stop("This implementation expects exactly 2 conditions")
  }
  # Create DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = colData,
                                design = ~ condition + batch)
  
  # Estimate size factors and dispersions
  dds <- estimateSizeFactors(dds)
  # dds <- estimateDispersions(dds)
  # Extract required values
  size_factors <- sizeFactors(dds)
  dispersions <- dispersions(dds)
  # For 2 conditions, we have 4 possible rhythmicity patterns in Step 1:
  rhythm_patterns <- c(
    "none",           # No rhythmicity in either condition
    "shared",         # Same rhythmicity in both conditions
    "independent",    # Different rhythmicity in each condition
    "condition1",     # Rhythmicity only in condition 1
    "condition2"      # Rhythmicity only in condition 2
  )
  # Calculate BIC and BICW
  # Step 1: Assess rhythmicity
  step1_fits <- list()
  print('Start fitting rhythmicity:')
  pb <- txtProgressBar(min = 0, max = nrow(counts), style = 3)
  for(gene in 1:nrow(counts)) {
    tryCatch({
      gene_fits <- list()
      n_param <- c()
      for(pattern in rhythm_patterns) {
        design <- create_design_matrix(
          colData$time,
          colData$condition,
          colData$batch,
          pattern
        )
        fit <- fit_nb_model(
          counts[gene,],
          design,
          size_factors)
        gene_fits[[pattern]] <- fit
        n_param <- c(n_param, ncol(design))
      }
      # Select best model using BICW
      weights <- calculate_model_weights(gene_fits, ncol(counts), n_param)
      step1_fits[[gene]] <- list(
        fits = gene_fits,
        best_model = names(which.max(weights$BICW)),
        weights = weights
      )
      setTxtProgressBar(pb, gene)
    }, error = function(e) {
      message(sprintf("Error fitting gene %s: %s", rownames(counts)[gene], e$message))
    })
  }
  # Remove genes that could not be used:
  names(step1_fits) <- rownames(counts)
  step1_fits <- Filter(Negate(is.null), step1_fits)
  step2_fits <- list()
  print('Start fitting mean:')
  pb <- txtProgressBar(min = 0, max = length(step1_fits), style = 3)
  
  for(gene in 1:length(step1_fits)) {
    n_param <- c()
    geneName <- names(step1_fits)[gene]
    best_rhythm <- step1_fits[[gene]]$best_model
    gene_fits <- list()
    
    # Fit models with different and shared means
    for(mean_pattern in c("different", "shared")) {
      if(mean_pattern == "shared") {
        design <- model.matrix(~ 1 + batch, data = colData)
      } else {
        design <- create_design_matrix(
          colData$time,
          colData$condition,
          colData$batch,
          best_rhythm
        )
      }
      fit <- fit_nb_model(
        counts[geneName,],
        design,
        size_factors)
      n_param <- c(n_param, ncol(design))
      gene_fits[[mean_pattern]] <- fit
    }
    weights <- calculate_model_weights(gene_fits, ncol(counts), n_param)
    step2_fits[[gene]] <- list(
      fits = gene_fits,
      best_model = names(which.max(weights$BICW)),
      weights = weights
    )
    setTxtProgressBar(pb, gene)
  }
  names(step2_fits) <- names(step1_fits)
  return(list(
    step1_fits = step1_fits,
    step2_fits = step2_fits,
    size_factors = size_factors,
    gene = names(step1_fits)
  ))
}

fit_temporal_expression_mixed <- function(counts, colData, period = 24, n_cores) {
  
  print(paste0('Condition1: ', unique(colData$condition)[1]))
  print(paste0('Condition2: ', unique(colData$condition)[2]))
  # Validate input
  if(length(unique(colData$condition)) != 2) {
    stop("This implementation expects exactly 2 conditions")
  }
  # Create DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = colData,
                                design = ~ condition + batch)
  
  # Estimate size factors and dispersions
  dds <- estimateSizeFactors(dds)
  # dds <- estimateDispersions(dds)
  # Extract required values
  size_factors <- sizeFactors(dds)
  dispersions <- dispersions(dds)
  # For 2 conditions, we have 4 possible rhythmicity patterns in Step 1:
  rhythm_patterns <- c(
    "none",           # No rhythmicity in either condition
    "shared",         # Same rhythmicity in both conditions
    "independent",    # Different rhythmicity in each condition
    "condition1",     # Rhythmicity only in condition 1
    "condition2"      # Rhythmicity only in condition 2
  )
  # Calculate BIC and BICW
  # Step 1: Assess rhythmicity
  process_gene <- function(gene_idx) {
    tryCatch({
      gene_fits <- list()
      n_param <- c()
      
      for(pattern in rhythm_patterns) {
        design <- create_design_matrix_random(
          colData$time,
          colData$condition,
          colData$batch,
          sample = mergedInfoInput$sample,
          pattern
        )
        fit <- fit_nb_mixed(
          counts[gene_idx,],
          design,
          size_factors)
        
        gene_fits[[pattern]] <- fit
        n_param <- c(n_param, ncol(design))
      }
      # Select best model using BICW
      weights <- calculate_model_weights(gene_fits, ncol(counts), n_param)
  
      result <- list(
        fits = gene_fits,
        best_model = names(which.max(weights$BICW)),
        weights = weights
      )
      
      names(result) <- rownames(counts)[gene_idx]
      return(result)
      
    }, error = function(e) {
      message(sprintf("Error fitting gene %s: %s", rownames(counts)[gene_idx], e$message))
      return(NULL)
    })
  }
  # Run the parallel processing
  print('Start fitting rhythmicity in parallel:')
  step1_fits <- mclapply(
    1:nrow(counts),
    process_gene,
    mc.cores = n_cores
  )
  # Remove genes that could not be used:
  names(step1_fits) <- rownames(counts)
  # Add names:
  for (i in names(step1_fits)) {
    names(step1_fits[i][[1]]) <-   (c(i, 'best_model', 'weights'))
  }
  # Remove genes that could not be used:
  step1_fits <- Filter(Negate(is.null), step1_fits)
  step2_fits <- list()
  saveRDS(step1_fits, '/scratch/project_mnt/S0007/solal/NMP/Results/Circadian/intermediary/step1_fits.rds')
  print('Start fitting mean:')
  pb <- txtProgressBar(min = 0, max = length(step1_fits), style = 3)
  
  for(gene in 1:length(step1_fits)) {
    n_param <- c()
    geneName <- names(step1_fits)[gene]
    best_rhythm <- step1_fits[[gene]]$best_model
    gene_fits <- list()
    
    # Fit models with different and shared means
    for(mean_pattern in c("different", "shared")) {
      if(mean_pattern == "shared") {
        design <- model.matrix(~ 1 + batch, data = colData)
      } else {
        design <- create_design_matrix_random(
          colData$time,
          colData$condition,
          colData$batch,
          sample = mergedInfoInput$sample,
          best_rhythm
        )
      }
      fit <- fit_nb_mixed(
        counts[geneName,],
        design,
        size_factors)
      n_param <- c(n_param, ncol(design))
      gene_fits[[mean_pattern]] <- fit
    }
    weights <- calculate_model_weights(gene_fits, ncol(counts), n_param)
    step2_fits[[gene]] <- list(
      fits = gene_fits,
      best_model = names(which.max(weights$BICW)),
      weights = weights
    )
    setTxtProgressBar(pb, gene)
  }
  names(step2_fits) <- names(step1_fits)
  return(list(
    step1_fits = step1_fits,
    step2_fits = step2_fits,
    size_factors = size_factors,
    gene = names(step1_fits)
  ))
}





