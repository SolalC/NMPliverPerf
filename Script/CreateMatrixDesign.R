# Fixed
create_design_matrix <- function(times, conditions, batches, rhythm_pattern, period = 24) {
  # Basic mean terms for two conditions
  X <- model.matrix(~ 0 + conditions + batches)
  colnames(X) <- c("mean_cond1", "mean_cond2", paste0("batch", 1:(length(unique(batches))-1)))
  # Add rhythmic terms based on pattern
  # rhythm_pattern can be: "none", "shared", "independent"
  if(rhythm_pattern != "none") {
    if(rhythm_pattern == "shared") {
      X <- cbind(X,
                 cos = cos((2 * pi * times)/period),
                 sin = sin((2 * pi * times)/period))
    } else if(rhythm_pattern == "independent") {
      cond1_mask <- conditions == unique(conditions)[1]
      cond2_mask <- conditions == unique(conditions)[2]
      X <- cbind(X,
                 cos1 = cos((2 * pi * times)/period) * cond1_mask,
                 sin1 = sin((2 * pi * times)/period) * cond1_mask,
                 cos2 = cos((2 * pi * times)/period) * cond2_mask,
                 sin2 = sin((2 * pi * times)/period) * cond2_mask)
    }
    else if(rhythm_pattern == "condition1") {
      cond1_mask <- conditions == unique(conditions)[1]
      X <- cbind(X,
                 cos = cos((2 * pi * times)/period) * cond1_mask,
                 sin = sin((2 * pi * times)/period) * cond1_mask)
    }
    else if(rhythm_pattern == "condition2") {
      cond2_mask <- conditions == unique(conditions)[2]
      X <- cbind(X,
                 cos = cos((2 * pi * times)/period) * cond2_mask,
                 sin = sin((2 * pi * times)/period) * cond2_mask)
    }
  }
  return(X)
}
# Design matrix including a random effect:
create_design_matrix_random <- function(times, conditions, batches, sample, rhythm_pattern, period = 24) {
  # Basic mean terms for two conditions
  X <- model.matrix(~ 0 + conditions + batches)
  colnames(X) <- c("mean_cond1", "mean_cond2", paste0("batch", 1:(length(unique(batches))-1)))
  # Add rhythmic terms based on pattern
  # rhythm_pattern can be: "none", "shared", "independent"
  if(rhythm_pattern != "none") {
    if(rhythm_pattern == "shared") {
      X <- cbind(X,
                 cos = cos((2 * pi * times)/period),
                 sin = sin((2 * pi * times)/period))
    } else if(rhythm_pattern == "independent") {
      cond1_mask <- conditions == unique(conditions)[1]
      cond2_mask <- conditions == unique(conditions)[2]
      X <- cbind(X,
                 cos1 = cos((2 * pi * times)/period) * cond1_mask,
                 sin1 = sin((2 * pi * times)/period) * cond1_mask,
                 cos2 = cos((2 * pi * times)/period) * cond2_mask,
                 sin2 = sin((2 * pi * times)/period) * cond2_mask)
    }
    else if(rhythm_pattern == "condition1") {
      cond1_mask <- conditions == unique(conditions)[1]
      X <- cbind(X,
                 cos = cos((2 * pi * times)/period) * cond1_mask,
                 sin = sin((2 * pi * times)/period) * cond1_mask)
    }
    else if(rhythm_pattern == "condition2") {
      cond2_mask <- conditions == unique(conditions)[2]
      X <- cbind(X,
                 cos = cos((2 * pi * times)/period) * cond2_mask,
                 sin = sin((2 * pi * times)/period) * cond2_mask)
    }
  }
  # Add the random effect to the design matrix:
  X <- cbind(X, as.integer(as.factor(sample)))
  colnames(X)[c(ncol(X))] <- c('sample')
  return(X)
}