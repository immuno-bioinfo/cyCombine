#' Compute Local Inverse Simpson's Index (LISI)
#'
#' NOTE: reimplemented from https://github.com/immunogenomics/LISI
#'
#' @param df A data.frame with cells (rows) and markers in columns.
#' @param markers Character vector or integer indices of columns containing features.
#'        Default is all numeric columns excluding metadata specified in `label_cols`.
#' @param label_cols Character vector of column names in `df` to compute LISI for.
#' @param perplexity The effective number of neighbors.
#' @param nn_eps Error bound for nearest neighbor search.
#'
#' @return A tibble with the original metadata plus new columns `lisi_<label>` for each label.
#' @export
compute_lisi <- function(
    df,
    markers = get_markers(df),
    label_cols = "batch",
    perplexity = 30,
    nn_eps = 0
) {

  N <- nrow(df)

  if (perplexity >= N) {
    perplexity <- max(1, N - 1)
    warning("Cannot calculate more neighbors than there are points. Setting perplexity to ", perplexity)
  }

  # Compute KNN
  dknn <- RANN::nn2(df[, markers], k = perplexity + 1, eps = nn_eps)

  # Remove self from neighbors
  nn_idx <- dknn$nn.idx[, -1]
  nn_dists <- dknn$nn.dists[, -1, drop = FALSE]

  # Internal function to compute LISI for a single label column
  compute_single_lisi <- function(label_colname) {
    labels_vec <- df[[label_colname]]

    if (any(is.na(labels_vec))) {
      message(sprintf("LISI: Cannot compute LISI on missing values for '%s'. Returning NA.", label_colname))
      return(rep(NA_real_, N))
    }

    labels_fact <- as.integer(factor(labels_vec))
    n_batches <- length(unique(labels_fact))

    if (n_batches <= 1) return(rep(1, N)) # No diversity possible

    # Vectorized Simpson Index Calculation
    # For each cell i, calculate sum(p_j^2) over neighbors j, where p_j is proportion of label k in neighborhood
    simpson_vals <- apply(nn_idx, 1, function(neighbors) {
      neighbor_labels <- labels_fact[neighbors]
      # Count frequencies of each label in the neighborhood
      counts <- tabulate(neighbor_labels, nbins = n_batches)
      probs <- counts / perplexity
      sum(probs^2)
    })

    # Handle division by zero if perplexity is 0
    simpson_vals[simpson_vals == 0] <- .Machine$double.eps

    return(1 / simpson_vals)
  }

  # Apply to all label columns
  lisi_results <- lapply(label_cols, function(x) {
    tibble::tibble(!!paste0("lisi_", x) := compute_single_lisi(x))
  })
  lisi_results <- do.call(dplyr::bind_cols, lisi_results)

  # Bind back to original metadata
  # lisi <- df |>
  #   dplyr::bind_cols(lisi_results)
  return(lisi_results)
}



#' Compute LISI split by groups
#'
#' NOTE: reimplemented from https://github.com/carmonalab/scIntegrationMetrics
#'
#' @param df A data.frame containing features and metadata.
#' @param markers Character vector or indices of feature columns.
#' @param label_cols Character vector of columns to compute LISI for.
#' @param split_by Character column name to split the data by (e.g., "cell_type").
#' @param normalize Logical. Normalize LISI to \[0, 1\]?
#' @param metricsLabels Character vector of levels in `split_by` to process. Default: all.
#' @param min.cells.split Minimum rows required in a split to process.
#' @param min.vars.label Minimum unique values in `label_cols` required in a split.
#' @param ... Additional arguments passed to `compute_lisi`.
#'
#' @return A named list of tibbles, one per group in `split_by`.
#' @export
compute_cilisi <- function(
    df,
    markers = NULL,
    label_cols = "batch",
    split_by = "label",
    downsample = TRUE,
    metricsLabels = NULL,
    return_mean = TRUE,
    normalize = TRUE,
    min.cells.split = 10,
    min.vars.label = 2,
    ...
) {
  # Identify groups to process
  if (is.null(metricsLabels)) {
    metricsLabels <- unique(na.omit(df[[split_by]]))
  }

  # Split data into a list of data.frames
  if (downsample) {
    sample_size <- min(10000, table(df[,label_cols[[1]]]))
    df <- df |>
      dplyr::group_by(dplyr::all_of(label_cols)) |>
      dplyr::slice_sample(n = sample_size)
    }
  df_split <- split(df, df[[split_by]])

  # names(data_split) <- unique(metricsLabels)


  # Compute LISI for each valid group
  results <- lapply(df_split, function(df_label) {

    # Run LISI
    res <- compute_lisi(
      df = df_label,
      markers = markers,
      label_cols = label_cols
    )

    # Normalize
    if (normalize) {
      # Calculate max possible LISI (number of unique levels per label)
      max_lisi <- sapply(label_cols, function(lc) length(unique(df_label[[lc]])))

      # Normalize: (value - 1) / (max - 1)
      res <- res %>%
        dplyr::mutate(dplyr::across(dplyr::starts_with("lisi_"), ~ {
          label_name <- sub("lisi_", "", dplyr::cur_column())
          max_val <- max_lisi[label_name]
          if (max_val <= 1) return(.x)
          (.x - 1) / (max_val - 1)
        }, .names = "{.col}"))
    }

    res
  })

  if(return_mean) results <- mean(unlist(results))

  return(results)
}
