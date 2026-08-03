



rm(list = ls())

  cat("\014")  
# ============================================================
# Real-data normalization workflow
# Input: merged_data(1).csv
# Last column: Groups
#
# Output:
#   10 aligned data blocks:
#   1 Raw data
#   2 CLR+1
#   3 CLR-BMR
#   4 TSS
#   5 Rarefaction
#   6 CSS
#   7 edgeR-TMM
#   8 DESeq2
#   9 ALDEx2
#  10 ANCOM-style ALR output
#
# Final column order:
# taxa, Methods, Groups, observation
# ============================================================

set.seed(1)

# -------------------------
# File settings
# -------------------------
input_file <- "D:\\merged_data.csv"
output_file <- "D:\\real_data_10_methods_combined.csv"

# -------------------------
# Read the real dataset
# -------------------------
real_data <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(real_data) < 2) {
  stop("The input file must contain taxa columns and a final Groups column.")
}

if (tail(colnames(real_data), 1) != "Groups") {
  stop(
    "The last column must be named 'Groups'. Detected last column: ",
    tail(colnames(real_data), 1)
  )
}

# Save the group column before removing it from the taxa matrix.
Groups_original <- real_data$Groups

if (anyNA(Groups_original)) {
  stop("The Groups column contains missing values.")
}

# Remove Groups before applying normalization.
taxa_data <- real_data[, setdiff(colnames(real_data), "Groups"), drop = FALSE]
taxa_cols <- colnames(taxa_data)

# Convert and validate taxa columns.
for (j in seq_along(taxa_cols)) {
  taxa_data[[j]] <- suppressWarnings(as.numeric(taxa_data[[j]]))
}

X_check <- as.matrix(taxa_data)

if (anyNA(X_check)) {
  stop("At least one taxa value is missing or non-numeric after conversion.")
}

if (any(X_check < 0)) {
  stop("Taxa counts must be non-negative.")
}

if (any(rowSums(X_check) == 0)) {
  stop("At least one sample has total count zero.")
}

n_samples <- nrow(taxa_data)

# Preserve the original sample ordering in every method block.
observation_original <- seq_len(n_samples)

message("Input samples: ", n_samples)
message("Taxa columns: ", length(taxa_cols))
message(
  "Group counts: ",
  paste(names(table(Groups_original)), table(Groups_original), sep = "=", collapse = ", ")
)

# Working data frame used by methods that require the group labels.
analysis_data <- data.frame(
  taxa_data,
  Groups = Groups_original,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ============================================================
# Normalization functions
# ============================================================

# -------------------------
# Rarefaction
# -------------------------
rarefy_only_taxa <- function(df, taxa_cols, depth = NULL) {
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required.")
  }

  X <- as.matrix(df[, taxa_cols, drop = FALSE])

  if (anyNA(X)) stop("Taxa matrix contains NA values.")
  if (any(X < 0)) stop("Counts must be non-negative.")

  libsize <- rowSums(X)

  if (any(libsize == 0)) {
    stop("Rarefaction cannot be applied to samples with zero total counts.")
  }

  if (is.null(depth)) {
    depth <- min(libsize)
  }

  if (depth <= 0 || depth > min(libsize)) {
    stop("Rarefaction depth must be positive and no greater than the minimum library size.")
  }

  X_rar <- vegan::rrarefy(X, sample = depth)

  out <- as.data.frame(X_rar, check.names = FALSE)
  colnames(out) <- taxa_cols
  rownames(out) <- rownames(df)
  out
}

# -------------------------
# CLR with pseudocount 1
# -------------------------
clr_only_taxa <- function(df, taxa_cols, pseudocount = 1) {
  X <- as.matrix(df[, taxa_cols, drop = FALSE])

  if (anyNA(X)) stop("Taxa matrix contains NA values.")
  if (any(X < 0)) stop("Counts must be non-negative.")

  X <- X + pseudocount
  logX <- log(X)
  clrX <- logX - rowMeans(logX)

  out <- as.data.frame(clrX, check.names = FALSE)
  colnames(out) <- taxa_cols
  rownames(out) <- rownames(df)
  out
}

# -------------------------
# CLR helper
# -------------------------
clr_from_positive_proportions <- function(X_prop, taxa_cols, row_ids = NULL) {
  X_prop <- as.matrix(X_prop)

  if (anyNA(X_prop)) stop("CLR input contains NA values.")
  if (any(X_prop <= 0)) stop("CLR input must be strictly positive.")

  logX <- log(X_prop)
  clrX <- logX - rowMeans(logX)

  out <- as.data.frame(clrX, check.names = FALSE)
  colnames(out) <- taxa_cols

  if (!is.null(row_ids)) {
    rownames(out) <- row_ids
  }

  out
}

# -------------------------
# Bayesian-multiplicative/CZM zero replacement
# -------------------------
czm_replace <- function(X_prop) {
  if (!requireNamespace("zCompositions", quietly = TRUE)) {
    stop("Package 'zCompositions' is required.")
  }

  X_repl <- zCompositions::cmultRepl(
    X_prop,
    label = 0,
    method = "CZM",
    output = "prop",
    z.delete = FALSE,
    z.warning = FALSE
  )

  X_repl <- as.matrix(X_repl)

  if (anyNA(X_repl)) stop("CZM replacement produced NA values.")
  if (any(X_repl <= 0)) stop("CZM replacement produced non-positive values.")

  sweep(X_repl, 1, rowSums(X_repl), "/")
}

# -------------------------
# CLR after multiplicative replacement
# -------------------------
clr_bmr_only_taxa <- function(df, taxa_cols) {
  X <- as.matrix(df[, taxa_cols, drop = FALSE])

  if (anyNA(X)) stop("Taxa matrix contains NA values.")
  if (any(X < 0)) stop("Counts must be non-negative.")

  libsize <- rowSums(X)

  if (any(libsize == 0)) {
    stop("CLR-BMR cannot be applied to samples with zero total counts.")
  }

  X_prop <- sweep(X, 1, libsize, "/")

  if (any(X_prop == 0)) {
    X_prop <- czm_replace(X_prop)
  }

  clr_from_positive_proportions(
    X_prop = X_prop,
    taxa_cols = taxa_cols,
    row_ids = rownames(df)
  )
}

# -------------------------
# TSS
# -------------------------
tss_only_taxa <- function(df, taxa_cols) {
  X <- as.matrix(df[, taxa_cols, drop = FALSE])

  if (anyNA(X)) stop("Taxa matrix contains NA values.")
  if (any(X < 0)) stop("Counts must be non-negative.")

  libsize <- rowSums(X)

  if (any(libsize == 0)) {
    stop("TSS cannot be applied to samples with zero total counts.")
  }

  X_tss <- sweep(X, 1, libsize, "/")

  out <- as.data.frame(X_tss, check.names = FALSE)
  colnames(out) <- taxa_cols
  rownames(out) <- rownames(df)
  out
}

# -------------------------
# CSS
# -------------------------
css_only_taxa <- function(df, taxa_cols) {
  if (!requireNamespace("metagenomeSeq", quietly = TRUE)) {
    stop("Package 'metagenomeSeq' is required.")
  }

  X <- as.matrix(df[, taxa_cols, drop = FALSE])

  if (anyNA(X)) stop("Taxa matrix contains NA values.")
  if (any(X < 0)) stop("Counts must be non-negative.")

  obj <- metagenomeSeq::newMRexperiment(t(X))
  p <- metagenomeSeq::cumNormStatFast(obj)
  obj <- metagenomeSeq::cumNorm(obj, p = p)

  X_css <- metagenomeSeq::MRcounts(obj, norm = TRUE, log = FALSE)

  out <- as.data.frame(t(X_css), check.names = FALSE)
  colnames(out) <- taxa_cols
  rownames(out) <- rownames(df)
  out
}

# -------------------------
# edgeR TMM-CPM
# -------------------------
edger_tmm_cpm_only_taxa <- function(df, taxa_cols) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' is required.")
  }

  X <- as.matrix(df[, taxa_cols, drop = FALSE])

  if (anyNA(X)) stop("Taxa matrix contains NA values.")
  if (any(X < 0)) stop("Counts must be non-negative.")

  dge <- edgeR::DGEList(counts = t(X))
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  tmm_cpm <- edgeR::cpm(dge, log = FALSE)

  out <- as.data.frame(t(tmm_cpm), check.names = FALSE)
  colnames(out) <- taxa_cols
  rownames(out) <- rownames(df)
  out
}

# -------------------------
# DESeq2 normalized counts
# -------------------------
deseq_norm_only_taxa <- function(df, taxa_cols, sf_type = "poscounts") {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required.")
  }

  X <- as.matrix(df[, taxa_cols, drop = FALSE])

  if (anyNA(X)) stop("Taxa matrix contains NA values.")
  if (any(X < 0)) stop("Counts must be non-negative.")

  countData <- t(round(X))
  storage.mode(countData) <- "integer"

  sample_ids <- paste0("S", seq_len(nrow(X)))
  rownames(countData) <- taxa_cols
  colnames(countData) <- sample_ids

  colData <- data.frame(
    row.names = sample_ids,
    dummy = rep(1, nrow(X))
  )

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = countData,
    colData = colData,
    design = ~ 1
  )

  dds <- DESeq2::estimateSizeFactors(dds, type = sf_type)
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)

  out <- as.data.frame(t(norm_counts), check.names = FALSE)
  colnames(out) <- taxa_cols
  rownames(out) <- rownames(df)
  out
}

# -------------------------
# ALDEx2 median CLR output
# -------------------------
aldex_clr_median_only_taxa <- function(
    df,
    taxa_cols,
    group_col = "Groups",
    mc.samples = 128,
    denom = "all"
) {
  if (!requireNamespace("ALDEx2", quietly = TRUE)) {
    stop("Package 'ALDEx2' is required.")
  }

  X <- as.matrix(df[, taxa_cols, drop = FALSE])

  if (anyNA(X)) stop("Taxa matrix contains NA values.")
  if (any(X < 0)) stop("Counts must be non-negative.")
  if (any(rowSums(X) == 0)) {
    stop("ALDEx2 cannot be applied to samples with zero total counts.")
  }

  if (!group_col %in% colnames(df)) {
    stop("Column '", group_col, "' was not found.")
  }

  keep_taxa <- taxa_cols[colSums(X) > 0]
  dropped_taxa <- setdiff(taxa_cols, keep_taxa)

  if (length(keep_taxa) < 2) {
    stop("Fewer than two non-zero taxa are available for ALDEx2.")
  }

  X_keep <- X[, keep_taxa, drop = FALSE]
  sample_ids <- paste0("S", seq_len(nrow(X_keep)))

  X_t <- t(round(X_keep))
  rownames(X_t) <- keep_taxa
  colnames(X_t) <- sample_ids

  conds <- as.character(df[[group_col]])

  clr_obj <- ALDEx2::aldex.clr(
    reads = X_t,
    conds = conds,
    mc.samples = mc.samples,
    denom = denom,
    verbose = FALSE,
    useMC = FALSE
  )

  if ("getMonteCarloInstances" %in% getNamespaceExports("ALDEx2")) {
    mc <- ALDEx2::getMonteCarloInstances(clr_obj)
  } else if (methods::is(clr_obj, "S4") &&
             "analysisData" %in% methods::slotNames(clr_obj)) {
    mc <- methods::slot(clr_obj, "analysisData")
  } else {
    stop("Could not extract ALDEx2 Monte Carlo instances.")
  }

  median_taxa_from_matrix <- function(m, keep_taxa) {
    m <- as.matrix(m)

    if (!is.null(rownames(m))) {
      common <- intersect(keep_taxa, rownames(m))
      if (length(common) == length(keep_taxa)) {
        vals <- apply(m[keep_taxa, , drop = FALSE], 1, median, na.rm = TRUE)
        return(vals)
      }
    }

    if (!is.null(colnames(m))) {
      common <- intersect(keep_taxa, colnames(m))
      if (length(common) == length(keep_taxa)) {
        vals <- apply(m[, keep_taxa, drop = FALSE], 2, median, na.rm = TRUE)
        return(vals)
      }
    }

    if (nrow(m) == length(keep_taxa)) {
      vals <- apply(m, 1, median, na.rm = TRUE)
      names(vals) <- keep_taxa
      return(vals)
    }

    if (ncol(m) == length(keep_taxa)) {
      vals <- apply(m, 2, median, na.rm = TRUE)
      names(vals) <- keep_taxa
      return(vals)
    }

    stop(
      "Cannot identify the taxa dimension in an ALDEx2 Monte Carlo matrix: ",
      paste(dim(m), collapse = " x ")
    )
  }

  if (is.list(mc) && length(mc) == nrow(X_keep)) {
    clr_med_keep <- sapply(
      mc,
      median_taxa_from_matrix,
      keep_taxa = keep_taxa
    )
  } else if (is.list(mc) && length(mc) == mc.samples) {
    arr_list <- lapply(mc, as.matrix)
    first <- arr_list[[1]]

    if (nrow(first) == length(keep_taxa) &&
        ncol(first) == nrow(X_keep)) {
      arr <- simplify2array(arr_list)
      clr_med_keep <- apply(arr, c(1, 2), median, na.rm = TRUE)
    } else if (nrow(first) == nrow(X_keep) &&
               ncol(first) == length(keep_taxa)) {
      arr <- simplify2array(arr_list)
      clr_med_keep <- t(apply(arr, c(1, 2), median, na.rm = TRUE))
    } else {
      stop("Cannot identify ALDEx2 list-by-instance dimensions.")
    }
  } else if (is.array(mc) && length(dim(mc)) == 3) {
    d <- dim(mc)

    if (d[1] == length(keep_taxa) && d[2] == nrow(X_keep)) {
      clr_med_keep <- apply(mc, c(1, 2), median, na.rm = TRUE)
    } else if (d[1] == nrow(X_keep) && d[2] == length(keep_taxa)) {
      clr_med_keep <- t(apply(mc, c(1, 2), median, na.rm = TRUE))
    } else if (d[2] == length(keep_taxa) && d[3] == nrow(X_keep)) {
      clr_med_keep <- apply(mc, c(2, 3), median, na.rm = TRUE)
    } else {
      stop("Unexpected ALDEx2 array dimensions: ", paste(d, collapse = " x "))
    }
  } else {
    stop("Unexpected ALDEx2 Monte Carlo object structure.")
  }

  clr_med_keep <- as.matrix(clr_med_keep)

  if (nrow(clr_med_keep) != length(keep_taxa) ||
      ncol(clr_med_keep) != nrow(X_keep)) {
    stop(
      "ALDEx2 output dimensions are inconsistent. Expected ",
      length(keep_taxa), " taxa x ", nrow(X_keep), " samples; obtained ",
      nrow(clr_med_keep), " x ", ncol(clr_med_keep), "."
    )
  }

  rownames(clr_med_keep) <- keep_taxa

  clr_med_full <- matrix(
    0,
    nrow = length(taxa_cols),
    ncol = nrow(X),
    dimnames = list(taxa_cols, rownames(df))
  )

  clr_med_full[keep_taxa, ] <- clr_med_keep[keep_taxa, , drop = FALSE]

  if (length(dropped_taxa) > 0) {
    message(
      "ALDEx2: all-zero taxa were reinserted as zero-valued columns: ",
      paste(dropped_taxa, collapse = ", ")
    )
  }

  out <- as.data.frame(t(clr_med_full), check.names = FALSE)
  colnames(out) <- taxa_cols
  rownames(out) <- rownames(df)
  out
}

# -------------------------
# ANCOM-style ALR output
#
# This reproduces the matrix-producing ANCOM/ALR procedure used
# in the supplied normalization script. It is not an ANCOM-BC
# differential-abundance result.
# -------------------------
ancom_alr_only_taxa <- function(
    df,
    taxa_cols,
    group_col = "Groups",
    pseudocount = 1,
    reference_taxon = NULL
) {
  X <- as.matrix(df[, taxa_cols, drop = FALSE])

  if (anyNA(X)) stop("Taxa matrix contains NA values.")
  if (any(X < 0)) stop("Counts must be non-negative.")
  if (any(rowSums(X) == 0)) {
    stop("ANCOM-style ALR cannot be applied to zero-total samples.")
  }

  if (!group_col %in% colnames(df)) {
    stop("Column '", group_col, "' was not found.")
  }

  if (is.null(reference_taxon)) {
    reference_taxon <- taxa_cols[which.max(colMeans(X))]
  }

  if (!reference_taxon %in% taxa_cols) {
    stop("reference_taxon must be one of the taxa columns.")
  }

  message("ANCOM-style ALR reference taxon: ", reference_taxon)

  X_pc <- X + pseudocount
  ref <- X_pc[, reference_taxon]

  X_alr <- log(sweep(X_pc, 1, ref, "/"))
  X_alr[, reference_taxon] <- 0

  out <- as.data.frame(X_alr, check.names = FALSE)
  colnames(out) <- taxa_cols
  rownames(out) <- rownames(df)
  out
}

# ============================================================
# Prepare one aligned method block
# ============================================================
prepare_method_block <- function(
    normalized_df,
    taxa_cols,
    method_number,
    groups,
    observations
) {
  normalized_df <- as.data.frame(normalized_df, check.names = FALSE)
  normalized_df <- normalized_df[, taxa_cols, drop = FALSE]

  if (nrow(normalized_df) != length(groups)) {
    stop(
      "Method ", method_number, " contains ", nrow(normalized_df),
      " rows, but ", length(groups), " group labels were supplied."
    )
  }

  if (nrow(normalized_df) != length(observations)) {
    stop(
      "Method ", method_number, " contains ", nrow(normalized_df),
      " rows, but ", length(observations), " observation IDs were supplied."
    )
  }

  out <- data.frame(
    normalized_df,
    Methods = rep(method_number, nrow(normalized_df)),
    Groups = groups,
    observation = observations,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  rownames(out) <- NULL
  out
}

# ============================================================
# Apply the ten procedures once to the real dataset
# ============================================================

raw_data <- taxa_data

clr_plus1 <- clr_only_taxa(
  analysis_data,
  taxa_cols,
  pseudocount = 1
)

clr_bmr <- clr_bmr_only_taxa(
  analysis_data,
  taxa_cols
)

tss_data <- tss_only_taxa(
  analysis_data,
  taxa_cols
)

rarefied_data <- rarefy_only_taxa(
  analysis_data,
  taxa_cols
)

css_data <- css_only_taxa(
  analysis_data,
  taxa_cols
)

tmm_data <- edger_tmm_cpm_only_taxa(
  analysis_data,
  taxa_cols
)

deseq2_data <- deseq_norm_only_taxa(
  analysis_data,
  taxa_cols,
  sf_type = "poscounts"
)

aldex2_data <- aldex_clr_median_only_taxa(
  analysis_data,
  taxa_cols,
  group_col = "Groups",
  mc.samples = 128,
  denom = "all"
)

ancom_data <- ancom_alr_only_taxa(
  analysis_data,
  taxa_cols,
  group_col = "Groups",
  pseudocount = 1
)

# ============================================================
# Reattach the saved Groups column and combine all blocks
#
# Method mapping:
#  1 Raw data
#  2 CLR+1
#  3 CLR-BMR
#  4 TSS
#  5 Rarefaction
#  6 CSS
#  7 edgeR-TMM
#  8 DESeq2
#  9 ALDEx2
# 10 ANCOM-style ALR
# ============================================================

raw_block <- prepare_method_block(
  raw_data, taxa_cols, 1, Groups_original, observation_original
)

clr_plus1_block <- prepare_method_block(
  clr_plus1, taxa_cols, 2, Groups_original, observation_original
)

clr_bmr_block <- prepare_method_block(
  clr_bmr, taxa_cols, 3, Groups_original, observation_original
)

tss_block <- prepare_method_block(
  tss_data, taxa_cols, 4, Groups_original, observation_original
)

rarefaction_block <- prepare_method_block(
  rarefied_data, taxa_cols, 5, Groups_original, observation_original
)

css_block <- prepare_method_block(
  css_data, taxa_cols, 6, Groups_original, observation_original
)

tmm_block <- prepare_method_block(
  tmm_data, taxa_cols, 7, Groups_original, observation_original
)

deseq2_block <- prepare_method_block(
  deseq2_data, taxa_cols, 8, Groups_original, observation_original
)

aldex2_block <- prepare_method_block(
  aldex2_data, taxa_cols, 9, Groups_original, observation_original
)

ancom_block <- prepare_method_block(
  ancom_data, taxa_cols, 10, Groups_original, observation_original
)

combined_10_methods <- rbind(
  raw_block,
  clr_plus1_block,
  clr_bmr_block,
  tss_block,
  rarefaction_block,
  css_block,
  tmm_block,
  deseq2_block,
  aldex2_block,
  ancom_block
)

rownames(combined_10_methods) <- NULL

# Enforce final order:
# taxa, Methods, Groups, observation
combined_10_methods <- combined_10_methods[
  ,
  c(taxa_cols, "Methods", "Groups", "observation"),
  drop = FALSE
]

# ============================================================
# Final validation
# ============================================================

expected_total_rows <- 10 * n_samples

if (nrow(combined_10_methods) != expected_total_rows) {
  stop(
    "Combined output contains ", nrow(combined_10_methods),
    " rows; expected ", expected_total_rows, "."
  )
}

rows_per_method <- table(combined_10_methods$Methods)

if (length(rows_per_method) != 10 ||
    any(rows_per_method != n_samples)) {
  stop("Each method must contain exactly ", n_samples, " samples.")
}

for (method_id in 1:10) {
  method_rows <- combined_10_methods[
    combined_10_methods$Methods == method_id,
    ,
    drop = FALSE
  ]

  if (!identical(as.character(method_rows$Groups), as.character(Groups_original))) {
    stop("Group order changed in method block ", method_id, ".")
  }

  if (!identical(method_rows$observation, observation_original)) {
    stop("Observation order changed in method block ", method_id, ".")
  }

  if (!identical(colnames(method_rows)[seq_along(taxa_cols)], taxa_cols)) {
    stop("Taxa-column order changed in method block ", method_id, ".")
  }
}

write.csv(
  combined_10_methods,
  file = output_file,
  row.names = FALSE,
  quote = TRUE
)

message("Combined real-data file saved to: ", normalizePath(output_file, mustWork = FALSE))
message("Final dimensions: ", nrow(combined_10_methods), " rows x ", ncol(combined_10_methods), " columns")

