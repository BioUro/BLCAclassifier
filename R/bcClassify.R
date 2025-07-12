
# ============================================================================
# BLCAclassifier: Bladder Cancer Immune Subtyping R Package
# Bug-Fixed Version
# ============================================================================

#' Bladder Cancer Immune Subtype Classification
#'
#' Classifies bladder cancer samples into immune subtypes using gene expression signatures
#'
#' @param expr_data Expression matrix with genes as rows and samples as columns
#' @return A character vector containing classification results for each sample
#' @export
#' @importFrom stats cor
#' @examples
#' \dontrun{
#' # Load expression data
#' expr_data <- read.csv("expression_data.csv", row.names = 1)
#' # Perform classification
#' results <- bcClassify(expr_data)
#' print(table(results))
#' }
bcClassify <- function(expr_data) {

  # Data validation
  if (!is.matrix(expr_data) && !is.data.frame(expr_data)) {
    stop("Expression data must be a matrix or data frame")
  }

  expr_data <- as.matrix(expr_data)
  cat("Number of samples:", ncol(expr_data), "\n")

  # Load template
  template_df <- load_ntp_template()

  # Gene matching with multiple strategies
  expr_genes <- rownames(expr_data)
  template_genes <- template_df$probe

  # Try direct matching first
  common_genes <- intersect(expr_genes, template_genes)

  # If few matches, try removing version numbers
  if (length(common_genes) < 50) {
    expr_genes_clean <- gsub("\\.[0-9]+$", "", expr_genes)
    common_genes <- intersect(expr_genes_clean, template_genes)
    if (length(common_genes) >= 50) {
      rownames(expr_data) <- expr_genes_clean
    }
  }

  # If still few matches, try case-insensitive
  if (length(common_genes) < 50) {
    expr_genes_upper <- toupper(rownames(expr_data))
    template_genes_upper <- toupper(template_genes)
    common_genes <- intersect(expr_genes_upper, template_genes_upper)
    if (length(common_genes) >= 50) {
      rownames(expr_data) <- toupper(rownames(expr_data))
      template_df$probe <- toupper(template_df$probe)
    }
  }

  cat("Matched genes:", length(common_genes), "/", length(template_genes), "\n")

  if (length(common_genes) < 5) {
    stop("Too few matched genes for reliable classification. Please check gene name format.")
  }

  # Filter data to matched genes
  expr_filtered <- expr_data[common_genes, , drop = FALSE]
  template_filtered <- template_df[template_df$probe %in% common_genes, ]

  # Check if all subtypes have at least one gene
  subtypes <- unique(template_filtered$class)
  cat("Available subtypes:", paste(subtypes, collapse = ", "), "\n")

  for (subtype in subtypes) {
    subtype_gene_count <- sum(template_filtered$class == subtype)
    cat("  -", subtype, ":", subtype_gene_count, "genes\n")
  }

  # Standardize expression data
  expr_scaled <- t(scale(t(expr_filtered), center = TRUE, scale = TRUE))

  # Handle potential NA values
  expr_scaled[is.na(expr_scaled)] <- 0

  # Perform classification
  predictions <- ntp_classify(expr_scaled, template_filtered)

  cat("Classification completed!\n")
  print(table(predictions))

  return(predictions)
}

#' Load NTP Template
#'
#' Internal function to load NTP template from package data
#'
#' @return Template data frame with probe and class columns
load_ntp_template <- function() {

  template_file <- system.file("extdata", "templates.qs", package = "BLCAclassifier")

  if (!file.exists(template_file)) {
    stop("Template file not found. Please ensure templates.qs is in inst/extdata/ directory.")
  }

  if (!requireNamespace("qs", quietly = TRUE)) {
    stop("qs package is required. Install with: install.packages('qs')")
  }

  template_data <- qs::qread(template_file)

  if (!is.data.frame(template_data) || !all(c("probe", "class") %in% colnames(template_data))) {
    stop("Invalid template format. Must be a data frame with 'probe' and 'class' columns.")
  }

  return(template_data)
}

#' NTP Classification Algorithm (Bug-Fixed)
#'
#' Internal function for gene signature-based classification
#'
#' @param expr_data Standardized expression matrix
#' @param template_df Template data frame
#' @return Vector of predicted subtypes
ntp_classify <- function(expr_data, template_df) {

  n_samples <- ncol(expr_data)
  sample_names <- colnames(expr_data)
  subtypes <- unique(template_df$class)

  predictions <- character(n_samples)
  names(predictions) <- sample_names

  # Calculate subtype scores for each sample
  for (i in seq_len(n_samples)) {
    sample_expr <- expr_data[, i]

    subtype_scores <- numeric(length(subtypes))
    names(subtype_scores) <- subtypes

    for (j in seq_along(subtypes)) {
      subtype <- subtypes[j]
      subtype_genes <- template_df$probe[template_df$class == subtype]
      available_genes <- intersect(subtype_genes, rownames(expr_data))

      if (length(available_genes) > 0) {
        # Calculate mean expression of signature genes
        score <- mean(sample_expr[available_genes], na.rm = TRUE)
        subtype_scores[j] <- if (is.na(score)) -Inf else score
      } else {
        subtype_scores[j] <- -Inf
      }
    }

    # Check if all scores are invalid
    if (all(is.infinite(subtype_scores)) || all(is.na(subtype_scores))) {
      # Assign to the first subtype as fallback
      predictions[i] <- subtypes[1]
      warning("Sample ", sample_names[i], " could not be reliably classified, assigned to ", subtypes[1])
    } else {
      # Find the best subtype
      valid_scores <- subtype_scores[!is.infinite(subtype_scores) & !is.na(subtype_scores)]
      if (length(valid_scores) > 0) {
        best_idx <- which.max(subtype_scores)
        predictions[i] <- names(subtype_scores)[best_idx]
      } else {
        predictions[i] <- subtypes[1]
      }
    }
  }

  return(predictions)
}

#' BLCAclassifier Package
#'
#' Bladder Cancer Immune Subtype Classification
#'
#' @section Usage:
#' \preformatted{
#' library(BLCAclassifier)
#' results <- bcClassify(expression_matrix)
#' }
#'
#' @keywords internal
"_PACKAGE"
