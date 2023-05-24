# Utility functions used in wranging raw data into shape for FacileDataSet
# ingestion.

#' Splits a SingleCellExperiment in a list of DGELists, split by the levels
#' one of the colData columns
#' 
#' This function is called from the dataset-wrangling.Rmd vignette, and
#' relies on the SingleCellExperiment package (which is in Suggests, not
#' Imports)
#' 
#' @export
#' @param x A SingleCellExperiment or SummarizedExperiment
#' @param split.by a column name in `colData(x)` that has categorical covariate
#'   (factor or character)
pb_split <- function(x, split.by = "cond") {
  reqpkg("SingleCellExperiment")
  split.levels <- x[[split.by]]
  stopifnot(is.character(split.levels) || is.factor(split.levels))
  sapply(unique(split.levels), function(dsname) {
    pbx <- x[, split.levels == dsname]
    y <- edgeR::DGEList(
      counts = SingleCellExperiment::counts(pbx),
      genes = as.data.frame(SummarizedExperiment::rowData(pbx)),
      # group was already defined in sce.all
      # group = paste0(pbx$cond, "__", pbx$cell_abbrev)
      samples = as.data.frame(SingleCellExperiment::colData(pbx)))
    y <- edgeR::calcNormFactors(y)
    des <- model.matrix(~ group, y$samples)
    keep <- edgeR::filterByExpr(y, design = des, min.count = 5, 
                                min.total.count = 10)
    edgeR::calcNormFactors(y[keep,,keep.lib.sizes = FALSE])
  }, simplify = FALSE)
}

#' Description of available assays in this packge
#' @export
#' @return a tibble of assay information
available_assays <- function() {
  dplyr::tribble(
    ~assay_name, ~assay_type, ~feature_type, ~storage_mode,
    "scRNAseq",  "pseudobulk", "ensgid",     "integer",
    "snRNAseq",  "pseudobulk", "ensgid",     "integer")
}

#' Returns the filepath for the raw assay data for a given assay(name)
#' @export
#' @param assay_name One of the (case insensitive) options in
#'   `available_assays()$assay_name`
#' @return file path to the serialized assaylist
assay_list_file_path <- function(assay_name) {
  assert_string(assay_name)
  ainfo <- available_assays() |> 
    mutate(name = tolower(assay_name)) |> 
    filter(.data$name == tolower(.env$assay_name))
  if (nrow(ainfo) == 0L) {
    stop("no information found for assay_name: ", assay_name)
  }
  outdir <- system.file("extdata", "assay-data",
                        package = "TestMultiModalFacileDataSet")
  file.path(outdir, paste0(ainfo$assay_name, "-assay-list.qs"))
}

#' Loads an assay list object
#' @export
#' @inheritParams assay_list_file_path
assay_lists_load <- function(assay_name) {
  fn <- assay_list_file_path(assay_name)
  qs::qread(fn)
}
