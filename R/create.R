#' Creata a FacileDataSet from the individual assay lists
#' 
#' @export
#' @param name The name of the datasets. Defaults to
#'   `"TestMultiModalFacileDataSet"`. Users may create different versions to
#'   test (ie. different assays to include, etc) and may want a more
#'   recognizable name for it.
#' @param path Local path to create the FacileDataSet directory. If `NULL`
#'   (default) a temporary directory will be created via`tempfile()` and will
#'   include `name` as the directory's prefix.
#' @param assays the nums of the assays to include. If `NULL` (default), all
#'   assays will be included.
create <- function(name = "TestMultiModalFacileDataSet",
                   path = NULL, assays = NULL) {
  if (FALSE) {
    name <- "TestMultiModalFacileDataSet"
    path <- NULL
    assays <- NULL
  }
  assert_string(name)
  if (is.null(path)) path <- tempfile(paste0(name, "_"))
  assert_directory(dirname(path), "w")
  
  all_assays <- available_assays()
  if (is.null(assays)) assays <- all_assays$assay_name
  assays <- assert_subset(tolower(assays), tolower(all_assays$assay_name))
  
  adata <- all_assays |> 
    mutate(assay_name = tolower(assay_name)) |> 
    semi_join(tibble(assay_name = assays), by = "assay_name")
  
  ainfo <- adata[1,]
  adat <- assay_lists_load(ainfo$assay_name)
  
  fds <- FacileData::as.FacileDataSet(
    adat,
    path = path,
    dataset_name = name,
    assay_name = ainfo$assay_name,
    assay_type = ainfo$assay_type,
    assay_description = ainfo$description,
    organism = "Homo sapiens")
  
  if (nrow(adata) > 1L) {
    for (i in 2:nrow(adata)) {
      ainfo <- adata[i,]
      adat <- assay_lists_load(ainfo$assay_name)
      
      if (is(adat[[1L]], "DGEList")) {
        assay_name <- "counts"
      } else if (is(adat[[1L]], "EList")) {
        assay_name <- "E"
      } else {
        stop("No handler for class yet: ", class(adat[[1]])[1L])
      }
      message("Adding assay: ", ainfo$assay_name)
      FacileData::addFacileAssaySet(
        fds,
        adat,
        facile_assay_name = ainfo$assay_name,
        facile_assay_type = ainfo$assay_type,
        facile_feature_type = ainfo$feature_type,
        facile_assay_description = ainfo$description,
        facile_feature_info = adat[[1]]$genes,
        storage_mode = ainfo$storage_mode,
        assay_name = assay_name)
    }
  }
}