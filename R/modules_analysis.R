


## cluster methods






## kegg analysis / go analysis

gprofiler2::gost(cls[[23]],sources = "KEGG",organism = "hsapiens")$result

bio_enrich <- function(module, custom_gmt = NULL, ...) {
  # Checks
  if (is.list(module)) {
    if (any(!unlist(lapply(module, is.vector, "character")))) {
      stop("If module is a list of modules, all elements of the list must be ",
           "vectors of gene names") }
    if (is.null(names(module)))
      warning("No name provided for the list of modules.")
  } else if (!is.vector(module, "character")) {
    stop("module must be either a list of modules or a single module ",
         "represented by a vector of gene names")
  } else if (!is.list(module) & length(module) < 2) {
    warning("module represent a single module and only contains one gene name")
  }
  if (!is.null(custom_gmt)) {
    if (is.list(custom_gmt)) {
      if (any(lapply(custom_gmt, function(x) !is.character(x) |
                     length(x) != 1 ) %>% unlist)) {
        stop("all element of the list must be string reprensenting paths")}
      if (!all(lapply(custom_gmt, file.exists)))
        stop("all custom_gmt path provided must exist")
    } else if (is.character(custom_gmt) & length(custom_gmt) == 1) {
      if (!file.exists(custom_gmt))
        stop("custom_gmt path provided does not exists")
    } else stop("custom_gmt must be a path or a list of path to gmt file(s)")
  }

  # Enrichment with gprofiler internal datasets
  enriched_modules <- gprofiler2::gost(query = module, ...)

  # Enrichment with custom gmt if provided
  if (!is.null(custom_gmt)) {
    if (!is.list(custom_gmt)) custom_gmt <- list(custom_gmt)
    list_res_custom_gmts <- lapply(custom_gmt, function(gmt) {
      gmt_id <- quiet(gprofiler2::upload_GMT_file(gmt))
      enriched_modules_gmt <- quiet(gprofiler2::gost(query = module,
                                                     organism = gmt_id, ...))
      if (is.null(enriched_modules_gmt))
        warning("No enrichment found on gmt ", gmt)
      return(enriched_modules_gmt)
    })

    # If no enrichment with custom_gmt, returning only classic gost enrichment
    if (all(lapply(list_res_custom_gmts, is.null) %>% unlist)) {
      warning("None of the custom_gmt file provided returned an enrichement")
      return(enriched_modules)
    }

    # Removing NULL output from the list to merge
    # (exist when at least one of the gmt provided return no enrichment for
    # any module)
    list_res_custom_gmts <- purrr::compact(list_res_custom_gmts)

    # Joining results
    enriched_modules <- join_gost(c(list(enriched_modules),
                                    list_res_custom_gmts))
  }

  return(enriched_modules)
}

## return a list，代表多个模块的某一类富集分析的结果，一个元素代表一个模块

























