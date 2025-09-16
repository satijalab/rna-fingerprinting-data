.registry <- data.frame(
  name = c('gwps', 'B_cell_immune_dictionary', 'CD4+_T_cell_immune_dictionary',
           'CD8+_T_cell_immune_dictionary', 'cDC1_immune_dictionary',
           'cDC2_immune_dictionary', 'eTAC_immune_dictionary',
           'ILC_immune_dictionary', 'Ki67+_T_cell_immune_dictionary',
           'Langerhans_immune_dictionary', 'LEC_immune_dictionary',
           'Macrophage_immune_dictionary', 'Mast_cell_immune_dictionary',
           'MigDC_immune_dictionary', 'Monocyte_immune_dictionary',
           'NK_cell_immune_dictionary', 'pDC_immune_dictionary',
           'Treg_immune_dictionary', 'γδ_T_cell_immune_dictionary'),
  description = c('K562 Genome-Wide Perturb-Seq (Replogle et al., 2022)',
                  'B cell Immune Dictionary (Cui et al., 2024)',
                  'CD4+ T cell Immune Dictionary (Cui et al., 2024)',
                  'CD8+ T cell Immune Dictionary (Cui et al., 2024)',
                  'cDC1 Immune Dictionary (Cui et al., 2024)',
                  'cDC2 Immune Dictionary (Cui et al., 2024)',
                  'eTAC Immune Dictionary (Cui et al., 2024)',
                  'ILC Immune Dictionary (Cui et al., 2024)',
                  'Ki67+ T cell Immune Dictionary (Cui et al., 2024)',
                  'Langerhans Immune Dictionary (Cui et al., 2024)',
                  'LEC Immune Dictionary (Cui et al., 2024)',
                  'Macrophage Immune Dictionary (Cui et al., 2024)',
                  'Mast cell Immune Dictionary (Cui et al., 2024)',
                  'MigDC Immune Dictionary (Cui et al., 2024)',
                  'Monocyte Immune Dictionary (Cui et al., 2024)',
                  'NK cell Immune Dictionary (Cui et al., 2024)',
                  'pDC Immune Dictionary (Cui et al., 2024)',
                  'Treg Immune Dictionary (Cui et al., 2024)',
                  'γδ T cell Immune Dictionary (Cui et al., 2024)'),
  package = c('gwps', 'B.cell.immune.dictionary', 'CD4.T.cell.immune.dictionary',
              'CD8.T.cell.immune.dictionary', 'cDC1.immune.dictionary',
              'cDC2.immune.dictionary', 'eTAC.immune.dictionary', 'ILC.immune.dictionary',
              'Ki67.T.cell.immune.dictionary', 'Langerhans.immune.dictionary',
              'LEC.immune.dictionary', 'Macrophage.immune.dictionary', 
              'Mast.cell.immune.dictionary', 'MigDC.immune.dictinary', 
              'Monocyte.immune.dictionary', 'NK.cell.immune.dictionary',
              'pDC.immune.dictionary', 'Treg.immune.dictionary',
              'gd.T.cell.immune.dictionary'),
  tarball = c('https://github.com/igrabski/fingerprinting.data.internal/releases/download/v1.0.0/gwps_0.0.0.9000.tar.gz',
              'https://github.com/igrabski/fingerprinting.data.internal/releases/download/v1.0.0/B.cell.immune.dictionary_0.0.0.9000.tar.gz',
              'https://github.com/igrabski/fingerprinting.data.internal/releases/download/v1.0.0/CD4.T.cell.immune.dictionary_0.0.0.9000.tar.gz',
              'https://github.com/igrabski/fingerprinting.data.internal/releases/download/v1.0.0/CD8.T.cell.immune.dictionary_0.0.0.9000.tar.gz',
              'https://github.com/igrabski/fingerprinting.data.internal/releases/download/v1.0.0/cDC1.immune.dictionary_0.0.0.9000.tar.gz',
              'https://github.com/igrabski/fingerprinting.data.internal/releases/download/v1.0.0/cDC2.immune.dictionary_0.0.0.9000.tar.gz',
              'https://github.com/igrabski/fingerprinting.data.internal/releases/download/v1.0.0/eTAC.immune.dictionary_0.0.0.9000.tar.gz',
              'https://github.com/igrabski/fingerprinting.data.internal/releases/download/v1.0.0/ILC.immune.dictionary_0.0.0.9000.tar.gz',
              'https://github.com/igrabski/fingerprinting.data.internal/releases/download/v1.0.0/Ki67.T.cell.immune.dictionary_0.0.0.9000.tar.gz',
              'https://github.com/igrabski/fingerprinting.data.internal/releases/download/v1.0.0/Langerhans.immune.dictionary_0.0.0.9000.tar.gz',
              'https://github.com/igrabski/fingerprinting.data.internal/releases/download/v1.0.0/LEC.immune.dictionary_0.0.0.9000.tar.gz',
              'https://github.com/igrabski/fingerprinting.data.internal/releases/download/v1.0.0/Macrophage.immune.dictionary_0.0.0.9000.tar.gz',
              'https://github.com/igrabski/fingerprinting.data.internal/releases/download/v1.0.0/Mast.cell.immune.dictionary_0.0.0.9000.tar.gz',
              'https://github.com/igrabski/fingerprinting.data.internal/releases/download/v1.0.0/MigDC.immune.dictionary_0.0.0.9000.tar.gz',
              'https://github.com/igrabski/fingerprinting.data.internal/releases/download/v1.0.0/Monocyte.immune.dictionary_0.0.0.9000.tar.gz',
              'https://github.com/igrabski/fingerprinting.data.internal/releases/download/v1.0.0/NK.cell.immune.dictionary_0.0.0.9000.tar.gz',
              'https://github.com/igrabski/fingerprinting.data.internal/releases/download/v1.0.0/pDC.immune.dictionary_0.0.0.9000.tar.gz',
              'https://github.com/igrabski/fingerprinting.data.internal/releases/download/v1.0.0/Treg.immune.dictionary_0.0.0.9000.tar.gz',
              'https://github.com/igrabski/fingerprinting.data.internal/releases/download/v1.0.0/gd.T.cell.immune.dictionary_0.0.0.9000.tar.gz'),
  stringsAsFactors = FALSE
)

#' List available datasets
#'
#' Returns a table of dataset names and descriptions
#'
#' @return A data frame with columns \code{name} and \code{description}
#' @export
AvailableData <- function() {
  data_info <- .registry[,c('name','description')]

  print(data_info, row.names = FALSE)
}

#' Install precomputed fingerprints
#' 
#' Installs requested precomputed fingerprints
#' 
#' @param name Name of dataset; see available datasets with AvailableData()
#' 
#' @return No return value; installs requested fingerprints
#' @export
InstallPrecomputedFingerprints <- function(name) {
  entry <- .registry[.registry$name == name, ]
  if (nrow(entry) == 0) stop("Unknown fingerprints: ", name)
  if (!requireNamespace(entry$package, quietly = TRUE)) {
    message("Installing ", name, " from GitHub release...")
    remotes::install_url(entry$tarball)
  } else {
    message(name, " is already installed.")
  }
}

#' @keywords internal
count_numbered_datasets <- function(name, pkg, max_check = 100) {
  count <- 0
  for (i in 1:max_check) {
    obj_name <- paste0(name, "_fingerprints", i)
    tmp_env <- new.env()
    loaded <- suppressWarnings(utils::data(list = obj_name,
                                           package = pkg,
                                           envir = tmp_env))
    
    if (obj_name %in% loaded && exists(obj_name, envir = tmp_env)) {
      count <- count + 1
    } else {
      break
    }
  }
  
  if (count == 0) {
    obj_name <- paste0(name, "_fingerprints")
    tmp_env <- new.env()
    loaded <- suppressWarnings(utils::data(list = obj_name,
                                           package = pkg,
                                           envir = tmp_env))
    if (obj_name %in% loaded && exists(obj_name, envir = tmp_env)) {
      count <- 1
    }
  }
  
  return(count)
}

#' Load pre-computed fingerprints
#'
#' Load in pre-computed fingerprints.
#'
#' @param name Name of pre-computed fingerprints
#' @param gene_key Formatting of genes (default ENSEMBL); ignored for mouse datasets
#'
#' @return A set of pre-computed fingerprints
#'
#' @importFrom fingerprinting.internal ConvertGenes
#' @export
#'
#' @examples
#' \dontrun{
#' gwps_fingerprints <- LoadPrecomputedFingerprints(
#' name = 'gwps'
#' )
#'
#' gwps_fingerprints <- LoadPrecomputedFingerprints(
#' name = 'gwps',
#' gene_key = 'SYMBOL'
#' )
#' }
#'
LoadPrecomputedFingerprints <- function(
    name,
    gene_key = 'ENSEMBL'
) {
  # Check for package
  pkg <- .registry$package[.registry$name == name]
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Fingerprints not installed. Use InstallPrecomputedFingerprints('", name, "') first.")
  }
  num <- count_numbered_datasets(name, pkg)
  
  # Get fingerprints
  fingerprints_comps <- lapply(1:num, function(i) {
    comp <- if (num == 1) paste0(name, "_fingerprints") else paste0(name, "_fingerprints", i)
    suppressWarnings(data(comp, package = pkg))
    get(comp)
  })
  fingerprints_Lambdas <- do.call(
    cbind, lapply(1:num,function(i)
      fingerprints_comps[[i]][[1]])
  )
  fingerprints_dists <- do.call(
    c, lapply(1:num,function(i)
      fingerprints_comps[[i]][[2]])
  )
  fingerprints <- list(
    Lambdas = fingerprints_Lambdas,
    dists = fingerprints_dists
  )
  
  # Get order
  hc_name <- paste0(name,"_hc")
  suppressWarnings(data(hc_name, package = pkg))
  hc <- get(hc_name)
  
  # Get pseudobulked data
  data_comps <- lapply(1:num,function(i) {
    comp <- if (num == 1) paste0(name, "_data") else paste0(name, "_data", i)
    suppressWarnings(data(comp, package = pkg))
    get(comp)
  })
  data <- do.call(
    cbind, lapply(1:num,function(i)
      data_comps[[i]])
  )
  
  data <- CreateSeuratObject(counts = data)
  data <- SetAssayData(data, layer = 'scale.data',
                       new.data = GetAssayData(data, layer = 'counts'))
  
  # Put into dictionary object
  dictionary <- structure(
    list(
      fingerprints = fingerprints,
      order = hc,
      data = data
    ),
    class = "dictionary"
  )
  
  if (gene_key != 'ENSEMBL') {
    dictionary <- fingerprinting.internal::ConvertGenes(
      dictionary,
      from = 'ENSEMBL',
      to = gene_key
    )
  }
  
  return(dictionary)
}