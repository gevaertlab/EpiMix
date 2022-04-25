#' The TCGA_GetData function
#' @description This function wraps the functions for downloading, pre-processing and analysis of the DNA methylation and gene expression data from the TCGA project.
#' @param CancerSite character string indicating the TCGA cancer code. The information can be found at: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
#' @param mode character string indicating the analytic mode to model DNA methylation. Should be one of the followings: "Regular", "Enhancer", "miRNA" or "lncRNA". Default: "Regular". See details for more information.
#' @param outputDirectory character string indicating the file path to save the output.
#' @param doBatchCorrection logical indicating whether to do batch effect correction during preprocessing. Default: False.
#' @param batch.correction.method character string indicating the method to perform batch effect correction. The value should be either "Seurat" or "Combat". Seurat is much fatster than the Combat. Default: "Seurat".
#' @param roadmap.epigenome.ids character vector indicating the epigenome ID(s) to be used for selecting enhancers. See details for more information. Default: NULL.
#' @param roadmap.epigenome.groups character vector indicating the tissue group(s) to be used for selecting enhancers. See details for more information. Default: NULL.
#' @param cores Number of CPU cores to be used for computation.
#' @export
#' @return The results from EpiMix is a list with the following components:
#' \item{MethylationDrivers}{CpG probes identified as differentially methylated by EpiMix.}
#' \item{NrComponents}{The number of methylation states found for each driver probe.}
#' \item{MixtureStates}{A list with the DM-values for each driver probe.
#' Differential Methylation values (DM-values) are defined as the difference between
#' the methylation mean of samples in one mixture component from the experiment group and the methylation mean
#' in samples from the control group, for a given probe.}
#' \item{MethylationStates}{Matrix with DM-values for all driver probes (rows) and all samples (columns).}
#' \item{Classifications}{Matrix with integers indicating to which mixture component each sample in the experiment group was assigned to, for each probe.}
#' \item{Models}{Beta mixture model parameters for each driver probe.}
#' \item{group.1}{sample names in group.1 (experimental group).}
#' \item{group.2}{sample names in group.2 (control group).}
#' \item{FunctionalPairs}{Dataframe with the prevalence of differential methyaltion for each CpG probe in the sample population, and fold change of mRNA expression and P values for each signifcant probe-gene pair.}
#' @details
#' mode:
#' EpiMix incorporates four alternative analytic modes for modeling DNA methylation: “Regular,” “Enhancer”, “miRNA” and “lncRNA”.
#' The four analytic modes target DNA methylation analysis on different genetic elements.
#' The Regular mode aims to model DNA methylation at proximal cis-regulatory elements of protein-coding genes.
#' The Enhancer mode targets DNA methylation analysis on distal enhancers.
#' The miRNA or lncRNA mode focuses on methylation analysis of miRNA- or lncRNA-coding genes.
#'
#' roadmap.epigenome.groups & roadmap.epigenome.ids:
#'
#' Since enhancers are cell-type or tissue-type specific, EpiMix needs to know the reference tissues or cell types in order to select proper enhancers.
#' EpiMix identifies enhancers from the RoadmapEpigenomic project (Nature, PMID: 25693563), in which enhancers were identified by ChromHMM in over 100 tissue and cell types.
#' Available epigenome groups (a group of relevant cell types) or epigenome ids (individual cell types) can be obtained from the original publication (Nature, PMID: 25693563, figure 2).
#' They can also be retrieved from the list.epigenomes() function. If both roadmap.epigenome.groups and roadmap.epigenome.ids are specified, EpiMix will select all the epigenomes from the combination of the inputs.
#' @examples
#' \dontrun{
#' # Example #1 - Regular mode
#' CancerSite <- "OV"
#' mode <- "Regular"
#' outputDirectory <- paste0(getwd(), "/")
#' EpiMixResults <- TCGA_GetData(CancerSite = CancerSite,
#'                               mode = mode,
#'                               outputDirectory = outputDirectory)
#'
#' # Example #2 - Enhancer mode
#' CancerSite <- "OV"
#' mode <- "Enhancer"
#' outputDirectory <- paste0(getwd(), "/")
#' roadmap.epigenome.ids = "E097"
#' EpiMixResults <- TCGA_GetData(CancerSite = CancerSite,
#'                               mode = mode,
#'                               roadmap.epigenome.ids = roadmap.epigenome.ids,
#'                               outputDirectory = outputDirectory)
#'
#' Example #3 - miRNA mode
#' CancerSite <- "OV"
#' mode <- "miRNA"
#' outputDirectory <- paste0(getwd(), "/")
#' EpiMixResults <- TCGA_GetData(CancerSite = CancerSite,
#'                               mode = mode,
#'                               outputDirectory = outputDirectory)
#'
#' #' Example #4 - lncRNA mode
#' CancerSite <- "OV"
#' mode <- "lncRNA"
#' outputDirectory <- paste0(getwd(), "/")
#' EpiMixResults <- TCGA_GetData(CancerSite = CancerSite,
#'                               mode = mode,
#'                               outputDirectory = outputDirectory)
#'
#' }
#'
TCGA_GetData <- function(CancerSite,
                         mode = "Regular",
                         outputDirectory = ".",
                         doBatchCorrection = FALSE,
                         batch.correction.method = "Seurat",
                         roadmap.epigenome.ids = NULL,
                         roadmap.epigenome.groups = NULL,
                         cores = 1) {
  # ---------------------------------------------------------------------------------------------
  # Step 1: Download and preprocess DNA methylation data
  # ---------------------------------------------------------------------------------------------
  # Downloading methylation data
  cat("Downloading methylation data for:", CancerSite, "\n")
  METdirectories <- TCGA_Download_DNAmethylation(CancerSite, outputDirectory)

  # Preprocess methylation data
  cat("Processing methylation data for:", CancerSite, "\n")
  METProcessedData <- TCGA_Preprocess_DNAmethylation(CancerSite,
                                                     METdirectories,
                                                     doBatchCorrection = doBatchCorrection,
                                                     batch.correction.method = batch.correction.method)

  cat("Saving methylation processed data for:", CancerSite, "\n")
  saveRDS(METProcessedData, paste0(outputDirectory, "/", "MET_", CancerSite, "_Processed.rds" ))

  # ---------------------------------------------------------------------------------------------
  # Step 2: Download and preprocess gene expression data
  # ---------------------------------------------------------------------------------------------
  # Downloading gene expression data
  GEdirectories <- TCGA_Download_GeneExpression(CancerSite,
                                                outputDirectory,
                                                mode = mode
                                                )

  # Processing gene expression data
  GEProcessedData <- TCGA_Preprocess_GeneExpression(CancerSite,
                                                    GEdirectories,
                                                    mode = mode,
                                                    doBatchCorrection = doBatchCorrection,
                                                    batch.correction.method =  batch.correction.method,
                                                    cores = cores
                                                    )
  # Saving gene expression processed data
  saveRDS(GEProcessedData, file = paste0(outputDirectory, "/", "GE_", CancerSite, "_Processed_", mode, ".rds"))

  # ---------------------------------------------------------------------------------------------
  # Step 3: Generate sample information
  # ---------------------------------------------------------------------------------------------
  sample.info = TCGA_GetSampleInfo(METProcessedData, CancerSite = CancerSite, TargetDirectory = outputDirectory)

  # ---------------------------------------------------------------------------------------------
  # Step 4: Run EpiMix
  # ---------------------------------------------------------------------------------------------
  EpiMixResults <- EpiMix(methylation.data = METProcessedData,
                          gene.expression.data = GEProcessedData,
                          mode = mode,
                          sample.info = sample.info,
                          group.1 = "Cancer",
                          group.2 = "Normal",
                          met.platform = "HM450",
                          genome = "hg38",
                          cores = cores,
                          roadmap.epigenome.groups = roadmap.epigenome.groups,
                          roadmap.epigenome.ids = roadmap.epigenome.ids,
                          OutputRoot = outputDirectory
  )
  return(EpiMixResults)
}
