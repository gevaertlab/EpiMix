#' @importFrom foreach foreach %dopar%
NULL

#' @importFrom parallel makeCluster stopCluster
NULL

#' The extractPriMiRNA function
#' @description Utility function to convert mature miRNA names to pri-miRNA names
#' @param str a character string for a mature miRNA name (e.g. "hsa-miR-34a-3p")
#' @return a character string for the corresponding pri-miRNA name (e.g. "hsa-mir-34a")
#'
.extractPriMiRNA <- function(str){
  # Split the string into elements separated by "-" character
  elements <- strsplit(str, "-")[[1]]
  # If the string contains less than three elements, return the original string
  if(length(elements) < 3){
    return(str)
  }

  # Replace "miR" to 'mir'
  if(elements[2] == "miR"){
    elements[2] <- "mir"
  }

  # Otherwise, paste together the first three elements with "-"
  new_str <- paste(elements[1:3], collapse = "-")
  # Return the new string
  return(new_str)
}

#' The find_miRNA_targets function
#' @description Detection potential target protein-coding genes for the differentially methylated miRNAs using messenger RNA expression data
#' @param EpiMixResults List of the result objects returned from the EpiMix function.
#' @param geneExprData Matrix of the messenger RNA expression data with genes in rows and samples in columns.
#' @param database character string indicating the database for retrieving miRNA targets. Default: "mirtarbase".
#' @param raw.pvalue.threshold Numeric value indicating the threshold of the raw P value for selecting the miRNA targets based on gene expression. Default: 0.05.
#' @param adjusted.pvalue.threshold Numeric value indicating the threshold of the adjusted P value for selecting the miRNA targets based on gene expression. Default: 0.2.
#' @return Matrix indicating the miRNA-target pairs, with fold changes of target gene expression and P values.
#' @param cores Number of CPU cores to be used for computation. Default: 1.
#' @export
#' @examples
#' \donttest{
#' library(multiMiR)
#' library(miRBaseConverter)
#'
#' data(mRNA.data)
#' data(Sample_EpiMixResults_miRNA)
#'
#' miRNA_targets <- find_miRNA_targets(
#'  EpiMixResults = Sample_EpiMixResults_miRNA,
#'  geneExprData = mRNA.data
#' )
#' }
#'

find_miRNA_targets <- function(EpiMixResults,
                               geneExprData,
                               database = 'mirtarbase',
                               raw.pvalue.threshold = 0.05,
                               adjusted.pvalue.threshold = 0.2,
                               cores = 1){

  if (!requireNamespace("multiMiR")) {
    message("This function requires the 'multiMiR' package.")
    return(invisible())
  }

  if (!requireNamespace("miRBaseConverter")) {
    message("This function requires the 'miRBaseConverter' package.")
    return(invisible())
  }

  DMvalues <- EpiMixResults$MethylationStates
  functionalPairs <- EpiMixResults$FunctionalPairs
  targetMiRNAs <- unique(functionalPairs$Gene)
  matureMiRNAs <- miRBaseConverter :: miRNA_PrecursorToMature(targetMiRNAs)
  matureMiRNAs <- c(matureMiRNAs$Mature1, matureMiRNAs$Mature2)
  matureMiRNAs <- matureMiRNAs[!is.na(matureMiRNAs)]
  cat(paste0("Looking for target genes of ", length(matureMiRNAs), " miRNAs"))

  # Retrieve miRNA target genes
  multiR <- multiMiR :: get_multimir(org = 'hsa', mirna = matureMiRNAs, table = database, summary = TRUE)
  targetMap <- multiR@summary
  targetGenes <- unique(targetMap$target_symbol)
  cat(paste0("Found ", length(targetGenes), " target genes"))
  targetMap$pri_mirna_id <- sapply(targetMap$mature_mirna_id, function(str) .extractPriMiRNA(str))
  targetMap <- targetMap[, c("pri_mirna_id", "target_symbol")]
  targetMap <- unique(targetMap)

  # Filter the miRNAs with the target gene expression data available
  targetMap <- targetMap[targetMap$target_symbol %in% rownames(geneExprData), ]

  # Compare the gene expression for miRNA target genes
  MET_matrix  <- EpiMixResults$MethylationStates
  MET_matrix <- filterMethMatrix(MET_matrix = MET_matrix, control.names = EpiMixResults$group.2,
                                 gene.expression.data = geneExprData)

  functionalPairs <- functionalPairs[functionalPairs$Probe %in% rownames(MET_matrix), ]
  functionalPairs <- functionalPairs[functionalPairs$Gene %in% targetMap$pri_mirna_id, ]

  iterations <- nrow(functionalPairs)
  pb <- utils::txtProgressBar(max = iterations, style = 3)

  cat("Registering sockets on multiple CPU cores...\n")
  cl <- parallel::makeCluster(cores)
  registerDoSNOW(cl)

  finalResults <- data.frame()
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  finalResults <- foreach::foreach(i = seq_len(iterations), .combine = rbind, .options.snow = opts,
                                   .verbose = FALSE) %dopar% {
                                     probeResults <- data.frame()
                                     probe <- functionalPairs$Probe[i]
                                     miRNA <- functionalPairs$Gene[i]
                                     targetGenes <- targetMap$target_symbol[which(targetMap$pri_mirna_id == miRNA)]
                                     dmValues <- MET_matrix[probe, ]
                                     for(gene in targetGenes){
                                       geneExprVals <- geneExprData[gene,]
                                       geneResults <- test_gene_expr(gene, probe, dmValues, geneExprVals, correlation = "positive")
                                       probeResults <- rbind(probeResults, geneResults)
                                     }
                                     probeResults["Adjusted.p"] <- p.adjust(as.vector(unlist(probeResults["Raw.p"])), method = "fdr")
                                     probeResults['miRNA'] = miRNA
                                     probeResults
                                   }

  orderColumns <- c("Probe", "State", "miRNA", "Gene", "Fold change of gene expression", "Comparators", "Raw.p", "Adjusted.p")
  finalResults <- finalResults[, orderColumns]
  sigResults <- finalResults[finalResults['Raw.p'] < raw.pvalue.threshold & finalResults['Adjusted.p'] < adjusted.pvalue.threshold, ]
  rownames(sigResults) <- NULL

  close(pb)
  parallel::stopCluster(cl)
  return(sigResults)
}














