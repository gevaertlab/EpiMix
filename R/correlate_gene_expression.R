# Utility functions to select differential DNA methylation that correlate with gene expression

#' @importFrom dplyr summarize group_by distinct mutate %>% if_else arrange
NULL

#' @importFrom tidyr separate_rows
NULL

#' @importFrom GenomicRanges makeGRangesFromDataFrame seqnames start end mcols
NULL

#' @importFrom S4Vectors queryHits subjectHits
NULL

#' @importFrom IRanges findOverlaps
NULL

#' @importFrom utils read.table
NULL

#' @importFrom stats pchisq coef confint
NULL

#' @importFrom downloader download
NULL

#' @importFrom rlang .data
NULL

#' @importFrom ExperimentHub ExperimentHub
NULL

#' @importFrom AnnotationHub query
NULL


#' The generateFunctionalPairs function
#' @description Wrapper function to get functional CpG-gene pairs
#' @param MET_matrix matrix of methylation states
#' @param MET_Control beta values of control groups
#' @param gene.expression.data matrix of gene expression data
#' @param ProbeAnnotation dataframe of probe annotation
#' @param raw.pvalue.threshold raw p value threshold
#' @param adjusted.pvalue.threshold  adjusted p value threshold
#' @param cores number of computational cores
#' @param mode character string indicating the analytic mode
#' @param correlation the expected relationship between DNAme and gene expression
#'
#' @return a dataframe of functional CpG-gene matrix

generateFunctionalPairs <- function(MET_matrix, MET_Control, gene.expression.data,
                                    ProbeAnnotation, raw.pvalue.threshold, adjusted.pvalue.threshold, cores, mode = "Regular", correlation = "negative") {
  MET_matrix <- filterMethMatrix(MET_matrix = MET_matrix, MET_Control = MET_Control,
                                 gene.expression.data = gene.expression.data)
  if (length(MET_matrix) == 0 | nrow(MET_matrix) == 0 | ncol(MET_matrix) == 0) {
    return(NULL)
  }
  ProbeAnnotation <- ProbeAnnotation[which(ProbeAnnotation$probe %in% rownames(MET_matrix)),
  ]

  uniqueGenes <- unique(ProbeAnnotation$gene)

  iterations <- length(uniqueGenes)
  pb <- utils::txtProgressBar(max = iterations, style = 3)

  FunctionalPairs <- data.frame()
  if (cores == "" | cores == 1) {
    for (i in seq_len(iterations)) {
      gene <- uniqueGenes[i]
      probes <- ProbeAnnotation$probe[which(ProbeAnnotation$gene == gene)]
      pairs <- getFunctionalProbes(gene, probes, MET_matrix, gene.expression.data,
                                   raw.pvalue.threshold = raw.pvalue.threshold, adjusted.pvalue.threshold = adjusted.pvalue.threshold, correlation = correlation)
      FunctionalPairs <- rbind(FunctionalPairs, pairs)
      utils::setTxtProgressBar(pb, i)
    }
  } else {
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    FunctionalPairs <- foreach::foreach(i = seq_len(iterations), .combine = rbind, .options.snow = opts,
                                        .verbose = FALSE) %dopar% {
                                          gene <- uniqueGenes[i]
                                          probes <- ProbeAnnotation$probe[which(ProbeAnnotation$gene == gene)]
                                          getFunctionalProbes(gene, probes, MET_matrix, gene.expression.data, raw.pvalue.threshold = raw.pvalue.threshold,
                                                              adjusted.pvalue.threshold = adjusted.pvalue.threshold, correlation = correlation)
                                        }
  }
  close(pb)
  return(FunctionalPairs)
}

#' The getFunctionalProbes function
#' @description Helper function to assess if the expression of a specific gene is reversely correlated with the methylation of a list of probes mapped to it.
#' @details This function is gene-centered, which is used in the regular mode and the lncRNA mode of EpiMix.
#' @param gene character string indicating the target gene to be modeled.
#' @param probes character vector indicating the probes mapped to the target gene.
#' @param MET_matrix methylation data matrix for CpGs from group.1 and group.2.
#' @param gene.expression.data gene expression data matrix.
#' @param correlation character indicating the direction of correlation between the methylation state of the CpG site and the gene expression levels. Can be either 'negative' or 'positive'.
#' @param raw.pvalue.threshold raw p value from testing DNA methylation and gene expression
#' @param adjusted.pvalue.threshold adjusted p value from testing DNA methylation and gene expression
#' @keywords internal
#'
#' @return dataframe with functional probe-gene pairs and corresponding p values obtained from the Wilcoxon test for gene expression and methylation.
#'
getFunctionalProbes <- function(gene, probes, MET_matrix, gene.expression.data, correlation = "negative",
                                raw.pvalue.threshold = 0.05, adjusted.pvalue.threshold = 0.01) {

  valid.probes <- c()
  mRNA.fold.change <- c()
  comparisons <- c()
  p <- c()

  for (probe in probes) {
    state <- fold_change <- cmp <- p_value <- NULL
    DM_values <- MET_matrix[probe, ]  # a single-row vector with the names as sample names, and the values as DM values
    state <- getMethStates_Helper(DM_values)
    if (state == "Hyper") {
      high.met.samples <- names(DM_values[DM_values > 0])
      low.met.samples <- names(DM_values[DM_values == 0])
      expr.low.values <- as.numeric(gene.expression.data[gene, high.met.samples])  # high methylation, low gene expression
      expr.high.values <- as.numeric(gene.expression.data[gene, low.met.samples])  # low methylation, high gene expression
      if (length(expr.high.values) < 3 || length(expr.low.values)< 3) next
      fold_change <- round(mean(expr.low.values)/mean(expr.high.values), 3)
      cmp <- "hyper vs normal"
      p_value <- wilcox.test(expr.high.values,
                                      expr.low.values,
                                      alternative = ifelse(correlation == "negative", "greater", "less"), exact = FALSE)$p.value

    } else if (state == "Dual") {
      high.met.samples <- names(DM_values[DM_values > 0])
      low.met.samples <- names(DM_values[DM_values < 0])
      expr.low.values <- as.numeric(gene.expression.data[gene, high.met.samples])  # high methylation, low gene expression
      expr.high.values <- as.numeric(gene.expression.data[gene, low.met.samples])  # low methylation, high gene expression
      if (length(expr.high.values) < 3 || length(expr.low.values)< 3) next
      fold_change <- round(mean(expr.high.values)/mean(expr.low.values), 3)
      cmp <- "hypo vs hyper"
      p_value <- wilcox.test(expr.high.values, expr.low.values, alternative = ifelse(correlation ==
                                                                                       "negative", "greater", "less"), exact = FALSE)$p.value
    } else {
      high.met.samples <- names(DM_values[DM_values == 0])
      low.met.samples <- names(DM_values[DM_values < 0])
      expr.low.values <- as.numeric(gene.expression.data[gene, high.met.samples])  # high methylation, low gene expression
      expr.high.values <- as.numeric(gene.expression.data[gene, low.met.samples])  # low methylation, high gene expression
      if (length(expr.high.values) < 3 || length(expr.low.values)< 3) next
      fold_change <- round(mean(expr.high.values)/mean(expr.low.values), 3)
      cmp <- "hypo vs normal"
      p_value <- wilcox.test(expr.high.values, expr.low.values, alternative = ifelse(correlation ==
                                                                                       "negative", "greater", "less"), exact = FALSE)$p.value
    }

    valid.probes <- append(valid.probes, probe)
    mRNA.fold.change <- append(mRNA.fold.change, fold_change)
    comparisons <- append(comparisons, cmp)
    p <- append(p, p_value)
  }

  # produce a dataframe for gene expression with methylation state and
  # prevalence information
  dataDEGs <- data.frame(Gene = rep(gene, length(valid.probes)), Probe = valid.probes)
  if(nrow(dataDEGs) == 0) return(NULL)
  dataDEGs["Fold change of gene expression"] <- mRNA.fold.change
  dataDEGs["Comparators"] <- comparisons
  dataDEGs["Raw.p"] <- p
  dataDEGs["Adjusted.p"] <- p.adjust(as.vector(unlist(dataDEGs["Raw.p"])), method = "fdr")
  dataDEGs <- dataDEGs[which(dataDEGs$Raw.p < raw.pvalue.threshold & dataDEGs$Adjusted.p <
                               adjusted.pvalue.threshold), ]
  return(dataDEGs)
}


#' The EpiMix_ModelGeneExpression function
#' @description use a linear regression filter to screen for probes that were negatively associated with gene expression.
#' @param methylation.data methylation data matrix.
#' @param gene.expression.data gene expression data matrix.
#' @param ProbeAnnotation dataframe of probe annotation
#' @param cores number of CPU cores used for computation
#' @param filter logical indicating whether to perform a linear regression to select functional probes
#' @param cluster logical indicating whether the CpGs were clustered using hierarchical clustering
#' @param correlation Character vector indicating the expected correlation between DNA methylation and gene expression. Can be either 'negative' or 'positive'. Default: 'negative'.
#' @return a character vector of probe names.
#' @import doSNOW
#' @keywords internal
#'
#'
EpiMix_ModelGeneExpression <- function(methylation.data, gene.expression.data, ProbeAnnotation,
                                       cores, filter, cluster, correlation = "negative") {
  FunctionalProbes <- character(0)
  # overlapping samples to select only samples with both methylation data and
  # gene expression data
  OverlapSamples <- intersect(colnames(methylation.data), colnames(gene.expression.data))
  cat("Found", length(OverlapSamples), "samples with both methylation and gene expression data.\n")
  gene.expression.data <- gene.expression.data[, OverlapSamples, drop = FALSE]
  methylation.data <- methylation.data[, OverlapSamples, drop = FALSE]

  # select the probes for those genes with the gene expression data available
  if(cluster){
    genes <- rownames(methylation.data)
    genes <- sapply(genes, function(x) unlist(strsplit(x, "---"))[1])
    methylation.data <- methylation.data[which(genes %in% rownames(gene.expression.data)), ]
  }else{
    ProbeAnnotation <- ProbeAnnotation[which(ProbeAnnotation$gene %in% rownames(gene.expression.data)),
    ]
    OverlapProbes <- intersect(rownames(methylation.data), ProbeAnnotation$probe)
    methylation.data <- methylation.data[OverlapProbes, , drop = FALSE]
  }

  uniqueProbes <- rownames(methylation.data)

  if (filter) {
    cat("Modeling gene expression...\n")
    PvalueThreshold <- 0.001
    RsquareThreshold <- 0.1

    iterations <- length(uniqueProbes)
    pb <- txtProgressBar(max = iterations, style = 3)

    if (cores == "" | cores == 1) {
      for (i in 1:iterations) {
        if(cluster){
          tmpGenes <- unlist(strsplit(uniqueProbes[i], "---"))[1]
        }else{
          tmpGenes <- ProbeAnnotation$gene[which(ProbeAnnotation$probe == uniqueProbes[i])]
        }
        for (gene in tmpGenes) {
          res <- lm(gene.expression.data[gene, ] ~ methylation.data[uniqueProbes[i],
          ])
          res.summary <- summary(res)
          if (correlation == "negative" & res$coefficients[2] < 0 & res.summary$coefficients[2, 4] <
              PvalueThreshold & res.summary$r.squared > RsquareThreshold) {
            FunctionalProbes <- c(FunctionalProbes, uniqueProbes[i])
            break
          }else if(correlation == "positive" & res$coefficients[2] > 0 & res.summary$coefficients[2, 4] <
                   PvalueThreshold & res.summary$r.squared > RsquareThreshold){
            FunctionalProbes <- c(FunctionalProbes, uniqueProbes[i])
            break
          }
        }
        setTxtProgressBar(pb, i)
      }
    } else {
      # set up a progress bar
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      i <- NULL  # to avoid 'no visible binding for global variable' in R CMD check
      FunctionalProbes <- foreach::foreach(i = 1:length(uniqueProbes), .combine = "c",
                                           .options.snow = opts) %dopar% {
                                             probe <- NULL
                                             if(cluster){
                                               tmpGenes <- unlist(strsplit(uniqueProbes[i], "---"))[1]
                                             }else{
                                               tmpGenes <- ProbeAnnotation$gene[which(ProbeAnnotation$probe == uniqueProbes[i])]
                                             }
                                             for (gene in tmpGenes) {
                                               res <- lm(gene.expression.data[gene, ] ~ methylation.data[uniqueProbes[i],
                                               ])
                                               res.summary <- summary(res)
                                               if (correlation == "negative" & res$coefficients[2] < 0 & res.summary$coefficients[2, 4] <
                                                   PvalueThreshold & res.summary$r.squared > RsquareThreshold) {
                                                 probe <- uniqueProbes[i]
                                                 break
                                               }else if(correlation == "positive" & res$coefficients[2] > 0 & res.summary$coefficients[2, 4] <
                                                        PvalueThreshold & res.summary$r.squared > RsquareThreshold){
                                                 probe <- uniqueProbes[i]
                                                 break
                                               }
                                             }
                                             probe
                                           }
      close(pb)
    }
    FunctionalProbes <- unique(FunctionalProbes)
    cat("\nFound", length(FunctionalProbes), "transcriptionally predictive probes.\n")
  } else {
    # No regression filter, but keep genes present in both methylation and
    # gene expression data, as the ones returned by the regression filter
    FunctionalProbes <- uniqueProbes
  }
  return(FunctionalProbes)
}





