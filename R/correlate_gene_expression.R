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

#' The .getComp function
#' @description Helper function to get a string indicating the comparison made for gene expression
#' @param state character string indicating the methylation state, can be either "Hyper", "Hypo", "Dual"
#' @keywords internal
#'
#' @return a list of sample names split by methylation group
#'

.getComp <- function(state){
  cmp <-  NULL
  if (state == "Hyper"){
    cmp <- "hyper vs normal"
  }else if (state == "Hypo"){
    cmp <- "hypo vs normal"
  }else if (state == "Dual"){
    cmp <- "hypo/normal vs hyper"
  }else{
    stop("Methylation state must be either 'Hyper', 'Hypo' or 'Dual'")
  }
  return(cmp)
}

#' The .getMetGroup function
#' @description Helper function to get sample names split by methylation group based on DM values
#' @param state character string indicating the methylation state, can be either "Hyper", "Hypo", "Dual"
#' @param DM_values a vector of DM values for the probe. The names of the vector are sample names.
#' @keywords internal
#'
#' @return a list of sample names split by methylation group
#'
.getMetGroup <- function(state, DM_values){
  high.met.samples <- low.met.samples <-  NULL
  if (state == "Hyper"){
    high.met.samples <- names(DM_values[DM_values > 0])
    low.met.samples <- names(DM_values[DM_values == 0])
  }else if (state == "Hypo"){
    high.met.samples <- names(DM_values[DM_values == 0])
    low.met.samples <- names(DM_values[DM_values < 0])
  }else if (state == "Dual"){
    high.met.samples <- names(DM_values[DM_values > 0])
    low.met.samples <- names(DM_values[DM_values <= 0])
  }else{
    stop("Methylation state must be either 'Hyper', 'Hypo' or 'Dual'")
  }
  return(list(high.met.samples = high.met.samples, low.met.samples = low.met.samples))
}


#' The test_gene_expr function
#' @description Helper function to test whether the expression levels of a gene is reversely correlated with the methylation state of a probe.
#' @param gene character string indicating a target gene to be modeled.
#' @param probe character string indicating a probe mapped to the target gene.
#' @param DM_values a vector of DM values for the probe. The names of the element should be sample names.
#' @param gene.expr.values a vector of gene expression values for the tested gene. The names of the vector are sample names.
#' @param correlation character indicating the direction of correlation between the methylation state of the CpG site and the gene expression levels. Can be either 'negative' or 'positive'.
#' @param raw.pvalue.threshold raw p value from testing DNA methylation and gene expression
#' @param adjusted.pvalue.threshold adjusted p value from testing DNA methylation and gene expression
#' @keywords internal
#'
#' @return dataframe with functional probe-gene pairs and corresponding p values obtained from the Wilcoxon test for gene expression and methylation.
#'
test_gene_expr <- function(gene, probe, DM_values, gene.expr.values, correlation = "negative"){
  state <- fold_change <- cmp <- p_value <- NULL
  state <- getMethStates_Helper(DM_values)
  cmp <- .getComp(state)

  metGroup <- .getMetGroup(state, DM_values)
  high.met.samples <- metGroup$high.met.samples
  low.met.samples <- metGroup$low.met.samples

  expr.high.values <- as.numeric(gene.expr.values[low.met.samples])  # low methylation, high gene expression
  expr.low.values <- as.numeric(gene.expr.values[high.met.samples])  # high methylation, low gene expression
  if (length(expr.high.values) < 3 || length(expr.low.values)< 3) return(NULL)
  p_value <- wilcox.test(expr.high.values,
                         expr.low.values,
                         alternative = ifelse(correlation == "negative", "greater", "less"), exact = FALSE)$p.value

  if(state == "Hyper"){
    fold_change <- round(mean(expr.low.values, na.rm = TRUE)/mean(expr.high.values, na.rm = TRUE), 3)
  }else{
    fold_change <- round(mean(expr.high.values, na.rm = TRUE)/mean(expr.low.values, na.rm = TRUE), 3)
  }

  # produce a dataframe for gene expression with methylation state and
  # prevalence information
  dataDEGs <- data.frame(Gene = gene, Probe = probe)
  if(nrow(dataDEGs) == 0) return(NULL)
  dataDEGs["State"] <- state
  dataDEGs["Fold change of gene expression"] <- fold_change
  dataDEGs["Comparators"] <- cmp
  dataDEGs["Raw.p"] <- p_value
  return(dataDEGs)
}

#' The generateFunctionalPairs function
#' @description Wrapper function to get functional CpG-gene pairs, used for Regular, miRNA and lncRNA modes
#' @param MET_matrix matrix of methylation states
#' @param control.names character vector indicating the samples names in the control group
#' @param gene.expression.data matrix of gene expression data
#' @param ProbeAnnotation dataframe of probe annotation
#' @param raw.pvalue.threshold raw p value threshold
#' @param adjusted.pvalue.threshold  adjusted p value threshold
#' @param cores number of computational cores
#' @param mode character string indicating the analytic mode
#' @param correlation the expected relationship between DNAme and gene expression
#'
#' @return a dataframe of functional CpG-gene matrix

generateFunctionalPairs <- function(MET_matrix, control.names, gene.expression.data,
                                    ProbeAnnotation, raw.pvalue.threshold, adjusted.pvalue.threshold,
                                    cores, mode = "Regular", correlation = "negative") {

  MET_matrix <- filterMethMatrix(MET_matrix = MET_matrix, control.names = control.names,
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
      gene.expr.values <- gene.expression.data[gene, ]
      probes <- ProbeAnnotation$probe[which(ProbeAnnotation$gene == gene)]
      for(probe in probes){
        DM_values <- MET_matrix[probe, ]
        pairs <- test_gene_expr(gene, probe, DM_values, gene.expr.values, correlation = correlation)
        FunctionalPairs <- rbind(FunctionalPairs, pairs)
      }
      utils::setTxtProgressBar(pb, i)
    }
  } else {
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    FunctionalPairs <- foreach::foreach(i = seq_len(iterations), .combine = rbind, .options.snow = opts,
                                        .verbose = FALSE) %dopar% {
                                          gene <- uniqueGenes[i]
                                          probes <- ProbeAnnotation$probe[which(ProbeAnnotation$gene == gene)]
                                          for(probe in probes){
                                            DM_values <- MET_matrix[probe, ]
                                            gene.expr.values <- gene.expression.data[gene, ]
                                            pairs <- test_gene_expr(gene, probe, DM_values, gene.expr.values, correlation = correlation)
                                            FunctionalPairs <- rbind(FunctionalPairs, pairs)
                                          }
                                        }
  }
  close(pb)
  FunctionalPairs["Adjusted.p"] <- p.adjust(as.vector(unlist(FunctionalPairs["Raw.p"])), method = "fdr")
  FunctionalPairs <- FunctionalPairs[which(FunctionalPairs$Raw.p < raw.pvalue.threshold & FunctionalPairs$Adjusted.p <
                               adjusted.pvalue.threshold), ]
  return(FunctionalPairs)
}


#' The getFunctionalGenes function
#' @description Helper function to assess if the methylation of a probe is reversely correlated with the expression of its nearby genes.
#' @details This function is probe-centered, which is used in the enhancer mode and the miRNA mode of EpiMix.
#' @param target.probe character string indicating the probe to be evaluated.
#' @param target.genes character vector indicating the nearby genes of the target probe.
#' @param MET_matrix methylation data matrix for CpGs from group.1 and group.2.
#' @param gene.expression.data gene expression data matrix.
#' @param ProbeAnnotation GRange object of CpG probe annotation.
#' @param raw.pvalue.threshold raw p value from testing DNA methylation and gene expression
#' @param adjusted.pvalue.threshold adjusted p value from testing DNA methylation and gene expression
#' @keywords internal
#' @return dataframe with functional probe-gene pair and p values from the Wilcoxon test for methylation and gene expression.
#'
#' @examples
#' \donttest{
#' data(Sample_EpiMixResults_Enhancer)
#' data(mRNA.data)
#' EpiMixResults <- Sample_EpiMixResults_Enhancer
#' target.probe <- EpiMixResults$FunctionalPairs$Probe[1]
#' target.genes <- EpiMixResults$FunctionalPairs$Gene
#' MET_matrix <- EpiMixResults$MethylationStates
#' ProbeAnnotation <- ExperimentHub::ExperimentHub()[["EH3675"]]
#' res <- getFunctionalGenes(target.probe, target.genes, MET_matrix, mRNA.data, ProbeAnnotation)
#' }
#'
getFunctionalGenes <- function(target.probe, target.genes, MET_matrix, gene.expression.data,
                               ProbeAnnotation, correlation = "negative", raw.pvalue.threshold = 0.05, adjusted.pvalue.threshold = 0.01) {


  DM_values <- MET_matrix[target.probe, ]
  # get fold changes and raw p values
  raw.pvals <- data.frame()
  for(gene in target.genes){
    gene.expr.values <- gene.expression.data[gene, ]
    pairs <- test_gene_expr(gene, target.probe, DM_values, gene.expr.values, correlation = correlation)
    raw.pvals  <- rbind(raw.pvals, pairs)
  }

  if (is.null(raw.pvals))
    return(NULL)

  # get permutation p values
  perm.pvals <- data.frame()
  random.genes <- getRandomGenes(target.probe = target.probe, gene.expression.data = gene.expression.data,
                                 ProbeAnnotation, genome = "hg38", perm = 1000)

  for(gene in random.genes){
    gene.expr.values <- gene.expression.data[gene, ]
    pairs <- test_gene_expr(gene, target.probe, DM_values, gene.expr.values, correlation = correlation)
    perm.pvals  <- rbind(perm.pvals, pairs)
  }

  dataDEGs <- Get.Pvalue.p(raw.pvals, perm.pvals)
  dataDEGs <- dataDEGs[which(dataDEGs$Raw.p < raw.pvalue.threshold & dataDEGs$Adjusted.p <
                               adjusted.pvalue.threshold), ]
  dataDEGs <- dplyr::distinct(dataDEGs)
  return(dataDEGs)
}

#' Calculate empirical Pvalue
#' @param U.matrix A data.frame of raw pvalue from U test. Output from .Stat.nonpara
#' @param permu data frame of permutation. Output from .Stat.nonpara.permu
#' @return A data frame with empirical Pvalue.
Get.Pvalue.p <- function(U.matrix, permu) {
  .Pvalue <- function(x, permu) {
    Gene <- as.character(x["Gene"])
    Raw.p <- as.numeric(x["Raw.p"])
    if (is.na(Raw.p)) {
      out <- NA
    } else {
      # num( Pp <= Pr) + 1 Pe = --------------------- x + 1 Pp = pvalue
      # probe (Raw.p) Pr = pvalue random probe (permu matrix) We have to
      # consider that floating Point Numbers are Inaccurate
      out <- (sum(permu$Raw.p - Raw.p < 10^-100, na.rm = TRUE) + 1)/(sum(!is.na(permu$Raw.p)) + 1)

    }
    return(out)
  }
  # message('Calculating empirical P value.\n')
  Pvalue <- unlist(apply(U.matrix, 1, .Pvalue, permu = permu))
  U.matrix$Adjusted.p <- Pvalue
  return(U.matrix)
}

#' The filterLinearProbes function
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
filterLinearProbes <- function(methylation.data, gene.expression.data, ProbeAnnotation,
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





