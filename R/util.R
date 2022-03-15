###################################################################################################################
#                                   Utility functions used by EpiMix
###################################################################################################################

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

utils::globalVariables(c("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene", "i", "."))

#' The EpiMix_getInfiniumAnnotation function
#' @description fetch the Infinium probe annotation from the seasameData library
#' @param plat character string indicating the methylation platform
#' @param genome character string indicating the version of genome build
#'
#' @return a GRange object of probe annotation
#' @keywords internal

EpiMix_getInfiniumAnnotation <- function(plat = "EPIC", genome = "hg38"){
  ProbeAnnotation = sesameData :: sesameDataGetAnno(paste0(plat, "/", plat, "." ,genome, ".manifest.rds"))
  # if(tolower(genome) == "hg19" & toupper(plat) == "HM27" ) ProbeAnnotation = sesameData::sesameDataGet("HM27.hg19.manifest")
  # if(tolower(genome) == "hg19" & toupper(plat) == "HM450" ) ProbeAnnotation = sesameData::sesameDataGet("HM450.hg19.manifest")
  # if(tolower(genome) == "hg19" & toupper(plat) == "EPIC" ) ProbeAnnotation = sesameData::sesameDataGet("EPIC.hg19.manifest")
  # if(tolower(genome) == "hg38" & toupper(plat) == "HM27" ) ProbeAnnotation = sesameData::sesameDataGet("HM27.hg38.manifest")
  # if(tolower(genome) == "hg38" & toupper(plat) == "HM450" ) ProbeAnnotation = sesameData::sesameDataGet("HM450.hg38.manifest")
  # if(tolower(genome) == "hg38" & toupper(plat) == "EPIC" ) ProbeAnnotation = sesameData::sesameDataGet("EPIC.hg38.manifest")
  return(ProbeAnnotation)
}

#' The convertAnnotToDF function
#' @description convert the probe annotation from the GRange object to a dataframe
#' @param annot a GRange object of probe annotation, can be the object returned from the getInfiniumAnnotation function.
#' @return a dataframe with chromosome, beginning and end position, mapped gene information for each CpG probe
#' @keywords internal
#'
convertAnnotToDF <- function(annot){
  df.annot = data.frame(
                        "CpG_chrm" = GenomicRanges :: seqnames(annot),
                        "CpG_beg" = GenomicRanges :: start(ranges(annot)),
                        "CpG_end" = GenomicRanges :: end(ranges(annot)),
                        "probeID" = names(annot),
                        "gene" = GenomicRanges :: mcols(annot)$gene
                        )
  return(df.annot)
}

#' The split.met.data function
#' @description Helper function to split the methylation data matrix into the experimental group and the control group
#' @param methylation.data methylation data matrix
#' @param sample.info sample information matrix
#' @param group.1 name of group.1
#' @param group.2 name of group.2
#' @keywords internal
#' @return a list with methylation data of group.1 and group.2

split.met.data <- function(methylation.data, sample.info, group.1, group.2){
  MET_Experiment <-  MET_Control <- NULL
  if(!is.null(sample.info) & !is.null(group.1) & !is.null(group.2)){
    Samples_Experiment <- intersect(colnames(methylation.data), sample.info$primary[which(sample.info$sample.type %in% group.1)])
    Samples_Control <- intersect(colnames(methylation.data), sample.info$primary[which(sample.info$sample.type %in% group.2)])
    cat("Found", length(Samples_Experiment), "samples in group.1 and", length(Samples_Control), "samples in group.2\n")
    if(length(Samples_Experiment) == 0) stop("Cannot find methylation data for samples in the group.1! Please check sample information.")
    if(length(Samples_Control) == 0) stop("Cannot find methylation data for samples in the group.2! Please check sample information.")
    MET_Experiment <- methylation.data[,Samples_Experiment,drop = FALSE]
    MET_Control <- methylation.data[,Samples_Control,drop = FALSE]
  }else{
    MET_Experiment = methylation.data
  }
  return(list(MET_Experiment = MET_Experiment, MET_Control = MET_Control))
}

#' The mapProbeGene function
#' @description since in the original probe annotation, a specific probe can be mapped to multiple genes, this function splits the rows and maps each probe to a signle gene in a row.
#' @param df.annot a dataframe with probe annotation, can be the object returned from the convertAnnotToDF function.
#' @return a dataframe with 1:1 mapping of probe and gene
#' @keywords internal
#'
mapProbeGene <- function(df.annot){
  df.annot <- df.annot[!is.na(df.annot$probeID) & !is.na(df.annot$gene), ]
  df.annot <- tidyr :: separate_rows(df.annot, .data$gene, sep = ";")
  df.annot <- dplyr :: distinct(df.annot)
  return(df.annot)
}


#' The getMethStates function
#' @description Helper function that adds a methyaltion state label to each driver probe
#' @param MethylMixResults the list object returned from the EpiMix function
#' @param DM.probes character vector of differentially methylated probes.
#' @keywords internal
#' @return a character vector with the methylation state ("Hypo", "Hyper" or "Dual") for each probe. The names for the vector are the probe names and the values are the methylation state.
#'
getMethStates <- function(MethylMixResults, DM.probes){
  mixture.states = MethylMixResults$MixtureStates
  mixture.states = mixture.states[DM.probes]
  meth.states = c()
  for(probe in DM.probes){
    state = NULL
    DMValues = mixture.states[[probe]]
    # after applying filter 1 in the main code, we only have probes with two states now (normal-hyper, normal-hypo or hypo-hyper)
    if(max(DMValues) > 0 & min(DMValues) < 0){
      state = "Dual"
    }else if(max(DMValues) > 0){
      state = "Hyper"
    }else{
      state = "Hypo"
    }
    meth.states = append(meth.states, state)
  }
  names(meth.states) = DM.probes
  return(meth.states)
}

#' The getFunctionalProbes function
#' @description Helper function to assess if the expression of a specific gene is reversely correlated with the methylation of a list of probes mapped to it.
#' @details This function is gene-centered, which is used in the regular mode and the lncRNA mode of EpiMix.
#' @param gene character string indicating the target gene to be modeled.
#' @param probes character vector indicating the probes mapped to the target gene.
#' @param MET_matrix methylation data matrix for CpGs from group.1 and group.2.
#' @param MET_Control methylation data matrix for CpGs from group.2.
#' @param exp gene expression data matrix.
#' @param methylation.states character vector indicating the methylation states of the target probes
#' @param raw.pvalue.threshold raw p value from testing DNA methylation and gene expression
#' @param adjusted.pvalue.threshold adjusted p value from testing DNA methylation and gene expression
#' @keywords internal
#'
#' @return dataframe with functional probe-gene pairs and corresponding p values obtained from the Wilcoxon test for gene expression and methylation.
#'
getFunctionalProbes <- function(gene,
                                probes,
                                MET_matrix,
                                MET_Control,
                                exp,
                                methylation.states,
                                raw.pvalue.threshold = 0.05,
                                adjusted.pvalue.threshold = 0.01){
  sample.number = ncol(MET_matrix)
  if(!is.null(MET_Control)){
    sample.number = ncol(MET_matrix) - length(intersect(colnames(MET_matrix), colnames(MET_Control)))
  }
  p = c()
  mRNA.fold.change = c()
  comparisons = c()
  hypo_prev = c()
  hyper_prev = c()
  states = c()
  for(probe in probes){
    prev.hypo <- prev.hyper <-  state <- NULL
    state = methylation.states[probe]
    DM_values = MET_matrix[probe, ] # a single-row vector with the names as sample names, and the values as DM values
    high.met.samples = names(DM_values[DM_values == max(DM_values)])
    low.met.samples = names(DM_values[DM_values == min(DM_values)])
    expr.low.values = as.numeric(exp[gene,high.met.samples])  # high methylation, low gene expression
    expr.high.values = as.numeric(exp[gene,low.met.samples])  # low methylation, high gene expression
    p_value = wilcox.test(expr.high.values, expr.low.values, alternative = "greater")$p.value
    if(state == "Hyper") {
      prev.hypo = 0
      prev.hyper = round(length(high.met.samples) / sample.number * 100, 3)
      fold_change = round(mean(expr.low.values) / mean(expr.high.values),3)
      cmp = "hyper vs normal"
    }else if(state == "Dual") {
      prev.hypo = round(length(low.met.samples) / sample.number * 100, 3)
      prev.hyper = round(length(high.met.samples) / sample.number * 100, 3)
      fold_change = round(mean(expr.high.values) / mean(expr.low.values), 3)
      cmp = "hypo vs hyper"
    }else{
      prev.hypo = round(length(low.met.samples) / sample.number * 100, 3)
      prev.hyper = 0
      fold_change = round(mean(expr.high.values) / mean(expr.low.values),3)
      cmp = "hypo vs normal"
    }
    states = append(states, state)
    p = append(p, p_value)
    hypo_prev= append(hypo_prev, prev.hypo)
    hyper_prev= append(hyper_prev, prev.hyper)
    comparisons = append(comparisons, cmp)
    mRNA.fold.change = append(mRNA.fold.change, fold_change)
  }

  # produce a dataframe for gene expression with methylation state and prevalence inforamtion
  dataDEGs <- data.frame(Gene = rep(gene,length(probes)),
                         Probe = probes,
                         State = states
  )
  dataDEGs['Proportion of hypo (%)'] <- hypo_prev
  dataDEGs['Proportion of hyper (%)'] <- hyper_prev
  dataDEGs['Fold change of gene expression'] <- mRNA.fold.change
  dataDEGs['Comparators'] <- comparisons
  dataDEGs['Raw.p'] <- p
  dataDEGs["Adjusted.p"] <- p.adjust(as.vector(unlist(dataDEGs['Raw.p'])), method='fdr')
  dataDEGs <- dataDEGs[which(dataDEGs$Raw.p<raw.pvalue.threshold & dataDEGs$Adjusted.p < adjusted.pvalue.threshold),]
  return(dataDEGs)
}

#' The splitmatix function
#'
#' @param x A matrix
#' @param by A character specify if split the matrix by row or column.
#' @keywords internal
#'
#' @return A list each of which is the value of each row/column in the matrix.
#'
splitmatrix <- function(x,by="row") {
  if(by %in% "row"){
    out <- split(x, rownames(x))
  }else if (by %in% "col"){
    out <- split(x, colnames(x))
  }
  return(out)
}

#' The get.prevalence function
#' @description Helper function to get the prevalence of the differential methylation of a CpG site in the study population.
#' @param target.probe target CpG site.
#' @param MET_matrix combined methylation state matrix for both case and control groups.
#' @param MET_control methylation state matrix only for the control group.
#' @param state methylation state for the target CpG, can be "Hyper", "Hypo" or "Dual".
#' @keywords internal
#'
#' @return a list of prevalence for the abnormal methylation
#'
get.prevalence <- function(target.probe, MET_matrix, MET_control, state){
  if(!is.null(MET_control)){
    control.samples <- colnames(MET_control)
    exp.samples <- setdiff(colnames(MET_matrix), control.samples)
  }else{
    exp.samples <- colnames(MET_matrix)
  }
  MET_exp <- MET_matrix[,exp.samples, drop = F]
  DM_values = MET_exp[target.probe, ]
  high.met.samples = names(DM_values[DM_values == max(DM_values)])
  low.met.samples = names(DM_values[DM_values == min(DM_values)])
  sample.number = ncol(MET_exp)
  if(state == "Hyper"){
    hypo_prev = 0
    hyper_prev = round(length(high.met.samples) / sample.number * 100, 3)
  }else if(state == "Hypo"){
    hypo_prev =  round(length(low.met.samples) / sample.number * 100, 3)
    hyper_prev = 0
  }else{
    hypo_prev = round(length(low.met.samples) / sample.number * 100, 3)
    hyper_prev = round(length(high.met.samples) / sample.number * 100, 3)
  }
  return(list(hypo = hypo_prev, hyper = hyper_prev))
}


#' The model.gene.expression function
#'
#' @param target.probe character string indicating the target CpG site.
#' @param MET_matrix combined methylation state matrix for both case and control groups.
#' @param gene.expression.data a matrix with gene expression data.
#' @param target.genes character vector with the names of the genes associated with the target CpGs.
#' @param state methylation state for the target CpG, can be "Hyper", "Hypo" or "Dual".
#' @param correlation the direction of correlation between the methylation state of the CpG site and the gene expression levels. Can be either "negative" or "positive".
#' @keywords internal
#'
#' @return a dataframe with fold change of gene expression and p values.
#'
model.gene.expression <- function(target.probe, MET_matrix, gene.expression.data, target.genes, state, correlation = "negative"){
  DM_values = MET_matrix[target.probe, ]
  high.met.samples = which(DM_values == max(DM_values))
  low.met.samples = which(DM_values == min(DM_values))
  Exps = gene.expression.data[target.genes, ,drop = F]
  if(length(Exps) == 0) return(NULL)
  fold.change <- unlist(lapply(splitmatrix(Exps),
                          function(x) {
                            if(state == "Hyper"){
                              round(mean(x[high.met.samples])/mean(x[low.met.samples]), 3)
                            }else if(state == "Dual"){
                              round(mean(x[low.met.samples])/mean(x[high.met.samples]), 3)
                            }else{
                              round(mean(x[low.met.samples])/mean(x[high.met.samples]), 3)
                            }
                          }
                       ))
  test.p <- unlist(lapply(splitmatrix(Exps),
                          function(x) {
                            wilcox.test(x[low.met.samples],
                                        x[high.met.samples],
                                        alternative = ifelse(correlation == "negative","greater","less"),
                                        exact = FALSE)$p.value
                          }
                   ))
  res <- data.frame(Gene = target.genes,
                    fold_change = fold.change[match(target.genes, names(fold.change))],
                    Raw.p = test.p[match(target.genes, names(test.p))],
                    stringsAsFactors = FALSE)
  return(res)
}


#' The get.chromosome function
#' @description given a list of genes, get the chromosomes of these genes.
#' @param genes character vector with the gene names
#' @param genome character string indicating the genome build version, can be either "hg19" or "hg38"
#' @keywords internal
#'
#' @return a dataframe for the mapping between genes and their chromosomes.
#'
get.chromosome <- function(genes, genome){
  gene.annotation = as.data.frame(getTSS(genome))
  gene.annotation = distinct(gene.annotation[,c("external_gene_name", "seqnames")])
  gene.annotation = gene.annotation[which(gene.annotation$external_gene_name %in% genes),]
  colnames(gene.annotation) = c("Gene", "Chr")
  return(gene.annotation)
}


#' The get.perm.pVals function
#' @description Helper function to perform permutation test for a target enhancer CpG. Used in the Enhancer mode.
#' @param target.probe character string indicating the target CpG for generating the permutation p values.
#' @param MET_matrix a matrix of methylaiton states for the patients.
#' @param gene.expression.data a matrix of gene expression data.
#' @param ProbeAnnotation GRange object of probe annotation.
#' @param genome character string indicating the genome build version, can be either "hg19" or "hg38".
#' @param correlation the direction of correlation between the methylation state of the CpG site and the gene expression levels. Can be either "negative" or "positive".
#' @param perm the number of permutation tests.
#' @keywords internal
#'
#' @return a dataframe for the permutation genes and p values for the target CpG site.
#'
get.perm.pVals <- function(target.probe,
                           MET_matrix,
                           gene.expression.data,
                           ProbeAnnotation,
                           genome = "hg38",
                           correlation = "negative",
                           perm = 1000){

  target.chr = as.character(seqnames(ProbeAnnotation)[which(names(ProbeAnnotation) == target.probe)])
  all.genes = rownames(gene.expression.data)
  gene.chr = get.chromosome(all.genes, genome)
  random.genes = gene.chr$Gene[which(!gene.chr$Chr %in% c(target.chr, "chrX", "chrY"))]
  if(length(unique(random.genes)) < perm){
    warning("There is no gene expression data for enough genes to generate ", perm, " permutations, using ", length(unique(random.genes)), " random genes instead.", immediate. = TRUE)
    perm = length(unique(random.genes))
  }
  random.genes = sample(gene.chr$Gene[which(!gene.chr$Chr %in% c(target.chr, "chrX", "chrY"))], perm)
  Exps = gene.expression.data[random.genes, ]
  DM_values = MET_matrix[target.probe, ] # a single-row vector with the names as sample names, and the values as DM values
  high.met.samples = which(DM_values == max(DM_values))
  low.met.samples = which(DM_values == min(DM_values))
  test.p <- unlist(lapply(splitmatrix(Exps),
                          function(x) {
                            wilcox.test(x[low.met.samples],
                                        x[high.met.samples],
                                        alternative = ifelse(correlation == "negative","greater","less"),
                                        exact = FALSE)$p.value
                          }

  ))
  test.p <- data.frame(Genes=random.genes,
                       Perm.p=test.p[match(random.genes, names(test.p))],
                       stringsAsFactors = FALSE)
}

#' The get.comparators function
#' @description Helper function to generate the string indicating the comparisons of gene expression. Used for generating the functional CpG-gene pair matrix.
#' @param state the methylation state of a CpG site. Can be either "Hyper", "Hypo" or "Dual".
#' @keywords internal
#' @return a string indicating the comparisions for gene expression.

get.comparators <- function(state){
  if(state == "Hyper") {
    cmp = "hyper vs normal"
  }else if(state == "Dual") {
    cmp = "hypo vs hyper"
  }else{
    cmp = "hypo vs normal"
  }
  return(cmp)
}

#' The getFunctionalGenes function
#' @description  Helper function to assess if the methylation of a probe is reversely correlated with the expression of its nearby genes.
#' @details This function is probe-centered, which is used in the enhancer mode and the miRNA mode of EpiMix.
#' @param target.probe character string indicating the probe to be evaluated.
#' @param state character string indicating the methylation state of the probe, should be either "Hyper", "Hypo" or "Dual".
#' @param target.genes character vector indicating the nearby genes of the target probe.
#' @param MET_matrix methylation data matrix for CpGs from group.1 and group.2.
#' @param MET_Control methylation data matrix for CpGs from group.2.
#' @param gene.expression.data gene expression data matrix.
#' @param ProbeAnnotation GRange object of CpG probe annotation.
#' @param raw.pvalue.threshold raw p value from testing DNA methylation and gene expression
#' @param adjusted.pvalue.threshold adjusted p value from testing DNA methylation and gene expression
#' @keywords internal
#' @return dataframe with functional probe-gene pair and p values from the Wilcoxon test for methylation and gene expression.
#'
getFunctionalGenes <- function(target.probe,
                               state,
                               target.genes,
                               MET_matrix,
                               MET_Control,
                               gene.expression.data,
                               ProbeAnnotation,
                               raw.pvalue.threshold = 0.05,
                               adjusted.pvalue.threshold = 0.01){
  target.genes = intersect(target.genes, rownames(gene.expression.data))
  if(length(target.genes) == 0) return(NULL)
  # calculate prevalence
  prev = get.prevalence(target.probe, MET_matrix, MET_Control, state)
  hypo_prev = prev$hypo
  hyper_prev = prev$hyper
  # modeling gene expression
  gene.expr <- model.gene.expression(target.probe, MET_matrix, gene.expression.data, target.genes, state)
  if(is.null(gene.expr)) return(NULL)
  # find the permutation p values of this probe
  perm.pvals = get.perm.pVals(target.probe,
                              MET_matrix,
                              gene.expression.data,
                              ProbeAnnotation,
                              genome = "hg38",
                              correlation = "negative")
  cmp = get.comparators(state)
  dataDEGs <- data.frame(Gene = gene.expr$Gene,
                         Probe = rep(target.probe,length(gene.expr$Gene)),
                         State = rep(state,length(gene.expr$Gene)),
                         Prev_hypo = rep(hypo_prev, length(gene.expr$Gene)),
                         Prev_hypo =  rep(hyper_prev, length(gene.expr$Gene)),
                         fold_change = gene.expr$fold_change,
                         Comparators = rep(cmp, length(gene.expr$Gene)),
                         Raw.p = gene.expr$Raw.p)
  dataDEGs = Get.Pvalue.p(dataDEGs, perm.pvals)
  dataDEGs = dataDEGs[which(dataDEGs$Raw.p<raw.pvalue.threshold & dataDEGs$Adjusted.p < adjusted.pvalue.threshold),]
  colnames(dataDEGs) <- c("Gene", "Probe", "State", "Proportion of hypo (%)", "Proportion of hyper (%)", "Fold change of gene expression", "Comparators", "Raw.p", "Adjusted.p")
  dataDEGs = dplyr :: distinct(dataDEGs)
  return(dataDEGs)
}

#' Calculate empirical Pvalue
#' @param U.matrix A data.frame of raw pvalue from U test. Output from .Stat.nonpara
#' @param permu data frame of permutation. Output from .Stat.nonpara.permu
#' @return A data frame with empirical Pvalue.
Get.Pvalue.p <- function(U.matrix, permu){
  .Pvalue <- function(x,permu){
    Gene <- as.character(x["Gene"])
    Raw.p <- as.numeric(x["Raw.p"])
    if(is.na(Raw.p)){
      out <- NA
    } else {
      #       num( Pp <= Pr) + 1
      # Pe = ---------------------
      #            x + 1
      # Pp = pvalue probe (Raw.p)
      # Pr = pvalue random probe (permu matrix)
      # We have to consider that floating Point Numbers are Inaccurate
      out <- (sum(permu$Perm.p - Raw.p < 10^-100, na.rm=TRUE) + 1) / (sum(!is.na(permu$Perm.p)) + 1)
    }
    return(out)
  }
  #message("Calculating empirical P value.\n")
  Pvalue <- unlist(apply(U.matrix,1,.Pvalue,permu=permu))
  U.matrix$Adjusted.p <- Pvalue
  return(U.matrix)
}

#' The filterMethMatrix function
#' @details This function filters methylation states from the beta mixture modeling for each probe. The filtered probes can be used to model gene expression by Wilcoxon test.
#' @param MET_matrix a matrix of methylation states from the EpiMix results
#' @param gene.expression.data a matrix with gene expression data
#' @keywords internal
#' @return a matrix of methylation states for each differentially methylated probe with probes in rows and patient in columns.
#'
filterMethMatrix <- function(MET_matrix, gene.expression.data){
  # Filter 1: excluding probes that present only 1 methylation state, in which case we can not compare the gene expression for one methylation component versus another
  METcounts <- apply(MET_matrix, 1, function(x)length(unique(x)))
  MET_matrix <- MET_matrix[METcounts!=1,,drop = FALSE]

  # Filter 2: filter on samples which have gene expression data available
  if(nrow(MET_matrix) > 0){
    common.samps <- intersect(colnames(MET_matrix),colnames(gene.expression.data))
    MET_matrix <- MET_matrix[,common.samps]
  }
  # Filter 3: filtering out probes with the minority group having sample number < 3. Otherwise, the Wilcoxon test won't work.
  if(ncol(MET_matrix) > 0){
    METminCounts <- apply(MET_matrix, 1, function(x)min(as.numeric(table(x))))
    MET_matrix <- MET_matrix[METminCounts>=3,]
  }
  return(MET_matrix)
}


#' The EpiMix_ModelGeneExpression function
#' @description pre-select funcitonal probes.
#' @param methylation.data methylation data matrix.
#' @param gene.expression.data gene expression data matrix.
#' @param ProbeAnnotation dataframe of probe annotation
#' @param cores number of CPU cores used for computation
#' @param filter logical indicating whether to perform a linear regression to select functional probes
#' @return a character vector of probe names.
#' @import doSNOW
#' @keywords internal
#'
#'
EpiMix_ModelGeneExpression <- function(methylation.data, gene.expression.data, ProbeAnnotation, cores, filter){
  FunctionalProbes = character(0)

  # overlapping samples to select only samples with both methylation data and gene expression data
  OverlapSamples = intersect(colnames(methylation.data), colnames(gene.expression.data))
  cat("Found", length(OverlapSamples), "samples with both methylation and gene expression data.\n")
  gene.expression.data = gene.expression.data[, OverlapSamples, drop = FALSE]
  methylation.data = methylation.data[, OverlapSamples, drop = FALSE]

  # select the probes for those genes with the gene expression data available
  ProbeAnnotation = ProbeAnnotation[which(ProbeAnnotation$gene %in% rownames(gene.expression.data)), ]
  OverlapProbes = intersect(rownames(methylation.data), ProbeAnnotation$probeID)
  methylation.data = methylation.data[OverlapProbes, , drop = FALSE]

  uniqueProbes = rownames(methylation.data)

  if(filter){
    cat("Modeling gene expression...\n")
    PvalueThreshold = 0.001
    RsquareThreshold = 0.1

    iterations = length(uniqueProbes)
    pb <- txtProgressBar(max = iterations, style = 3)

    if(cores == "" | cores == 1){
      for (i in 1:iterations){
        tmpGenes = ProbeAnnotation$gene[which(ProbeAnnotation$probeID == uniqueProbes[i])]
        for(gene in tmpGenes){
          res = lm(gene.expression.data[gene, ] ~ methylation.data[uniqueProbes[i], ])
          res.summary = summary(res)
          if (res$coefficients[2] < 0 & res.summary$coefficients[2, 4] < PvalueThreshold & res.summary$r.squared  > RsquareThreshold) {
            FunctionalProbes = c(FunctionalProbes, uniqueProbes[i])
            break
          }
        }
        setTxtProgressBar(pb,i)
      }
    } else{
      # set up a progress bar
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      i <- NULL # to avoid "no visible binding for global variable" in R CMD check
      FunctionalProbes = foreach::foreach(i = 1 : length(uniqueProbes), .combine = 'c', .options.snow= opts) %dopar% {
        probe = NULL
        tmpGenes = ProbeAnnotation$gene[which(ProbeAnnotation$probeID == uniqueProbes[i])]
        for(gene in tmpGenes){
          res = lm(gene.expression.data[gene, ] ~ methylation.data[uniqueProbes[i], ])
          res.summary = summary(res)
          if (res$coefficients[2] < 0 & res.summary$coefficients[2, 4] < PvalueThreshold & res.summary$r.squared  > RsquareThreshold) {
            probe = uniqueProbes[i]
            break
          }
        }
        probe
      }
      close(pb)
    }
    FunctionalProbes = unique(FunctionalProbes)
    cat("\nFound", length(FunctionalProbes), "transcriptionally predictive probes.\n")
  } else{
    # No regression filter, but keep genes present in both methylation and gene expression data, as the ones returned by the regression filter
    FunctionalProbes =  uniqueProbes
  }
  return(FunctionalProbes)
}

#' getRoadMapEnhancerProbes
#' @details get the CpG probes that locate at the enhancer regions identified by the Roadmap epigenomics project
#' @param met.platform character string indicating the methylation platform, can be either "EPIC" or "HM450"
#' @param genome character string indicating the genome build version, can be either "hg19" or "hg38"
#' @param functional.regions character vector indicating the MNEMONIC chromatin states that will be retrieved from the Roadmap epigenomics. Default values are the active enhancers:"EnhA1", "EnhA2".
#' @param listOfEpigenomes character vector indicting which epigenome(s) to use for finding enhancers.
#' @return a dataframe with enhancer probes and their chromosome coordinates
#' @keywords internal
#' @examples
#'\dontrun{
#' met.platform = "EPIC"
#' genome = "hg38"
#' listOfEpigenomes = c("E034", "E045", "E047")
#' functional.regions = c("EnhA1", "EnhA2", "EnhG1", "EnhG2")
#' df.enhancer.probes <-  getEnhancerProbes(met.platform = met.platform,
#'                                         genome = genome,
#'                                         functional.regions = functional.regions,
#'                                         listOfEpigenomes = listOfEpigenomes)
#'
#' }

getRoadMapEnhancerProbes <- function(met.platform = "EPIC", genome = "hg38", functional.regions=c("EnhA1", "EnhA2"), listOfEpigenomes = NULL){

  # Step 1. Find all the filenames ending with "_18_core_K27ac_hg38lift_mnemonics.bed.gz" from the Roadmap Epigenomics web portal
  K27Ac_url <- "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/"
  urlData = RCurl :: getURL(K27Ac_url)
  urlData2 = unlist(strsplit(urlData,"\\n"))
  filenames = as.matrix(urlData2[grep("_18_core_K27ac_hg38lift_mnemonics.bed.gz",urlData2)])
  filenames = unlist(strsplit(filenames, ">|<"))
  filenames = filenames[grep("_18_core_K27ac_hg38lift_mnemonics.bed.gz",filenames)]
  filenames = filenames[-grep("a href",filenames)]
  if(!is.null(listOfEpigenomes)){
    listOfEpigenomes = toupper(listOfEpigenomes)
    filenames = filenames[grep(paste(listOfEpigenomes,collapse="|"), filenames)]
  }

  # Step 2.Set the destination file directory and start to downloading files
  dir = file.path(getwd(), "ReferenceEpigenomes")
  dir.create(dir, showWarnings = FALSE)

  cat("Downloading chromatin states from the Roadmap Epigenomics...\n")
  for(file in filenames){
    file.source = paste0(K27Ac_url,file)
    destfile = paste0(dir, "/", file)
    if(!file.exists(destfile)) {
      if(Sys.info()["sysname"] == "Windows") mode <- "wb" else  mode <- "w"
      downloader::download(file.source, destfile = destfile, mode = mode)
    }
  }

  # Step 3. Get probe annotation as a GRange object
  cat("Getting probe annotation...\n")
  ProbeAnnotation <- EpiMix_getInfiniumAnnotation(plat = met.platform, genome = genome)

  # Step 4. Get overlaps between probes and enhancers
  enhancerProbes = character(0)
  for(file in filenames){
    destfile = paste0(dir, "/", file)
    genome.state <- as.data.frame(read.table(file = destfile,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
    colnames(genome.state)=c("chr","start","end", "state")
    genome.state$state = sapply(genome.state$state, function(x) unlist(strsplit(x, "_"))[2])
    genome.state = genome.state[genome.state$state %in% functional.regions, ]
    gr.enhancer = GenomicRanges :: makeGRangesFromDataFrame(genome.state, ignore.strand = TRUE, keep.extra.columns = TRUE)
    target.probes = ProbeAnnotation[unique(S4Vectors:: queryHits(IRanges :: findOverlaps(ProbeAnnotation,gr.enhancer, ignore.strand = TRUE)))]
    this.epigenome = unlist(strsplit(file, "_"))[1]
    cat("\tIdentifed", length(target.probes), "enhancer CpGs from the epigenome", this.epigenome, "\n")
    enhancerProbes = c(enhancerProbes, names(target.probes))
    enhancerProbes = unique(enhancerProbes)
  }
  # Step 5. Generate a dataframe for enhancer probes with their coordinates
  probe.ID = names(ProbeAnnotation)
  probe.chr = GenomicRanges :: seqnames(ProbeAnnotation)
  probe.start.pos = GenomicRanges :: start(ranges(ProbeAnnotation))
  probe.start.end = GenomicRanges :: end(ranges(ProbeAnnotation))
  df.ProbeAnnotation <- data.frame(probeID = probe.ID, chr=probe.chr,start =probe.start.pos, end = probe.start.end)
  df.enhancer.probes <-  df.ProbeAnnotation[df.ProbeAnnotation$probeID %in% enhancerProbes, ]
  return(df.enhancer.probes)
}


#' mapTranscriptToGene
#' @description map the miRNA precursor names to HGNC
#' @param transcripts vector with the name of miRNA precursors
#' @keywords internal
#'
#' @return a dataframe with two columns: "Transcript" indicating the miRNA precursor names, "Gene_name" indicating the actual human gene names (HGNC)
#'

mapTranscriptToGene <- function(transcripts){

  # first, we need to get the precursor names for the mature miRNA transcripts
  #df.miRNA.mapping = data.frame(miRBaseConverter :: miRNA_MatureToPrecursor(transcripts))
  #df.miRNA.mapping = df.miRNA.mapping[!is.na(df.miRNA.mapping$Precursor), ]

  df.miRNA.mapping = data.frame(Transcript = transcripts)
  df.miRNA.mapping = df.miRNA.mapping[!is.na(df.miRNA.mapping$Transcript), , drop = F]
  gene_name = c()
  for (tr in transcripts){
    elements <- family <- number <- name <- NULL
    elements = unlist(strsplit(tr, '-'))
    family = toupper(elements[2])
    # sort out family names for the let-7 miRNAs (hsa-let-7a -> MIRLET7) from the MIR miRNAs ("hsa-miR-16" -> MIR16)
    if(family == "LET"){
      family = "MIRLET"
    }
    # sort out the different cases for "hsa-mir-25 -> MIR25" vs. "hsa-mir-24-2 -> MIR24-2"
    if(length(elements) == 3){
      number = elements[3]
    } else{
      # when the gene ID has an additional number after the slash, the actual gene number will depend on whether the last digit is a letter or a number
      # if the last digit is a number, we will need a slash in between, otherwise, we don't need a slash
      # e.g. hsa-mir-26a-1 => MIR26A1, hsa-mir-24-1 = > MIR24-1
      if(is.na(as.numeric(elements[3][length(elements[3])]))){
        number = paste0(elements[3], elements[4])
      }else{
        number = paste0(elements[3], "-", elements[4])
      }
    }
    # we need upper case for the letter in the number ("hsa-miR-301b" -> MIR301B)
    number = toupper(number)
    name = paste(family, number, sep = "")
    gene_name = c(gene_name, name)
  }
  df.miRNA.mapping$Gene_name = gene_name
  df.miRNA.mapping = df.miRNA.mapping[!is.na(df.miRNA.mapping$Gene_name), , drop = F]
  return(df.miRNA.mapping)
}


#' The get.survival.probe function
#' @description Get probes whose methylation state is predictive of patient survival
#' @param EpiMixResults List of objects returned from the EpiMix function
#' @param TCGA_CancerSite TCGA cancer code (e.g. "LUAD")
#' @param clinical.data (If the TCGA_CancerSite parameter has been specified, this parameter is optional) Dataframe with survival information. Must contain at least three columns: "sample.id", "days_to_death", "days_to_last_follow_up".
#' @param pval.threshold numeric value indicting the p value threshold for selecting the survival predictive probes. Survival time is compared by log-rank test. Default: 0.05
#' @param OutputRoot path to save the output. If not null, the return value will be saved as "Survival)Probes.csv".
#' @return a dataframe with probes whose methylation state is predictive of patient survival and the p value.
#' @export
#' @examples
#' \dontrun{
#' library(survival)
#'
#' data("Sample_EpiMixResults_miRNA")
#'
#' # Set the TCGA cancer site.
#' CancerSite = "LUAD"
#'
#' # Find survival-associated CpGs/genes
#' survival.CpGs <- get.survival.probe (EpiMixResults = Sample_EpiMixResults_miRNA,
#'                                      TCGA_CancerSite = CancerSite)
#' }
#'

get.survival.probe <- function(EpiMixResults,
                               TCGA_CancerSite = NULL,
                               clinical.data = NULL,
                               pval.threshold = 0.05,
                               OutputRoot = ""){

  if(!requireNamespace("survival")){
    message("This function requires the 'survival' package.")
    return(invisible())
  }

  if(is.null(TCGA_CancerSite) & is.null(clinical.data)){
    stop("Please provide the value for either TCGA_CancerSite or clinical.data")
  }

  FunctionalPairs = EpiMixResults$FunctionalPairs

  # Group the genes associated with the same CpG together
  ProbeAnnotation =  FunctionalPairs %>%
                     dplyr :: group_by(.data$Probe) %>%
                     dplyr :: mutate(Genes = paste(.data$Gene, collapse = ";")) %>%
                     dplyr :: select(.data$Probe, .data$Genes)

  # Download clinical data
  if(!is.null(TCGA_CancerSite)){
    cat(paste0("Downloading patient clinical data for ", TCGA_CancerSite, "\n"))
    clinical.data.directory = get_firehoseData(TCGA_acronym_uppercase = TCGA_CancerSite,
                                               saveDir = tempdir(),
                                               dataFileTag = "Merge_Clinical.Level_1",
                                               )
    clinical.data = data.table::fread(paste0(clinical.data.directory,TCGA_CancerSite,".merged_only_clinical_clin_format.txt"))
    clinical.data=as.matrix(clinical.data)
    rownames(clinical.data)=clinical.data[,1]
    clinical.data = clinical.data[,-1]
    clinical.data = clinical.data[-1,]
    survival.info = data.frame(sample.id = toupper(clinical.data["patient.bcr_patient_barcode", ]),
                               days_to_death = as.numeric(clinical.data["patient.days_to_death", ]),
                               days_to_last_follow_up = as.numeric(clinical.data["patient.days_to_last_followup", ]))
  }else{
    survival.info = clinical.data
  }

  # Construct a dataframe for survival information. Status: 1 = censored, 2 = dead
  survival.info <- survival.info %>%
    dplyr :: mutate(
                    time =  dplyr :: if_else(is.na(.data$days_to_death), .data$days_to_last_follow_up, .data$days_to_death),
                    status =  dplyr :: if_else(is.na(.data$days_to_death), 1, 2)
                    )
  survival.info <- survival.info[!is.na(survival.info$time), ]
  survival.info <- survival.info[survival.info$time != 0, ]

  # Survival analysis
  meth.state <- getMethStates(EpiMixResults, EpiMixResults$MethylationDrivers)
  DMvalues <- EpiMixResults$MethylationStates
  DMvalues <-TCGA_GENERIC_CleanUpSampleNames(DMvalues, 12)

  cat("Finding survival associated genes\n")
  target.probes <- rownames(DMvalues)
  iterations <- length(target.probes)
  Probe <- State <- pval <- hazard.ratio <- low.conf <- high.conf <- character(0)
  for(i in seq(1:iterations)){
    target.probe <- state <- meth.values <- abnormal <- normal <- mixture.group <- target.survival <- sdf <- P.value <- sde <- HR <- low.cl <- high.cl <-  NULL
    target.probe <- target.probes[i]
    state <- meth.state[target.probe]
    DM.value <- DMvalues[target.probe,]
    if(state == "Hypo"){
      abnormal <- names(DM.value)[which(DM.value < 0)]
      normal <-  names(DM.value)[which(DM.value == 0)]
    }else if(state == "Hyper"){
      abnormal <- names(DM.value)[which(DM.value > 0)]
      normal <-  names(DM.value)[which(DM.value == 0)]
    }else{
      abnormal <- names(DM.value)[which(DM.value < 0)]
      normal <-  names(DM.value)[which(DM.value > 0)]
    }

    if(length(abnormal) < 20 | length(normal) < 20){
      next()
    }
    normal <- data.frame(sample.id = normal, State = 1)
    abnormal <- data.frame(sample.id = abnormal, State = 2)
    mixture.group <- rbind(normal, abnormal)
    target.survival <-merge(survival.info, mixture.group)
    sdf <- survival :: survdiff(Surv(time, status) ~ State, data = target.survival)
    P.value <- 1 - stats :: pchisq(sdf$chisq, length(sdf$n) - 1)
    sde <- survival :: coxph(Surv(time, status) ~ State, data = target.survival)
    HR <- stats :: coef(summary(sde))[, 2]
    low.cl <- exp(stats :: confint(sde))[,1]
    high.cl <- exp(stats :: confint(sde))[,2]
    Probe <- c(Probe, target.probe)
    State <- c(State, state)
    pval <- c(pval, P.value)
    hazard.ratio <- c(hazard.ratio, HR)
    low.conf <- c(low.conf, low.cl)
    high.conf <- c(high.conf, high.cl)
  }
  survival.results <- data.frame(Probe = Probe,
                                 State = State,
                                 HR = hazard.ratio,
                                 lower.cl = low.conf,
                                 higher.cl = high.conf,
                                 p.value = pval)
  survival.results <- merge(survival.results, ProbeAnnotation)
  survival.results <- survival.results[order(survival.results$p.value), ]
  survival.results <- survival.results[survival.results$p.value < pval.threshold, , drop = FALSE]
  survival.results <- survival.results %>%
                      dplyr :: select(.data$Probe, .data$Genes, .data$State, .data$HR, .data$lower.cl, .data$higher.cl,.data$p.value)
  rownames(survival.results) = NULL
  cat("Found", nrow(survival.results), "survival predictive CpGs\n")

  if(OutputRoot!="" & length(survival.results) > 0){
    cat(paste0("Saving the result to ", OutputRoot, "/", "Survival_CpGs.csv"))
    utils :: write.csv(survival.results, paste0(OutputRoot, "/", "Survival_CpGs.csv"), row.names = FALSE)
  }
  return(survival.results)
}


#' EpiMix_PlotSurvival function
#' @description function to plot Kaplan-meier survival curves for patients with different methylation state of a specific probe.
#' @param EpiMixResults List of objects returned from the EpiMix function
#' @param plot.probe Character string with the name of the probe
#' @param TCGA_CancerSite TCGA cancer code (e.g. "LUAD")
#' @param clinical.df (If the TCGA_CancerSite parameter has been specified, this parameter is optional) Dataframe with survival information. Must contain at least three columns: "sample.id", "days_to_death", "days_to_last_follow_up".
#' @param font.legend numeric value indicating the font size of the figure legend. Default: 16
#' @param font.x numeric value indicating the font size of the x axis label. Default: 16
#' @param font.y numeric value indicating the font size of the y axis label. Default: 16
#' @param font.tickslab numeric value indicating the font size of the axis tick label. Default: 14
#' @param legend numeric vector indicating the x,y coordinate for positioning the figure legend. c(0,0) indicates bottom left, while c(1,1) indicates top right. Default: c(0.8,0.9). If "None", legend will be removed.
#' @param width numeric value with graph width in pixels to save. Default: 700
#' @param height numeric values with graph height in pixels to save. Default: 480
#' @param OutputRoot file path to save the graph.
#'
#' @return Kaplan-meier survival curve showing the survival time for patients with different methylation states of the probe.
#' @export
#' @examples
#' \dontrun{
#' library(survival)
#' library(survminer)
#'
#' # Select the target CpG site whose methylation states will be evaluated
#' Probe = "cg00909706"
#'
#' # Generate the graph
#' EpiMix_PlotSurvival(EpiMixResults = Sample_EpiMixResults_miRNA,
#'                     plot.probe = Probe,
#'                     TCGA_CancerSite = CancerSite)
#' }
#'

EpiMix_PlotSurvival <- function(EpiMixResults,
                                plot.probe,
                                TCGA_CancerSite = NULL,
                                clinical.df = NULL,
                                font.legend = 16,
                                font.x = 16,
                                font.y = 16,
                                font.tickslab = 14,
                                legend= c(0.8,0.9),
                                width = 700,
                                height = 480,
                                OutputRoot = ""){

  if(!requireNamespace("survival")){
    message("This function requires the 'survival' package.")
    return(invisible())
  }

  if(!requireNamespace("survminer")){
    message("This function requires the 'survival' package.")
    return(invisible())
  }

  Classifications <- EpiMixResults$Classifications
  Classifications <-TCGA_GENERIC_CleanUpSampleNames(Classifications, 12)

  mixture.group <- Classifications[plot.probe,]
  mixture.group <- data.frame(sample.id = names(mixture.group), State = mixture.group)
  rownames(mixture.group) = NULL

  # Download clinical data
  if(!is.null(TCGA_CancerSite)){
    file_path = paste0("gdac_20160128/gdac.broadinstitute.org_", TCGA_CancerSite, ".Merge_Clinical.Level_1.2016012800.0.0/", TCGA_CancerSite, ".merged_only_clinical_clin_format.txt")
    if(!file.exists(file_path)){
      cat(paste0("Downloading patient clinical data for ", TCGA_CancerSite))
      clinical.data.directory = get_firehoseData(TCGA_acronym_uppercase = TCGA_CancerSite,
                                                 dataFileTag = "Merge_Clinical.Level_1",
                                                 saveDir = tempdir())
      clinical.data = data.table::fread(paste0(clinical.data.directory, TCGA_CancerSite, ".merged_only_clinical_clin_format.txt"))
    }else{
      clinical.data = data.table::fread(file_path)
    }
    clinical.data=as.matrix(clinical.data)
    rownames(clinical.data)=clinical.data[,1]
    clinical.data = clinical.data[,-1]
    clinical.data = clinical.data[-1,]
    survival.info = data.frame(sample.id = toupper(clinical.data["patient.bcr_patient_barcode", ]),
                               days_to_death = as.numeric(clinical.data["patient.days_to_death", ]),
                               days_to_last_follow_up = as.numeric(clinical.data["patient.days_to_last_followup", ]))
  }else{
    survival.info = clinical.df
  }

  # Construct a dataframe for survival information. Status: 1 = censored, 2 = dead
  survival.info <- survival.info %>%
    mutate(
      time = dplyr :: if_else(is.na(.data$days_to_death), .data$days_to_last_follow_up, .data$days_to_death),
      status = dplyr :: if_else(is.na(.data$days_to_death), 1, 2)
    )
  survival.info <- survival.info[!is.na(survival.info$time), ]
  survival.info <- survival.info[survival.info$time != 0, ]
  target.survival <-merge(survival.info, mixture.group)
  target.survival <- target.survival[order(target.survival$State), ]
  survival <- survminer :: ggsurvplot(
    survminer :: surv_fit(survival :: Surv(time, status) ~ State, data = target.survival),
    pval = TRUE,
    legend.title = "mixture component",
    legend.labs = unique(target.survival$State),
    xlab = "Days",
    ylab = "Overall survival probability",
    palette = RColorBrewer::brewer.pal(8, "Set1")[1:length(unique(target.survival$State))],
    font.legend = c(font.legend, "bold"),
    font.x = c(font.x, "bold"),
    font.y = c(font.y, "bold"),
    font.tickslab = c(font.tickslab),
    legend = legend
  )
  if(OutputRoot != ""){
    png(filename = paste0(OutputRoot, "/", plot.probe, "_survival.png"),
        width=700, height=480)
    print(survival)
    dev.off()
  }
  return(survival)
}


#' The function.enrich function
#' @description Perform functional enrichment analysis for the differentially methylated genes occurring in the significant CpG-gene pairs.
#' @param EpiMixResults List of the result objects returned from the EpiMix function.
#' @param methylation.state character string indicating whether to use all the differentially methylated genes or only use the hypo- or hyper-methylated genes for enrichment analysis. Can be either "all", "Hyper" or "Hypo".
#' @param enrich.method character string indicating the method to perform enrichment analysis, can be either "GO" or "KEGG".
#' @param ont character string indicating the aspect for GO analysis. Can be one of "BP" (i.e., biological process), "MF" (i.e., molecular function), and "CC" (i.e., cellular component) subontologies, or "ALL" for all three.
#' @param simplify boolean value indicating whether to remove redundancy of enriched GO terms.
#' @param cutoff if simplify is TRUE, this is the threshold for similarity cutoff of the ajusted p value.
#' @param pvalueCutoff adjusted pvalue cutoff on enrichment tests to report
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param qvalueCutoff qvalue cutoff on enrichment tests to report as significant. Tests must pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on adjusted pvalues and iii) qvalueCutoff on qvalues to be reported.
#' @param save.dir path to save the enrichment table.
#' @return a clusterProfiler enrichResult instance
#' @export
#' @examples
#' \dontrun{
#' library(clusterProfiler)
#' library(org.Hs.eg.db)
#'
#' data(Sample_EpiMixResults_Regular)
#'
#' # Check the functions of both the hypo- and hypermethylated genes.
#' methylation.state = "all"
#'
#' # Use the gene ontology for functional analysis.
#' enrich.method = "GO"
#'
#' # Use the "biological process" sub-term
#' selected.pathways = "BP"
#'
#' # Perform enrichment analysis
#' enrich.results <- function.enrich(
#'   EpiMixResults = Sample_EpiMixResults_Regular,
#'   methylation.state = methylation.state,
#'   enrich.method = enrich.method,
#'   ont = selected.pathways,
#'   simplify = TRUE,  # get rid of overlapping pathways
#'   save.dir = ""
#' )
#'
#' }
#'
function.enrich <- function(EpiMixResults,
                            methylation.state = "all",
                            enrich.method = "GO",
                            ont = "BP",
                            simplify = TRUE,
                            cutoff = 0.7,
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.2,
                            save.dir = "."){

  if(!requireNamespace("clusterProfiler")){
    message("This function requires the 'clusterProfiler' package.")
    return(invisible())
  }

  if(!requireNamespace("org.Hs.eg.db")){
    message("This function requires the 'org.Hs.eg.db' package.")
    return(invisible())
  }

  # input check
  if(enrich.method != "GO" & enrich.method != "KEGG"){
    stop("enrich.method must be 'GO' or 'KEGG' !")
  }
  if(methylation.state != "all" & methylation.state != "Hypo" & methylation.state != "Hyper"){
    stop("methylaiton.state must be 'all', 'Hyper', or 'Hypo' !")
  }

  FunctionalPairs <- EpiMixResults$FunctionalPairs
  if(methylation.state == "Hypo"){
    FunctionalPairs <- FunctionalPairs[FunctionalPairs$State == "Hypo",]
  }else if(methylation.state == "Hyper"){
    FunctionalPairs <- FunctionalPairs[FunctionalPairs$State == "Hyper",]
  }
  df.gene = FunctionalPairs %>%
            dplyr :: select(.data$Gene, .data$`Fold change of gene expression`) %>%
            dplyr :: group_by(.data$Gene) %>%
            dplyr :: summarize(Avg.expr = mean(.data$`Fold change of gene expression`, na.rm = TRUE))
  gene_id_map = AnnotationDbi::select(org.Hs.eg.db, keys = df.gene$Gene, columns = c("ENTREZID"), keytype = "SYMBOL")
  if(sum(is.na(gene_id_map$ENTREZID)) > 0){
    absentGenes = gene_id_map$SYMBOL[which(is.na(gene_id_map$ENTREZID))]
    warning(paste0("Can not find ENTREZID for ", length(absentGenes)," genes.", "These genes can not be included in the enrichment analysis:\n", paste0(absentGenes, collapse = ",")))
  }
  colnames(gene_id_map) = c("Gene", "ENTREZID")
  df.gene = merge(df.gene, gene_id_map)
  gene.vector = df.gene$Avg.expr
  names(gene.vector) = df.gene$ENTREZID
  if(enrich.method == "GO"){
    ego <- clusterProfiler :: enrichGO(names(gene.vector), OrgDb = "org.Hs.eg.db", ont = ont, pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff, readable = TRUE)
  }else if(enrich.method == "KEGG"){
    ego <- clusterProfiler :: enrichKEGG(names(gene.vector), pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff, readable = TRUE)
  }
  if(simplify & enrich.method == "GO"){
    ego <- clusterProfiler :: simplify(ego, cutoff = cutoff)
  }

  if(length(save.dir) > 0){
    save.file.name <- paste0(save.dir, "/", "FunctionEnrichment_", enrich.method, ".csv")
    utils :: write.csv(ego@result, save.file.name, row.names = FALSE)
  }
  return(ego)
}

# Unregister sockets, in case the socket was not completely removed from last run
# unregister <- function() {
#   env <- foreach::.foreachGlobals
#   rm(list=ls(name=env), pos=env)
# }

# Function to retrieve data from the data package
EpiMix_GetData <- function(...)
{
  e <- new.env()
  name <- utils :: data(..., package = "EpiMix.data",envir = e)[1]
  e[[ls(envir = e)[1]]]
}


mapGeneToEntrezID <- function(gene_names)


  list.epigenomes<- function(){
    EpigenomeMap = EpiMix_GetData("EpigenomeMap")
    print(EpigenomeMap)
  }

list.chromatin.states <- function(){
  chromatin.states = c("EnhG1"	= "Genic enhancer1",
                       "EnhG2" = "Genic enhancer2",
                       "EnhA1" =	"Active Enhancer 1",
                       "EnhA2" =	"Active Enhancer 2",
                       "EnhWk"	= "Weak Enhancer",
                       "EnhBiv" = "Bivalent Enhancer")
  get("chromatin.states")
}

#' The validEpigenomes function
#' @description check user input for roadmap epigenome groups or ids
#' @param roadmap.epigenome.groups epigenome groups
#' @param roadmap.epigenome.ids epigenome ids
#'
#' @return a character vector of selected epigenome ids
#'
validEpigenomes <- function(roadmap.epigenome.groups, roadmap.epigenome.ids){
  EpigenomeMap = EpiMix_GetData("EpigenomeMap")
  selectedEpigenomes = character(0)
  if(!is.null(roadmap.epigenome.groups)){
    if("ES_deriv." %in% roadmap.epigenome.groups) roadmap.epigenome.groups[which(roadmap.epigenome.groups == "ES_deriv.")] = "ES_deriv"
    if("Mesench." %in% roadmap.epigenome.groups) roadmap.epigenome.groups[which(roadmap.epigenome.groups == "Mesench.")] = "Mesench"
    if("Myosat." %in% roadmap.epigenome.groups) roadmap.epigenome.groups[which(roadmap.epigenome.groups == "Myosat.")] = "Myosat"
    if("Neurosph." %in% roadmap.epigenome.groups) roadmap.epigenome.groups[which(roadmap.epigenome.groups == "Neurosph.")] = "Neurosph"

    overlapGroups = intersect(names(EpigenomeMap), roadmap.epigenome.groups)
    if(length(overlapGroups) !=  length(roadmap.epigenome.groups)){
      warning(paste0("Can not find epigenome group(s): ",
                     roadmap.epigenome.groups[which(!roadmap.epigenome.groups %in% overlapGroups)]),
              immediate. = TRUE)
      cat("Available epigenome groups: \n", paste0(names(EpigenomeMap), collapse = ",  "), "\n")
    }
    for(gr in overlapGroups){
      selectedEpigenomes = c(selectedEpigenomes,  EpigenomeMap[[gr]])
    }
  }
  if(!is.null(roadmap.epigenome.ids)){
    all.ids = character(0)
    for (i in 1:length(EpigenomeMap)){
      all.ids = c(all.ids, EpigenomeMap[[i]])
    }
    overlapIDs = intersect(roadmap.epigenome.ids, all.ids)
    if(length(overlapIDs) != length(roadmap.epigenome.ids)){
      warning(paste0("Can not find some of the epigenome ids: ",roadmap.epigenome.ids[which(!roadmap.epigenome.ids %in% overlapIDs)]), immediate. = TRUE)
      cat("Available epigenome ids: \n", paste0(all.ids[order(all.ids)], collapse = ", "), "\n")
    }
    selectedEpigenomes = c(selectedEpigenomes, overlapIDs)
  }
  selectedEpigenomes = unique(selectedEpigenomes)
  return(selectedEpigenomes)
}











