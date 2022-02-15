###################################################################################################################
#                                   Utility functions used by EpiMix
###################################################################################################################
#
#' The EpiMix_getInfiniumAnnotation function
#' @description fetch the Infinium probe annotation from the seasameData as a GRange object
#' @param plat character string indicating the methylation platform
#' @param genome character string indicating the version of genome build
#'
#' @return a GRange object of probe annotation
#' @keywords internal

EpiMix_getInfiniumAnnotation <- function(plat = "EPIC", genome = "hg38"){
  ProbeAnnotation = NULL
  if(tolower(genome) == "hg19" & toupper(plat) == "HM27" ) ProbeAnnotation = sesameData::sesameDataGet("HM27.hg19.manifest")
  if(tolower(genome) == "hg19" & toupper(plat) == "HM450" ) ProbeAnnotation = sesameData::sesameDataGet("HM450.hg19.manifest")
  if(tolower(genome) == "hg19" & toupper(plat) == "EPIC" ) ProbeAnnotation = sesameData::sesameDataGet("EPIC.hg19.manifest")
  if(tolower(genome) == "hg38" & toupper(plat) == "HM27" ) ProbeAnnotation = sesameData::sesameDataGet("HM27.hg38.manifest")
  if(tolower(genome) == "hg38" & toupper(plat) == "HM450" ) ProbeAnnotation = sesameData::sesameDataGet("HM450.hg38.manifest")
  if(tolower(genome) == "hg38" & toupper(plat) == "EPIC" ) ProbeAnnotation = sesameData::sesameDataGet("EPIC.hg38.manifest")
  return(ProbeAnnotation)
}

#' The convertAnnotToDF function
#' @description convert the probe annotation from the GRange object to a dataframe and extract chromosome information
#' @param annot a GRange object of probe annotation, can be the object returned from the getInfiniumAnnotation function.
#' @return a dataframe with chromosome, beginning and end position, mapped gene information for each CpG probe
#' @import GenomicRanges
#' @keywords internal
#'
convertAnnotToDF <- function(annot){
  df.annot = data.frame("CpG_chrm" = seqnames(annot),
                        "CpG_beg" = start(ranges(annot)),
                        "CpG_end" = end(ranges(annot)),
                        "probeID" = names(annot),
                        "gene" = mcols(annot)$gene)
 return(df.annot)
}


#' The mapProbeGene function
#' @param df.annot a dataframe with probe annotation, can be the object returned from the convertAnnotToDF function.
#' @description since in the original probe annotation, a specific probe can be mapped to multiple genes, this function splits the rows and maps each probe to a signle gene in a row.
#' @return a dataframe with 1:1 mapping of probe and gene
#' @import tidyr
#' @import dplyr
#' @keywords internal
#'

mapProbeGene <- function(df.annot){
  df.annot <- df.annot[!is.na(df.annot$probeID) & !is.na(df.annot$gene), ]
  df.annot <- tidyr :: separate_rows(df.annot, gene, sep = ";")
  df.annot <- dplyr :: distinct(df.annot)
  return(df.annot)
}


#' The getMethStates function
#' @description Helper function that adds a methyaltion state label to each driver probe
#' @param MethylMixResults the list object returned from the EpiMix function
#' @param DM.probes character vector of differentially methylated probes.
#' @return a character vector with the methylation state ("Hypo", "Hyper" or "Dual") for each probe. The names for the vector are the probe names and the values are the methylation state.
#' @keywords internal
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
#' @description  Helper function to assess if the expression of a specific gene is reversely correlated with the methylation of a list of probes mapped to it.
#' @details This function is gene-centered, which is used in the regular mode and the lncRNA mode of EpiMix.
#' @param gene character string indicating the target gene.
#' @param probes character vector indicating the probes mapped to the gene.
#' @param MET_matrix methylation data matrix.
#' @param exp gene expression data matrix,
#' @return dataframe with functional probe-gene pairs and corresponding p values obtained from the Wilcoxon test for gene expression and methylation.
#'
getFunctionalProbes <- function(gene, probes, MET_matrix, MET_Control, exp, methylation.states,raw.pvalue.threshold = 0.05, adjusted.pvalue.threshold = 0.01){
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

# splitmatix
# @param x A matrix
# @param by A character specify if split the matrix by row or column.
# @return A list each of which is the value of each row/column in the matrix.
splitmatrix <- function(x,by="row") {
  if(by %in% "row"){
    out <- split(x, rownames(x))
  }else if (by %in% "col"){
    out <- split(x, colnames(x))
  }
  return(out)
}

get.chromosome <- function(genes, genome){
  gene.annotation = as.data.frame(getTSS(genome))
  gene.annotation = distinct(gene.annotation[,c("external_gene_name", "seqnames")])
  gene.annotation = gene.annotation[which(gene.annotation$external_gene_name %in% genes),]
  colnames(gene.annotation) = c("Gene", "Chr")
  return(gene.annotation)
}

get.perm.pVals <- function(target.probe,
                           MET_matrix,
                           gene.expression.data,
                           ProbeAnnotation,
                           genome = "hg38",
                           correlation = "negative"){

  target.chr = as.character(seqnames(ProbeAnnotation)[which(names(ProbeAnnotation) == target.probe)])
  all.genes = rownames(gene.expression.data)
  gene.chr = get.chromosome(all.genes, genome)
  random.genes = sample(gene.chr$Gene[which(!gene.chr$Chr %in% c(target.chr, "chrX", "chrY"))], 1000)
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
  test.p <- data.frame(GeneID=random.genes,
                       Perm.p=test.p[match(random.genes, names(test.p))],
                       stringsAsFactors = FALSE)
}




#' The getFunctionalGenes function
#' @description  Helper function to assess if the methylation of a probe is reversely correlated with the expression of its nearby genes.
#' @details This function is probe-centered, which is used in the enhancer mode and the miRNA mode of EpiMix.
#' @param probe character string indicating the probe to be evaluated.
#' @param genes character vector indicating the nearby genes for the target probe.
#' @param MET_matrix methylation data matrix.
#' @param exp gene expression data matrix.
#' @param state character string indicating the methylation state of the probe, should be either "Hyper", "Hypo" or "Dual".
#' @return dataframe with functional probe-gene pair and p values from the Wilcoxon test for methylation and gene expression.
#'
getFunctionalGenes <- function(target.probe,
                               target.genes,
                               MET_matrix,
                               MET_Control,
                               gene.expression.data,
                               ProbeAnnotation,
                               state,raw.pvalue.threshold = 0.05,
                               adjusted.pvalue.threshold = 0.01){
  sample.number = ncol(MET_matrix)
  if(!is.null(MET_Control)){
    sample.number = ncol(MET_matrix) - length(intersect(colnames(MET_matrix), colnames(MET_Control)))
  }
  DM_values = MET_matrix[target.probe, ] # a single-row vector with the names as sample names, and the values as DM values
  high.met.samples = names(DM_values[DM_values == max(DM_values)])
  low.met.samples = names(DM_values[DM_values == min(DM_values)])

  target.genes = intersect(target.genes, rownames(gene.expression.data))
  if(length(target.genes) == 0) return(NULL)
  gene.expression.data = gene.expression.data[target.genes, ,drop = F]
  expr.low = gene.expression.data[target.genes,high.met.samples, drop = F]  # high methylation, low gene expression
  expr.high = gene.expression.data[target.genes,low.met.samples, drop = F]  # low methylation, high gene expression

  # calculate prevalence
  hypo_prev <-  hyper_prev <- NULL
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

   # modeling gene expression
  p = c()
  mRNA.fold.change = c()
  comparisons = c()
  for (gene in target.genes){
    expr.high.values = as.numeric(expr.high[gene,])
    expr.low.values = as.numeric(expr.low[gene,])
    p_value = wilcox.test(expr.high.values, expr.low.values, alternative = "greater")$p.value
    if(state == "Hyper"){
      fold_change = round(mean(expr.low.values)/mean(expr.high.values), 3)
      cmp = "hyper vs normal"
    }else if(state == "Dual"){
      fold_change = round(mean(expr.high.values)/mean(expr.low.values), 3)
      cmp ="hypo vs hyper"
    }else{
      fold_change = round(mean(expr.high.values)/mean(expr.low.values), 3)
      cmp ="hypo vs normal"
    }
    p = append(p, p_value)
    mRNA.fold.change = append(mRNA.fold.change, fold_change)
    comparisons = append(comparisons, cmp)
  }

  # produce a dataframe for gene expression with methylation state and prevalence inforamtion
  dataDEGs <- data.frame(Gene = target.genes,
                         Probe = rep(target.probe,length(target.genes)),
                         State = rep(state,length(target.genes))
                         )

  # find the permutation p values of this probe
  perm.pvals = get.perm.pVals(target.probe,
                              MET_matrix,
                              gene.expression.data,
                              ProbeAnnotation,
                              genome = "hg38",
                              correlation = "negative")

  dataDEGs['Proportion of hypo (%)'] = rep(hypo_prev, length(target.genes))
  dataDEGs['Proportion of hyper (%)'] = rep(hyper_prev, length(target.genes))
  dataDEGs['Fold change of gene expression'] = mRNA.fold.change
  dataDEGs['Comparators'] = comparisons
  dataDEGs['Raw.p'] = p
  dataDEGs = Get.Pvalue.p(dataDEGs, perm.pvals)
  #dataDEGs["Adjusted.p"] <- p.adjust(as.vector(unlist(dataDEGs['Raw.p'])), method='fdr')
  dataDEGs <- dataDEGs[which(dataDEGs$Raw.p<raw.pvalue.threshold & dataDEGs$Adjusted.p < adjusted.pvalue.threshold),]
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
  message("Calculating empirical P value.\n")
  Pvalue <- unlist(apply(U.matrix,1,.Pvalue,permu=permu))
  U.matrix$Adjusted.p <- Pvalue
  return(U.matrix)
}








#' The filterMethMatrix function
#' @details This function filters methylation states from the beta mixture modeling for each probe. The filtered probes can be used to model gene expression by Wilcoxon test.
#' @param MethylMixResults
#'
#' @return A matrix of methylation states for each differentially methylated probe with probes in rows and patient in columns.
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
#' @import GenomicRanges
#' @importFrom utils read.table
#' @keywords internal
#' @examples
#'\dontrun{
# met.platform = "EPIC"
# genome = "hg38"
# listOfEpigenomes = c("E034", "E045", "E047")
# functional.regions = c("EnhA1", "EnhA2", "EnhG1", "EnhG2")
# target.Dir <- paste0(getwd(),"/Annotation")
# df.enhancer.probes <-  getEnhancerProbes(met.platform = met.platform,
#                                          genome = genome,
#                                          functional.regions = functional.regions,
#                                          listOfEpigenomes = listOfEpigenomes) # total 551,335 probes for EPIC array, including 481,961("EnhA1", "EnhA2") probes
#
# #write.csv(df.enhancer.probes,paste0(target.Dir, "/Roadmap_", met.platform,"_Enhancer_Probes.csv"), row.names = FALSE)
# saveRDS(df.enhancer.probes,paste0(target.Dir, "/Roadmap_", met.platform,"_Enhancer_Probes.rds"))
#' }

getRoadMapEnhancerProbes <- function(met.platform = "EPIC", genome = "hg38", functional.regions=c("EnhA1", "EnhA2"), listOfEpigenomes = NULL){

  # Step 1. Find all the filenames ending with "_18_core_K27ac_hg38lift_mnemonics.bed.gz" from the Roadmap Epigenomics web portal
  K27Ac_url <- "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/"
  urlData=getURL(K27Ac_url)
  urlData2=unlist(strsplit(urlData,"\\n"))
  filenames=as.matrix(urlData2[grep("_18_core_K27ac_hg38lift_mnemonics.bed.gz",urlData2)])
  filenames=unlist(strsplit(filenames, ">|<"))
  filenames=filenames[grep("_18_core_K27ac_hg38lift_mnemonics.bed.gz",filenames)]
  filenames=filenames[-grep("a href",filenames)]
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
    gr.enhancer = makeGRangesFromDataFrame(genome.state, ignore.strand = TRUE, keep.extra.columns = TRUE)
    target.probes = ProbeAnnotation[unique(queryHits(findOverlaps(ProbeAnnotation,gr.enhancer, ignore.strand = TRUE)))]
    this.epigenome = unlist(strsplit(file, "_"))[1]
    cat("\tIdentifed", length(target.probes), "enhancer probes from the epigenome", this.epigenome, "\n")
    enhancerProbes = c(enhancerProbes, names(target.probes))
    enhancerProbes = unique(enhancerProbes)
  }
  # Step 5. Generate a dataframe for enhancer probes with their coordinates
  probe.ID = names(ProbeAnnotation)
  probe.chr = seqnames(ProbeAnnotation)
  probe.start.pos = start(ranges(ProbeAnnotation))
  probe.start.end = end(ranges(ProbeAnnotation))
  df.ProbeAnnotation <- data.frame(probeID = probe.ID, chr=probe.chr,start =probe.start.pos, end = probe.start.end)
  df.enhancer.probes <-  df.ProbeAnnotation[df.ProbeAnnotation$probeID %in% enhancerProbes, ]
  return(df.enhancer.probes)
}


#' mapTranscriptToGene
#' @description map the miRNA precursor names to HGNC
#' @param transcripts vector with the name of miRNA precursors
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

#'
#' The getMiRNetwork function
#' @param EpiMixResults list object returned from EpiMix.
#' @param tissue character string indicating the target tissue type
#'
#' @return a list include three matrices: (1) target: map of miRNA genes to their respective target genes; (2) family: family enrichment analysis of the methylaiton-driven miRNAs; (3) pathway: pathway enrichment analysis of miRNA target genes.

getMiRNetwork <- function(EpiMixResults, tissue = "Lung", enrich.method = ""){

  if(!requireNamespace("miRNetR")){
    message("This function requires the 'miRNetR' package.")
    return(invisible())
  }

  if(!requireNamespace("jsonlite")){
    message("This function requires the 'jsonlite' package.")
    return(invisible())
  }

  Genes = unique(EpiMixResults$FunctionalPairs$Gene)
  #### Step 1. Initiate the dataSet object
  Init.Data("mir", "mirlist")
  #### Step 2. Set up the user input data
  if(tissue != ""){
    SetupMirListData(mirs = Genes, orgType = "hsa", idType = "mir_id", tissue = tissue)
  }else{
    SetupMirListData(mirs = Genes, orgType = "hsa", idType = "mir_id")
  }

  #### Step 3. Set up targets
  nms.vec = c("gene")
  SetCurrentDataMulti()

  #### Step 4. Perform miRNAs to target genes mapping,
  #### results are downloaded in the working directory ("mirnet_mir_target.csv")
  QueryMultiListMir()

  #### Step 5. Generate miRNA-gene network files
  CreateMirNets(net.type = "mir2gene")

  #### Step 6. Prepare network files,
  #### results are downloaded in the working directory ("node_table_mirnet_0.csv", "mirnet_0.json" and "mirnet.graphml")
  PrepareMirNet(mir.nm = "mirnet1", file.nm = "mirnet_0.json")

  #### Step 7. Perform miRNA family enrichment analysis,
  #### results are downloaded in the working directory ("network_enrichment_mirfamily_1.json" and "mirnet_enrichment.csv")
  PerformMirTargetEnrichAnalysis(
    adjust.type = "NA",
    fun.type = "mirfamily",
    file.nm = "network_enrichment_mirfamily_1",
    IDs = Genes,
    algo = "hyp"
  )
  resTable = read.csv("mirnet_enrichment.csv", header=T, as.is=T)

  # parse json to add the gene names to gene family
  miRFamily = fromJSON("network_enrichment_mirfamily_1.json")
  flatFamily = sapply(miRFamily$fun.anot, function(x) paste(x, collapse = ";"))
  mapFamilyToGenes = data.frame(Pathway = names(flatFamily), `Genes` = flatFamily)
  rownames(mapFamilyToGenes) <-  NULL
  resTable = merge(resTable, mapFamilyToGenes)
  resTable = resTable[, c("Pathway", "Total", "Expected", "Hits", "Genes", "Pval")]
  resTable = resTable[order(resTable$Pval),]
  colnames(resTable)[1] = "miR.family"

  #### Step 8. Perform pathway enrichment analysis
  PerformMirTargetEnrichAnalysis(
    adjust.type = "NA",
    fun.type = "kegg",
    file.nm = "network_enrichment_kegg_1",
    IDs = c(dataSet$mir.res$Target, unique(dataSet$mir.res$ID)),
    algo = "hyp"
  )

  pathway = read.csv("mirnet_enrichment.csv")
  HittedGenes = fromJSON("network_enrichment_kegg_1.json")
  flatGenes = sapply(HittedGenes$fun.anot, function(x) paste(x, collapse = ";"))
  mapPathwayToGenes = data.frame(Pathway = names(flatGenes), `Genes` =  flatGenes)
  rownames(mapPathwayToGenes) <-  NULL
  pathway = merge(pathway, mapPathwayToGenes)
  pathway = pathway[, c("Pathway", "Total", "Expected", "Hits", "Genes", "Pval")]
  pathway = pathway[order(pathway$Pval),]
  res = list(target = dataSet$mir.res, family = resTable, pathway = pathway)
  write.csv(res$target, "EpiMix_target_gene.csv", row.names = FALSE)
  write.csv(res$family, "EpiMix_miRfamily_enrichment.csv", row.names = FALSE)
  write.csv(res$pathway, "EpiMix_targetGene_pathway.csv", row.names = FALSE)
  return(res)
}

#' The get.survival.probe function
#' @description Get probes whose methylation state is predictive of patient survival
#' @param EpiMixResults List of objects returned from the EpiMix function
#' @param TCGA_CancerSite TCGA cancer code (e.g. "LUAD")
#' @param clinical.data (If the TCGA_CancerSite parameter has been specified, this parameter is optional) Dataframe with survival information. Must contain at least three columns: "sample.id", "days_to_death", "days_to_last_follow_up".
#' @param pval.threshold numeric value indicting the p value threshold for selecting the survival predictive probes. Survival time is compared by log-rank test. Default: 0.05
#' @param met.platform character string indicating the methylation platform. Default: "HM450"
#' @param mode character string indicating the EpiMix mode. Must be "Regular", "Enhancer", "miRNA" or "lncRNA"
#' @param OutputRoot path to save the output. If not null, the return value will be saved as "Survival)Probes.csv".
#' @return a dataframe with probes whose methylation state is predictive of patient survival and the p value.
#' @import GenomicRanges
#' @import dplyr
#' @importFrom stats pchisq
#' @export
#'

get.survival.probe <- function(EpiMixResults,
                               TCGA_CancerSite = NULL,
                               clinical.data = NULL,
                               pval.threshold = 0.05,
                               met.platform = "HM450",
                               mode = "Regular",
                               OutputRoot = "."){

  if(!requireNamespace("survival")){
    message("This function requires the 'survival' package.")
    return(invisible())
  }

  if(is.null(TCGA_CancerSite) & is.null(clinical.data)){
    stop("Please provide the value for either TCGA_CancerSite or clinical.data")
  }
  if(!mode %in% c("Regular", "Enhancer", "miRNA", "lncRNA")){
    stop("'mode' must be one of the following values: 'Regular', 'Enhancer', 'miRNA', 'lncRNA'")
  }

  # Retrieve probe annotation
  cat("Retriving probe annotation...\n")
  ProbeAnnotation = NULL
  if(mode == "Regular"){
    ProbeAnnotation = EpiMix_getInfiniumAnnotation(met.platform)
    ProbeAnnotation = convertAnnotToDF(ProbeAnnotation)
    ProbeAnnotation = ProbeAnnotation[, c("probeID", "gene")]
  }else if(mode == "Enhancer"){
    ProbeAnnotation = EpiMix_getInfiniumAnnotation(met.platform)
    ProbeAnnotation = convertAnnotToDF(ProbeAnnotation)
    geneAnnot <- getTSS(genome = "hg38") #ELMER function to retrieve a GRange object that contains coordinates of promoters for human genome.
    DM.probes = EpiMixResults$MethylationDrivers
    DM.probes.annotation = ProbeAnnotation[match(DM.probes,ProbeAnnotation$probeID), ,drop = F]
    probe.gr = GenomicRanges::GRanges(seqnames =  DM.probes.annotation$CpG_chrm,
                                      range = IRanges(start = DM.probes.annotation$CpG_beg, end = DM.probes.annotation$CpG_end),
                                      name = DM.probes.annotation$probeID
    )
    names(probe.gr) = DM.probes.annotation$probeID
    NearbyGenes <- GetNearGenes(geneAnnot = geneAnnot,
                                TRange = probe.gr,
                                numFlankingGenes = 20)
    ProbeAnnotation  <- NearbyGenes %>% group_by(ID) %>% mutate(gene = paste(Symbol, collapse = ";")) %>% select(ID, gene)
    ProbeAnnotation <-  distinct(ProbeAnnotation)
  }else if(mode == "miRNA"){
    ProbeAnnotation <- readRDS(paste0("Annotation/", met.platform ,"_miRNA_probes.rds"))
    ProbeAnnotation <- ProbeAnnotation %>% group_by(probe) %>% mutate(Genes = paste(gene, collapse = ";"))
    ProbeAnnotation <- distinct(ProbeAnnotation[, c("probe", "Genes")])
  }else{
    ProbeAnnotation <- readRDS(paste0("Annotation/", met.platform ,"_lncRNA_probes.rds"))
    ProbeAnnotation <- data.frame(probe = ProbeAnnotation, gene = names(ProbeAnnotation))
    ProbeAnnotation <- ProbeAnnotation %>% group_by(probe) %>% mutate(Genes = paste(gene, collapse = ";"))
    ProbeAnnotation <- distinct(ProbeAnnotation[, c("probe", "Genes")])
  }

  colnames(ProbeAnnotation) = c("Probe", "Genes")

  # Download clinical data
  if(!is.null(TCGA_CancerSite)){
    cat(paste0("Downloading patient clinical data for ", TCGA_CancerSite))
    clinical.data.directory = get_firehoseData(TCGA_acronym_uppercase = TCGA_CancerSite, dataFileTag = "Merge_Clinical.Level_1")
    clinical.data = data.table::fread(paste0(clinical.data.directory, "LUAD.merged_only_clinical_clin_format.txt"))
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
    mutate(
      time = if_else(is.na(days_to_death), days_to_last_follow_up, days_to_death),
      status = if_else(is.na(days_to_death), 1, 2)
    )
  survival.info <- survival.info[!is.na(survival.info$time), ]
  survival.info <- survival.info[survival.info$time != 0, ]

  # Survival analysis
  Classifications <- EpiMixResults$Classifications
  Classifications <-TCGA_GENERIC_CleanUpSampleNames(Classifications, 12)

  METcounts <- apply(Classifications, 1, function(x)length(unique(x)))
  Classifications <- Classifications[METcounts!=1,,drop = FALSE]
  if(ncol(Classifications) > 0){
    METminCounts <- apply(Classifications, 1, function(x)min(as.numeric(table(x))))
    Classifications <- Classifications[METminCounts>=20,]
  }

  cat("Performing survival analysis\n")
  target.probes <- rownames(Classifications)
  iterations <- length(target.probes)
  Probe <- State <- pval <- character(0)
  for(i in seq(1:iterations)){
    target.probe <- state <- meth.values <- hyper.samples <- hypo.samples <- normal.samples <- mixture.group <- target.survival <- sdf <- P.value <- NULL
    target.probe <- target.probes[i]
    mixture.group <- Classifications[target.probe,]
    mixture.group <- data.frame(sample.id = names(mixture.group), State = mixture.group)
    rownames(mixture.group) = NULL
    target.survival <-merge(survival.info, mixture.group)
    sdf <- survival :: survdiff(Surv(time, status) ~ State, data = target.survival)
    P.value <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    Probe <- c(Probe, target.probe)
    State <- c(State, state)
    pval <- c(pval, P.value)
  }
  survival.results <- data.frame(Probe = Probe, p.value = pval)
  survival.results <- merge(survival.results, ProbeAnnotation)
  survival.results <- survival.results[order(survival.results$p.value), ]
  survival.results <- survival.results[survival.results$p.value < pval.threshold, ]
  rownames(survival.results) = NULL
  cat("Found ", nrow(survival.results), "survival predictive probes\n")
  if(!is.null(OutputRoot)){
    cat(paste0("Saving the result to ", OutputRoot, "/Survival_", mode, "_Probes.csv"))
    write.csv(survival.results, paste0(OutputRoot, "/Survival_", mode, "_Probes.csv"))
  }
  return(survival.results)
}


#' The plot.survival.probe function
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
#' @import dplyr
#' @export
#'
#'
plot.survival.probe <- function(EpiMixResults,
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
                                OutputRoot = "."){

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
      clinical.data.directory = get_firehoseData(TCGA_acronym_uppercase = TCGA_CancerSite, dataFileTag = "Merge_Clinical.Level_1")
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
      time = if_else(is.na(days_to_death), days_to_last_follow_up, days_to_death),
      status = if_else(is.na(days_to_death), 1, 2)
    )
  survival.info <- survival.info[!is.na(survival.info$time), ]
  survival.info <- survival.info[survival.info$time != 0, ]
  target.survival <-merge(survival.info, mixture.group)
  target.survival <- target.survival[order(target.survival$State), ]
  survival <- survminer :: ggsurvplot(
    survminer :: surv_fit(Surv(time, status) ~ State, data = target.survival),
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
    png(file= paste0(OutputRoot, "/", plot.probe, "_survival.png"),
        width=700, height=480)
    print(survival)
    #dev.off()
  }
  return(survival)
}



#' The function.enrich function
#' @description Perform functional enrichment analysis for the differentially methylated genes occurring in the significant probe-gene pairs.
#' @param EpiMixResults List of objects returned from the EpiMix function
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
#' @import dplyr
#' @export

#'
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
    dplyr::select(Gene, `Fold change of gene expression`) %>%
    group_by(Gene) %>% summarize(Avg.expr = mean(`Fold change of gene expression`, na.rm = TRUE))
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
   ego <- enrichGO(names(gene.vector), OrgDb = "org.Hs.eg.db", ont = ont, pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff, readable = TRUE)
  }else if(enrich.method == "KEGG"){
    ego <- enrichKEGG(names(gene.vector), pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff, readable = TRUE)
  }
  if(simplify & enrich.method == "GO"){
    ego <- simplify(ego, cutoff = cutoff)
  }

  if(length(save.dir) > 0){
    save.file.name <- paste0(save.dir, "/", "FunctionEnrichment_", enrich.method, ".csv")
    write.csv(ego@result, save.file.name, row.names = FALSE)
  }
  return(ego)
}


# Unregister sockets, in case the socket was not completely removed from last run
unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

# Function to retrieve data from the data package
EpiMix_GetData <- function(...)
{
  e <- new.env()
  name <- data(..., package = "EpiMix.data",envir = e)[1]
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
  selectedEpigenomes = character(0)
  if(!is.null(roadmap.epigenome.groups)){
    if("ES_deriv." %in% roadmap.epigenome.groups) roadmap.epigenome.groups[which(roadmap.epigenome.groups == "ES_deriv.")] = "ES_deriv"
    if("Mesench." %in% roadmap.epigenome.groups) roadmap.epigenome.groups[which(roadmap.epigenome.groups == "Mesench.")] = "Mesench"
    if("Myosat." %in% roadmap.epigenome.groups) roadmap.epigenome.groups[which(roadmap.epigenome.groups == "Myosat.")] = "Myosat"
    if("Neurosph." %in% roadmap.epigenome.groups) roadmap.epigenome.groups[which(roadmap.epigenome.groups == "Neurosph.")] = "Neurosph"

    EpigenomeMap = EpiMix_GetData("EpigenomeMap")
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



#' The addTSS function
#' @param FunctionalPairs Functional pairs from EpiMix Results
#' @param met.platform methylation platform
#' @param genome genome build version
#'
#' @return
#'
#'EpiMixResults <- readRDS("Results/TCGA_LUAD/EpiMix_Results_regular.rds")
#FunctionalPairs <- EpiMixResults$FunctionalPairs

addTSS <- function(FunctionalPairs, met.platform = "HM450", genome = "hg38"){
  ProbeAnnotation = as.data.frame(EpiMix_getInfiniumAnnotation(plat = met.platform, genome = genome))
  new_df <- FunctionalPairs[, c("Gene", "Probe")]
  new_df$CpG_start <- ProbeAnnotation$start[match(new_df$Probe, rownames(ProbeAnnotation))]
  TSS <- as.data.frame(getTSS(genome = genome))
  TSS <- TSS[which(TSS$external_gene_name %in% new_df$Gene), ]
  TSS <- TSS %>%
    group_by(external_gene_name) %>%
    arrange(transcription_start_site, .by_group = TRUE)
  new_df$TSS <- TSS$transcription_start_site[match(new_df$Gene, TSS$external_gene_name)]
  new_df$distTSS <- new_df$CpG_start - new_df$TSS
  FunctionalPairs <- merge(x = FunctionalPairs, y = new_df, by = c("Gene", "Probe"), all.x = TRUE)
  cols <- c("Gene", "Probe", "State", "distTSS", "Proportion of hypo (%)", "Proportion of hyper (%)", "Fold change of gene expression",
            "Comparators",  "Raw.p", "Adjusted.p")
  FunctionalPairs <- FunctionalPairs[, cols]
  FunctionalPairs <- FunctionalPairs %>%
    group_by(Gene) %>%
    arrange(distTSS, .by_group = TRUE)
  # new_df <- FunctionalPairs[, c("Gene", "Probe")]
  # new_df$probe_start <- ProbeAnnotation$start[match(new_df$Probe, rownames(ProbeAnnotation))]
  # TSS <- TSS[which(TSS$external_gene_name %in% new_df$Gene), ]
  # TSS <- TSS[, c("external_gene_name", "transcription_start_site")]
  # TSS <- distinct(TSS)
  # colnames(TSS) <- c("Gene", "tss")
  # new_df <- merge(x = new_df, y = TSS, by = c("Gene"), all.x = TRUE)
  # new_df$distance <- new_df$probe_start - new_df$tss
  # new_df %>%
  #   group_by(Gene, Probe) %>%
  #   summarise(disTSS = paste0(disTSS, collapse = ";"))
  return(FunctionalPairs)
}

# library("GenomicFeatures")
# library("TxDb.Hsapiens.UCSC.hg38.knownGene")
# genome <- TxDb.Hsapiens.UCSC.hg38.knownGene
# # the plec gene
# plec_gene = genes(genome)[which(genes(genome)$gene_id == 5339),]
# # get the exons with the gene coordinates
# plec_exons = subsetByOverlaps(exons(genome), plec_gene)
# # same for transcripts
# plec_t = subsetByOverlaps(transcripts(genome), plec_gene)











