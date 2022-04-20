# Fix for NOTEs from R CMD check (no visible binding for global variable)

#' @importFrom grDevices dev.off pdf png
NULL

#' @importFrom graphics lines par plot title
NULL

#' @importFrom stats anova aov as.dist cor cutree dbeta density dnorm hclust lm optim prcomp qqline qqnorm qqplot quantile rgamma t.test var wilcox.test p.adjust
NULL

#' @importFrom utils download.file read.csv tail untar write.table write.csv data setTxtProgressBar txtProgressBar
NULL

#' @importFrom methods as is
NULL

#' @importFrom foreach foreach %dopar%
NULL

#' @importFrom parallel makeCluster stopCluster
NULL

#' The EpiMix function
#' @description EpiMix uses a model-based approach to identify functional changes DNA methylation that affect gene expression.
#' @param methylation.data Matrix of the DNA methylation data with CpGs in rows and samples in columns.
#' @param gene.expression.data  Matrix of the gene expression data with genes in rows and samples in columns.
#' @param mode Character string indicating the analytic mode to model DNA methylation.
#' Should be one of the followings: "Regular", "Enhancer", "miRNA" or "lncRNA". Default: "Regular". See details for more information.
#' @param sample.info Dataframe that maps each sample to a study group.
#' Should contain two columns: the first column (named "primary") indicates the sample names, and the second column (named "sample.type") indicating which study group each sample belongs to (e.g.,“Cancer” vs. “Normal”,  “Experiment” vs. “Control”). Sample names in the "primary" column must coincide with the column names of the methylation.data.
#' @param group.1 Character vector indicating the name(s) for the experiment group.
#' @param group.2 Character vector indicating the names(s) for the control group.
#' @param promoters Logic indicating whether to focus the analysis on CpGs associated with promoters (2000 bp upstream and 1000 bp downstream of the transcription start site). This parameter is only used for the Regular mode.
#' @param met.platform Character string indicating the microarray type for collecting the DNA methylation data. The value should be either "HM27", "HM450" or "EPIC". Default: "HM450"
#' @param genome Character string indicating the genome build version to be used for CpG annotation. Should be either "hg19" or "hg38". Default: "hg38".
#' @param listOfGenes Character vector used for filtering the genes to be evaluated.
#' @param raw.pvalue.threshold Numeric value indicating the threshold of the raw P value for selecting the functional CpG-gene pairs. Default: 0.05.
#' @param adjusted.pvalue.threshold Numeric value indicating the threshold of the adjusted P value for selecting the function CpG-gene pairs. Default: 0.05.
#' @param numFlankingGenes Numeric value indicating the number of flanking genes whose expression is to be evaluated for selecting the functional enhancers. Default: 20.
#' @param roadmap.epigenome.groups (parameter used for the "Enhancer" mode) Character vector indicating the tissue group(s) to be used for selecting the enhancers. See details for more information. Default: NULL.
#' @param roadmap.epigenome.ids (parameter used for the "Enhancer" mode) Character vector indicating the epigenome ID(s) to be used for selecting the enhancers. See details for more information. Default: NULL.
#' @param chromatin.states (parameter used for the "Enhancer" mode) Character vector indicating the chromatin states to be used for selecting the enhancers. To get the available chromatin states, please run the list.chromatin.states() function. Default: c("EnhA1", "EnhA2", "EnhG1", "EnhG2").
#' @param NoNormalMode Logical indicating if the methylation states found in the experiment group should be compared to the control group. Default: FALSE.
#' @param cores Number of CPU cores to be used for computation. Default: 1.
#' @param OutputRoot File path to store the EpiMix result object. Default: "." (current directory)
#' @param MixtureModelResults Pre-computed EpiMix results, used for generating functional probe-gene pair matrix. Default: NULL
#' @import EpiMix.data
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
#' Since enhancers are cell-type or tissue-type specific, EpiMix needs to know the reference tissues or cell types in order to select the proper enhancers.
#' EpiMix identifies enhancers from the RoadmapEpigenomic project (Nature, PMID: 25693563), which enhancers were identified by ChromHMM in over 100 tissue and cell types.
#' Available epigenome groups (a group of relevant cell types) or epigenome ids (individual cell types) can be obtained from the original publication (Nature, PMID: 25693563, figure 2).
#' They can also be retrieved from the list.epigenomes() function. If both roadmap.epigenome.groups and roadmap.epigenome.ids are specified, EpiMix will select all the epigenomes from the combination of the inputs.
#'
#' @export
#' @examples
#' \dontrun{
#' # Regular mode
#' library(EpiMix)
#' data(MET.data)
#' data(mRNA.data)
#' data(sample.info)
#' EpiMixResults_Regular <- EpiMix(methylation.data = MET.data,
#'                                 gene.expression.data = mRNA.data,
#'                                 sample.info = sample.info,
#'                                 group.1 = "Cancer",
#'                                 group.2 = "Normal",
#'                                 met.platform = "HM450")
#'}
#' \dontrun{
#' # Enhancer mode
#' library(EpiMix)
#' data(MET.data)
#' data(mRNA.data)
#' data(sample.info)
#' mode = "Enhancer"
#' roadmap.epigenome.ids = "E096"
#' EpiMixResults_Enhancer <- EpiMix(methylation.data = MET.data,
#'                                 gene.expression.data = mRNA.data,
#'                                 sample.info = sample.info,
#'                                 mode = mode,
#'                                 group.1 = "Cancer",
#'                                 group.2 = "Normal",
#'                                 roadmap.epigenome.ids = roadmap.epigenome.ids,
#'                                 met.platform = "HM450")
#'}

EpiMix <- function(methylation.data,
                   gene.expression.data,
                   mode = "Regular",
                   sample.info,
                   group.1,
                   group.2,
                   promoters = FALSE,
                   met.platform = "HM450",
                   genome = "hg38",
                   listOfGenes = NULL,
                   raw.pvalue.threshold = 0.05,
                   adjusted.pvalue.threshold = 0.05,
                   numFlankingGenes = 20,
                   roadmap.epigenome.groups = NULL,
                   roadmap.epigenome.ids = NULL,
                   chromatin.states = c("EnhA1", "EnhA2", "EnhG1", "EnhG2"),
                   NoNormalMode = FALSE,
                   cores = 1,
                   MixtureModelResults = NULL,
                   OutputRoot = "."
) {

  ### Process 1: Check user input
  if (missing(methylation.data)) stop("Need to provide DNA methylation matrix\n")
  #if (missing(gene.expression.data)) stop("Need to provide gene expression matrix\n")
  stopifnot(
    is.matrix(methylation.data),
    is.null(gene.expression.data) | is.matrix(gene.expression.data) | is.data.frame(gene.expression.data),
    is(mode, "character"),
    class(sample.info) %in% c("data.frame", "matrix"),
    class(group.1) %in% c("character", "NULL"),
    class(group.2) %in% c("character", "NULL"),
    is(promoters, "logical"),
    is(met.platform, "character"),
    is(genome, "character"),
    class(listOfGenes) %in% c("character", "NULL"),
    is(raw.pvalue.threshold, "numeric"),
    is(adjusted.pvalue.threshold, "numeric"),
    is(numFlankingGenes, "numeric"),
    class(roadmap.epigenome.groups) %in% c("character", "NULL"),
    class(roadmap.epigenome.ids) %in% c("character", "NULL"),
    class(chromatin.states) %in% c("character", "NULL"),
    is(NoNormalMode, "logical"),
    is(cores, "numeric"),
    is(OutputRoot, "character")
  )

  if(nrow(methylation.data) == 0){
    stop("methylation.data matrix is empty, please check the methyaltion.data matrix\n")
  }
  if(!is.null(gene.expression.data) & length(gene.expression.data) == 0){
    stop("gene.expression.data matrix is empty, please check the gene.expression.data matrix\n")
  }
  if (!mode %in% c("Regular", "Enhancer", "miRNA", "lncRNA")) stop ("'mode' must be one of the followings: 'Regular', 'Enhancer', 'miRNA', 'lncRNA'")
  if (is.null(sample.info)){
    warning("'sample.info' is not provided. EpiMix will treat all samples as one group (i.e., no comparison will be made between the experiment and the control group) !!!\n")
  }
  if (!is.null(sample.info) & (is.null(group.1) | is.null(group.2))){
    warning("Only zero or one group is sepecfied. EpiMix will treat all samples as one group (i.e., no comparison will be made between the experiment and the control group) !!!\n")
  }

  if (length(met.platform)!=1 & !toupper(met.platform) %in% c("EPIC", "HM450", "HM27")) stop("'met.platform' must be either 'EPIC', 'HM450', 'HM27'\n")
  if (length(genome)!=1 & tolower(genome) %in% c("hg19", "hg38")) stop("'genome' must to be either 'hg19' or 'hg38'\n")

  met.platform = toupper(met.platform)
  genome = tolower(genome)

  if(!is.null(OutputRoot) & length(OutputRoot) > 0){
    dir.create(OutputRoot,showWarnings=FALSE)
  }

  ### Process 2: Set up the parallel back-end
  if(cores>1){
    #unregister()
    cat("Registering sockets on multiple CPU cores...\n")
    cl <- parallel :: makeCluster(cores)
    registerDoSNOW(cl)
  }

  ### Process 3: filter CpG probes and samples
  if(!is.null(sample.info) & !is.null(group.1) & !is.null(group.2)){
    target.samples = sample.info$primary[which(sample.info$sample.type %in% c(group.1, group.2))]
    sample.info = sample.info[sample.info$primary %in% target.samples,]
    overlapSamples.met = intersect(colnames(methylation.data), target.samples)
    methylation.data = methylation.data[,overlapSamples.met, drop = FALSE]
    if(!is.null(gene.expression.data)){
      overlapSamples.exp = intersect(colnames(gene.expression.data), target.samples)
      gene.expression.data = gene.expression.data[,overlapSamples.exp, drop = FALSE]
    }
  }

  ### Process 5: run MethylMix
  #--------------------------------------------Regular mode------------------------------------------------------
  if(mode == "Regular"){
    cat("Running", mode, "mode...\n")

    ### Step 1: filter CpGs based on user-specified conditions
    ProbeAnnotation <- filterProbes(mode = mode,
                                     gene.expression.data = gene.expression.data,
                                     listOfGenes = listOfGenes,
                                     promoters = promoters,
                                     met.platform = met.platform,
                                     genome = genome)
    overlapProbes = unique(intersect(ProbeAnnotation$probe, rownames(methylation.data)))
    methylation.data = methylation.data[overlapProbes,,drop = FALSE]

    ### Step 2: modeling the gene expression using the methylation data (beta values scale) to select functional probes
    FunctionalProbes = NULL
    if(is.null(MixtureModelResults) & !is.null(gene.expression.data)){
      FunctionalProbes =  EpiMix_ModelGeneExpression(methylation.data, gene.expression.data, ProbeAnnotation, cores = cores, filter = TRUE)
      if(length(FunctionalProbes) == 0){
        stop("No transcriptionally predicitve CpGs were found.")
      }
      if(OutputRoot != ""){
        saveRDS(FunctionalProbes, paste0(OutputRoot, "/", "FunctionalProbes_", mode, ".rds"))
      }
    }

    ### Step 3: split methylation data into group.1 and group.2
    MET_Experiment <-  MET_Control <- NULL
    methylation.data <- split.met.data(methylation.data, sample.info, group.1, group.2)
    MET_Experiment <- methylation.data$MET_Experiment
    MET_Control <- methylation.data$MET_Control
    rm(methylation.data); gc()

    ### Step 4: modeling the methylation data as a mixture of beta distributions
    if(is.null(FunctionalProbes)){
      FunctionalProbes = rownames(MET_Experiment)
    }
    if(!is.null(MixtureModelResults)){
      MethylMixResults = readRDS(MixtureModelResults)
    }else{
      MethylMixResults <- MethylMix_MixtureModel(MET_Experiment, MET_Control, FunctionalProbes, NoNormalMode)
    }

    if(!is.null(MethylMixResults$MethylationDrivers)){
      cat("Found", length(MethylMixResults$MethylationDrivers), "differentially methylated CpGs\n")
    }
    ### Step 5: optionally, write the intermediate output to file (test purpose only!!!)
    if (OutputRoot != "") {
      saveRDS(MethylMixResults, file = paste0(OutputRoot, "/", "EpiMix_Results_",mode,".rds"))
    }

    ### Step 6: select functional CpGs and calculate prevalence and fold change
    MET_matrix <- MethylMixResults$MethylationStates
    prev.data <- get.prevalence(MethylMixResults)
    if(is.null(gene.expression.data)){
      # No gene expression data, just report the gene names and prevalence
      prev.data <- addGeneNames(prev.data, ProbeAnnotation)
      MethylMixResults$FunctionalPairs = prev.data %>% dplyr :: select(.data$Gene, .data$Probe, .data$'Prevalence of hypo (%)', data$'Prevalence of hyper (%)')
      return(MethylMixResults)
    }

    ### Step 7:  modeling the gene expression and select the functional enhancer probes
    cat("Identifying functional CpG-gene pairs...\n")
    FunctionalPairs <- generateFunctionalPairs(MET_matrix,
                                               MET_Control,
                                               gene.expression.data,
                                               ProbeAnnotation,
                                               raw.pvalue.threshold,
                                               adjusted.pvalue.threshold,
                                               cores)
    if(is.null(FunctionalPairs)){
      cat("Not enough differentially methylated genes or not sufficient gene expression data, returning EpiMix results...\n")
      prev.data <- addGeneNames(prev.data, ProbeAnnotation)
      MethylMixResults$FunctionalPairs = prev.data %>% dplyr :: select(.data$Probe, .data$Gene, .data$'Prevalence of hypo (%)', data$'Prevalence of hyper (%)')
      return(MethylMixResults)
    }

    # Add in the prevalence information
    FunctionalPairs <- merge(x = prev.data, y = FunctionalPairs, by = "Probe")
    col.order <- c("Gene", "Probe", "State", "Prevalence of hypo (%)", "Prevalence of hyper (%)", "Fold change of gene expression",
                   "Comparators", "Raw.p", "Adjusted.p")
    FunctionalPairs <- FunctionalPairs[, col.order]
    FunctionalPairs <- FunctionalPairs[order(FunctionalPairs$Gene), ]
    rownames(FunctionalPairs) = c()
    MethylMixResults$FunctionalPairs = FunctionalPairs
  }

  #--------------------------------------------miRNA mode------------------------------------------------------
  if(mode == "miRNA"){
    cat("Running", mode, "mode...\n")
    cat("Please be mindful that the gene expression data are expected to be data obtained from microRNA-seq.\n")
    ### Step 1: filter CpGs based on user-specified conditions
    ProbeAnnotation <- filterProbes(mode = mode,
                                    gene.expression.data = gene.expression.data,
                                    listOfGenes = listOfGenes,
                                    promoters = promoters,
                                    met.platform = met.platform,
                                    genome = genome)
    overlapProbes = unique(intersect(ProbeAnnotation$probe, rownames(methylation.data)))
    methylation.data = methylation.data[overlapProbes,,drop = FALSE]

    ### Step 2: split methylation data into group.1 and group.2
    MET_Experiment <-  MET_Control <- NULL
    methylation.data <- split.met.data(methylation.data, sample.info, group.1, group.2)
    MET_Experiment <- methylation.data$MET_Experiment
    MET_Control <- methylation.data$MET_Control
    rm(methylation.data); gc()

    ### Step 3: modeling the methylation data as a mixture of beta distributions of miRNAs
    if(!is.null(MixtureModelResults)){
      MethylMixResults = readRDS(MixtureModelResults)
    }else{
      FunctionalGenes = rownames(MET_Experiment)
      MethylMixResults <- MethylMix_MixtureModel(MET_Experiment, MET_Control, FunctionalGenes, NoNormalMode = FALSE)
      if (OutputRoot != "") { # Save the intermediate results for test purpose
        saveRDS(MethylMixResults, file = paste0(OutputRoot, "/", "EpiMix_Results_",mode,".rds"))
      }
    }
    cat("Found", length(MethylMixResults$MethylationDrivers), "differentially methylated probes\n")

    ### Step 4: optionally, write the intermediate output to file
    if (OutputRoot != "") {
      saveRDS(MethylMixResults, file = paste0(OutputRoot, "/", "EpiMix_Results_",mode,".rds"))
    }

    ### Step 5: select functional CpGs and calculate prevalence and fold change
    MET_matrix <- MethylMixResults$MethylationStates
    prev.data <- get.prevalence(MethylMixResults)
    if(is.null(gene.expression.data)){
      # No gene expression data, just report the gene names and prevalence
      prev.data <- addGeneNames(prev.data, ProbeAnnotation)
      MethylMixResults$FunctionalPairs = prev.data %>% dplyr :: select(.data$Gene, .data$Probe, .data$'Prevalence of hypo (%)', data$'Prevalence of hyper (%)')
      return(MethylMixResults)
    }

    ### Step 6: identify transcriptionally predictive probes
    cat("Identifying functional CpG-gene pairs...\n")
    FunctionalPairs <- generateFunctionalPairs(MET_matrix,
                                               MET_Control,
                                               gene.expression.data,
                                               ProbeAnnotation,
                                               raw.pvalue.threshold,
                                               adjusted.pvalue.threshold,
                                               cores)
    if(is.null(FunctionalPairs)){
      cat("Not enough differentially methylated genes or not sufficient gene expression data, returning EpiMix results...\n")
      prev.data <- addGeneNames(prev.data, ProbeAnnotation)
      MethylMixResults$FunctionalPairs = prev.data %>% dplyr :: select(.data$Probe, .data$Gene, .data$'Prevalence of hypo (%)', data$'Prevalence of hyper (%)')
      return(MethylMixResults)
    }

    # Add in the prevalence information
    FunctionalPairs <- merge(x = prev.data, y = FunctionalPairs, by = "Probe")
    col.order <- c("Gene", "Probe", "State", "Prevalence of hypo (%)", "Prevalence of hyper (%)", "Fold change of gene expression",
                   "Comparators", "Raw.p", "Adjusted.p")
    FunctionalPairs <- FunctionalPairs[, col.order]
    FunctionalPairs <- FunctionalPairs[order(FunctionalPairs$Gene), ]
    rownames(FunctionalPairs) = c()
    MethylMixResults$FunctionalPairs = FunctionalPairs
  }

  #--------------------------------------------LncRNA Mode--------------------------------------------------------------------------------------------
  if(mode == "lncRNA"){
    cat("Running", mode, "mode...\n")
    cat("We recommend using the kallisto-sleuth pipline to process the RNA-seq data in order to detect more lncRNAs.
        Please see our publication for details: PMID: 31808800\n")

    ### Step 1: filter CpGs based on user-specified conditions
    ProbeAnnotation <- filterProbes(mode = mode,
                                    gene.expression.data = gene.expression.data,
                                    listOfGenes = listOfGenes,
                                    promoters = promoters,
                                    met.platform = met.platform,
                                    genome = genome)
    overlapProbes = unique(intersect(ProbeAnnotation$probe, rownames(methylation.data)))
    methylation.data = methylation.data[overlapProbes,,drop = FALSE]


    ### Step 2: split methylation data into group.1 and group.2
    MET_Experiment <-  MET_Control <- NULL
    methylation.data <- split.met.data(methylation.data, sample.info, group.1, group.2)
    MET_Experiment <- methylation.data$MET_Experiment
    MET_Control <- methylation.data$MET_Control
    rm(methylation.data); gc()

    ### Step 4: modeling the methylation data as a mixture of beta distributions of lncRNA probes
    if(!is.null(MixtureModelResults)){
      MethylMixResults = readRDS(MixtureModelResults)
    }else{
      FunctionalGenes = rownames(MET_Experiment)
      MethylMixResults <- MethylMix_MixtureModel(MET_Experiment, MET_Control, FunctionalGenes, NoNormalMode = FALSE)
      if (OutputRoot != "") { # Save the intermedidate results for test purpose
        saveRDS(MethylMixResults, file = paste0(OutputRoot, "/", "EpiMix_Results_",mode,".rds"))
      }
    }
    cat("Found", length(MethylMixResults$MethylationDrivers), "differentially methylated probes\n")

    ### Step 5: select functional CpGs and calculate prevalence and fold change
    MET_matrix <- MethylMixResults$MethylationStates
    prev.data <- get.prevalence(MethylMixResults)
    if(is.null(gene.expression.data)){
      # No gene expression data, just report the gene names and prevalence
      prev.data <- addGeneNames(prev.data, ProbeAnnotation)
      MethylMixResults$FunctionalPairs = prev.data %>% dplyr :: select(.data$Gene, .data$Probe, .data$'Prevalence of hypo (%)', data$'Prevalence of hyper (%)')
      return(MethylMixResults)
    }

    ### Step 6: identify transcriptionally predictive probes
    cat("Identifying functional CpG-gene pairs...\n")
    FunctionalPairs <- generateFunctionalPairs(MET_matrix,
                                               MET_Control,
                                               gene.expression.data,
                                               ProbeAnnotation,
                                               raw.pvalue.threshold,
                                               adjusted.pvalue.threshold,
                                               cores)
    if(is.null(FunctionalPairs)){
      cat("Not enough differentially methylated genes or not sufficient gene expression data, returning EpiMix results...\n")
      prev.data <- addGeneNames(prev.data, ProbeAnnotation)
      MethylMixResults$FunctionalPairs = prev.data %>% dplyr :: select(.data$Probe, .data$Gene, .data$'Prevalence of hypo (%)', data$'Prevalence of hyper (%)')
      return(MethylMixResults)
    }

    # Add in the prevalence information
    FunctionalPairs <- merge(x = prev.data, y = FunctionalPairs, by = "Probe")
    col.order <- c("Gene", "Probe", "State", "Prevalence of hypo (%)", "Prevalence of hyper (%)", "Fold change of gene expression",
                   "Comparators", "Raw.p", "Adjusted.p")
    FunctionalPairs <- FunctionalPairs[, col.order]
    FunctionalPairs <- FunctionalPairs[order(FunctionalPairs$Gene), ]
    rownames(FunctionalPairs) = c()
    MethylMixResults$FunctionalPairs = FunctionalPairs
  }

  #--------------------------------------------Enhancer Mode--------------------------------------------------------------------------------------------
  if(mode == "Enhancer"){
    cat("Running", mode, "mode...\n")
    ### Step 1: filter enhancer probes
    cat("Fetching probe annotation...\n")
    suppressMessages({
      ProbeAnnotation = EpiMix_getInfiniumAnnotation(plat = met.platform, genome = genome)
    })
    selectedEpigenomes = validEpigenomes(roadmap.epigenome.groups, roadmap.epigenome.ids)
    if(length(selectedEpigenomes) == 0){
      stop("Must input valid Roadmap Epigenome IDs or epigenome groups")
    }

    ### Step 2: split methylation data into group.1 and group.2
    MET_Experiment <-  MET_Control <- NULL
    methylation.data <- split.met.data(methylation.data, sample.info, group.1, group.2)
    MET_Experiment <- methylation.data$MET_Experiment
    MET_Control <- methylation.data$MET_Control
    #rm(methylation.data); gc()

    ### Step 3: modeling the methylation data as a mixture of beta distributions of enhancer probes
    if(!is.null(MixtureModelResults)){
      MethylMixResults = readRDS(MixtureModelResults)
    }else{
      cat("Fetching enhancer CpGs from Roadmap Epigenomics...\n")
      RoadMap.enhancer.probes <- getRoadMapEnhancerProbes(met.platform = met.platform,
                                                          genome = genome,
                                                          functional.regions = chromatin.states,
                                                          listOfEpigenomes = selectedEpigenomes,
                                                          ProbeAnnotation = ProbeAnnotation)
      # distal enhancer probes = intersect(distal probes, enhancer probes)
      distal.probes <- names(get.feature.probe(met.platform = met.platform, genome = genome))
      enhancer.probes = intersect(distal.probes, RoadMap.enhancer.probes$probeID)
      presentEnhancerProbes = intersect(rownames(MET_Experiment), enhancer.probes)
      cat("Found", length(presentEnhancerProbes), "CpGs associated with distal enhancers in the methylation dataset\n")
      MET_Experiment <- MET_Experiment[presentEnhancerProbes, ,drop = FALSE]
      MET_Control <- MET_Control[presentEnhancerProbes, ,drop = FALSE]

      FunctionalGenes = rownames(MET_Experiment)
      MethylMixResults <- MethylMix_MixtureModel(MET_Experiment, MET_Control, FunctionalGenes, NoNormalMode = FALSE)
      if (OutputRoot != "") { # Save the intermediate MethylMix results
        saveRDS(MethylMixResults, file = paste0(OutputRoot, "/", "EpiMix_Results_",mode,".rds"))
      }
    }
    cat("Found", length(MethylMixResults$MethylationDrivers), "differentially methylated CpGs\n")

    # Calculate the prevalence for differential DNAme
    MET_matrix <- MethylMixResults$MethylationStates
    prev.data <- get.prevalence(MethylMixResults)
    if(is.null(gene.expression.data)){
      # No gene expression data, just report the gene names and prevalence
      prev.data['Gene'] <- ProbeAnnotation$gene_HGNC[which(names(ProbeAnnotation) %in% rownames(MET_matrix))]
      MethylMixResults$FunctionalPairs = prev.data %>% dplyr :: select(.data$Gene, .data$Probe, .data$'Prevalence of hypo (%)', data$'Prevalence of hyper (%)')
      return(MethylMixResults)
    }

    ### Step 4:  modeling the gene expression and select the functional enhancer probes
    cat("Modeling the gene expression for enhancers...\n")
    MET_matrix<- filterMethMatrix(MET_matrix = MET_matrix, MET_Control = MET_Control, gene.expression.data = gene.expression.data)

    if(length(MET_matrix) == 0){
      cat("Not enough differentially methylated genes or not sufficient gene expression data, returning EpiMix results...\n")
      prev.data['Gene'] <-  ProbeAnnotation$gene_HGNC[which(names(ProbeAnnotation) %in% rownames(MET_matrix))]
      MethylMixResults$FunctionalPairs = prev.data %>% dplyr :: select(.data$Probe, .data$Gene, .data$'Prevalence of hypo (%)', data$'Prevalence of hyper (%)')
      return(MethylMixResults)
    }

    # Get nearby genes for the differentially methylated CpGs
    DM.probes = rownames(MET_matrix)
    geneAnnot <- getTSS(genome = genome) #ELMER function to retrieve a GRange object that contains coordinates of promoters for human genome.
    DM.probes.annotation = ProbeAnnotation[which(names(ProbeAnnotation) %in% DM.probes), ,drop = FALSE]
    NearbyGenes <- GetNearGenes(geneAnnot = geneAnnot,
                                TRange = DM.probes.annotation,
                                numFlankingGenes = numFlankingGenes)

    # Find functional CpGs

    cat("Looking for differentially methylated enhancers associated with gene expression\n")
    iterations = length(DM.probes)
    pb <- utils :: txtProgressBar(max = iterations, style = 3)

    FunctionalPairs = data.frame()
    if(cores == "" | cores == 1){
      for(i in 1:iterations){
        target.probe = DM.probes[i]
        target.genes =  intersect(NearbyGenes$Symbol[which(NearbyGenes$ID == target.probe)],rownames(gene.expression.data))
        if (length(target.genes)==0) next()
        pairs = getFunctionalGenes(target.probe = target.probe,
                                   target.genes = target.genes,
                                   MET_matrix = MET_matrix,
                                   gene.expression.data = gene.expression.data,
                                   ProbeAnnotation = ProbeAnnotation,
                                   raw.pvalue.threshold = raw.pvalue.threshold,
                                   adjusted.pvalue.threshold = adjusted.pvalue.threshold)
        FunctionalPairs = rbind(FunctionalPairs, pairs)
        utils :: setTxtProgressBar(pb,i)
      }
    }
    else{
      progress <- function(n) utils :: setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      FunctionalPairs <- foreach :: foreach(i = 1:iterations, .combine = rbind, .options.snow= opts, .verbose = FALSE)  %dopar% {
        target.probe = DM.probes[i]
        target.genes =  NearbyGenes$Symbol[which(NearbyGenes$ID == target.probe)]
        target.genes = intersect(target.genes, rownames(gene.expression.data))
        if(length(target.genes)>0) {
          getFunctionalGenes(target.probe = target.probe,
                             target.genes = target.genes,
                             MET_matrix = MET_matrix,
                             gene.expression.data = gene.expression.data,
                             ProbeAnnotation = ProbeAnnotation,
                             raw.pvalue.threshold = raw.pvalue.threshold,
                             adjusted.pvalue.threshold = adjusted.pvalue.threshold)
        }
      }
    }
    close(pb)
    # Add in the prevalence information
    FunctionalPairs <- merge(x = prev.data, y = FunctionalPairs, by = "Probe")
    rownames(FunctionalPairs) = c()
    col.order <- c("Gene", "Probe", "State", "Prevalence of hypo (%)", "Prevalence of hyper (%)", "Fold change of gene expression",
                   "Comparators", "Raw.p", "Adjusted.p")
    FunctionalPairs <- FunctionalPairs[, col.order]
    FunctionalPairs <- FunctionalPairs[order(FunctionalPairs$Gene), ]
    MethylMixResults$FunctionalPairs = FunctionalPairs
  }

  if(!is.null(MethylMixResults$FunctionalPairs)){
    cat("Found", nrow(MethylMixResults$FunctionalPairs), "functional probe-gene pairs.\n")
    MethylMixResults$FunctionalPairs =  MethylMixResults$FunctionalPairs[order(MethylMixResults$FunctionalPairs$Gene), ]
  }

  # Save the sample names for the experiment and the control groups (used for the EpiMix_plotModel function)
  MethylMixResults$group.1 = colnames(MET_Experiment)
  MethylMixResults$group.2 = colnames(MET_Control)

  # Save the output
  if (!is.null(OutputRoot) & OutputRoot != "") {
    cat("Saving the EpiMix results to the output directory...\n")
    saveRDS(MethylMixResults, file = paste0(OutputRoot, "/", "EpiMix_Results_",mode,".rds"))
    utils :: write.csv(MethylMixResults$FunctionalPairs, paste0(OutputRoot, "/", "FunctionalPairs_", mode,".csv"), row.names = FALSE )
  }

  # Clean up the environment
  if(cores > 1){
    parallel :: stopCluster(cl)
  }
  return(MethylMixResults)
}



