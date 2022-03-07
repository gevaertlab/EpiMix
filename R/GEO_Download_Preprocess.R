###################################################################################################################
#                                  Functions to download methylation and gene expression data from GEO
###################################################################################################################
#'@importFrom Biobase exprs pData phenoData
NULL

#' @importFrom GEOquery getGEO
NULL

#' @importFrom biomaRt useDataset getBM
NULL


#' The GEO_Download_DNAmethylation function
#' @description Download the methylation data and the associated sample phenotypic data from the GEO database.
#' @param AccessionID character string indicating GEO accession number. Currently support the GEO series (GSE) data type.
#' @param targetDirectory character string indicting the file path to save the data. Default: '.' (current directory).
#' @param DownloadData logical indicating whether the actual data should be downloaded (Default: TRUE). If False, the desired directory where the downloaded data should have been saved is returned.
#' @return a list with two elements. The first element ("$MethylationData") indicating the file path to the downloaded methylation data. The second element ("$PhenotypicData") indicating the file path to the sample phenotypic data.
#' @export
#' @keywords download

GEO_Download_DNAMethylation <- function(AccessionID,targetDirectory = '.', DownloadData = TRUE) {

  # set the environmental variable to deal with the long file name issue
  Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)
  dir_MET_data = paste0(targetDirectory,"/",AccessionID, "_MET_data.rds")
  dir_Pheno_data = paste0(targetDirectory,"/", AccessionID,"_Pheno_data.csv.gz")
  if(DownloadData){
    cat("Downloading methylation data from GEO", AccessionID, "\n")
    dir.create(targetDirectory,showWarnings=FALSE)
    all.data <- GEOquery :: getGEO(GEO = AccessionID, destdir = targetDirectory)
    expression.data <- Biobase :: exprs(all.data[[1]])
    if(nrow(expression.data) <=1){
      warning("The methylation data matrix has zero feature! Please check whether the actual data are only saved as a supplmentary file in the GEO database.\n In this case, you will have to manually download the data from the supplmentary files from GEO.\n")
    }
    sample.metadata <- Biobase :: pData(Biobase :: phenoData(all.data[[1]]))
    cat("Saving methylation data into the target dirctory...\n")
    saveRDS(expression.data, dir_MET_data)
    cat("Saving sample phenotypic data into the target dirctory...\n")
    utils :: write.csv(sample.metadata, file = gzfile(dir_Pheno_data))
    cat("Methylation data are saved at:",dir_MET_data, "\n")
    cat("Phenotypic data are saved at:", dir_Pheno_data, "\n")
  }
  METdirectories = list(MethylationData = dir_MET_data, PhenotypicData=dir_Pheno_data)
  return(METdirectories)
}


#' The GEO_Download_GeneExpression function
#' @description Download the gene expression data and the associated sample phenotypic data from the GEO database.
#' @param AccessionID character string indicating the GEO accession number. Currently support the GEO series (GSE) data type.
#' @param targetDirectory character string indicting the file path to save the data. Default: '.' (current directory)
#' @param DownloadData logical indicating whether the actual data should be downloaded (Default: TRUE). If False, the desired directory where the downloaded data should have been saved is returned.
#' @return a list with two elements. The first element ("$GeneExpressionData") indicating the file path to the downloaded methylation data. The second element ("$PhenotypicData") indicating the file path to the sample phenotypic data.
#' @export
#' @keywords download

GEO_Download_GeneExpression <- function(AccessionID,targetDirectory = '.', DownloadData = TRUE) {

  # set the environmental variable to deal with the long file name issue
  Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)
  dir_expression_data = paste0(targetDirectory,"/",AccessionID, "_expression_data.rds")
  dir_Pheno_data = paste0(targetDirectory,"/", AccessionID,"_Pheno_data.csv.gz")

  if(DownloadData){
    cat("Downloading gene expression data from GEO", AccessionID, "\n")
    dir.create(targetDirectory,showWarnings=FALSE)
    all.data <- GEOquery :: getGEO(GEO = AccessionID, destdir = targetDirectory)
    expression.data <- Biobase :: exprs(all.data[[1]])
    if(nrow(expression.data) <=1){
      warning("The gene expression data matrix has zero feature! Please check whether the actual data are only saved as a supplmentary file in the GEO database.\n In this case, you will have to manually download the data from the supplmentary files from GEO.\n")
    }
    sample.metadata <- Biobase :: pData(Biobase :: phenoData(all.data[[1]]))

    cat("Saving gene expression data into the target dirctory...\n")
    saveRDS(expression.data, dir_expression_data)
    cat("Saving sample phenotypic data into the target dirctory...\n")
    utils :: write.csv(sample.metadata, file = gzfile(dir_Pheno_data))
    cat("Gene expression data are saved at:", dir_expression_data, "\n")
    cat("Phenotypic data are saved at:", dir_Pheno_data, "\n")
  }
  GEdirectories = list(GeneExpressionData = dir_expression_data, PhenotypicData=dir_Pheno_data)
  return(GEdirectories)
}

#' the GEO_getSampleMap function
#' @description auxiliary function to generate a sample map for DNA methylation data and gene expression data
#' @param METdirectories list of the file paths to the downloaded DNA methylation datasets, which can be the output from the GEO_Download_DNAMethylation function.
#' @param GEdirectories list of the file paths to the downloaded gene expression datasets, which can be the output from the GEO_Download_GeneExpression function.
#' @param targetDirectory file path to save the output. Default: '.' (current directory)
#' @return dataframe with three columns: $assay (character string indicating the type of the experiment, can be either "DNA methylation" or "Gene expression"), $primary(character string indicating the actual sample names), $colnames (character string indicating the actual column names for each samples in DNA methylation data and gene expression data)
#' @export

 GEO_getSampleMap <- function(METdirectories, GEdirectories, targetDirectory='.'){

  pData.met = METdirectories[[2]]
  pData.exp = GEdirectories[[2]]
  pData.met = read.csv(pData.met, row.names = 1)
  pData.exp = read.csv(pData.exp, row.names = 1)
  sample_met = data.frame(assay = rep("DNA methylation", nrow(pData.met)), primary = pData.met$title, colnames = pData.met$geo_accession)
  cat("Found", nrow(pData.met), "samples in DNA methyatlion data.\n")
  sample_exp = data.frame(assay = rep("Gene expression", nrow(pData.exp)), primary = pData.exp$title, colnames = pData.exp$geo_accession)
  cat("Found", nrow(pData.exp), "samples in gene expression data.\n")
  sample.map = rbind(sample_met, sample_exp)
  if(targetDirectory!='' | !is.null(targetDirectory)) {
    cat("Saving the sample map to files...\n")
    utils :: write.csv(sample.map, paste0(targetDirectory,"/", "sample.map.csv"), row.names = FALSE)
    cat("Sample map has been saved to:", paste0(targetDirectory,"/", "sample.map.csv\n"))
  }
  overlapSamples = intersect(sample.map$primary[which(sample.map$assay == "DNA methylation")], sample.map$primary[which(sample.map$assay == "Gene expression")])
  if(length(overlapSamples) == 0){
    warning("No overlap samples were found between the DNA methyation data and the gene expression data !!!\nPlease manually check the actual sample names ('primary' column) in the output.\n")
  }
  return(sample.map)
}


#' The GEO_GetSampleInfo function
#' @description auxiliary function to generate a sample information dataframe that indicates which study group each sample belongs to.
#' @param METdirectories  list of the file paths to the downloaded DNA methylation data, which can be the output from the GEO_Download_DNAMethylation function.
#' @param group.column character string indicating the column in the phenotypic data that defines the study group of each sample. The values in this column will be used to split the experiment and the control group.
#' @param targetDirectory file path to save the output. Default: '.' (current directory)
#' @return a dataframe with two columns: a "primary" column indicating the actual sample names, a "sample.type" column indicating the study group for each sample.
#' @export

GEO_GetSampleInfo <- function(METdirectories, group.column, targetDirectory = '.'){
  pData = METdirectories[[2]]
  pData = read.csv(pData, row.names = 1)
  if(!group.column %in% colnames(pData)){
    message("Can not find group.column in the columns of the phenotypicData (pData) !!!\n")
    cat("Avalible columns include:\n", paste0(colnames(pData), ","), "\n")
    stop("Please re-specifiy a correct group.column name.\n")
  }
  sample.info <- data.frame(primary = pData$title, sample.type = pData[,group.column])
  if(targetDirectory!='' | !is.null(targetDirectory)){
    write.csv(sample.info, paste0(targetDirectory, "/", "sample.info.csv"), row.names = FALSE)
  }
  return (sample.info)
}

#' The GEO_Preprocess_DNAMethylation function
#'
#' @description Preprocess DNA methylation data from the GEO database.
#' @details
#' The data preprocessing pipeline includes:
#' (1) eliminating samples and genes with too many NAs, imputing NAs.
#' (2) (optional) mapping the column names of the DNA methylation data to the actual sample names based on the information from "sample.map".
#' (3) (optional) removing CpG probes on the sex chromosomes or the user-defined chromosomes.
#' (4) (optional) doing Batch correction.
#' If both sample.info and group.1 and group.2 information are provided, the function will perform missing value estimation and batch correction on group.1 and group.2 separately. This will ensure that the true difference between group.1 and group.2 will not be obscured by missing value estimation and batch correction.
#' @param methylation.data matrix of DNA methylation data with CpG in rows and sample names in columns.
#' @param met.platform character string indicating the type of the Illumina Infinium BeadChip for collecting the methylation data. Should be either "HM450" or "EPIC". Default: "EPIC"
#' @param genome character string indicating the genome build version for retrieving the probe annotation. Should be either "hg19" or "hg38". Default: "hg38".
#' @param sample.info dataframe that maps each sample to a study group. Should contain two columns: the first column (named: "primary") indicating the sample names, and the second column (named: "sample.type") indicating which study group each sample belongs to (e.g., “Experiment” vs. “Control”,  “Cancer” vs. “Normal”). Sample names in the "primary" column must coincide with the column names of the methylation.data. Please see details for more information. Default: NULL.
#' @param group.1 character vector indicating the name(s) for the experiment group. The values must coincide with the values in the "sample.type" of the sample.info dataframe.Please see details for more information. Default: NULL.
#' @param group.2 character vector indicating the names(s) for the control group. The values must coincide with the values in the "sample.type" of the sample.info dataframe. Please see details for more information. Default: NULL.
#' @param sample.map dataframe for mapping the GEO accession ID (column names) to the actual sample names. Can be the output from the GEO_getSampleMap function. Default: NULL.
#' @param rm.chr character vector indicating the probes on which chromosomes to be removed. Default: "chrX", "chrY".
#' @param MissingValueThresholdGene threshold for missing values per gene. Genes with a percentage of NAs greater than this threshold are removed. Default: 0.3.
#' @param MissingValueThresholdSample threshold for missing values per sample. Samples with a percentage of NAs greater than this threshold are removed. Default: 0.1.
#' @param doBatchCorrection logical indicating whether to perform batch correction. If TRUE, the batch data need to be provided.
#' @param BatchData dataframe with batch information. Should contain two columns: the first column indicating the actual sample names, the second column indicating the batch. Users are expected to retrieve the batch information from the GEO on their own, but this can also be done using the GEO_getSampleInfo function with the "group.column" as the column indicating the batch for each sample. Defualt": NULL.
#' @param batch.correction.method character string indicating the method that will be used for batch correction. Should be either "Seurat" or "Combat". Default: "Seurat".
#' @param cores number of CPU cores to be used for batch effect correction. Defaut: 1.
#' @return DNA methylation data matrix with probes in rows and samples in columns.
#' @import doSNOW
#' @import doParallel
#' @export
#' @keywords preprocess

GEO_Preprocess_DNAMethylation <- function(methylation.data,
                                          met.platform = "EPIC",
                                          genome = "hg38",
                                          sample.info = NULL,
                                          group.1 = NULL,
                                          group.2 = NULL,
                                          sample.map = NULL,
                                          rm.chr = c("chrX", "chrY"),
                                          MissingValueThresholdGene = 0.2,
                                          MissingValueThresholdSample = 0.2,
                                          doBatchCorrection = FALSE,
                                          BatchData = NULL,
                                          batch.correction.method = "Seurat",
                                          cores = 1){

  # check data
  if(nrow(methylation.data) <=1 | ncol(methylation.data) <=1){
    stop("methylation.data is empty!\n Please check whether the actual data are saved as supplementary files in GEO.\n")
  }

  ### Step 1: convert column names to the actual patient names
  if(!is.null(sample.map)){
    cat("Mapping column names of the DNA methylation data to the actual sample names...\n")
    sample.map = sample.map[sample.map$assay == "DNA methylation", ]
    overlapSamples = intersect(sample.map$colnames, colnames(methylation.data))
    if(length(overlapSamples) > 0){
      methylation.data = methylation.data[,overlapSamples, drop = F]
      presentSamples = match(colnames(methylation.data), sample.map$colnames)
      sampleNames = sample.map$primary
      sampleNames = sampleNames[presentSamples]
      colnames(methylation.data) = sampleNames
    }else{
      warning("No overlap samples were found between the sample map and the column names of the methylation data! No mapping of column names was performed!")
    }
  }

  ### Step 2: Filter out SNP probes
  cat("\tFiltering out SNP probes...\n")
  GoodProbes= rownames(methylation.data)[startsWith(rownames(methylation.data), "cg")]
  NrProbesToRemove=length(rownames(methylation.data))-length(GoodProbes)
  methylation.data=methylation.data[GoodProbes,]
  cat("\tRemoved",NrProbesToRemove,"probes with SNPs.\n")

  ### Step 3:Filter out CpG probes in the user-specified chromosomes
  cat("\tFetching probe annotation for", met.platform, "\n")
  ProbeAnnotation = EpiMix_getInfiniumAnnotation(plat = met.platform, genome = genome)
  ProbeAnnotation = convertAnnotToDF(ProbeAnnotation)
  ProbesToRemove = ProbeAnnotation$probeID[ProbeAnnotation$CpG_chrm %in% rm.chr]
  methylation.data <- methylation.data[-which(rownames(methylation.data) %in% ProbesToRemove),  ,drop = FALSE]
  cat("Removing", length(which(rownames(methylation.data) %in% ProbesToRemove)), "CpG probes on", paste0(rm.chr, ','), "\n")

  ### Step 4: Split the experiment and the control group
  MET_Experiment <-  MET_Control <- NULL
  if(!is.null(sample.info) & !is.null(group.1) & !is.null(group.2)){
    overlapSamples = intersect(sample.info$primary, colnames(methylation.data))
    methylation.data = methylation.data[,overlapSamples,drop = F]
    sample.info = sample.info[sample.info$primary %in% overlapSamples,]
    cat("Found", length(overlapSamples), "samples with sample information.\n")
    Samples_Experiment <- sample.info[sample.info$sample.type %in% group.1,"primary"]
    Samples_Control <- sample.info[sample.info$sample.type %in% group.2,"primary"]
    MET_Experiment <- methylation.data[,Samples_Experiment,drop = FALSE]
    MET_Control <- methylation.data[,Samples_Control,drop = FALSE]
    if(ncol(MET_Experiment) == 0) stop("Cannot find methylation data with samples in the group.1 ! The sample names must overlap with the column names of the methylation data.")
    if(ncol(MET_Control) == 0) stop("Cannot find methylation data with samples in the group.2 ! The sample names must overlap with the column names of the methylation data.")
    cat("There are", ncol(MET_Experiment), "samples in the", paste0(group.1, collapse = " and "), "group, and", ncol(MET_Control), "samples in the", paste0(group.2, collapse = " and "), "group.\n")
  }else{
    warning("sample.info or group.1 and group.2 information is not provided, the funciton will perform the missing value estimation and the batch correction on the entire dataset.")
    MET_Experiment = methylation.data
  }
  rm(methylation.data); gc()

  ### Step 4: Missing value estimation
  if(sum(is.na(MET_Experiment)) > 0){
    cat("\tMissing value estimation on group.1...\n")
    MET_Experiment <- GEO_EstimateMissingValues_Methylation(MET_Experiment, MissingValueThresholdGene, MissingValueThresholdSample)
  }
  if(!is.null(MET_Control) & ncol(MET_Control) > 0 & sum(is.na(MET_Control)) > 0){
    cat("\tMissing value estimation on group.2...\n")
    MET_Control <- GEO_EstimateMissingValues_Methylation(MET_Control, MissingValueThresholdGene, MissingValueThresholdSample)
  }

  ### Step 5: Batch correction
  if(doBatchCorrection){
    MinInBatch = 0
    if(batch.correction.method == "Combat"){
      MinInBatch = 5
      if(cores > 1){
        #unregister()
        cat("Registering sockets on multiple CPU cores...\n")
        cl <- parallel :: makeCluster(cores)
      }
    }
    cat("Performing batch correction on group.1...\n")
    MET_Experiment <- CorrectBatchEffect(MET_Experiment, BatchData, batch.correction.method, MinInBatch = MinInBatch, featurePerSet = 50000)
    if(!is.null(MET_Control) & ncol(MET_Control) > 0){
      cat("Performing batch correction on group.2...\n")
      MET_Control <- CorrectBatchEffect(MET_Control, BatchData, batch.correction.method, MinInBatch = MinInBatch, featurePerSet = 50000)
    }

    # Set values <0 to 0 and >1 to 1, because of batch correction
    MET_Experiment[MET_Experiment<0]=0
    MET_Experiment[MET_Experiment>1]=1
    if (!is.null(MET_Control) & ncol(MET_Control) > 0) {
      MET_Control[MET_Control<0]=0
      MET_Control[MET_Control>1]=1
    }
    if(cores > 1 & batch.correction.method == "Combat") parallel :: stopCluster(cl)
  }


  ### Step 6: combine MET_Experiment and MET_Control into one matrix
  if(!is.null(MET_Control) & ncol(MET_Control) > 0){
    overlapProbes <- intersect(rownames(MET_Experiment),rownames(MET_Control))
    cat("Found", length(overlapProbes), "overlapping probes between group.1 and group.2 after preprocessing...\n")
    MET_Experiment = MET_Experiment[overlapProbes,,drop = F]
    MET_Control = MET_Control[overlapProbes,,drop = F]
    MET_Data = cbind(MET_Experiment, MET_Control)
    return(MET_Data)
  } else{
    return(MET_Experiment)
  }
}


#' The GEO_Preprocess_GeneExpression function
#' @description Preprocess the gene expression data from the GEO database.
#' @details
#' The preprocessing pipeline includes:
#' (1) eliminating samples and genes with too many NAs and imputing NAs.
#' (2) if the gene names (rownames) in the gene expression data are ensembl_gene_ids or ensembl_transcript_ids, translate the gene names or the transcript names to human gene symbols (HGNC).
#' (3) mapping the column names of the gene expression data to the actual sample names based on the information from "sample.map".
#' (4) doing batch correction.
#'
#' @param gene.expression.data a matrix of gene expression data with gene in rows and samples in columns.
#' @param sample.info dataframe that maps each sample to a study group. Should contain two columns: the first column (named: "primary") indicating the sample names, and the second column (named: "sample.type") indicating which study group each sample belongs to (e.g., “Experiment” vs. “Control”,  “Cancer” vs. “Normal”). Sample names in the "primary" column must coincide with the column names of the methylation.data. Please see details for more information. Default: NULL.
#' @param group.1 character vector indicating the name(s) for the experiment group. The values must coincide with the values in the "sample.type" of the sample.info dataframe.Please see details for more information. Default: NULL.
#' @param group.2 character vector indicating the names(s) for the control group. The values must coincide with the values in the "sample.type" of the sample.info dataframe. Please see details for more information. Default: NULL.
#' @param sample.map dataframe for mapping the GEO accession ID (column names) to the actual sample names. Can be the output from the GEO_getSampleMap function. Default: NULL.
#' @param MissingValueThresholdGene threshold for missing values per gene. Genes with a percentage of NAs greater than this threshold are removed. Default is 0.3.
#' @param MissingValueThresholdSample threshold for missing values per sample. Samples with a percentage of NAs greater than this threshold are removed. Default is 0.1.
#' @param doBatchCorrection logical indicating whether to perform batch correction. If TRUE, the batch data need to be provided.
#' @param BatchData dataframe with batch information. Should contain two columns: the first column indicating the actual sample names, the second column indicating the batch. Users are expected to retrieve the batch information from GEO on their own, but this can also be done using the GEO_getSampleInfo function with the "group.column" as the column indicating the batch for each sample. Defualt": NULL.
#' @param batch.correction.method character string indicating the method that be used for batch correction. Should be either "Seurat" or "Combat". Default: "Seurat".
#' @param cores number of CPU cores to be used for batch effect correction. Default: 1
#' @return gene expression data matrix with genes in rows and samples in columns.
#' @export
#' @keywords preprocess

GEO_Preprocess_GeneExpression <- function(gene.expression.data,
                                          sample.info = NULL,
                                          group.1 = NULL,
                                          group.2 = NULL,
                                          sample.map = NULL,
                                          MissingValueThresholdGene = 0.3,
                                          MissingValueThresholdSample = 0.1,
                                          doBatchCorrection = FALSE,
                                          BatchData = NULL,
                                          batch.correction.method = "Seurat",
                                          cores = 1){
  # set up some parameters
  if(cores > 1 & batch.correction.method == "Combat"){
    #unregister()
    cat("Registering sockets on multiple CPU cores...\n")
    cl <- parallel :: makeCluster(cores)
  }

   ### Step 1: check data
  gene.expression.data = as.matrix(gene.expression.data)
  if(nrow(gene.expression.data) <=1 | ncol(gene.expression.data) <=1){
    stop("the gene.expression.data matrix is empty!\n Please check whether the actual data are saved in supplementary files in GEO \n")
  }

  ### Step 2: convert column names to the actual patient names
  if(!is.null(sample.map)){
    cat("Mapping column names to the actual sample names...\n")
    sample.map = sample.map[sample.map$assay == "Gene expression", ]
    overlapSamples = intersect(sample.map$colnames, colnames(gene.expression.data))
    if(length(overlapSamples) > 0){
      gene.expression.data = gene.expression.data[, overlapSamples]
      presentSamples = match(colnames(gene.expression.data), sample.map$colnames)
      sampleNames = sample.map$primary
      sampleNames = sampleNames[presentSamples]
      colnames(gene.expression.data) = sampleNames
    }else{
      warning("No overlap samples were found between the sample map and the column names of the gene expression data! No mapping of the column names to the actual sample names was performed!")
    }
  }

  # Step 3: Mapping ENSG gene names to external gene names and remove duplicated genes
  if(any(grepl("ENSG",rownames(gene.expression.data))) | any(grepl("ENST",rownames(gene.expression.data)))){
    cat("Mapping transcripts to genes...\n")
    gene.expression.data = convertGeneNames(gene.expression.data)
  }

  if(length(rownames(gene.expression.data)) != length(unique(rownames(gene.expression.data)))){
    gene.expression.data = removeDuplicatedGenes(gene.expression.data)
  }

  ### Step 4: Split the experiment and the control group
  MET_Experiment <-  MET_Control <- NULL
  if(!is.null(sample.info) & !is.null(group.1) & !is.null(group.2)){
    overlapSamples = intersect(sample.info$primary, colnames(gene.expression.data))
    gene.expression.data = gene.expression.data[,overlapSamples,drop = F]
    sample.info = sample.info[sample.info$primary %in% overlapSamples,]
    cat("Found", length(overlapSamples), "samples with sample information.\n")
    Samples_Experiment <- sample.info[sample.info$sample.type %in% group.1,"primary"]
    Samples_Control <- sample.info[sample.info$sample.type %in% group.2,"primary"]
    MET_Experiment <- gene.expression.data[,Samples_Experiment,drop = FALSE]
    MET_Control <- gene.expression.data[,Samples_Control,drop = FALSE]
    if(ncol(MET_Experiment) == 0) stop("Cannot find methylation data with samples in the group.1 ! The sample names must overlap with the column names of the gene expression data.")
    if(ncol(MET_Control) == 0) stop("Cannot find methylation data with samples in the group.2 ! The sample names must overlap with the column names of the gene expression data.")
    cat("There are", ncol(MET_Experiment), "samples in the", paste0(group.1, collapse = " and "), "group, and", ncol(MET_Control), "samples in the", paste0(group.2, collapse = " and "), "group.\n")
  }else{
    warning("sample.info or group.1 and group.2 information is not provided, the funciton will perform the missing value estimation and the batch correction on the entire dataset.")
    MET_Experiment = gene.expression.data
  }
  rm(gene.expression.data); gc()

  ### Step 5: Missing value estimation
  if(sum(is.na(MET_Experiment)) > 0){
    cat("\tMissing value estimation on group.1...\n")
    MET_Experiment <- GEO_EstimateMissingValues_Methylation(MET_Experiment, MissingValueThresholdGene, MissingValueThresholdSample)
  }
  if(!is.null(MET_Control)){
    if(sum(is.na(MET_Control)) > 0){
      cat("\tMissing value estimation on group.2...\n")
      MET_Control <- GEO_EstimateMissingValues_Methylation(MET_Control, MissingValueThresholdGene, MissingValueThresholdSample)
    }
  }

  ### Step 6: Batch correction
  if(doBatchCorrection){
    MinInBatch = 0
    if(batch.correction.method == "Combat"){
      MinInBatch = 5
    }
    cat("Performing batch correction on group.1...\n")
    MET_Experiment <- CorrectBatchEffect(MET_Experiment, BatchData, batch.correction.method, MinInBatch = MinInBatch, featurePerSet = 50000)
    if(!is.null(MET_Control)){
      cat("Performing batch correction on group.2...\n")
      MET_Control <- CorrectBatchEffect(MET_Control, BatchData, batch.correction.method, MinInBatch = MinInBatch, featurePerSet = 50000)
    }

    # Set values <0 to 0 and >1 to 1, because of batch correction
    MET_Experiment[MET_Experiment<0]=0
    MET_Experiment[MET_Experiment>1]=1
    if (!is.null(MET_Control) & ncol(MET_Control) > 0) {
      MET_Control[MET_Control<0]=0
      MET_Control[MET_Control>1]=1
    }
  }

  if(cores > 1) parallel :: stopCluster(cl)

  ### Step 6: combine MET_Experiment and MET_Control into one matrix
  if(!is.null(MET_Control)){
    overlapProbes <- intersect(rownames(MET_Experiment),rownames(MET_Control))
    MET_Experiment = MET_Experiment[overlapProbes,,drop = F]
    MET_Control = MET_Control[overlapProbes,,drop = F]
    MET_Data = cbind(MET_Experiment, MET_Control)
    return(as.matrix(MET_Data))
  } else{
    return(as.matrix(MET_Experiment))
  }
}

#' The convertGeneNames function
#' @description auxiliary function to translate ensembl_gene_ids or ensembl_transcript_ids to human gene symbols (HGNC)
#' @param gene.expression.data gene expression data matrix with the rownames to be the ensembl_gene_ids or ensembl_transcript_ids
#' @return gene expression matrix with rownames translated to human gene symbols (HGNC)
#' @keywords internal

convertGeneNames <- function(gene.expression.data){
    filter = NULL
    if(any(grepl("ENSG",rownames(gene.expression.data)))){
      cat("Mapping ENSG gene names to human gene symbols..\n")
      if(all(grepl("\\.", rownames(gene.expression.data)))){
        filter = "ensembl_gene_id_version"
      }else{
        filter = "ensembl_gene_id"
      }
    }else{
      if(all(grepl("\\.", rownames(gene.expression.data)))){
        cat("Mapping ENST transcript names to human gene symbols..\n")
        filter = "ensembl_transcript_id_version"
      }else{
        filter = "ensembl_transcript_id"
      }
    }
    cat("Retrieving the transcript annotation from the Ensembl server ...\n")
    mart <- biomaRt :: useDataset("hsapiens_gene_ensembl", biomaRt ::  useMart("ensembl"))
    ensemblID <-  rownames(gene.expression.data)
    ensembl_gene_map <- biomaRt :: getBM(filters= filter,
                                         attributes= c(filter,"hgnc_symbol", "chromosome_name"),
                                         values = ensemblID, mart= mart)
    cat("Removing", sum(ensembl_gene_map$hgnc_symbol==""), "transcripts that can not be mapped to human gene symbols.\n")
    ensembl_gene_map <- ensembl_gene_map[ensembl_gene_map$hgnc_symbol!="",]

    #get rid of those genes on the haplotypic region (chromosome names starting with "CHR")
    ensembl_gene_map <- ensembl_gene_map[!startsWith(ensembl_gene_map$chromosome_name, "CHR"), ] #16,033
    cat("Removed", sum(startsWith(ensembl_gene_map$chromosome_name, "CHR")), "trnascripts mapped to the haplotypic region\n")

    #filter only genes with ensembl annotation
    overlapGenes = intersect(ensembl_gene_map[[filter]], rownames(gene.expression.data)) #16,033
    #cat("Removing", nrow(gene.expression.data) - length(overlapGenes), "transcripts that can not be mapped to human genes.\n")
    gene.expression.data = gene.expression.data[overlapGenes,,drop = F]
    presentEnsemblID = match(rownames(gene.expression.data), ensembl_gene_map$ensembl_gene_id)
    geneNames <- ensembl_gene_map$hgnc_symbol
    geneNames <- geneNames[presentEnsemblID]
    rownames(gene.expression.data) <- geneNames
    return(gene.expression.data)
}


#' The removeDuplicatedGenes function
#' @description sum up the transcript expression values if a gene has multiple transcripts
#' @param GEN_data gene expression data matrix
#' @return gene expression data matrix with duplicated genes removed
#'
removeDuplicatedGenes <- function(GEN_data){
  cat("Removing duplicated transcripts...\n")
  duplicatedGenes = rownames(GEN_data)[duplicated(rownames(GEN_data))]
  singleGenes = GEN_data[!rownames(GEN_data) %in% duplicatedGenes,,drop = F]

  uniqueGenes = unique(duplicatedGenes)
  exp = data.frame()
  iterations = length(uniqueGenes)
  pb <- txtProgressBar(max = iterations, style = 3)
  for(i in 1:iterations){
    gene = uniqueGenes[i]
    temp = GEN_data[which(rownames(GEN_data) == gene),,drop = F]
    if(nrow(temp) > 1){
      temp = matrix(apply(temp,2, sum), nrow =1)
      colnames(temp) = colnames(data)
      rownames(temp) = gene
    }
    exp = rbind(exp, temp)
    setTxtProgressBar(pb,i)
  }

  cleaned_GEN_Data = rbind(singleGenes, exp)
  cat("\nRemoved", nrow(GEN_data) - nrow( cleaned_GEN_Data), "duplicated transcripts\n")
  return(cleaned_GEN_Data)
}

#' The GEO_EstimateMissingValues_Molecular function
#'
#' Internal. Removes samples and genes with more missing values than the MissingValueThreshold, and imputes remaining missing values using Tibshirani's KNN method.
#' @param MET_Data methylation data or gene expression data matrix.
#' @param MissingValueThresholdGene threshold for missing values per gene. Genes with a percentage of NAs greater than this threshold are removed. Default is 0.3.
#' @param MissingValueThresholdSample threshold for missing values per sample. Samples with a percentage of NAs greater than this threshold are removed. Default is 0.1.
#' @return the dataset with imputed values and possibly some genes or samples deleted.
#' @keywords internal
#'
GEO_EstimateMissingValues_Molecular <- function(MET_Data, MissingValueThresholdGene = 0.3, MissingValueThresholdSample = 0.1) {

  # FIRST REMOVING BAD SAMPLEs

  # removing clones with too many missing values
  NrMissingsPerGene=apply(MET_Data,1,function(x) sum(is.na(x))) /ncol(MET_Data)
  cat("Removing",sum(NrMissingsPerGene>MissingValueThresholdGene),"genes with more than",MissingValueThresholdGene*100,"% missing values.\n")
  if (sum(NrMissingsPerGene>MissingValueThresholdGene)>0) MET_Data=MET_Data[NrMissingsPerGene<MissingValueThresholdGene,]

  # removing patients with too many missings values
  NrMissingsPerSample=apply(MET_Data,2,function(x) sum(is.na(x))) /nrow(MET_Data)
  cat("Removing",sum(NrMissingsPerSample>MissingValueThresholdSample),"samples with more than",MissingValueThresholdSample*100,"% missing values.\n")
  if (sum(NrMissingsPerSample>MissingValueThresholdSample)>0) MET_Data=MET_Data[,NrMissingsPerSample<MissingValueThresholdSample]

  # knn impute using Tibshirani's method
  if (length(colnames(MET_Data))>1) {
    k=15
    KNNresults=impute::impute.knn(as.matrix(MET_Data),k)
    MET_Data_KNN=KNNresults$data
    # cleaning up sample names
    return(MET_Data_KNN)

  } else {
    # when only 1 sample,need to make a matrix again
    #MET_Data=as.matrix(MET_Data)
    return(MET_Data)
  }
}

#' The GEO_EstimateMissingValues_Methylation function
#'
#' Internal. Removes samples and probes with more missing values than the MissingValueThreshold, and imputes remaining missing values using Tibshirani's KNN method.
#' @param MET_Data methylation data or gene expression data matrix.
#' @param MissingValueThresholdGene threshold for missing values per gene. Genes with a percentage of NAs greater than this threshold are removed. Default is 0.3.
#' @param MissingValueThresholdSample threshold for missing values per sample. Samples with a percentage of NAs greater than this threshold are removed. Default is 0.1.
#' @return the dataset with imputed values and possibly some genes or samples deleted.
#' @keywords internal
#'
GEO_EstimateMissingValues_Methylation <- function(MET_Data, MissingValueThresholdGene = 0.3, MissingValueThresholdSample = 0.3) {

  # removing probes with too many missing values
  NrMissingsPerGene=apply(MET_Data,1,function(x) sum(is.na(x))) /ncol(MET_Data)
  cat("Removing",sum(NrMissingsPerGene>MissingValueThresholdGene),"probes with more than",MissingValueThresholdGene*100,"% missing values.\n")
  if (sum(NrMissingsPerGene>MissingValueThresholdGene)>0) MET_Data=MET_Data[NrMissingsPerGene<MissingValueThresholdGene,]

  # removing samples with too many missings values
  NrMissingsPerSample=apply(MET_Data,2,function(x) sum(is.na(x))) /nrow(MET_Data)
  cat("Removing",sum(NrMissingsPerSample>MissingValueThresholdSample),"samples with more than",MissingValueThresholdSample*100,"% missing values.\n")
  if (sum(NrMissingsPerSample>MissingValueThresholdSample)>0) MET_Data=MET_Data[,NrMissingsPerSample<MissingValueThresholdSample]

  # knn impute using Tibshirani's method
  if (length(colnames(MET_Data))>1) {
    k=15
    KNNresults=impute::impute.knn(as.matrix(MET_Data),k)
    MET_Data_KNN=KNNresults$data
    # cleaning up sample names
    return(MET_Data_KNN)

  } else {
    # when only 1 sample,need to make a matrix again
    #MET_Data=as.matrix(MET_Data)
    return(MET_Data)
  }
}



