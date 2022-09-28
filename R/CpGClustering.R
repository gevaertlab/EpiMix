#' The ClusterProbes function
#'
#' This function uses the annotation for Illumina methylation arrays to map each probe to a gene. Then, for each gene,
#' it clusters all its CpG sites using hierchical clustering and Pearson correlation as distance and complete linkage.
#' If data for normal samples is provided, only overlapping probes between cancer and normal samples are used.
#' Probes with SNPs are removed.
#' This function is prepared to run in parallel if the user registers a parallel structure, otherwise it runs sequentially.
#' This function also cleans up the sample names, converting them to the 12 digit format.
#' @param MET_data data matrix for methylation.
#' @param ProbeAnnotation GRange object for probe annoation.
#' @param CorThreshold correlation threshold for cutting the clusters.
#' @return List with the clustered data sets and the mapping between probes and genes.
#' @export
#' @keywords cluter_probes
#' @importFrom foreach %dopar%

ClusterProbes <- function(MET_data, ProbeAnnotation, CorThreshold = 0.4) {

  ###### only iterating over genes that have probes present
  # Getting the positions relative to probe annotation of the probes present in this data set.
  ProbeAnnotation = data.frame(ProbeAnnotation)
  PresentProbes=intersect(ProbeAnnotation$probe, rownames(MET_data))
  UniqueGenes = unique(sort(ProbeAnnotation$gene[which(ProbeAnnotation$probe %in% PresentProbes)]))
  UniqueGenes=UniqueGenes[which(UniqueGenes != "")]

  # Cluster CpG sites at gene-level
  cat("Clustering",length(rownames(MET_data)),"probes in CpG site clusters.\n")
  iterations <- length(UniqueGenes)
  pb <- utils::txtProgressBar(max = iterations, style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  tmpClusterResults = foreach::foreach(i=seq_len(iterations), .export='TCGA_GENERIC_MET_ClusterProbes_Helper_ClusterGenes_with_hclust',
                                       .options.snow = opts, .verbose = FALSE) %dopar% {
    TCGA_GENERIC_MET_ClusterProbes_Helper_ClusterGenes_with_hclust(UniqueGenes[i], ProbeAnnotation,MET_Cancer = MET_data,  MET_Normal = NULL, CorThreshold)
                                       }
  close(pb)

  # Check the number of total clusters
  ProbeSum = 0
  METmatrixSum = 0
  for ( i in 1:length(UniqueGenes)){
    ProbeSum = ProbeSum + nrow(tmpClusterResults[[i]][[3]])
    METmatrixSum = METmatrixSum + nrow(tmpClusterResults[[i]][[1]])
  }

  # Initiate new matrices to save clustered methylation data and probe mapping
  MET_data_C=matrix(0,METmatrixSum,length(colnames(MET_data)))
  colnames(MET_data_C)=colnames(MET_data)
  ProbeMapping=matrix(0,ProbeSum,2)

  # Map clustered methylation data and probe mapping to new matrices
  ProbeCounter=1
  METmatrixCounter=1
  ClusteredRownames=c()

  for ( i in 1:length(UniqueGenes)){
    MET_data_C[METmatrixCounter:(METmatrixCounter-1+nrow(tmpClusterResults[[i]][[1]])),]=tmpClusterResults[[i]][[1]]
    ProbeMapping[ProbeCounter:(ProbeCounter-1+nrow(tmpClusterResults[[i]][[3]])),]=tmpClusterResults[[i]][[3]]
    ClusteredRownames=c(ClusteredRownames,rownames(tmpClusterResults[[i]][[1]]))
    METmatrixCounter=METmatrixCounter+nrow(tmpClusterResults[[i]][[1]])
    ProbeCounter=ProbeCounter+nrow(tmpClusterResults[[i]][[3]])
  }

  # remove excessively large matrix
  MET_data_C=MET_data_C[1:length(ClusteredRownames),]
  rownames(MET_data_C)=ClusteredRownames

  cat("\nFound",length(rownames(MET_data_C)),"CpG site clusters.\n")
  return(list(MET_data_clustered = MET_data_C, ProbeMapping=ProbeMapping))
}

#' The TCGA_GENERIC_MET_ClusterProbes_Helper_ClusterGenes_with_hclust function
#'
#' Internal. Cluster probes into genes.
#' @param Gene gene.
#' @param ProbeAnnotation data set matching probes to genes.
#' @param MET_Cancer data matrix for cancer samples.
#' @param MET_Normal data matrix for normal samples.
#' @param CorThreshold correlation threshold for cutting the clusters.
#' @return List with the clustered data sets and the mapping between probes and genes.
#' @keywords internal
#'
TCGA_GENERIC_MET_ClusterProbes_Helper_ClusterGenes_with_hclust <- function(Gene,
                                                                           ProbeAnnotation,
                                                                           MET_Cancer,
                                                                           MET_Normal = NULL,
                                                                           CorThreshold = 0.4) {

  # first lookup the probes matching a single gene. DO NOT USE grep, it does
  # not do exact matching, but looks for the pattern anywhere !!!
  Probes <- ProbeAnnotation[which(ProbeAnnotation[, 2] == Gene), 1]
  Probes <- Probes[which(Probes %in% rownames(MET_Cancer))]
  Probes <- unique(Probes)

  METcancer_Clustered <- matrix(0, 0, length(colnames(MET_Cancer)))
  if (!is.null(MET_Normal))
    METnormal_Clustered <- matrix(0, 0, length(colnames(MET_Normal)))
  Clusternames <- c()
  InverseCorrelationThreshold <- 1 - CorThreshold
  GeneClustersForProbeMapping <- array(dim = length(Probes))
  if (length(Probes) > 1) {
    ProbeCorrelation <- cor(t(MET_Cancer[Probes, ]), method = "pearson")
    ClusterResults <- hclust(as.dist(1 - ProbeCorrelation), method = "complete",
                             members = NULL)
    # plot(ClusterResults)
    Clusters <- cutree(ClusterResults, h = InverseCorrelationThreshold)
    for (i in 1:length(unique(Clusters))) {
      tmpGeneProbes <- Probes[Clusters == i]
      if (length(tmpGeneProbes) > 1) {

        tmpAveragedProfile <- colMeans(MET_Cancer[tmpGeneProbes, ])
        METcancer_Clustered <- rbind(METcancer_Clustered, tmpAveragedProfile)
        # Same for normal
        if (!is.null(MET_Normal))
          tmpAveragedProfile <- colMeans(MET_Normal[tmpGeneProbes, ])
        if (!is.null(MET_Normal))
          METnormal_Clustered <- rbind(METnormal_Clustered, tmpAveragedProfile)
      } else {
        METcancer_Clustered <- rbind(METcancer_Clustered, MET_Cancer[tmpGeneProbes,
        ])
        if (!is.null(MET_Normal))
          METnormal_Clustered <- rbind(METnormal_Clustered, MET_Normal[tmpGeneProbes,
          ])
      }
      Clusternames <- c(Clusternames, paste(Gene, "---Cluster", i, sep = ""))
      pos <- which(Probes %in% tmpGeneProbes)
      GeneClustersForProbeMapping[pos] <- paste(Gene, "---Cluster", i, sep = "")
    }
    rownames(METcancer_Clustered) <- Clusternames
    if (!is.null(MET_Normal))
      rownames(METnormal_Clustered) <- Clusternames

  } else {
    METcancer_Clustered <- MET_Cancer[Probes, , drop = FALSE]
    if (!is.null(MET_Normal))
      METnormal_Clustered <- MET_Normal[Probes, , drop = FALSE]
    Clusternames <- Gene
    GeneClustersForProbeMapping <- Gene
    rownames(METcancer_Clustered) <- Clusternames
    if (!is.null(MET_Normal))
      rownames(METnormal_Clustered) <- Clusternames
  }

  ProbeMapping <- t(rbind(Probes, GeneClustersForProbeMapping))

  if (is.null(MET_Normal))
    METnormal_Clustered <- NULL
  return(list(METcancer_Clustered, METnormal_Clustered, ProbeMapping))
}


#' The translateMethylMixResults function
#' @description unfold clustered MethylMix results to single CpGs
#' @param MethylMixResults list of MethylMix output
#' @param probeMapping dataframe of probe to gene-cluster mapping
#' @return list of unfolded MethylMix results
#'
#'
translateMethylMixResults <- function(MethylMixResults, probeMapping){

  # Convert MethylationStates
  MethylationStates <- data.frame(MethylMixResults$MethylationStates, check.names = FALSE) # (573, 582)
  colnames(probeMapping) <- c("CpG", "Gene_Cluster") # (23835, 2)
  MethylationStates['Gene_Cluster'] <- rownames(MethylationStates)
  MethylationStates <- merge(MethylationStates, probeMapping, by = "Gene_Cluster")
  MethylationStates <- MethylationStates[!duplicated(MethylationStates$CpG),]
  rownames(MethylationStates) <- MethylationStates$CpG
  keeps <- setdiff(colnames(MethylationStates), c("CpG", "Gene_Cluster"))
  MethylationStates <- MethylationStates[, keeps] # (797,582)
  MethylationStates <- as.matrix(MethylationStates)

  # Convert Nr component
  NrComponents <- apply(MethylationStates, 1, function(x) length(unique(x)))

  # Convert Mixture States
  mixture.states <- apply(MethylationStates, 1, function(x) sort(unique(x)))

  # Convert Classifications
  Classfications <- data.frame(MethylMixResults$Classifications, check.names = FALSE) # (573, 582)
  colnames(probeMapping) <- c("CpG", "Gene_Cluster") # (23835, 2)
  Classfications ['Gene_Cluster'] <- rownames(Classfications )
  Classfications  <- merge(Classfications, probeMapping, by = "Gene_Cluster")
  Classfications  <- Classfications [!duplicated(Classfications$CpG),]
  rownames(Classfications) <- Classfications$CpG
  keeps <- setdiff(colnames(Classfications), c("CpG", "Gene_Cluster"))
  Classfications  <- Classfications [, keeps] # (797,582)
  Classfications  <- as.matrix(Classfications)

  # Convert Methylation Drivers
  MethylationDrivers <-  rownames(MethylationStates)

  MethylMixResults$NrComponents <- NrComponents
  MethylMixResults$MethylationStates <- MethylationStates
  MethylMixResults$MixtureStates <- mixture.states
  MethylMixResults$Classifications <- Classfications
  MethylMixResults$MethylationDrivers <- MethylationDrivers

  return(MethylMixResults)
}
