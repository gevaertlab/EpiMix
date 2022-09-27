# Helper function used by EpiMix

#' The filterProbes function
#' @description filter CpG sites based on user-specified conditions
#' @param mode analytic mode
#' @param gene.expression.data matrix of gene expression data
#' @param listOfGenes list of genes of interest
#' @param promoters logic indicating whether to filter CpGs on promoters
#' @param met.platform methylation platform
#' @param genome genome build version
#'
#' @return filtered ProbeAnnotation

filterProbes <- function(mode, gene.expression.data, listOfGenes, promoters, met.platform,
    genome) {
    cat("Fetching probe annotation...\n")
    ProbeAnnotation <- getProbeAnnotation(mode = mode, met.platform = met.platform,
        genome = genome)
    if (!is.null(listOfGenes)) {
        warning("Please input the selected gene names in the upper case format: e.g., IGF1,  CCND3, MIRLET7A1, MIR10226, LINC01409, TTLL10-AS1",
            immediate. = TRUE)
       ProbeAnnotation <- ProbeAnnotation[which(ProbeAnnotation$gene %in% listOfGenes),]
       cat("Found", length(unique(ProbeAnnotation$probe)), "CpGs associated with the user-specified genes.\n")
    }
    if (!is.null(gene.expression.data)) {
        # we only select the CpGs associated with genes with expression
        # data available
       ProbeAnnotation <- ProbeAnnotation[ProbeAnnotation$gene %in% rownames(gene.expression.data), ]
    }
    if (promoters) {
        cat("Selecting CpGs associated with gene promoters...\n")
        if(mode == "Regular" | mode == "lncRNA"){
          promoters <- getFeatureProbe(met.platform = met.platform,
                                       genome = genome,
                                       promoter = TRUE,
                                       TSS.range = list(upstream = 2000, downstream = 500))
          promoter.probes <- names(promoters)
          overlapProbes <- intersect(promoter.probes, ProbeAnnotation$probe)
        }else if(mode == "miRNA"){
          overlapProbes <- ProbeAnnotation$probe[which(ProbeAnnotation$promoter ==1)]
        }
        cat("Found", length(overlapProbes), "CpGs on promoters\n")
        ProbeAnnotation <- ProbeAnnotation[which(ProbeAnnotation$probe %in% overlapProbes), ]
    }
    return(ProbeAnnotation)
}

#' The getProbeAnnotation function
#' @description Helper function to get the probe annotation based on mode
#' @param mode analytic mode
#' @param met.platform methylation platform
#' @param genome genome build version
#'
#' @return a ProbeAnnotation dataframe consisting of two columns: probe, gene

getProbeAnnotation <- function(mode, met.platform, genome) {
    ProbeAnnotation <- NULL
    if (mode == "Regular") {
        manifest <- EpiMix_getInfiniumAnnotation(plat = met.platform, genome = genome)
        ProbeAnnotation <- data.frame(probe = names(manifest), gene = manifest$gene)
        ProbeAnnotation <- .mapProbeGene(ProbeAnnotation)
    } else if (mode == "miRNA") {
        if (met.platform == "HM27") {
            ProbeAnnotation <- EpiMix_GetData("HM27_miRNA_probes")
        } else if (met.platform == "HM450") {
            ProbeAnnotation <- EpiMix_GetData("HM450_miRNA_probes")
        } else if (met.platform == "EPIC") {
            ProbeAnnotation <- EpiMix_GetData("EPIC_miRNA_probes")
        }
    } else if (mode == "lncRNA") {
        if (met.platform == "HM27") {
            ProbeAnnotation <- EpiMix_GetData("HM27_lncRNA_probes")
        } else if (met.platform == "HM450") {
            ProbeAnnotation <- EpiMix_GetData("HM450_lncRNA_probes")
        } else if (met.platform == "EPIC") {
            ProbeAnnotation <- EpiMix_GetData("EPIC_lncRNA_probes")
        }
        ProbeAnnotation <- data.frame(probe = ProbeAnnotation, gene = names(ProbeAnnotation))
    }
    return(ProbeAnnotation)
}

#' The addGeneNames function
#' @description Given a dataframe with a column of probe names, add the gene names
#' @param df_data a dataframe with a column named Probe
#' @param ProbeAnnotation a dataframe with ProbeAnnotation, including one column named 'probe' and another column named 'gene'
#'
#' @return a dataframe with added gene names
#'
addGeneNames <- function(df_data, ProbeAnnotation) {
    ProbeAnnotation <- ProbeAnnotation %>%
        dplyr::group_by(.data$probe) %>%
        dplyr::mutate(Genes = paste(.data$gene, collapse = ";")) %>%
        dplyr::select(.data$probe, .data$Genes)
    colnames(ProbeAnnotation) <- c("Probe", "Gene")
    df_data <- merge(df_data, ProbeAnnotation, all.x = TRUE)
    df_data <- distinct(df_data)
    return(df_data)
}



