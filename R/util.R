# Utility function used by EpiMix

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

utils::globalVariables(c("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene", "i",
    "."))

#' The EpiMix_getInfiniumAnnotation function
#' @description fetch the Infinium probe annotation from the seasameData library
#' @param plat character string indicating the methylation platform
#' @param genome character string indicating the version of genome build
#'
#' @return a GRange object of probe annotation
#' @keywords internal

EpiMix_getInfiniumAnnotation <- function(plat = "EPIC", genome = "hg38") {
    hubID <- NULL
    if (tolower(genome) == "hg19" & toupper(plat) == "HM27")
        hubID <- "EH3672"
    if (tolower(genome) == "hg38" & toupper(plat) == "HM27")
        hubID <- "EH3673"
    if (tolower(genome) == "hg19" & toupper(plat) == "HM450")
        hubID <- "EH3674"
    if (tolower(genome) == "hg38" & toupper(plat) == "HM450")
        hubID <- "EH3675"
    if (tolower(genome) == "hg19" & toupper(plat) == "EPIC")
        hubID <- "EH3670"
    if (tolower(genome) == "hg38" & toupper(plat) == "EPIC")
        hubID <- "EH3671"
    ProbeAnnotation <- ExperimentHub::ExperimentHub()[[hubID]]
    return(ProbeAnnotation)
}

#' The convertAnnotToDF function
#' @description convert the probe annotation from the GRange object to a dataframe
#' @param annot a GRange object of probe annotation, can be the object returned from the getInfiniumAnnotation function.
#' @return a dataframe with chromosome, beginning and end position, mapped gene information for each CpG probe
#' @keywords internal
#'
convertAnnotToDF <- function(annot) {
    df.annot <- data.frame(CpG_chrm = GenomicRanges::seqnames(annot), CpG_beg = GenomicRanges::start(ranges(annot)),
        CpG_end = GenomicRanges::end(ranges(annot)), probeID = names(annot), gene = GenomicRanges::mcols(annot)$gene)
    return(df.annot)
}

#' The .splitMetData function
#' @description Helper function to split the methylation data matrix into the experimental group and the control group
#' @param methylation.data methylation data matrix
#' @param sample.info sample information matrix
#' @param group.1 name of group.1
#' @param group.2 name of group.2
#' @keywords internal
#' @return a list with methylation data of group.1 and group.2

.splitMetData <- function(methylation.data, sample.info, group.1, group.2) {
    MET_Experiment <- MET_Control <- NULL
    if (!is.null(sample.info) & !is.null(group.1) & !is.null(group.2)) {
        Samples_Experiment <- intersect(colnames(methylation.data), sample.info$primary[which(sample.info$sample.type %in%
            group.1)])
        Samples_Control <- intersect(colnames(methylation.data), sample.info$primary[which(sample.info$sample.type %in%
            group.2)])
        cat("Found", length(Samples_Experiment), "samples in group.1 and", length(Samples_Control),
            "samples in group.2\n")
        if (length(Samples_Experiment) == 0)
            stop("Cannot find methylation data for samples in the group.1! Please check sample information.")
        if (length(Samples_Control) == 0)
            stop("Cannot find methylation data for samples in the group.2! Please check sample information.")
        MET_Experiment <- methylation.data[, Samples_Experiment, drop = FALSE]
        MET_Control <- methylation.data[, Samples_Control, drop = FALSE]
    } else {
        MET_Experiment <- methylation.data
    }
    return(list(MET_Experiment = MET_Experiment, MET_Control = MET_Control))
}

#' The .mapProbeGene function
#' @description since in the original probe annotation, a specific probe can be mapped to multiple genes, this function splits the rows and maps each probe to a signle gene in a row.
#' @param df.annot a dataframe with probe annotation, can be the object returned from the convertAnnotToDF function.
#' @return a dataframe with 1:1 mapping of probe and gene
#' @keywords internal
#'
.mapProbeGene <- function(df.annot) {
    df.annot <- df.annot[!is.na(df.annot$gene), ]
    df.annot <- tidyr::separate_rows(df.annot, .data$gene, sep = ";")
    df.annot <- dplyr::distinct(df.annot)
    return(df.annot)
}


#' The getMethStates function
#' @description Helper function that adds a methyaltion state label to each driver probe
#' @param MethylMixResults the list object returned from the EpiMix function
#' @param DM.probes character vector of differentially methylated probes.
#' @keywords internal
#' @return a character vector with the methylation state ('Hypo', 'Hyper' or 'Dual') for each probe. The names for the vector are the probe names and the values are the methylation state.
#'
getMethStates <- function(MethylMixResults, DM.probes) {
    mixture.states <- MethylMixResults$MixtureStates
    mixture.states <- mixture.states[DM.probes]
    meth.states <- c()
    for (probe in DM.probes) {
        state <- NULL
        DMValues <- mixture.states[[probe]]
        state <- getMethStates_Helper(DMValues = DMValues)
        meth.states <- append(meth.states, state)
    }
    names(meth.states) <- DM.probes
    return(meth.states)
}

#' The getMethStates_Helper function
#' @description helper function to determine the methylation state based on DM values
#' @param DMValues a character vector indicating the DM values of a CpG site
#' @return a character string incdicating the methylation state of the CpG
#'
getMethStates_Helper <- function(DMValues) {
    if (max(DMValues) > 0 & min(DMValues) < 0) {
        state <- "Dual"
    } else if (max(DMValues) > 0) {
        state <- "Hyper"
    } else {
        state <- "Hypo"
    }
    return(state)
}

#' The splitmatix function
#'
#' @param x A matrix
#' @param by A character specify if split the matrix by row or column.
#' @keywords internal
#'
#' @return A list each of which is the value of each row/column in the matrix.
#'
splitmatrix <- function(x, by = "row") {
    if (by %in% "row") {
        out <- split(x, rownames(x))
    } else if (by %in% "col") {
        out <- split(x, colnames(x))
    }
    return(out)
}

#' The get.prevalence function
#' @description Helper function to get the methylation state and the
#' prevalence of the differential methylation of a CpG sites in the study population.
#' @param MET_matrix matrix of methylation states
#' @keywords internal
#'
#' @return a list of prevalence for the abnormal methylation
#'
get.prevalence <- function(MethylMixResults) {
    MET_matrix <- MethylMixResults$MethylationStates
    hypo_prev <- apply(MET_matrix, 1, function(x) round(length(x[x < 0])/length(x) *
        100, 3))
    hyper_prev <- apply(MET_matrix, 1, function(x) round(length(x[x > 0])/length(x) *
        100, 3))
    prev.data <- data.frame(Probe = rownames(MET_matrix), row.names = NULL)
    prev.data["Prevalence of hypo (%)"] <- hypo_prev
    prev.data["Prevalence of hyper (%)"] <- hyper_prev
    return(prev.data)
}



#' The get.chromosome function
#' @description given a list of genes, get the chromosomes of these genes.
#' @param genes character vector with the gene names
#' @param genome character string indicating the genome build version, can be either 'hg19' or 'hg38'
#' @keywords internal
#'
#' @return a dataframe for the mapping between genes and their chromosomes.
#'
get.chromosome <- function(genes, genome) {
    gene.annotation <- as.data.frame(getTSS(genome))
    gene.annotation <- distinct(gene.annotation[, c("external_gene_name", "seqnames")])
    gene.annotation <- gene.annotation[which(gene.annotation$external_gene_name %in%
        genes), ]
    colnames(gene.annotation) <- c("Gene", "Chr")
    return(gene.annotation)
}


#' The getRandomGenes function
#' @description Helper function to get a set of random genes located on different chromosomes of the target CpG.
#' @param target.probe character string indicating the target CpG for generating the permutation p values.
#' @param gene.expression.data a matrix of gene expression data.
#' @param ProbeAnnotation GRange object of probe annotation.
#' @param genome character string indicating the genome build version, can be either 'hg19' or 'hg38'.
#' @param perm the number of permutation tests. Default: 1000
#' @keywords internal
#'
#' @return a dataframe for the permutation genes and p values for the target CpG site.
#'
getRandomGenes <- function(target.probe, gene.expression.data, ProbeAnnotation,
    genome = "hg38", perm = 1000) {

    # get the chromosome of the current CpG
    target.chr <- as.character(seqnames(ProbeAnnotation)[which(names(ProbeAnnotation) ==
        target.probe)])
    # get the chromosomes of all the genes with gene expression data avaliable
    all.genes <- rownames(gene.expression.data)
    gene.chr <- get.chromosome(all.genes, genome)
    # get the genes that are not on the same chromosome of the current CpG
    random.genes <- gene.chr$Gene[which(!gene.chr$Chr %in% c(target.chr, "chrX",
        "chrY"))]
    if (length(unique(random.genes)) < perm) {
        warning("There is not enough genes to generate ", perm, " permutations. Using ",
            length(unique(random.genes)), " random genes instead.", immediate. = TRUE)
        perm <- length(unique(random.genes))
    }
    random.genes <- sample(gene.chr$Gene[which(!gene.chr$Chr %in% c(target.chr, "chrX",
        "chrY"))], perm)
    return(random.genes)
}

#' The filterMethMatrix function
#' @details This function filters methylation states from the beta mixture modeling for each probe. The filtered probes can be used to model gene expression by Wilcoxon test.
#' @param MET_matrix a matrix of methylation states from the EpiMix results
#' @param control.names a character vector of control sample names
#' @param gene.expression.data a matrix with gene expression data
#' @keywords internal
#' @return a matrix of methylation states for each differentially methylated probe with probes in rows and patient in columns.
#'
filterMethMatrix <- function(MET_matrix, control.names, gene.expression.data) {
    # Combine the methylation states for the experiment group and the control
    # group into one matrix
    MET_matrix_control <- matrix(0, nrow(MET_matrix), length(control.names))
    colnames(MET_matrix_control) <- control.names
    MET_matrix <- cbind(MET_matrix, MET_matrix_control)

    # Filter 1: filter on samples which have gene expression data available
    if (nrow(MET_matrix) > 0) {
        common.samps <- intersect(colnames(MET_matrix), colnames(gene.expression.data))
        MET_matrix <- MET_matrix[, common.samps, drop = FALSE]
    }

    # Filter 2: filtering out probes with only one methylation state
    if (ncol(MET_matrix) > 0) {
      METstatesCounts <- apply(MET_matrix, 1, function(x) length(unique(x)))
      MET_matrix <- MET_matrix[METstatesCounts > 1, ,drop = FALSE]
    }

    # Filter 3: filtering out probes with the minority group having sample
    # number < 3. Otherwise, the Wilcoxon test won't work.
    if (nrow(MET_matrix) > 0) {
        METminCounts <- apply(MET_matrix, 1, function(x) min(as.numeric(table(x))))
        MET_matrix <- MET_matrix[METminCounts >= 3, ,drop = FALSE]
    }
    return(MET_matrix)
}


#' getRoadMapEnhancerProbes
#' @details get the CpG probes that locate at the enhancer regions identified by the Roadmap epigenomics project
#' @param met.platform character string indicating the methylation platform, can be either 'EPIC' or 'HM450'
#' @param genome character string indicating the genome build version, can be either 'hg19' or 'hg38'
#' @param functional.regions character vector indicating the MNEMONIC chromatin states that will be retrieved from the Roadmap epigenomics. Default values are the active enhancers:'EnhA1', 'EnhA2'.
#' @param listOfEpigenomes character vector indicting which epigenome(s) to use for finding enhancers.
#' @param ProbeAnnotation GRange object of probe annotation.
#' @return a dataframe with enhancer probes and their chromosome coordinates
#' @keywords internal
#' @examples
#'\donttest{
#' met.platform = 'EPIC'
#' genome = 'hg38'
#' listOfEpigenomes = c('E034', 'E045', 'E047')
#' functional.regions = c('EnhA1', 'EnhA2', 'EnhG1', 'EnhG2')
#' df.enhancer.probes <-  getEnhancerProbes(met.platform = met.platform,
#'                                         genome = genome,
#'                                         functional.regions = functional.regions,
#'                                         listOfEpigenomes = listOfEpigenomes)
#'
#' }

getRoadMapEnhancerProbes <- function(met.platform = "EPIC", genome = "hg38", functional.regions = c("EnhA1",
    "EnhA2"), listOfEpigenomes = NULL, ProbeAnnotation) {

    # Step 1. Find all the filenames ending with
    # '_18_core_K27ac_hg38lift_mnemonics.bed.gz' from the Roadmap Epigenomics
    # web portal
    K27Ac_url <- "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/"
    urlData <- RCurl::getURL(K27Ac_url)
    urlData2 <- unlist(strsplit(urlData, "\\n"))
    filenames <- as.matrix(urlData2[grep("_18_core_K27ac_hg38lift_mnemonics.bed.gz",
        urlData2)])
    filenames <- unlist(strsplit(filenames, ">|<"))
    filenames <- filenames[grep("_18_core_K27ac_hg38lift_mnemonics.bed.gz", filenames)]
    filenames <- filenames[-grep("a href", filenames)]
    if (!is.null(listOfEpigenomes)) {
        listOfEpigenomes <- toupper(listOfEpigenomes)
        filenames <- filenames[grep(paste(listOfEpigenomes, collapse = "|"), filenames)]
    }

    # Step 2.Set the destination file directory and start to downloading files
    # dir = file.path(getwd(), 'ReferenceEpigenomes') dir.create(dir,
    # showWarnings = FALSE)
    dir <- tempdir()

    cat("Downloading chromatin states from the Roadmap Epigenomics...\n")
    for (file in filenames) {
        file.source <- paste0(K27Ac_url, file)
        destfile <- paste0(dir, "/", file)
        if (!file.exists(destfile)) {
            if (Sys.info()["sysname"] == "Windows")
                mode <- "wb" else mode <- "w"
            downloader::download(file.source, destfile = destfile, mode = mode)
        }
    }

    # Step 3. Get overlaps between probes and enhancers
    enhancerProbes <- character(0)
    for (file in filenames) {
        destfile <- paste0(dir, "/", file)
        genome.state <- as.data.frame(read.table(file = destfile, header = FALSE,
            sep = "\t", stringsAsFactors = FALSE, quote = ""))
        colnames(genome.state) <- c("chr", "start", "end", "state")
        genome.state$state <- vapply(genome.state$state, function(x) unlist(strsplit(x,
            "_"))[2], character(1))
        genome.state <- genome.state[genome.state$state %in% functional.regions,
            ]
        gr.enhancer <- GenomicRanges::makeGRangesFromDataFrame(genome.state, ignore.strand = TRUE,
            keep.extra.columns = TRUE)
        target.probes <- ProbeAnnotation[unique(S4Vectors::queryHits(IRanges::findOverlaps(ProbeAnnotation,
            gr.enhancer, ignore.strand = TRUE)))]
        this.epigenome <- unlist(strsplit(file, "_"))[1]
        cat("\tIdentifed", length(target.probes), "enhancer CpGs from the epigenome",
            this.epigenome, "\n")
        enhancerProbes <- c(enhancerProbes, names(target.probes))
        enhancerProbes <- unique(enhancerProbes)
    }
    # Step 4. Generate a dataframe for enhancer probes with their coordinates
    probe.ID <- names(ProbeAnnotation)
    probe.chr <- GenomicRanges::seqnames(ProbeAnnotation)
    probe.start.pos <- GenomicRanges::start(ranges(ProbeAnnotation))
    probe.start.end <- GenomicRanges::end(ranges(ProbeAnnotation))
    df.ProbeAnnotation <- data.frame(probeID = probe.ID, chr = probe.chr, start = probe.start.pos,
        end = probe.start.end)
    df.enhancer.probes <- df.ProbeAnnotation[df.ProbeAnnotation$probeID %in% enhancerProbes,
        ]
    return(df.enhancer.probes)
}


#' mapTranscriptToGene
#' @description map the miRNA precursor names to HGNC
#' @param transcripts vector with the name of miRNA precursors
#' @keywords internal
#'
#' @return a dataframe with two columns: 'Transcript' indicating the miRNA precursor names, 'Gene_name' indicating the actual human gene names (HGNC)
#'

mapTranscriptToGene <- function(transcripts) {

    # first, we need to get the precursor names for the mature miRNA
    # transcripts df.miRNA.mapping = data.frame(miRBaseConverter ::
    # miRNA_MatureToPrecursor(transcripts)) df.miRNA.mapping =
    # df.miRNA.mapping[!is.na(df.miRNA.mapping$Precursor), ]

    df.miRNA.mapping <- data.frame(Transcript = transcripts)
    df.miRNA.mapping <- df.miRNA.mapping[!is.na(df.miRNA.mapping$Transcript), , drop = FALSE]
    gene_name <- c()
    for (tr in transcripts) {
        elements <- family <- number <- name <- NULL
        elements <- unlist(strsplit(tr, "-"))
        family <- toupper(elements[2])
        # sort out family names for the let-7 miRNAs (hsa-let-7a -> MIRLET7)
        # from the MIR miRNAs ('hsa-miR-16' -> MIR16)
        if (family == "LET") {
            family <- "MIRLET"
        }
        # sort out the different cases for 'hsa-mir-25 -> MIR25' vs.
        # 'hsa-mir-24-2 -> MIR24-2'
        if (length(elements) == 3) {
            number <- elements[3]
        } else {
            # when the gene ID has an additional number after the slash, the
            # actual gene number will depend on whether the last digit is a
            # letter or a number if the last digit is a number, we will need a
            # slash in between, otherwise, we don't need a slash e.g.
            # hsa-mir-26a-1 => MIR26A1, hsa-mir-24-1 = > MIR24-1
            if (is.na(as.numeric(elements[3][length(elements[3])]))) {
                number <- paste0(elements[3], elements[4])
            } else {
                number <- paste0(elements[3], "-", elements[4])
            }
        }
        # we need upper case for the letter in the number ('hsa-miR-301b' ->
        # MIR301B)
        number <- toupper(number)
        name <- paste(family, number, sep = "")
        gene_name <- c(gene_name, name)
    }
    df.miRNA.mapping$Gene_name <- gene_name
    df.miRNA.mapping <- df.miRNA.mapping[!is.na(df.miRNA.mapping$Gene_name), , drop = FALSE]
    return(df.miRNA.mapping)
}


#' The functionEnrich function
#' @description Perform functional enrichment analysis for the differentially methylated genes occurring in the significant CpG-gene pairs.
#' @param EpiMixResults List of the result objects returned from the EpiMix function.
#' @param methylation.state character string indicating whether to use all the differentially methylated genes or only use the hypo- or hyper-methylated genes for enrichment analysis. Can be either 'all', 'Hyper' or 'Hypo'.
#' @param enrich.method character string indicating the method to perform enrichment analysis, can be either 'GO' or 'KEGG'.
#' @param ont character string indicating the aspect for GO analysis. Can be one of 'BP' (i.e., biological process), 'MF' (i.e., molecular function), and 'CC' (i.e., cellular component) subontologies, or 'ALL' for all three.
#' @param simplify boolean value indicating whether to remove redundancy of enriched GO terms.
#' @param cutoff if simplify is TRUE, this is the threshold for similarity cutoff of the ajusted p value.
#' @param pvalueCutoff adjusted pvalue cutoff on enrichment tests to report
#' @param pAdjustMethod one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none'
#' @param qvalueCutoff qvalue cutoff on enrichment tests to report as significant. Tests must pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on adjusted pvalues and iii) qvalueCutoff on qvalues to be reported.
#' @param save.dir path to save the enrichment table.
#' @return a clusterProfiler enrichResult instance
#' @export
#' @examples
#' \donttest{
#' library(clusterProfiler)
#' library(org.Hs.eg.db)
#'
#' data(Sample_EpiMixResults_Regular)
#'
#' enrich.results <- function.enrich(
#'  EpiMixResults = Sample_EpiMixResults_Regular,
#'  enrich.method = 'GO',
#'  ont = 'BP',
#'  simplify = TRUE,
#'  save.dir = ''
#' )
#' }
#'
functionEnrich <- function(EpiMixResults, methylation.state = "all", enrich.method = "GO",
    ont = "BP", simplify = TRUE, cutoff = 0.7, pvalueCutoff = 0.05, pAdjustMethod = "BH",
    qvalueCutoff = 0.2, save.dir = ".") {

    if (!requireNamespace("clusterProfiler")) {
        message("This function requires the 'clusterProfiler' package.")
        return(invisible())
    }

    if (!requireNamespace("org.Hs.eg.db")) {
        message("This function requires the 'org.Hs.eg.db' package.")
        return(invisible())
    }

    # input check
    if (enrich.method != "GO" & enrich.method != "KEGG") {
        stop("enrich.method must be 'GO' or 'KEGG' !")
    }
    if (methylation.state != "all" & methylation.state != "Hypo" & methylation.state !=
        "Hyper") {
        stop("methylaiton.state must be 'all', 'Hyper', or 'Hypo' !")
    }

    FunctionalPairs <- EpiMixResults$FunctionalPairs
    if (methylation.state == "Hypo") {
        FunctionalPairs <- FunctionalPairs[FunctionalPairs$State == "Hypo", ]
    } else if (methylation.state == "Hyper") {
        FunctionalPairs <- FunctionalPairs[FunctionalPairs$State == "Hyper", ]
    }
    df.gene <- FunctionalPairs %>%
        dplyr::select(.data$Gene, .data$`Fold change of gene expression`) %>%
        dplyr::group_by(.data$Gene) %>%
        dplyr::summarize(Avg.expr = mean(.data$`Fold change of gene expression`,
            na.rm = TRUE))

    gene_id_map <- AnnotationDbi::select(org.Hs.eg.db, keys = df.gene$Gene, columns = c("ENTREZID"),
        keytype = "SYMBOL")
    if (sum(is.na(gene_id_map$ENTREZID)) > 0) {
        absentGenes <- gene_id_map$SYMBOL[which(is.na(gene_id_map$ENTREZID))]
        warning(paste0("Can not find ENTREZID for ", length(absentGenes), " genes.",
            "These genes can not be included in the enrichment analysis:\n", paste0(absentGenes,
                collapse = ",")))
    }
    colnames(gene_id_map) <- c("Gene", "ENTREZID")
    df.gene <- merge(df.gene, gene_id_map)
    gene.vector <- df.gene$Avg.expr
    names(gene.vector) <- df.gene$ENTREZID
    if (enrich.method == "GO") {
        ego <- clusterProfiler::enrichGO(names(gene.vector), OrgDb = "org.Hs.eg.db",
            ont = ont, pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
            qvalueCutoff = qvalueCutoff, readable = TRUE)
    } else if (enrich.method == "KEGG") {
        ego <- clusterProfiler::enrichKEGG(names(gene.vector), pvalueCutoff = pvalueCutoff,
            pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff)
    }
    if (simplify & enrich.method == "GO") {
        ego <- clusterProfiler::simplify(ego, cutoff = cutoff)
    }

    if (save.dir != "" & length(save.dir) > 0) {
        save.file.name <- paste0(save.dir, "/", "FunctionEnrichment_", enrich.method,
            "_", methylation.state, ".csv")
        utils::write.csv(ego@result, save.file.name, row.names = FALSE)
    }
    return(ego)
}

# Unregister sockets, in case the socket was not completely removed from last
# run unregister <- function() { env <- foreach::.foreachGlobals
# rm(list=ls(name=env), pos=env) }

# Function to retrieve data from the data package
EpiMix_GetData <- function(...) {
    e <- new.env()
    name <- utils::data(..., package = "EpiMix.data", envir = e)[1]
    e[[ls(envir = e)[1]]]
}



listEpigenomces <- function() {
    EpigenomeMap <- EpiMix_GetData("EpigenomeMap")
    print(EpigenomeMap)
}


listChromatinStates<- function() {
    chromatin.states <- c(EnhG1 = "Genic enhancer1", EnhG2 = "Genic enhancer2", EnhA1 = "Active Enhancer 1",
        EnhA2 = "Active Enhancer 2", EnhWk = "Weak Enhancer", EnhBiv = "Bivalent Enhancer")
    get("chromatin.states")
}

#' The validEpigenomes function
#' @description check user input for roadmap epigenome groups or ids
#' @param roadmap.epigenome.groups epigenome groups
#' @param roadmap.epigenome.ids epigenome ids
#'
#' @return a character vector of selected epigenome ids
#'
validEpigenomes <- function(roadmap.epigenome.groups, roadmap.epigenome.ids) {
    EpigenomeMap <- EpiMix_GetData("EpigenomeMap")
    selectedEpigenomes <- character(0)
    if (!is.null(roadmap.epigenome.groups)) {
        if ("ES_deriv." %in% roadmap.epigenome.groups)
            roadmap.epigenome.groups[which(roadmap.epigenome.groups == "ES_deriv.")] <- "ES_deriv"
        if ("Mesench." %in% roadmap.epigenome.groups)
            roadmap.epigenome.groups[which(roadmap.epigenome.groups == "Mesench.")] <- "Mesench"
        if ("Myosat." %in% roadmap.epigenome.groups)
            roadmap.epigenome.groups[which(roadmap.epigenome.groups == "Myosat.")] <- "Myosat"
        if ("Neurosph." %in% roadmap.epigenome.groups)
            roadmap.epigenome.groups[which(roadmap.epigenome.groups == "Neurosph.")] <- "Neurosph"

        overlapGroups <- intersect(names(EpigenomeMap), roadmap.epigenome.groups)
        if (length(overlapGroups) != length(roadmap.epigenome.groups)) {
            warning(paste0("Can not find epigenome group(s): ", roadmap.epigenome.groups[which(!roadmap.epigenome.groups %in%
                overlapGroups)]), immediate. = TRUE)
            cat("Available epigenome groups: \n", paste0(names(EpigenomeMap), collapse = ",  "),
                "\n")
        }
        for (gr in overlapGroups) {
            selectedEpigenomes <- c(selectedEpigenomes, EpigenomeMap[[gr]])
        }
    }
    if (!is.null(roadmap.epigenome.ids)) {
        all.ids <- character(0)
        for (i in 1:length(EpigenomeMap)) {
            all.ids <- c(all.ids, EpigenomeMap[[i]])
        }
        overlapIDs <- intersect(roadmap.epigenome.ids, all.ids)
        if (length(overlapIDs) != length(roadmap.epigenome.ids)) {
            warning(paste0("Can not find some of the epigenome ids: ", roadmap.epigenome.ids[which(!roadmap.epigenome.ids %in%
                overlapIDs)]), immediate. = TRUE)
            cat("Available epigenome ids: \n", paste0(all.ids[order(all.ids)], collapse = ", "),
                "\n")
        }
        selectedEpigenomes <- c(selectedEpigenomes, overlapIDs)
    }
    selectedEpigenomes <- unique(selectedEpigenomes)
    return(selectedEpigenomes)
}












