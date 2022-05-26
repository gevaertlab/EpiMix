###################################################################################################################
#                                  Functions to get adjacent genes of CpGs
###################################################################################################################

#'@importFrom dplyr slice left_join group_by_ full_join filter do
NULL

#'@import ELMER.data
NULL

# NearGenes
# @param Target A charactor which is name of TRange or one of rownames of TBed.
# @param Gene A GRange object contains coordinates of promoters for human genome.
# @param geneNum A number determine how many gene will be collected from each
# side of target (number shoule be even).
# @param TRange A GRange object contains coordinate of targets.
# @return A data frame of nearby genes and information: genes' IDs, genes' symbols,
# distance with target and side to which the gene locate to the target.
#'@importFrom GenomicRanges strand<-
NearGenes <- function (Target = NULL,
                       Gene = NULL,
                       geneNum = 20,
                       TRange = NULL){
  # Algorithm:
  # 1) get the follow gene (overlapping genes are diconsidered) to be the first in L1 (index variable)
  #                 probe
  #     -------       O      ------
  #    |      |      ||     |     |
  #    -------       ||     ------
  #  follow gene          precede gene
  # 2) Sort genes (by start)
  # 3.1) Get 9 genes before index (L1)
  #      If we only have l genes (l < 9) due to end of genomic region, get the l ones and get more 10-l to the right
  # 3.2) Get 10 after index (L1)
  #      If we only have r genes (r < 10) due to end of genomic region, get the r ones and get more 10-r to the left
  # Where 10 is genum/2
  Gene$GENEID <- Gene$ensembl_gene_id
  if("external_gene_name" %in% colnames(S4Vectors::mcols(Gene))){
    Gene$SYMBOL <- Gene$external_gene_name
  } else  if("external_gene_id" %in% colnames(S4Vectors::mcols(Gene))){
    Gene$SYMBOL <- Gene$external_gene_id
  } else {
    stop("No gene symbol column found (expected external_gene_id or external_gene_name")
  }
  if(is.null(Gene) | is.null(Target)){
    stop ("Target and Genes should both be defined")
  }
  if(is.null(TRange)){
    stop( "TRange must be defined")
  }else{
    # Just to be sure we have only one probe. To be removed ?
    regionInfo <- TRange[names(TRange) %in% Target]
  }
  GeneIDs <- c()
  Distances <- c()
  strand(Gene) <- "*"
  # We will get only genes in the same same chromossome
  Gene <- Gene[as.character(seqnames(Gene)) %in% as.character(seqnames(regionInfo))]
  if(length(Gene)==0){
    warning(paste0(Target," don't have any nearby gene in the given gene list."))
    Final <- NA
  } else {
    Gene <- sort(Gene)
    index <- follow(regionInfo,Gene)
    #left side
    Leftlimit <- geneNum/2
    Rightlimit <- geneNum/2
    n <- 1
    if(is.na(index)){
      index<- 0
      Leftlimit <- 0
      Left <- c()
    }else if(index==1){
      Left <- index
      Leftlimit <- length(Left)
    }else{
      Left <- index
      while(length(Left) < Leftlimit){
        # If the gene is not in the list already add it, otherwise go to the next
        if(!as.character(Gene$GENEID[index-n]) %in% as.character(Gene$GENEID[Left])) Left <- c((index-n),Left)

        # Is it the first gene? If so there is nothing in the left anymore
        if((index-n)==1) Leftlimit <- length(Left)
        n <- n + 1
      }
    }

    Right <- c()
    n <- 1
    if(index==length(Gene) ||
       all(unique(Gene$GENEID[(index+1):length(Gene)]) %in% as.character(Gene$GENEID[index]))){
      Rightlimit <- length(Right)
    }else{
      while(length(Right) < Rightlimit){
        if(!as.character(Gene$GENEID[index+n])%in% as.character(Gene$GENEID[c(Right,Left)]))
          Right <- c(Right,(index+n))

        if(index+n==length(Gene)){
          Rightlimit <- length(Right)
        } else{
          n <- n+1
        }
      }
    }

    if(Rightlimit < geneNum/2){
      n <- 1
      if(Left[1]-n > 0){
        while((length(Left)+length(Right)) < geneNum){
          if(!as.character(Gene$GENEID[Left[1]-n])%in%as.character(Gene$GENEID[c(Left,Right)]))
            Left <- c((Left[1]-n),Left)
          n <- n+1
        }
      }
    }

    if(Leftlimit < geneNum/2){
      n <- 1
      m <- length(Right)
      if(Right[m]+n < length(Gene)+1)
        while((length(Left)+length(Right)) < geneNum){
          if(!as.character(Gene$GENEID[Right[m]+n])%in%as.character(Gene$GENEID[c(Left,Right)]))
            Right <- c(Right,(Right[m]+n))
          n <- n+1
        }
    }
    Whole <- c(Left,Right)
    GeneIDs <- Gene$GENEID[Whole]
    Symbols <- Gene$SYMBOL[Whole]
    #Distances <-  suppressWarnings(distance(Gene[Whole],regionInfo))
    Distances <-  distance(Gene[Whole],regionInfo)
    if(Rightlimit < 1){
      Sides <- paste0("L",length(Left):1)
    } else if( Leftlimit < 1){
      Sides <- paste0("R",1:length(Right))
    } else{
      Sides <- c(paste0("L",length(Left):1),paste0("R",1:length(Right)))
    }

    Final <- data.frame(Target=rep(Target,length(GeneIDs)),GeneID=GeneIDs,
                        Symbol=Symbols,Distance=Distances, Side=Sides,
                        stringsAsFactors = FALSE)
    Final <- Final[order(Final$Side,Final$Distance),]
  }
  return(Final)
}

#' GetNearGenes to collect nearby genes for one locus.
#' @description
#' GetNearGenes is a function to collect equal number of gene on each side of one locus.
#' It can receite either multi Assay Experiment with both DNA methylation and gene Expression matrix
#' and the names of probes to select nearby genes, or it can receive two granges objects TRange and geneAnnot.
#' @param data A multi Assay Experiment with both DNA methylation and gene Expression objects
#' @param probes Name of probes to get nearby genes (it should be rownames of the DNA methylation
#' object in the data argument object)
#' @param geneAnnot A GRange object  or Summarized Experiment object that contains coordinates of promoters for
#' human genome.
#' @param numFlankingGenes A number determines how many gene will be collected totally.
#' Then the number devided by 2 is the number of genes collected from
#' each side of targets (number shoule be even) Default to 20.
#' @param TRange A GRange object or Summarized Experiment object that contains coordinates of a list of targets loci.
#' @return A data frame of nearby genes and information: genes' IDs, genes' symbols,
#' distance with target and side to which the gene locate to the target.
#' @importFrom GenomicRanges strand follow distance
#' @importFrom plyr alply
#' @importFrom doParallel registerDoParallel
#' @importFrom SummarizedExperiment rowRanges
#' @references
#' Yao, Lijing, et al. "Inferring regulatory element landscapes and transcription
#' factor networks from cancer methylomes." Genome biology 16.1 (2015): 1.
GetNearGenes <- function(data = NULL,
                         probes = NULL,
                         geneAnnot = NULL,
                         TRange = NULL,
                         numFlankingGenes = 20){
  message("Searching for the ", numFlankingGenes, " near genes")
  # Comment out the original part from the ELMER function, since EpiMix does not use this condition

  # if(!is.null(data)){
  #   if(is.null(probes)) stop("Please set the probes argument (names of probes to select nearby genes)")
  #   TRange <- subset(getMet(data), rownames(getMet(data)) %in% probes)
  #   geneAnnot <- getExp(data)
  # }
  if(is.null(TRange)){
    stop("TRange must be defined")
  }
  tssAnnot <- NULL

  # Comment out the original portion of code from ELMER, since EpiMix does not use this

  # if(is.null(geneAnnot)){
  #   if("genome" %in% names(metadata(data))){
  #     genome <- metadata(data)$genome
  #     tssAnnot <- getTSS(genome = genome)
  #     geneAnnot <- get.GRCh(genome = genome,as.granges = TRUE)
  #   }
  # }

  # if(class(TRange) == class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
  #   TRange <- rowRanges(TRange)
  # }
  # if(class(geneAnnot) == class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
  #   geneAnnot <- rowRanges(geneAnnot)
  # }

  if(is.null(names(TRange))) {
    if(is.null(TRange$name)) stop("No probe names found in TRange")
    names(TRange) <- TRange$name
  }

  NearGenes <-
    getRegionNearGenes(
      numFlankingGenes = numFlankingGenes,
      geneAnnot = geneAnnot,
      tssAnnot = tssAnnot,
      TRange = TRange
    )
  return(NearGenes)
}


#' @title Calculate the distance between probe and gene TSS
#' @description Calculate the distance between probe and gene TSS
#' @param data A multi Assay Experiment with both DNA methylation and gene Expression objects
#' @param NearGenes A list or a data frame with the pairs gene probes
#' @param cores Number fo cores to be used. Deafult: 1
#' @param met.platform DNA methyaltion platform to retrieve data from: EPIC or 450K (default)
#' @param genome Which genome build will be used: hg38 (default) or hg19.
#' @return a dataframe of nearest genes with distance to TSS.
addDistNearestTSS <- function(
  data,
  NearGenes,
  genome,
  met.platform,
  cores = 1
) {

  if(missing(NearGenes)) stop("Please set NearGenes argument")

  # used to recover TSS information
  if(missing(data) & missing(genome)) {
    stop("Please set data argument or genome arguments")
  }

  # For a given probe/region and gene find nearest TSS distance
  if(!missing(data)){
    # Comment out the original portion of code from ELMER, since EpiMix does not use this
    #tss <- getTSS(metadata(data)$genome)
  } else {
    tss <- getTSS(genome = genome)
  }

  message("Update the distance to gene to distance to the nearest TSS of the gene")

  # If our input has the probe names we will have to recover the probe metadata to map.
  region <- FALSE
  if(!missing(data)){
    # Comment out the original portion from ELMER, since EpiMix does not use this function
    #met <- rowRanges(getMet(data))
  } else if(!missing(met.platform)){
    met <- EpiMix_getInfiniumAnnotation(plat = met.platform, genome = genome)
  } else {
    region <- TRUE
    met <- NearGenes %>% tidyr::separate("Target", c("seqnames","start","end"),
                                         sep = ":|-", remove = FALSE,
                                         convert = FALSE, extra = "warn", fill = "warn") %>%
      makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  }

  NearGenes <- calcDistNearestTSS(links = NearGenes,TRange = met,tssAnnot = tss)

  return(NearGenes)
}

#' @title Calculate distance from region to nearest TSS
#' @description
#' Idea
#' For a given region R linked to X genes G
#' merge R with nearest TSS for G (multiple)
#' this will increse nb of lines
#' i.e R1 - G1 - TSS1 - DIST1
#'    R1 - G1 - TSS2 - DIST2
#'    To vectorize the code:
#' make a granges from left and onde from right
#' and find distance
#' collapse the results keeping min distance for equals values
#' @param links Links to calculate the distance
#' @param TRange Genomic coordinates for Tartget region
#' @param tssAnnot TSS annotation
#' @importFrom dplyr slice left_join group_by_
#' @return dataframe of genomic distance from TSS
#' @author Tiago C. Silva
calcDistNearestTSS <- function(
  links,
  TRange,
  tssAnnot
){

  message("calculating Distance to nearest TSS")
  if(!is(tssAnnot,"GenomicRanges")){
    stop("tssAnnot is not a GenomicRanges")
  }

  if(!is(TRange,"GenomicRanges")){
    stop("tssAnnot is not a GenomicRanges")
  }

  if(!"ID" %in% colnames(values(TRange))){
    TRange$ID <- names(TRange)
  }

  if(!"ensembl_gene_id" %in% colnames(links)){
    colnames(links)[grep("GeneID", colnames(links))] <- "ensembl_gene_id"
  }

  merged <- dplyr::left_join(
    links,
    tibble::as_tibble(tssAnnot),
    #suppressWarnings(tibble::as_tibble(tssAnnot)),
    by = c("ensembl_gene_id")
  )

  merged <- dplyr::left_join(
    merged,
    #suppressWarnings(tibble::as_tibble(TRange)),
    tibble::as_tibble(TRange),
    by = c("ID")
  )

  # In case a gene was removed from newer versions and not mapped
  merged <- merged[!is.na(merged$transcription_start_site),]

  left <- makeGRangesFromDataFrame(
    merged,
    start.field = "transcription_start_site",
    end.field = "transcription_start_site",
    seqnames.field = "seqnames.x",
    strand.field = "strand.x",
    ignore.strand = FALSE
  )

  right <- makeGRangesFromDataFrame(
    merged,
    start.field = "start.y",
    end.field = "end.y",
    strand.field = "strand.y",
    seqnames.field = "seqnames.y",
    ignore.strand = FALSE
  )

  merged$DistanceTSS <- distance(left,right,ignore.strand = TRUE)
  merged <- unique(merged[,c("ID","ensembl_gene_id","DistanceTSS")])

  ret <- merged %>%
    dplyr::group_by(.data$ID,.data$ensembl_gene_id) %>%
    dplyr::slice(which.min(.data$DistanceTSS))

  #ret <- ret[match(links %>% tidyr::unite(ID,Target,GeneID) %>% pull(ID),
  #          ret %>% tidyr::unite(ID,Target,GeneID) %>% pull(ID)),]
  links <- dplyr::full_join(links,ret)
  colnames(links)[1:3] <- c("ID","GeneID","Symbol")
  return(links)
}

#' @title  Identifies nearest genes to a region
#' @description
#'  Auxiliary function for GetNearGenes
#'  This will get the closest genes (n=numFlankingGenes) for a target region (TRange)
#'  based on a genome of refenrece gene annotation (geneAnnot). If the
#'  transcript level annotation (tssAnnot) is provided the Distance will be updated to
#'  the distance to the nearest TSS.
#' @param geneAnnot A GRange object contains gene coordinates of for human genome.
#' @param tssAnnot A GRange object contains tss coordinates of for human genome.
#' @param numFlankingGenes A number determine how many gene will be collected from each
# side of target (number shoule be even).
#' @param TRange A GRange object contains coordinate of targets.
#' @return A data frame of nearby genes and information: genes' IDs, genes' symbols,
# distance with target and side to which the gene locate to the target.
#' @importFrom GenomicRanges nearest precede follow
#' @importFrom tibble as_tibble
#' @importFrom dplyr group_by do filter
#' @importFrom progress progress_bar
#' @author
#' Tiago C Silva (maintainer: tiagochst@usp.br)
getRegionNearGenes <- function(TRange = NULL,
                               numFlankingGenes = 20,
                               geneAnnot = NULL,
                               tssAnnot = NULL){


  pb <- progress::progress_bar$new(total = numFlankingGenes * 2)

  TRange$ID <- names(TRange)

  # We will consider the input at a gene level only
  if(! "ensembl_gene_id" %in%  colnames(S4Vectors::mcols(geneAnnot))){
    message("geneAnnot needs the following column ensembl_gene_id")
  }
  geneAnnot <- geneAnnot[!duplicated(geneAnnot$ensembl_gene_id)]


  # Optimized version
  # Idea: vectorize search
  # 1) For all regions, get nearest gene
  # 2) check follow and overlapping genes recursively
  # 3) check precede and overlapping genes recursively
  # 4) map the positions based on min distance (L1)
  # The input data has to be at gene level and not transcript which would broke
  # some of the optimizations for which we remove the genes already evaluated
  all <- 1:length(TRange)
  nearest.idx <-
    nearest(TRange,
            geneAnnot,
            select = "all",
            ignore.strand = TRUE)
  idx <- tibble::as_tibble(nearest.idx)
  evaluating <- idx$queryHits
  ret <-
    cbind(
      tibble::as_tibble(geneAnnot[idx$subjectHits]),
      #suppressWarnings(tibble::as_tibble(geneAnnot[idx$subjectHits])),
      tibble::tibble(
        "ID" = names(TRange)[idx$queryHits],
        "Distance" = distance(TRange[idx$queryHits],
                              geneAnnot[idx$subjectHits], select = "all", ignore.strand = TRUE) *
          ifelse(start(TRange[evaluating]) < start(geneAnnot[idx$subjectHits]), 1,-1)

      )
    )
  for (i in 1:(numFlankingGenes)) {
    idx <-
      unique(rbind(
      tibble::as_tibble(
          findOverlaps(
            geneAnnot[idx$subjectHits],
            geneAnnot,
            ignore.strand = TRUE,
            type = "any",
            select = "all"
          )
        ),
       tibble::as_tibble(
          precede(
            geneAnnot[idx$subjectHits],
            geneAnnot,
            select = "all",
            ignore.strand = TRUE
          )
        )
      ))
    idx$evaluating <-  evaluating[idx$queryHits]
    # remove same target gene and probe if counted twice
    idx <- idx[!duplicated(idx[, 2:3]), ]

    # todo remove already evaluated previously (we don't wanna do it again)
    idx <-
      idx[!paste0(geneAnnot[idx$subjectHits]$ensembl_gene_id, names(TRange)[idx$evaluating]) %in% paste0(ret$ensembl_gene_id, ret$ID), ]
    evaluating <- evaluating[idx$queryHits]
    ret <-
      rbind(ret, # keep old results
            cbind(
                tibble::as_tibble(geneAnnot[idx$subjectHits]),
              tibble::tibble(
                "ID" = names(TRange)[evaluating],
                "Distance" = ifelse(start(TRange[evaluating]) < start(geneAnnot[idx$subjectHits]), 1,-1) *
                  distance(TRange[evaluating],
                           geneAnnot[idx$subjectHits], select = "all",
                           ignore.strand = TRUE)
              )
            ))
    pb$tick()
  }
  ret <- ret[!duplicated(ret[,c("ensembl_gene_id","ID")]),]

  idx <- tibble::as_tibble(nearest.idx)
  evaluating <- idx$queryHits
  for (i in 1:(numFlankingGenes)) {
    idx <-
      unique(rbind(
      tibble::as_tibble(
          findOverlaps(
            geneAnnot[idx$subjectHits],
            geneAnnot,
            ignore.strand = TRUE,
            type = "any",
            select = "all"
          )
       ),
       tibble::as_tibble(
          follow(
            geneAnnot[idx$subjectHits],
            geneAnnot,
            select = "all",
            ignore.strand = TRUE
          )
        )
      ))
    idx$evaluating <-  evaluating[idx$queryHits]
    idx <- idx[!duplicated(idx[, 2:3]), ]
    idx <-
      idx[!paste0(geneAnnot[idx$subjectHits]$ensembl_gene_id, names(TRange)[idx$evaluating]) %in% paste0(ret$ensembl_gene_id, ret$ID), ]
    evaluating <- evaluating[idx$queryHits]
    ret <-
      rbind(ret, cbind(
        tibble::as_tibble(geneAnnot[idx$subjectHits]),
        tibble::tibble(
          "ID" = names(TRange)[evaluating],
          "Distance" = ifelse(start(TRange[evaluating]) < start(geneAnnot[idx$subjectHits]), 1,-1) *
            distance(TRange[evaluating],
                     geneAnnot[idx$subjectHits], select = "all",ignore.strand = TRUE)
        )
      ))
    pb$tick()
  }
  ret <- ret[!duplicated(ret[,c("ensembl_gene_id","ID")]),]
  ret <- ret[order(ret$Distance),]

  ret <- ret[, c("ID",
                 "ensembl_gene_id",
                 grep("external_gene_", colnames(ret), value = TRUE),
                 "Distance")]

  message("Identifying gene position for each probe")
  f <- function(x) {
    center <- which(abs(x$Distance) == min(abs(x$Distance)))[1]
    pos <- setdiff(-center:(nrow(x) - center), 0)
    x$Side <- ifelse(pos > 0, paste0("R", abs(pos)), paste0("L", abs(pos)))
    out <-  x %>% dplyr::filter(x$Side %in% c(paste0("R", 1:(numFlankingGenes / 2)),
                                              paste0("L", 1:(numFlankingGenes / 2))
    ))
    if (nrow(out) < numFlankingGenes) {
      if (paste0("R", floor(numFlankingGenes / 2)) %in% out$Side) {
        cts <- length(grep("L", sort(x$Side), value = TRUE))
        out <- x %>% dplyr::filter(x$Side %in% c(paste0("R", 1:(numFlankingGenes - cts)),
                                               grep("L", sort(out$Side), value = TRUE)))
      } else {
        cts <- length(grep("R", sort(x$Side), value = TRUE))
        out <- x %>% dplyr::filter(x$Side %in%
                                     c(paste0("L", 1:(numFlankingGenes - cts)),
                                       grep("R", sort(out$Side), value = TRUE))
        )
      }
    }
    out <- out[order(out$Distance), ]
    return(out)
  }
  ret <- ret %>% group_by(.data$ID) %>% do(f(.))

  if (!is.null(tssAnnot)) {
    message("Calculating distance to nearest TSS")
    ret <- calcDistNearestTSS(
      links = ret,
      TRange = TRange,
      tssAnnot = tssAnnot
    )
  }

  if(any(grepl("external_gene_", colnames(ret)))){
    colnames(ret)[1:3] <- c("ID", "GeneID", "Symbol")
  } else {
    colnames(ret)[1:2] <- c("ID", "GeneID")
  }
  pb$terminate()
  return(ret)
}


#' getTSS to fetch GENCODE gene annotation (transcripts level) from Bioconductor package biomaRt
#' If upstream and downstream are specified in TSS list, promoter regions of GENCODE gene will be generated.
#' @description
#' getTSS to fetch GENCODE gene annotation (transcripts level) from Bioconductor package biomaRt
#' If upstream and downstream are specified in TSS list, promoter regions of GENCODE gene will be generated.
#' @param TSS A list. Contains upstream and downstream like TSS=list(upstream, downstream).
#'  When upstream and downstream is specified, coordinates of promoter regions with gene annotation will be generated.
#' @param genome Which genome build will be used: hg38 (default) or hg19.
#' @return GENCODE gene annotation if TSS is not specified. Coordinates of GENCODE gene promoter regions if TSS is specified.
#' @author Lijing Yao (maintainer: lijingya@usc.edu)
#' @import GenomeInfoDb
#' @importFrom GenomicFeatures transcripts
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom biomaRt useEnsembl
getTSS <- function(genome = "hg38",
                   TSS = list(upstream = NULL, downstream = NULL)){

  if (tolower(genome) == "hg38") {
    tss <- getdata("Human_genes__GRCh38_p12__tss")
  } else {
    tss <- getdata("Human_genes__GRCh37_p13__tss")
  }

  tss$chromosome_name <-  paste0("chr", tss$chromosome_name)
  tss$strand[tss$strand == 1] <- "+"
  tss$strand[tss$strand == -1] <- "-"
  tss <- makeGRangesFromDataFrame(tss,
                                  strand.field = "strand",
                                  start.field = "transcript_start",
                                  end.field = "transcript_end",
                                  keep.extra.columns = TRUE)

  if (!is.null(TSS$upstream) & !is.null(TSS$downstream))
    tss <- promoters(tss, upstream = TSS$upstream, downstream = TSS$downstream)
  return(tss)
}


getdata <- function(...)
{
  e <- new.env()
  name <- data(..., package = "ELMER.data",envir = e)[1]
  e[[ls(envir = e)[1]]]
}


#' @title getFeatureProbe to select probes within promoter regions or distal regions.
#' @description
#' getFeatureProbe is a function to select the probes falling into
#' distal feature regions or promoter regions.
#' @importFrom GenomicRanges promoters
#' @description This function selects the probes on HM450K that either overlap
#' distal biofeatures or TSS promoter.
#' @param promoter A logical.If TRUE, function will ouput the promoter probes.
#' If FALSE, function will ouput the distal probes overlaping with features. The
#' default is FALSE.
#' @param met.platform DNA methyaltion platform to retrieve data from: EPIC or 450K (default)
#' @param genome Which genome build will be used: hg38 (default) or hg19.
#' @param feature A GRange object containing biofeature coordinate such as
#' enhancer coordinates.
#' If NULL only distal probes (2Kbp away from TSS will be selected)
#' feature option is only usable when promoter option is FALSE.
#' @param TSS A GRange object contains the transcription start sites. When promoter is FALSE, Union.TSS
#' in \pkg{ELMER.data} will be used for default. When promoter is TRUE, UCSC gene TSS will
#' be used as default (see detail). User can specify their own preference TSS annotation.
#' @param TSS.range A list specify how to define promoter regions.
#' Default is upstream =2000bp and downstream=2000bp.
#' @param rm.chr A vector of chromosome need to be remove from probes such as chrX chrY or chrM
#' @return A GRange object containing probes that satisfy selecting critiria.
#' @importFrom S4Vectors queryHits subjectHits
#' @import GenomicRanges
#' @details
#'  In order to get real distal probes, we use more comprehensive annotated TSS by both
#'  GENCODE and UCSC. However, to get probes within promoter regions need more
#'  accurate annotated TSS such as UCSC. Therefore, there are different settings for
#'  promoter and distal probe selection. But user can specify their own favorable
#'  TSS annotation. Then there won't be any difference between promoter and distal
#'  probe selection.
#'  @return A GRanges object contains the coordinate of probes which locate
#'  within promoter regions or distal feature regions such as union enhancer from REMC and FANTOM5.
#'  @usage getFeatureProbe(feature,
#'                           TSS,
#'                           TSS.range = list(upstream = 2000, downstream = 2000),
#'                           promoter = FALSE, rm.chr = NULL)
#'
getFeatureProbe <- function(feature = NULL,
                              TSS,
                              genome = "hg38",
                              met.platform = "HM450",
                              TSS.range = list(upstream = 2000, downstream = 2000),
                              promoter = FALSE,
                              rm.chr = NULL){
  probe <- EpiMix_getInfiniumAnnotation(toupper(met.platform),genome)
  # We will remove the rs probes, as they should not be used in the analysis
  probe <- probe[!grepl("rs",names(probe)),]
  probe <- probe[!probe$MASK_general,] # remove masked probes   #766,499
  if(!is.null(rm.chr)) probe <- probe[!as.character(seqnames(probe)) %in% rm.chr]  # 749312

  if(missing(TSS)){
    # The function getTSS gets the transcription coordinantes from Ensemble (GENCODE)
    TSS <- getTSS(genome = genome)   # 208,423
  }
    promoters <- promoters(TSS,
                           upstream = TSS.range[["upstream"]],
                           downstream = TSS.range[["downstream"]])

  if(!promoter){
    probe <- probe[setdiff(1:length(probe),unique(queryHits(findOverlaps(probe,promoters,ignore.strand=TRUE))))]  #598,225


    if(is.null(feature)) {
      message("Returning distal probes: ", length(probe))
      return(probe)
    }
    if(is(feature,"GRanges")) {
      probe <- probe[unique(queryHits(findOverlaps(probe,feature)))]
      message("Returning distal probes overlapping with features: ", length(probe))

    } else {
      stop("feature is not GRanges object.")
    }
  } else {
    probe <- probe[unique(queryHits(findOverlaps(probe,promoters,ignore.strand=TRUE)))]
    message("Returning promoter probes: ", length(probe))
  }
  return(probe)
}














