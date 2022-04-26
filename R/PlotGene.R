#' The EpiMix_PlotGene function
#' @description plot the genomic coordinate, DM values and chromatin state for each CpG probe of a specific gene.
#' @details this function requires R package dependencies: karyoploteR, TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.eg.db
#' @param gene.name character string indicating the name of the gene to be plotted.
#' @param EpiMixResults the resulting list object returned from the function of EpiMix.
#' @param met.platform character string indicating the type of the microarray where the DNA methylation data were collected. The value should be either "HM27", "HM450" or "EPIC". Default: "HM450"
#' @param roadmap.epigenome.id character string indicating the epigenome id (EID) for a reference tissue or cell type. Default: "E002"
#' @param left.gene.margin  numeric value indicating the number of extra nucleotide bases to be plotted on the left side of the target gene. Default: 10000.
#' @param right.gene.margin numeric value indicating the number of extra nucleotide bases to be plotted on the right side of the target gene. Default: 10000.
#' @param gene.name.font numeric value indicating the font size for the gene name. Default: 0.7.
#' @param plot.transcripts logic indicating whether to plot each individual transcript of the gene. Default: TRUE. If False, the gene will be plotted with a single rectangle, without showing the structure of individual transcripts.
#' @param plot.transcripts.structure logic indicating whether to plot the transcript structure (introns and exons). Non-coding exons are shown in green and the coding exons are shown in red. Default: TRUE.
#' @param show.probe.name logic indicating whether to show the name(s) for each differentially methylated CpG probe. Default: TRUE
#' @param probe.name.font numeric value indicating the font size of the name(s) for the differentially methylated probe(s) in pixels. Default: 0.6.
#' @param y.label.font font size of the y axis label
#' @param y.label.margin distance between y axis label and y axis
#' @param axis.number.font font size of axis ticks and numbers
#' @param chromatin.label.font font size of the labels of the histone proteins
#' @param chromatin.label.margin distance between the histone protein labels and axis
#' @importFrom biomaRt useDataset getBM useEnsembl
#' @importFrom GenomicRanges makeGRangesFromDataFrame seqnames start end mcols
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom RCurl getURL
#' @importFrom RColorBrewer brewer.pal
#' @return plot of the genomic coordinate, DM values and chromatin state for each CpG probe of a specific gene.
#' @details
#' roadmap.epigenome.id: since the chromatin state is tissue or cell-type specific,
#' EpiMix needs to know the reference tissue or cell type in order to retrieve the proper DNase-seq and histone ChIP-seq data.
#' Available epigenome ids can be obtained from the Roadmap Epigenomic study (Nature, PMID: 25693563, figure 2).
#' They can also be retrieved from the list.epigenomes() function.
#' @export
#' @examples
#' \dontrun{
#' library(karyoploteR)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(org.Hs.eg.db)
#' library(regioneR)
#'
#' data(Sample_EpiMixResults_Regular)
#'
#' gene.name = "CCND2"
#'
#' roadmap.epigenome.id = "E096"
#'
#' EpiMix_PlotGene(gene.name = gene.name,
#'                 EpiMixResults = Sample_EpiMixResults_Regular,
#'                 met.platform = "HM450",
#'                 roadmap.epigenome.id = roadmap.epigenome.id)
#' }
#'
EpiMix_PlotGene <- function(gene.name,
                      EpiMixResults,
                      met.platform = "HM450",
                      roadmap.epigenome.id = "E002",
                      left.gene.margin = 10000,
                      right.gene.margin = 10000,
                      gene.name.font = 0.7,
                      show.probe.name =  TRUE,
                      probe.name.font = 0.6,
                      plot.transcripts = TRUE,
                      plot.transcripts.structure = TRUE,
                      y.label.font = 0.8,
                      y.label.margin = 0.1,
                      axis.number.font = 0.5,
                      chromatin.label.font = 0.7,
                      chromatin.label.margin = 0.02
                      ){

  if(!requireNamespace("karyoploteR")){
    message("This function requires the 'karyoploteR' package.")
    return(invisible())
  }

  if(!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene")){
    message("This function requires the 'TxDb.Hsapiens.UCSC.hg19.knownGene' package.")
    return(invisible())
  }

  if(!requireNamespace("org.Hs.eg.db")){
    message("This function requires the 'org.Hs.eg.db' package.")
    return(invisible())
  }

  if(!requireNamespace("regioneR")){
    message("This function requires the 'org.Hs.eg.db' package.")
    return(invisible())
  }

  # Retrieve the genetic coordinates of the gene
  cat("Retrieving genetic information for", gene.name, "from Ensembl...\n")

  # the bigwig visulization of histone modifciations and DNA-binding proteins only support hg19.
  genome = "hg19"
  GRCh = NULL
  if(genome == "hg19"){
    GRCh = 37
  }else if(genome == "hg38"){
    GRCh = 38
  }
  mart = biomaRt :: useEnsembl(biomart="ensembl",GRCh=GRCh, dataset = "hsapiens_gene_ensembl")
  genePosition <- biomaRt :: getBM(filters= "hgnc_symbol", attributes= c("chromosome_name","start_position", "end_position"),values = gene.name, mart= mart)
  chr_ID <- genePosition$chromosome_name # "16"
  chr_name <- paste0("chr",chr_ID)   # "chr16"

  # Retrieve probe information
  cat("Retrieving probe annotation...\n")
  suppressMessages({
    ProbeAnnotation = EpiMix_getInfiniumAnnotation(plat = met.platform, genome = genome)
  })
  ProbeAnnotation = convertAnnotToDF(ProbeAnnotation)
  ProbeAnnotation = mapProbeGene(ProbeAnnotation)
  ProbeAnnotation = ProbeAnnotation[ProbeAnnotation$gene == gene.name, ]
  ProbeAnnotation = ProbeAnnotation[order(ProbeAnnotation$CpG_beg), ]

  # Find the start and the end position of the gene. Note: gene can be inversely positioned on the genome.
  start_position = genePosition$start_position
  end_position  = genePosition$end_position
  start_position <- start_position - left.gene.margin
  end_position <-  end_position + right.gene.margin
  gene.region <-  regioneR :: toGRanges(data.frame(chr_name, start_position,end_position))

 # Set up some parameters for graph
  pp <- karyoploteR :: getDefaultPlotParams(plot.type=1)
  pp$leftmargin <- 0.15
  pp$topmargin <- 15
  pp$bottommargin <- 15
  pp$ideogramheight <- 5
  pp$data1inmargin <- 10
  pp$data1outmargin <- 0

  kp <- karyoploteR :: plotKaryotype(genome = genome, zoom = gene.region, plot.type=1, plot.params = pp)
  karyoploteR :: kpAddBaseNumbers(kp, tick.dist = (end_position-start_position)/6, minor.tick.dist =(end_position-start_position)/12,
                   add.units = TRUE, cex=0.5, tick.len = 3)

  # Plot gene structure with transcripts
  cat("Retrieving transcript information from TxDB...\n")

  # if(genome == "hg19"){
  #   TxDB_genome<- TxDb.Hsapiens.UCSC.hg19.knownGene
  # }else{
  #   TxDB_genome<- TxDb.Hsapiens.UCSC.hg38.knownGene
  #   }

  TxDB_genome<- TxDb.Hsapiens.UCSC.hg19.knownGene
  suppressMessages(
    genes.data <- karyoploteR :: makeGenesDataFromTxDb(TxDB_genome, karyoplot=kp)
  )

  # gene.data contain the transcript information for all the genes in the designated plot region, we need to filter the information to the gene specified by the users
  suppressMessages({
    gene_id_map = AnnotationDbi::select(org.Hs.eg.db, keys = gene.name, columns = c("ENTREZID"), keytype = "SYMBOL")
  })

  gene.id = gene_id_map$ENTREZID
  genes.data$genes = genes.data$genes[gene.id]
  genes.data$transcripts = genes.data$transcripts[gene.id]
  tx.id = as.character(genes.data$transcripts[[1]]$tx_id)
  genes.data$coding.exons = genes.data$coding.exons[tx.id]
  genes.data$non.coding.exons = genes.data$non.coding.exons[tx.id]
  genes.data <- karyoploteR :: addGeneNames(genes.data)
  #genes.data <- mergeTranscripts(genes.data)

  r0_Gene = 0.02
  r1_Gene = 0.15
  karyoploteR :: kpDataBackground(kp, r0=r0_Gene - 0.02,r1=r1_Gene + 0.01, col="#AACBFF")
  karyoploteR :: kpPlotGenes(
                              kp, data=genes.data,data.panel = 1,
                              r0=r0_Gene, r1=r1_Gene,
                              gene.name.position = "left",
                              gene.name.cex = gene.name.font,
                              plot.transcripts =  plot.transcripts,
                              plot.transcripts.structure = plot.transcripts.structure,
                              non.coding.exons.col = "green",
                              non.coding.exons.border.col = "green",
                              coding.exons.col = "red",
                              coding.exons.border.col ="red"
                              )

  karyoploteR :: kpAddLabels(kp, labels="Transcripts",r0=r0_Gene,r1=r1_Gene, srt=90, pos=1, label.margin = y.label.margin, cex = y.label.font)

  # Plot DM values for CpGs
  # Find differentially methylated CpG probes for the target gene
  CpGProbes <- ProbeAnnotation$probeID
  MethylationDrivers <- EpiMixResults$MethylationDrivers
  MixtureStates <- EpiMixResults$MixtureStates

  # Obtain vectors recording x,y coordinates for each CpG probes
  x = character(0)
  y = character(0)

  for (i in c(1:length(CpGProbes))){
    CpG_Position =  ProbeAnnotation$CpG_beg[i]
    x = c(x, CpG_Position)
    DM_Value = 0
    if(CpGProbes[i] %in% MethylationDrivers){
      values = unlist(MixtureStates[CpGProbes[i]])
      if (min(values) < 0) DM_Value = min(values) # This will change in the future for dual methylated CpGs
      if (max(values) > 0) DM_Value = max(values)
    }
    y = c(y, DM_Value)
  }
  x <- as.numeric(x)
  y <- as.numeric(y)

   # Plot axis and labels
  ymax <- ymin <- 0
  if(max(y) > 0){
    ymax = round(max(y), digits = 1) + 0.1
  }
  if(min(y) < 0){
    ymin = round(min(y), digits = 1) - 0.1
  }

  r0_DMValue = r1_Gene + 0.03
  r1_DMValue = r0_DMValue + 0.20

  karyoploteR :: kpAxis(kp,r0=r0_DMValue,r1=r1_DMValue, ymin = ymin, ymax = ymax, tick.pos = c(ymin, 0, ymax), cex = axis.number.font)
  karyoploteR :: kpAddLabels(kp, labels="DM value",r0=r0_DMValue,r1=r1_DMValue, srt=90, pos=1, label.margin = y.label.margin, cex = y.label.font)

  # Plot points and lines
  karyoploteR :: kpPoints(kp, chr=chr_name, x=x, y=y, ymin = ymin, ymax = ymax,r0 = r0_DMValue,r1 = r1_DMValue, cex=0.7, col="blue")
  karyoploteR :: kpLines(kp, chr=chr_name, x=x, y=y, ymin = ymin, ymax = ymax,r0 = r0_DMValue,r1 = r1_DMValue, lwd=1.5)

  # Plot the name for differentially methylated CpG probes
  if(show.probe.name){
    pos = 1
    if (max(y)>0) pos = 3
    karyoploteR :: kpText(kp, chr=chr_name, x=x[which(y!=0)], y=y[which(y!=0)], ymin = ymin, ymax = ymax,r0 = r0_DMValue,r1 = r1_DMValue, labels= CpGProbes[which(y!=0)], col="red", pos = pos, offest = 0.3, cex = probe.name.font)
  }

  # Plot two lines encompassing the differentially methylated CpG probes
  #karyoploteR :: kpSegments(kp, chr=chr_name, x0=x[which(y!=0)]-20, x1=x[which(y!=0)]-20, y0=rep(r0_DMValue, length(x[which(y!=0)])), y1 =rep(1, length(x[which(y!=0)])))
  #karyoploteR :: kpSegments(kp, chr=chr_name, x0=x[which(y!=0)]+20, x1=x[which(y!=0)]+20, y0=rep(r0_DMValue, length(x[which(y!=0)])), y1 =rep(1, length(x[which(y!=0)])))

   # Plot signals for histone marks and DNA binding proteins
  # Construct a vector with urls of the bigwig files
  base_url <- "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/"
  urlData = RCurl :: getURL(base_url)
  urlData2 = unlist(strsplit(urlData,"\\n"))
  filenames = as.matrix(urlData2[grep(roadmap.epigenome.id,urlData2)])
  filenames = unlist(strsplit(filenames, ">|<"))
  filenames = filenames[grep( roadmap.epigenome.id,filenames)]
  filenames = filenames[-grep("a href",filenames)]
  DNA.binding = character(0)
  for(file in filenames){
    start = unlist(gregexpr("-", file)) + 1
    end = unlist(gregexpr("pval", file)) - 2
    protein.name = substring(file, start, end)
    DNA.binding[protein.name] = paste0(base_url, file)
  }

  r0_Histone = r1_DMValue + 0.03
  total.tracks <- length(DNA.binding)
  colors <- RColorBrewer :: brewer.pal(total.tracks, "Set3")
  out.at <- karyoploteR :: autotrack(1:total.tracks, total.tracks, margin = 0.3, r0=r0_Histone)

  karyoploteR :: kpAddLabels(kp, labels = "Protein signals", r0 = out.at$r0, r1=out.at$r1, cex=y.label.font,
                             srt=90, pos=1, label.margin = y.label.margin)

  for(i in seq_len(total.tracks)){
    at <-  karyoploteR :: autotrack(i, total.tracks, r0=out.at$r0, r1=out.at$r1, margin = 0.1)
    kp <-  karyoploteR :: kpPlotBigWig(kp, data=DNA.binding[i], ymax="visible.region",
                                       r0=at$r0, r1=at$r1, col = colors[i])
    computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
    karyoploteR :: kpAxis(kp, ymin=0, ymax=computed.ymax, tick.pos = computed.ymax,
                          r0=at$r0, r1=at$r1, cex=axis.number.font)
    karyoploteR :: kpAddLabels(kp, labels = names(DNA.binding)[i], r0=at$r0, r1=at$r1,
                              cex= chromatin.label.font, label.margin = chromatin.label.margin)
  }
}


