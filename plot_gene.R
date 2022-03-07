#' The plot.gene function
#' @description plot the genomic coordinate, DM values and chromatin state for each CpG probe of a specific gene.
#' @details this function requires R package dependencies: karyoploteR, TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.eg.db
#' @param gene.name character string indicating the name of the gene to be plotted
#' @param EpiMixResults the list object returned from EpiMix
#' @param met.platform  character string indicating the type of Illumina Infinium DNA methylation BeadChip for collecting the DNA methylation data. The value should be either "HM450" or "EPIC". Default: "HM450"
#' @param left.gene.margin  numeric value indicating the number of extra nucleotide bases to be plotted on the left side of the gene. Default: 10000.
#' @param right.gene.margin numeric value indicating the number of extra nucleotide bases to be plotted on the right side of the gene. Default: 10000.
#' @param gene.name.font numeric value indicating the font size for the gene name. Default: 0.7.
#' @param plot.transcripts logical indicating whether to plot each individual transcripts for the gene. Default: TRUE. If False, the gene will plotted with a single rectangle, without showing each individual transcripts.
#' @param plot.transcripts.structure logical indicating whether to plot the transcript structure (introns and exons). If TRUE, the non-coding exons will be shown in green and the coding exons will be shown in red. Default: TRUE.
#' @param show.probe.name logical indicating whether to show the name(s) for each differentially methylated CpG probe. Default: TRUE
#' @param probe.name.font numeric value indicating the font size of the name(s) for the differentially methylated probe(s) in pixels. Default: 0.6.
#' @param roadmap.epigenome.id character string indicating the epigenome id (EID) for the reference cell used for retrieving the DNase-seq and histone-ChIP signals. The EID can be found in the figure 2 of this paper: PMID: 25693563
#' @param axis.label.font numeric value indicating the font size for the axis label. Default: 0.5.
#' @param axis.label.margin numeric value indicating the margin between the label and the axis. Default: 0.1.
#' @param y.label.font font size of the y axis label
#' @param y.label.margin distance between y axis label and y axis
#' @param axis.number.font font size of axis ticks and numbers
#' @param chromatin.label.font font size of the labels of the histone proteins
#' @param chromatin.label.margin distance between the histone protein labels and axis
#' @import biomaRt
#' @import GenomicRanges
#' @import RColorBrewer
#' @return plot of the genomic coordinate, DM values and chromatin state for each CpG probe of a specific gene.
#' @export

plot.gene <- function(gene.name,
                      EpiMixResults,
                      met.platform = "HM450",
                      left.gene.margin = 10000,
                      right.gene.margin = 10000,
                      gene.name.font = 0.7,
                      show.probe.name =  TRUE,
                      probe.name.font = 0.6,
                      plot.transcripts = TRUE,
                      plot.transcripts.structure = TRUE,
                      roadmap.epigenome.id = "E002",
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
  mart = useEnsembl(biomart="ensembl",GRCh=GRCh, dataset = "hsapiens_gene_ensembl")
  genePosition <- getBM(filters= "hgnc_symbol", attributes= c("chromosome_name","start_position", "end_position"),values = gene.name, mart= mart)
  chr_ID <- genePosition$chromosome_name # "16"
  chr_name <- paste0("chr",chr_ID)   # "chr16"

  # Retrieve probe information
  cat("Retrieving probe annotation...\n")
  ProbeAnnotation = EpiMix_getInfiniumAnnotation(plat = met.platform, genome = genome)
  ProbeAnnotation = convertAnnotToDF(ProbeAnnotation)
  ProbeAnnotation = mapProbeGene(ProbeAnnotation)
  ProbeAnnotation = ProbeAnnotation[ProbeAnnotation$gene == gene.name, ]
  ProbeAnnotation = ProbeAnnotation[order(ProbeAnnotation$CpG_beg), ]

  # Find the start and the end position of the gene. Note: gene can be inversely positioned on the genome.
  start_position = genePosition$start_position
  end_position  = genePosition$end_position
  start_position <- start_position - left.gene.margin
  end_position <-  end_position + right.gene.margin
  gene.region <- toGRanges(data.frame(chr_name, start_position,end_position))

 # Set up some parameters for graph
  pp <- getDefaultPlotParams(plot.type=1)
  pp$leftmargin <- 0.15
  pp$topmargin <- 15
  pp$bottommargin <- 15
  pp$ideogramheight <- 5
  pp$data1inmargin <- 10
  pp$data1outmargin <- 0

  kp <- plotKaryotype(genome = genome, zoom = gene.region, plot.type=1, plot.params = pp)
  kpAddBaseNumbers(kp, tick.dist = (end_position-start_position)/6, minor.tick.dist =(end_position-start_position)/12,
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
    genes.data <- makeGenesDataFromTxDb(TxDB_genome, karyoplot=kp)
  )

  # gene.data contain the transcript information for all the genes in the designated plot region, we need to filter the information to the gene specified by the users
  gene_id_map = AnnotationDbi::select(org.Hs.eg.db, keys = gene.name, columns = c("ENTREZID"), keytype = "SYMBOL")
  gene.id = gene_id_map$ENTREZID
  genes.data$genes = genes.data$genes[gene.id]
  genes.data$transcripts = genes.data$transcripts[gene.id]
  tx.id = as.character(genes.data$transcripts[[1]]$tx_id)
  genes.data$coding.exons = genes.data$coding.exons[tx.id]
  genes.data$non.coding.exons = genes.data$non.coding.exons[tx.id]
  genes.data <- addGeneNames(genes.data)
  #genes.data <- mergeTranscripts(genes.data)

  r0_Gene = 0.02
  r1_Gene = 0.15
  kpDataBackground(kp, r0=r0_Gene - 0.02,r1=r1_Gene + 0.01, col="#AACBFF")

  kpPlotGenes(kp, data=genes.data,data.panel = 1,
              r0=r0_Gene, r1=r1_Gene,
              gene.name.position = "left",
              gene.name.cex = gene.name.font,
              plot.transcripts =  plot.transcripts,
              plot.transcripts.structure = plot.transcripts.structure,
              non.coding.exons.col = "green",
              non.coding.exons.border.col = "green",
              coding.exons.col = "red",
              coding.exons.border.col ="red")

  kpAddLabels(kp, labels="Transcripts",r0=r0_Gene,r1=r1_Gene, srt=90, pos=1, label.margin = y.label.margin, cex = y.label.font)

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

  kpAxis(kp,r0=r0_DMValue,r1=r1_DMValue, ymin = ymin, ymax = ymax, tick.pos = c(ymin, 0, ymax), cex = axis.number.font)
  kpAddLabels(kp, labels="DM value",r0=r0_DMValue,r1=r1_DMValue, srt=90, pos=1, label.margin = y.label.margin, cex = y.label.font)

  # Plot points and lines
  kpPoints(kp, chr=chr_name, x=x, y=y, ymin = ymin, ymax = ymax,r0 = r0_DMValue,r1 = r1_DMValue, cex=0.7, col="blue")
  kpLines(kp, chr=chr_name, x=x, y=y, ymin = ymin, ymax = ymax,r0 = r0_DMValue,r1 = r1_DMValue, lwd=1.5)

  # Plot the name for differentially methylated CpG probes
  if(show.probe.name){
    pos = 1
    if (max(y)>0) pos = 3
    kpText(kp, chr=chr_name, x=x[which(y!=0)], y=y[which(y!=0)], ymin = ymin, ymax = ymax,r0 = r0_DMValue,r1 = r1_DMValue, labels= CpGProbes[which(y!=0)], col="red", pos = pos, offest = 0.3, cex = probe.name.font)
  }

  # Plot two lines encompassing the differentially methylated CpG probes
  #kpSegments(kp, chr=chr_name, x0=x[which(y!=0)]-20, x1=x[which(y!=0)]-20, y0=rep(r0_DMValue, length(x[which(y!=0)])), y1 =rep(1, length(x[which(y!=0)])))
  #kpSegments(kp, chr=chr_name, x0=x[which(y!=0)]+20, x1=x[which(y!=0)]+20, y0=rep(r0_DMValue, length(x[which(y!=0)])), y1 =rep(1, length(x[which(y!=0)])))

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
  colors <-   RColorBrewer :: brewer.pal(total.tracks, "Set3")
  out.at <- autotrack(1:total.tracks, total.tracks, margin = 0.3, r0=r0_Histone)

  kpAddLabels(kp, labels = "Protein signals", r0 = out.at$r0, r1=out.at$r1, cex=y.label.font,
              srt=90, pos=1, label.margin = y.label.margin)

  for(i in seq_len(total.tracks)){
    at <- autotrack(i, total.tracks, r0=out.at$r0, r1=out.at$r1, margin = 0.1)
    kp <- kpPlotBigWig(kp, data=DNA.binding[i], ymax="visible.region",
                       r0=at$r0, r1=at$r1, col = colors[i])
    computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
    kpAxis(kp, ymin=0, ymax=computed.ymax, tick.pos = computed.ymax,
           r0=at$r0, r1=at$r1, cex=axis.number.font)
    kpAddLabels(kp, labels = names(DNA.binding)[i], r0=at$r0, r1=at$r1,
                cex= chromatin.label.font, label.margin = chromatin.label.margin)
  }
}


