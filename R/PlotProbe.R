#' The EpiMix_PlotProbe function
#' @description plot the genomic coordinate and the chromatin state of a specific CpG probe and the nearby genes.
#' @details this function requires additional dependencies: karyoploteR, TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.eg.db
#' @param probe.name character string indicating the CpG probe name.
#' @param EpiMixResults  resulting list object returned from EpiMix.
#' @param met.platform  character string indicating the type of micro-array where the DNA methylation data were collected.Can be either "HM27", "HM450" or "EPIC". Default: "HM450"
#' @param roadmap.epigenome.id character string indicating the epigenome id (EID) for a reference tissue or cell type. Default: "E002"
#' @param numFlankingGenes numeric value indicating the number of flanking genes to be plotted with the CpG probe. Default: 20 (10 gene upstream and 10 gene downstream).
#' @param left.gene.margin  numeric value indicating the number of extra nucleotide bases to be plotted on the left side of the image. Default: 10000.
#' @param right.gene.margin numeric value indicating the number of extra nucleotide bases to be plotted on the right side of the image. Default: 10000.
#' @param gene.name.pos  integer indicating the position for plotting the gene name relative to the gene structure. Should be 1 or 2 or 3 or 4, indicating bottom, left, top, and right, respectively.
#' @param gene.name.size numeric value indicating the font size of the gene names in pixels.
#' @param gene.arrow.length numeric value indicating the size of the arrow which indicates the positioning of the gene.
#' @param gene.line.width numeric value indicating the line width for the genes.
#' @param plot.chromatin.state logical indicating whether to plot the DNase-seq and histone ChIP-seq signals. Warnings: If the "numFlankingGenes" is a larger than 15, plotting the chromatin state may flood the internal memory.
#' @param y.label.font font size of the y axis label.
#' @param y.label.margin distance between y axis label and y axis.
#' @param axis.number.font font size of axis ticks and numbers.
#' @param chromatin.label.font font size of the labels of the histone proteins.
#' @param chromatin.label.margin distance between the histone protein labels and axis.
#' @importFrom biomaRt useDataset getBM useEnsembl
#' @importFrom GenomicRanges makeGRangesFromDataFrame seqnames start end mcols
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom RColorBrewer brewer.pal
#' @importFrom RCurl getURL
#' @return plot with CpG probe and nearby genes. Genes whose expression is significantly negatively associated with the methylation of the probe are shown in red, while the others are shown in black.
#' @details
#' roadmap.epigenome.id: since the chromatin state is tissue or cell-type specific,
#' EpiMix needs to know the reference tissue or cell type in order to retrieve the proper DNase-seq and histone ChIP-seq data.
#' Available epigenome ids can be obtained from the Roadmap Epigenomic study (Nature, PMID: 25693563, figure 2).
#' They can also be retrieved from the list.epigenomes() function.
#' @export
#' @examples
#' {
#' library(karyoploteR)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(org.Hs.eg.db)
#' library(regioneR)
#'
#' data(Sample_EpiMixResults_Regular)
#'
#' # The CpG site to plot
#' probe.name = "cg00374492"
#'
#' # The number of adjacent genes to be plotted
#' numFlankingGenes = 10
#'
#' # Set up the reference cell/tissue type
#' roadmap.epigenome.id = "E096"
#'
#' # Generate the plot
#' EpiMix_PlotProbe(probe.name = probe.name,
#'                  EpiMixResults = Sample_EpiMixResults_Regular,
#'                  met.platform = "HM450",
#'                  roadmap.epigenome.id = roadmap.epigenome.id,
#'                  numFlankingGenes = numFlankingGenes)
#'
#' }
#'
EpiMix_PlotProbe <- function(probe.name,
                             EpiMixResults,
                             met.platform = "HM450",
                             roadmap.epigenome.id = "E002",
                             numFlankingGenes = 20,
                             left.gene.margin = 10000,
                             right.gene.margin = 10000,
                             gene.name.pos = 2,
                             gene.name.size = 0.5,
                             gene.arrow.length = 0.05,
                             gene.line.width = 2,
                             plot.chromatin.state = TRUE,
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

  if(numFlankingGenes > 10){
    warnings("Setting the number of flanking genes > 10 may blow the internal memory.Try to set
             a smaller gene number as possible")
  }

  # get nearby genes of the probe
  genome = "hg19"
  cat("Retrieving probe annotation...\n")
  suppressMessages({
    ProbeAnnotation = EpiMix_getInfiniumAnnotation(plat = met.platform, genome = genome)
  })
  ProbeAnnotation = ProbeAnnotation[probe.name, ]

  geneAnnot <- getTSS(genome = genome)
  NearbyGenes <- GetNearGenes(geneAnnot = geneAnnot,
                              TRange =  ProbeAnnotation,
                              numFlankingGenes = numFlankingGenes)

  # Retrieve the genetic coordinates of each near genes
  cat("Retrieving gene information from Ensembl...\n")
  GRCh = NULL
  if(genome == "hg19"){
    GRCh = 37
  }else if(genome == "hg38"){
    GRCh = 38
  }
  mart = biomaRt :: useEnsembl(biomart="ensembl",GRCh=GRCh, dataset = "hsapiens_gene_ensembl")
  genePosition <- biomaRt :: getBM(filters= "hgnc_symbol",
                                    attributes= c("hgnc_symbol","chromosome_name","start_position", "end_position", "strand"),
                                    values = NearbyGenes$Symbol,
                                    mart= mart)

  # if a gene is postioned on the opposite strand, flip the start and end position
  for (i in 1:nrow(genePosition)){
    start = genePosition$start_position[i]
    end = genePosition$end_position[i]
    if(genePosition$strand[i] == -1){
      genePosition$start_position[i] = end
      genePosition$end_position[i] = start
    }
  }

  chr_ID <- genePosition$chromosome_name[1] # "16"
  chr_name <- paste0("chr",chr_ID)   # "chr16"

  # Find the range of the genomic region to be plotted
  positions = c(genePosition$start_position, genePosition$end_position)
  left_end = min(positions) - left.gene.margin
  right_end = max(positions) + right.gene.margin
  gene.region <- regioneR :: toGRanges(data.frame(chr_name, left_end,right_end))

  # Split probe annotation by gene expression state
  FunctionalProbes = EpiMixResults$FunctionalPairs
  FunctionalProbes = FunctionalProbes[FunctionalProbes$Probe == probe.name, ]
  DEGenes_genePosition = genePosition[genePosition$hgnc_symbol %in% FunctionalProbes$Gene, ]
  Other_genePosition = genePosition[!genePosition$hgnc_symbol %in% FunctionalProbes$Gene, ]

  # Set up some parameters for graph
  pp <- karyoploteR :: getDefaultPlotParams(plot.type=2)
  pp$leftmargin <- 0.15
  pp$topmargin <- 15
  pp$bottommargin <- 15
  pp$ideogramheight <- 5
  pp$data1inmargin <- 10
  pp$data1outmargin <- 0

  kp <- karyoploteR :: plotKaryotype(genome = genome, zoom = gene.region, plot.type = 2, plot.params = pp)
  karyoploteR :: kpAddBaseNumbers(kp,
                                 tick.dist = (right_end-left_end)/6,
                                 minor.tick.dist =(right_end-left_end)/12,
                                 add.units = TRUE,
                                 cex=0.5,
                                 tick.len = 3)

  # Plot the probe
  r0_probe = 0
  r1_probe = r0_probe + 0.20
  probe_pos = GenomicRanges :: start(ranges(ProbeAnnotation))
  karyoploteR :: kpLines(kp,
                        chr = chr_name,
                        x = c(probe_pos,  probe_pos),
                        data.panel = 2,
                        y = c(0.1,0.9),
                        r0 = r0_probe,
                        r1= r1_probe,
                        col = "purple"
                       )

  karyoploteR :: kpText(kp,
                         chr = chr_name,
                         x = c(probe_pos),
                         y = c(0.5),
                         data.panel = 2,
                         r0 = r0_probe,
                         r1= r1_probe,
                         col = "purple",
                         label = names(ProbeAnnotation),
                         pos = 4,
                         cex = 0.7
                         )


  karyoploteR :: kpAddLabels(kp, data.panel = 2, labels="Probe",r0=r0_probe,r1=r1_probe,
                             srt=90, pos=1, label.margin = y.label.margin, cex = y.label.font)


  # Plot the gene position
  r0_Gene = r1_probe + 0.04
  r1_Gene = 0.98
  karyoploteR :: kpDataBackground(kp, data.panel = 2, r0=r0_Gene - 0.03,r1=r1_Gene + 0.01, col="#AACBFF")

  # Plot the genes whose expression were not different
  if(nrow(Other_genePosition) > 0){
    karyoploteR ::  kpArrows(kp,
                             chr=chr_name,
                             data.panel = 2,
                             r0 = r0_Gene,
                             r1 = r1_Gene,
                             x0=Other_genePosition$start_position,
                             x1=Other_genePosition$end_position,
                             y0=seq(0.2,0.8,length.out =nrow(Other_genePosition)),
                             y1=seq(0.2,0.8,length.out =nrow(Other_genePosition)),
                             lwd = gene.line.width,
                             length = gene.arrow.length)

    karyoploteR :: kpText(kp,
                         chr=chr_name,
                         data.panel = 2,
                         r0 = r0_Gene,
                         r1 = r1_Gene,
                         x=Other_genePosition$start_position,
                         y=seq(0.2,0.8,length.out =nrow(Other_genePosition)) + 0.05,
                         labels = Other_genePosition$hgnc_symbol,
                         pos = gene.name.pos,
                         cex = gene.name.size
                  )
  }

  # Plot the genes whose expression was different and show this gene in red

  if(nrow(DEGenes_genePosition) > 0){
    karyoploteR :: kpArrows(kp,
                           chr=chr_name,
                           data.panel = 2,
                           r0 = r0_Gene,
                           r1 = r1_Gene,
                           x0=DEGenes_genePosition$start_position,
                           x1=DEGenes_genePosition$end_position,
                           y0=seq(0.2,0.7,length.out =nrow(DEGenes_genePosition)),
                           y1=seq(0.2,0.7,length.out =nrow(DEGenes_genePosition)),
                           lwd = gene.line.width,
                           col = "red",
                           length = gene.arrow.length)


    karyoploteR :: kpText(kp,
                         chr=chr_name,
                         data.panel = 2,
                         r0 = r0_Gene,
                         r1 = r1_Gene,
                         x=DEGenes_genePosition$start_position,
                         y=seq(0.2,0.7,length.out =nrow(DEGenes_genePosition)) + 0.05,
                         labels = DEGenes_genePosition$hgnc_symbol,
                         pos = gene.name.pos,
                         offset = 0.2,
                         cex = gene.name.size,
                         col = "red"
                  )
  }
  karyoploteR :: kpAddLabels(kp, labels="Genes",data.panel = 2, r0=r0_Gene,r1=r1_Gene, srt=90,
                             pos=1, label.margin = y.label.margin, cex = y.label.font)


# Plot histone marks and DNase1 sensitivity
  cat("Retrieving signal intensity of histone marks and DNase sensitivity from RoadmapEpigenomics...\n")
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

  r0_Histone = 0
  total.tracks <- length(DNA.binding)
  colors <-  RColorBrewer  :: brewer.pal(total.tracks, "Set3")
  out.at <- karyoploteR :: autotrack(1:total.tracks, total.tracks, margin = 0.3, r0=r0_Histone)

  karyoploteR :: kpAddLabels(kp, labels = "Protein signals", r0 = out.at$r0, r1=out.at$r1, cex=y.label.font,
                            srt=90, pos=1, label.margin = y.label.margin)

  for(i in seq_len(total.tracks)){
    at <- karyoploteR :: autotrack(i, total.tracks, r0=out.at$r0, r1=out.at$r1, margin = 0.1)
    kp <- karyoploteR :: kpPlotBigWig(kp, data=DNA.binding[i], ymax="visible.region",
                                     r0=at$r0, r1=at$r1, col = colors[i])
    computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
    karyoploteR :: kpAxis(kp, ymin=0, ymax=computed.ymax, tick.pos = computed.ymax,
                          r0=at$r0, r1=at$r1, cex=axis.number.font)
    karyoploteR :: kpAddLabels(kp, labels = names(DNA.binding)[i], r0=at$r0, r1=at$r1,
                               cex= chromatin.label.font, label.margin = chromatin.label.margin)
  }
}
