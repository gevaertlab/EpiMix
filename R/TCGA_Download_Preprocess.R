#' The TCGA_Download_DNAmethylation function
#'
#' @description Download DNA methylation data from TCGA.
#' @param CancerSite character of length 1 with TCGA cancer code.
#' @param TargetDirectory character with directory where a folder for downloaded files will be created.
#' @param downloadData logical indicating if data should be downloaded (default: TRUE). If false, the url of the desired data is returned.
#' @return list with paths to downloaded files for both 27k and 450k methylation data.
#' @export
#' @keywords download
#' @examples
#' \donttest{
#' METdirectories <- TCGA_Download_DNAmethylation(CancerSit = 'OV', TargetDirectory = tempdir())
#' }
#'
TCGA_Download_DNAmethylation <- function(CancerSite, TargetDirectory, downloadData = TRUE) {


    options(timeout = 1e+05)
    dir.create(TargetDirectory, showWarnings = FALSE)

    # download the 27k data
    dataType <- "stddata"
    dataFileTag <- "Merge_methylation__humanmethylation27"
    cat("Searching 27k MET data for:", CancerSite, "\n")
    METdirectory27k <- get_firehoseData(downloadData, TargetDirectory, CancerSite,
        dataType, dataFileTag)

    # download the 450k data
    dataFileTag <- "Merge_methylation__humanmethylation450"
    cat("Searching 450k MET data for:", CancerSite, "\n")
    METdirectory450k <- get_firehoseData(downloadData, TargetDirectory, CancerSite,
        dataType, dataFileTag)

    return(METdirectories = list(METdirectory27k = METdirectory27k, METdirectory450k = METdirectory450k))
}

#' The get_firehoseData function
#'
#' Gets data from TCGA's firehose.
#' @param downloadData logical indicating if data should be downloaded (default: TRUE). If false, the url of the desired data is returned.
#' @param saveDir path to directory to save downloaded files.
#' @param TCGA_acronym_uppercase TCGA's cancer site code.
#' @param dataType type of data in TCGA (default: 'stddata').
#' @param dataFileTag name of the file to be downloaded (the default is to download RNAseq data, but this can be changed to download other data).
#' @param FFPE logical indicating if FFPE data should be downloaded (default: FALSE).
#' @param fileType type of downloaded file (default: 'fileType', other type not admitted at the moment).
#' @param gdacURL gdac url.
#' @param untarUngzip logical indicating if the gzip file downloaded should be untarred (default: TRUE).
#' @param printDisease_abbr if TRUE data is not downloaded but all the possible cancer sites codes are shown (default: FALSE).
#' @return DownloadedFile path to directory with downloaded files.
#' @keywords internal
#'
get_firehoseData <- function(downloadData = TRUE, saveDir = "./", TCGA_acronym_uppercase = "LUAD",
    dataType = "stddata", dataFileTag = "mRNAseq_Preprocess.Level_3", FFPE = FALSE,
    fileType = "tar.gz", gdacURL = "http://gdac.broadinstitute.org/runs/", untarUngzip = TRUE,
    printDisease_abbr = FALSE) {

    # Cases Shipped by BCR # Cases with Data* Date Last Updated (mm/dd/yy)
    cancers <- c("Acute Myeloid Leukemia [LAML] \n", "Adrenocortical carcinoma [ACC] \n",
        "Bladder Urothelial Carcinoma [BLCA] \n", "Brain Lower Grade Glioma [LGG] \n",
        "Breast invasive carcinoma [BRCA] \n", "Cervical squamous cell carcinoma and endocervical adenocarcinoma [CESC] \n",
        "Cholangiocarcinoma [CHOL] \n", "Colon adenocarcinoma [COAD] \n", "Esophageal carcinoma [ESCA] \n",
        "Glioblastoma multiforme [GBM] \n", "Head and Neck squamous cell carcinoma [HNSC]   \n",
        "Kidney Chromophobe [KICH] \n", "Kidney renal clear cell carcinoma [KIRC]   \n",
        "Kidney renal papillary cell carcinoma [KIRP]  \n", "Liver hepatocellular carcinoma [LIHC]  \n",
        "Lung adenocarcinoma [LUAD]    \n", "Lung squamous cell carcinoma [LUSC] \n",
        "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma [DLBC]    \n", "Mesothelioma [MESO] \n",
        "Ovarian serous cystadenocarcinoma [OV]    \n", "Pancreatic adenocarcinoma [PAAD]   \n",
        "Pheochromocytoma and Paraganglioma [PCPG] \n", "Prostate adenocarcinoma [PRAD] \n",
        "Rectum adenocarcinoma [READ]  \n", "Sarcoma [SARC] \n", "Skin Cutaneous Melanoma [SKCM] \n",
        "Stomach adenocarcinoma [STAD] \n", "Testicular Germ Cell Tumors [TGCT] \n",
        "Thymoma [THYM] \n", "Thyroid carcinoma [THCA]  \n", "Uterine Carcinosarcoma [UCS]    \n",
        "Uterine Corpus Endometrial Carcinoma [UCEC]   \n", "Uveal Melanoma [UVM] \n",
        "Colorectal Adenocarcinoma [COADREAD] \n")

    if (printDisease_abbr) {
        return(cat("here are the possible TCGA database disease acronyms. \nRe-run this function with printDisease_abbr=FALSE to then run an actual query.\n\n",
            cancers))
    }
    gdacURL_orig <- gdacURL

    # New code to handle dates - Marcos
    gdacURLnew <- paste0(gdacURL_orig, dataType, "__latest/data/", TCGA_acronym_uppercase,
        "/")
    urlDataNew <- RCurl::getURL(gdacURLnew)
    urlDataNew <- limma::strsplit2(urlDataNew, "href=\\\"")  #regular expressions: need \ to have R recognize any ' or \ that's actually in our text
    getlatestdate <- urlDataNew[grep("^201[0-9][01][0-9][0123][0-9]", urlDataNew)]
    getlatestdate <- substring(getlatestdate, 1, 8)
    gdacURLnew <- paste0(gdacURLnew, getlatestdate, "/")
    urlData <- RCurl::getURL(gdacURLnew)
    urlData <- limma::strsplit2(urlData, "href=\\\"")  #regular expressions: need \ to have R recognize any ' or \ that's actually in our text
    lastDateCompress <- lastDate <- getlatestdate  # for compatibility with the rest of old code
    gdacURL <- gdacURLnew  # for compatibility with the rest of old code
    # end New code

    # remove any FFPE datasets, or only keep those depending on user inputs.
    if (FFPE) {
        urlData <- urlData[grep("FFPE", urlData)]
        if (length(urlData) == 0) {
            stop("\nNo FFPE data found for this query. Try FFPE=FALSE.\n")
        }
    } else {
        # we DON'T want FFPE data. but if no FFPE data to begin with: don't
        # subset on this.
        if (length(grep("FFPE", urlData)) > 0) {
            urlData <- urlData[-grep("FFPE", urlData)]
        }
        if (length(urlData) == 0) {
            stop("\nNo non-FFPE data found for this query. Try FFPE=TRUE.\n")
        }
    }
    # now get full dataset name.
    fileName <- urlData[grep(dataFileTag, urlData)]

    if (length(fileName) == 0) {
        # warnMessage <- paste0('\nNot returning any viable url data paths
        # after searching by date for disease ',TCGA_acronym_uppercase,' for
        # data type ',dataFileTag ,'.No data was downloaded.\n')
        # warning(warnMessage)
        cat("\tThere is no", dataFileTag, "data for", TCGA_acronym_uppercase, "\n")
        return(NA)
    }
    # some redundancy..but that' OK because we'll add back on the unique tar.gz
    # file tag. first file is one we want - not md5 file.
    fileName <- limma::strsplit2(fileName, "tar.gz")[1, 1]
    fileName <- paste(fileName, fileType, sep = "")

    # final download url
    gdacURL <- paste(gdacURL, fileName, sep = "")
    # Directory for downloads
    saveDir <- paste(saveDir, "/", "gdac_", lastDateCompress, "/", sep = "")

    if (!grepl("Windows", Sys.info()["sysname"])) {
        # Not Windows
        tarfile <- paste0(saveDir, fileName)
        finalDir <- strsplit(tarfile, paste0(".", fileType))[[1]][1]
        if (downloadData) {
            cat("\tDownloading", dataFileTag, "data, version:", lastDate, "\n")
            cat("\tThis may take 10-60 minutes depending on the size of the data set.\n")
            dir.create(saveDir, showWarnings = FALSE)
            # download file
            download.file(gdacURL, destfile = paste0(saveDir, fileName), quiet = FALSE,
                mode = "wb")
            # this assumes a tar.gz file.
            if (fileType == "tar.gz" && untarUngzip) {
                cat("\tUnpacking data.\n")
                tarfile <- paste0(saveDir, fileName)
                untar(tarfile, exdir = saveDir)
                # remove tarred file
                fileToRemove <- limma::strsplit2(gdacURL, "/")[, ncol(limma::strsplit2(gdacURL,
                  "/"))]
                removed <- file.remove(paste0(saveDir, fileToRemove))
            } else if (untarUngzip) {
                warning("File expansion/opening only built in for tar.gz files at the moment.\n")
            }
            cat("\tFinished downloading", dataFileTag, "data to", finalDir, "\n")
        } else {
            cat("\tdownload data url is :\n ", gdacURL, "\n")
        }
        DownloadedFile <- paste0(finalDir, "/")
        return(DownloadedFile)
    } else {
        # new code to handle long names in windows - Marcos WINDOWS: name of
        # file can be too long (and it's repeated in the folder and the file)
        # and windows internally will use another representation for the name
        # of the folder, so then when we want to load the file it says it
        # doesn't exist. So I'm changing the name of the folder to prevent
        # this. We can't change the name of the file as it's used in the
        # Preprocess functions
        idx <- which(vapply(c("methylation27", "methylation450", "mRNAseq", "transcriptome"),
            grepl, dataFileTag, FUN.VALUE = integer(1)))[1]
        newtag <- ifelse(idx == 1, "meth27", ifelse(idx == 2, "meth450", "geneexpr"))
        nameForFolder <- paste(TCGA_acronym_uppercase, dataType, newtag, sep = "_")
        nameForDownloadedFile <- paste0(nameForFolder, ".", fileType)
        nameForDownloadedFileFullPath <- paste0(saveDir, nameForDownloadedFile)
        finalDir <- paste0(saveDir, nameForFolder)
        if (downloadData) {
            cat("\tDownloading", dataFileTag, "data, version:", lastDate, "\n")
            cat("\tThis may take 10-60 minutes depending on the size of the data set.\n")
            dir.create(saveDir, showWarnings = FALSE)

            # Create a virtual drive to overcome long names issue in Windows
            saveDir2 <- gsub("\\\\", "/", saveDir)
            saveDir2 <- substr(saveDir2, 1, nchar(saveDir2) - 1)
            system(paste("subst x:", saveDir2))

            download.file(gdacURL, destfile = paste0("x://", nameForDownloadedFile),
                quiet = FALSE, mode = "wb")
            # this assumes a tar.gz file.
            if (fileType == "tar.gz" && untarUngzip) {
                cat("\tUnpacking data.\n")
                untar(nameForDownloadedFileFullPath, exdir = saveDir)
                # untar(paste0('x://', nameForDownloadedFile), exdir = 'x://')
                # # doesn't work because it calls system and system doesnt know
                # about the virtual drive
                removed <- file.remove(paste0("x://", nameForDownloadedFile))
                # Anyway I change folder name to make it shorter
                changed <- file.rename(from = paste0("x://", gsub(".tar.gz", "",
                  fileName)), to = paste0("x://", nameForFolder))
                system("subst x: /D")  # stop the virtual drive
            } else if (untarUngzip) {
                warning("File expansion/opening only built in for tar.gz files at the moment.\n")
            }
            cat("\tFinished downloading", dataFileTag, "data to", finalDir, "\n")
        } else {
            cat("\tdownload data url is :\n ", gdacURL, "\n")
        }
        DownloadedFile <- paste0(finalDir, "/")
        return(DownloadedFile)
        # end new code
    }
}

#' The TCGA_Select_Dataset function
#' @description internal function to select which MET dataset to use
#' @param CancerSite TCGA cancer code
#' @param MET_Data_27K matrix of MET_Data_27K
#' @param MET_Data_450K matrix of MET_Data_450K
#' @param use450K logic indicating whether to force use 450K data
#'
#' @return the selected MET data set
#'
TCGA_Select_Dataset <- function(CancerSite, MET_Data_27K, MET_Data_450K, use450K) {
    MET_Data <- NULL
    if (use450K || CancerSite == "LAML") {
        # LAML is a special case, only using 450k data
        MET_Data <- MET_Data_450K
    } else {
        if (length(MET_Data_27K) != 0 & length(MET_Data_450K) != 0) {
            if (ncol(MET_Data_27K) > ncol(MET_Data_450K) & ncol(MET_Data_27K) > 50) {
                cat("Not enough 450k samples, only using the 27k (need min 50 samples).\n")
                MET_Data <- MET_Data_27K
            } else {
                cat("Not enough 27k samples, only using the 450k.\n")
                MET_Data <- MET_Data_450K
            }
        } else if (length(MET_Data_27K) != 0) {
            cat("\tOnly 27k samples.\n")
            MET_Data <- MET_Data_27K
        } else {
            cat("\tOnly 450k samples.\n")
            MET_Data <- MET_Data_450K
        }
    }
    return(MET_Data)
}

#' The TCGA_Preprocess_DNAmethylation function
#'
#' @description Pre-processes DNA methylation data from TCGA.
#' @param CancerSite character string indicating the TCGA cancer code.
#' @param METdirectories character vector with directories with the downloaded data. It can be the object returned by the TCGA_Download_DNAmethylation function.
#' @param doBatchCorrection logical indicating whether to perform batch correction. Default: False.
#' @param batch.correction.method character string indicating the method to perform batch correction. The value should be either 'Seurat' or 'Combat'. Default: 'Seurat'. Note: Seurat is much faster than the Combat.
#' @param MissingValueThreshold numeric values indicating the threshold for removing samples or genes with missing values.Default: 0.2.
#' @param cores integer indicating the number of cores to be used for performing batch correction with Combat.
#' @param use450K logic indicating whether to force use 450K, instead of 27K data.
#' @return pre-processed methylation data matrix with CpG probe in rows and samples in columns.
#' @details
#' Pre-process includes eliminating samples and genes with too many NAs, imputing NAs, and doing Batch correction.
#' If there are samples with both 27k and 450k data, the 27k data will be used only if the sample number in the 27k data is greater than the 450k data and there is more than 50 samples in the 27k data. Otherwise, the 450k data is used and the 27k data is discarded.
#' @return Pre-processed methylation data matrix with CpG probe in rows and samples in columns.
#' @export
#' @keywords preprocess
#' @examples
#' \donttest{
#' METdirectories <- TCGA_Download_DNAmethylation(CancerSite = 'OV', TargetDirectory = tempdir())
#' METProcessedData <- TCGA_Preprocess_DNAmethylation(CancerSite = 'OV',
#'                                                    METdirectories = METdirectories)
#' }

TCGA_Preprocess_DNAmethylation <- function(CancerSite, METdirectories, doBatchCorrection = FALSE,
    batch.correction.method = "Seurat", MissingValueThreshold = 0.2, cores = 1, use450K = FALSE) {

    cat("\tProcessing data for", CancerSite, "\n")

    # set up some paramerters
    if (cores > 1 & batch.correction.method == "Combat") {
        # unregister()
        cat("Registering sockets on multiple CPU cores...\n")
        cl <- parallel::makeCluster(cores)
    }

    # ---------------------------------------------------------------------------------------------
    # Step 1: Load 27K and 450K data and select which dataset to use
    # ---------------------------------------------------------------------------------------------
    MET_Data_27K <- c()
    MET_Data_450K <- c()
    if (!is.na(METdirectories$METdirectory27k)) {
        cat("\tLoading data for 27k.\n")
        MET_Data_27K <- TCGA_Load_MethylationData(METdirectories$METdirectory27k,
            ArrayType = "27K")
    }
    if (!is.na(METdirectories$METdirectory450k)) {
        cat("\tLoading data for 450k.\n")
        MET_Data_450K <- TCGA_Load_MethylationData(METdirectories$METdirectory450k,
            ArrayType = "450K")
    }

    # check if we want to combine 27k and 450k
    MET_Data <- TCGA_Select_Dataset(CancerSite = CancerSite, MET_Data_27K = MET_Data_27K,
        MET_Data_450K = MET_Data_450K, use450K = use450K)
    rm(MET_Data_27K)
    gc()
    rm(MET_Data_450K)
    gc()

    # ---------------------------------------------------------------------------------------------
    # Step 2: Split up normal and cancer data
    # ---------------------------------------------------------------------------------------------
    Samplegroups <- TCGA_GENERIC_GetSampleGroups(colnames(MET_Data))
    if (CancerSite == "LAML") {
        MET_Data_Cancer <- MET_Data[, Samplegroups$PeripheralBloodCancer]
    } else {
        MET_Data_Cancer <- MET_Data[, Samplegroups$Primary]
    }
    MET_Data_Cancer <- TCGA_GENERIC_CleanUpSampleNames(MET_Data_Cancer, 15)
    if (CancerSite == "LAML") {
        MET_Data_Normal <- MET_Data[, Samplegroups$BloodNormal]
    } else {
        MET_Data_Normal <- MET_Data[, Samplegroups$SolidNormal]
    }
    if (length(MET_Data_Normal) > 0) {
        MET_Data_Normal <- TCGA_GENERIC_CleanUpSampleNames(MET_Data_Normal, 15)
    }
    cat("There are", length(colnames(MET_Data_Cancer)), "cancer samples and", length(colnames(MET_Data_Normal)),
        "normal samples with DNA methylation data in", CancerSite, "\n")

    # Clear space
    rm(MET_Data)
    gc()

    # ---------------------------------------------------------------------------------------------
    # Step 3: Perform Missing value estimation
    # ---------------------------------------------------------------------------------------------
    cat("\tMissing value estimation for the cancer samples.\n")
    MET_Data_Cancer <- TCGA_Process_EstimateMissingValues(MET_Data_Cancer, MissingValueThreshold)
    if (length(MET_Data_Normal) > 0) {
        cat("\tMissing value estimation for the normal samples.\n")
        MET_Data_Normal <- TCGA_Process_EstimateMissingValues(MET_Data_Normal, MissingValueThreshold)
    }

    # ---------------------------------------------------------------------------------------------
    # Step 4: Perform batch correction
    # ---------------------------------------------------------------------------------------------
    if (doBatchCorrection) {
        BatchData <- EpiMix_GetData("TCGA_BatchData")
        MinPerBatch <- 5
        cat("\tBatch correction for the cancer samples.\n")

        # p_init = TCGA_GENERIC_CheckBatchEffect(MET_Data_Cancer, BatchData) #
        # p = 0.054

        # Remove samples with batch number 0
        if (length(-which(BatchData[, 2] == 0)) > 0) {
            BatchData <- BatchData[-which(BatchData[, 2] == 0), ]
        }

        MET_Data_Cancer <- CorrectBatchEffect(GEN_Data = MET_Data_Cancer, BatchData = BatchData,
            batch.correction.method = batch.correction.method, MinInBatch = MinPerBatch,
            featurePerSet = 50000)

        p_after <- TCGA_GENERIC_CheckBatchEffect(MET_Data_Cancer, BatchData)  # p = 0.999

        if (length(MET_Data_Normal) > 0) {
            cat("Batch correction for the normal samples.\n")
            MET_Data_Normal <- CorrectBatchEffect(GEN_Data = MET_Data_Normal, BatchData = BatchData,
                batch.correction.method = batch.correction.method, MinInBatch = MinPerBatch,
                featurePerSet = 50000)
        } else {
            MET_Data_Normal <- c()
        }
        # Set values <0 to 0 and >1 to 1, because of batch correction
        MET_Data_Cancer[MET_Data_Cancer < 0] <- 0
        MET_Data_Cancer[MET_Data_Cancer > 1] <- 1
        if (length(MET_Data_Normal) > 0) {
            MET_Data_Normal[MET_Data_Normal < 0] <- 0
            MET_Data_Normal[MET_Data_Normal > 1] <- 1
        }
    }

    if (cores > 1 & batch.correction.method == "Combat") {
        parallel::stopCluster(cl)
    }

    # ---------------------------------------------------------------------------------------------
    # Step 5: Merge cancer and normal data into one matrix
    # ---------------------------------------------------------------------------------------------
    # The EpiMix only takes one methylation matrix, so we combine cancer and
    # normal into one matrix
    if (!is.null(MET_Data_Normal)) {
        overlapProbes <- intersect(rownames(MET_Data_Cancer), rownames(MET_Data_Normal))
        MET_Data_Cancer <- MET_Data_Cancer[overlapProbes, , drop = FALSE]
        MET_Data_Normal <- MET_Data_Normal[overlapProbes, , drop = FALSE]
        MET_Data <- cbind(MET_Data_Cancer, MET_Data_Normal)
        return(MET_Data)
    } else {
        return(MET_Data_Cancer)
    }
}

#' The TCGA_Load_MethylationData function
#' @details load 27K or 450K methyaltion data into memory
#' @param METdirectory path to the 27K or 450K data
#' @param ArrayType character string indicating the array type, can be either '27K' or '450K'
#' @return matrix of methylation data with probes in rows and patient in columns
#' @keywords internal
#'

TCGA_Load_MethylationData <- function(METdirectory, ArrayType) {

    prefix <- NULL
    if (ArrayType == "27K") {
        prefix <- "methylation__humanmethylation27"
    } else {
        prefix <- "methylation__humanmethylation450"
    }

    if (grepl("Windows", Sys.info()["sysname"])) {
        # If Windows I'll create a virtual drive to handle the long file names
        # issue Create a virtual drive to overcome long names issue in Windows
        virtualDir <- METdirectory
        virtualDir <- gsub("\\\\", "/", virtualDir)
        virtualDir <- substr(virtualDir, 1, nchar(virtualDir) - 1)
        system(paste("subst x:", virtualDir))

        METfiles <- dir("x:")
        MatchedFilePosition <- grep(prefix, METfiles)
        Filename <- paste0("x://", METfiles[MatchedFilePosition])
        MET_Data <- TCGA_GENERIC_LoadIlluminaMethylationData(Filename)

        system("subst x: /D")  #stop virtual drive
    } else {
        # Not windows Load data
        METfiles <- dir(METdirectory)
        MatchedFilePosition <- grep(prefix, METfiles)
        Filename <- paste0(METdirectory, METfiles[MatchedFilePosition])
        MET_Data <- TCGA_GENERIC_LoadIlluminaMethylationData(Filename)
    }
    return(MET_Data)
}


#' The TCGA_GENERIC_LoadIlluminaMethylationData function
#'
#' Internal. Read in an illumina methylation file with the following format: header row with sample labels,
#' 2nd header row with 4 columns per sample: beta-value, geneSymbol, chromosome and GenomicCoordinate.
#' The first column has the probe names.
#' @param Filename name of the file with the data.
#' @return methylation data.
#' @keywords internal
#'
TCGA_GENERIC_LoadIlluminaMethylationData <- function(Filename) {

    # read in an illumina methylation file with the following format: header
    # row with sample labels 2nd header row with 4 columns per sample:
    # beta-value, geneSymbol, chromosome and GenomicCoordinate The first column
    # has the probe names.
    MET_Data <- data.table::fread(Filename)
    MET_Data <- as.matrix(MET_Data)
    Probes <- MET_Data[, 1]
    rownames(MET_Data) <- Probes
    MET_Data <- MET_Data[, -1]
    MET_Data <- MET_Data[-1, ]
    MET_Data <- MET_Data[, seq(1, ncol(MET_Data), 4)]
    class(MET_Data) <- "numeric"

    return(MET_Data)
}

#' The TCGA_GENERIC_GetSampleGroups function
#'
#' Internal. Looks for the group of the samples (normal/cancer).
#' @param SampleNames vector with sample names.
#' @return a list.
#' @keywords internal
#'
TCGA_GENERIC_GetSampleGroups <- function(SampleNames) {

    # First replace any . with - so the sample groups are uniform.
    SampleGroups <- list()

    # 1: Primary Tumor
    Matches <- regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]01[.|-]*", SampleNames,
        perl = FALSE, useBytes = FALSE)
    SampleGroups$Primary <- SampleNames[Matches == 1]

    # 2: Recurrent tumor
    Matches <- regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]02[.|-]*", SampleNames,
        perl = FALSE, useBytes = FALSE)
    SampleGroups$Recurrent <- SampleNames[Matches == 1]

    # 3: Primary blood derived cancer - peripheral blood
    Matches <- regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]03[.|-]*", SampleNames,
        perl = FALSE, useBytes = FALSE)
    SampleGroups$PeripheralBloodCancer <- SampleNames[Matches == 1]

    # 10: Blood derived normal
    Matches <- regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]10[.|-]*", SampleNames,
        perl = FALSE, useBytes = FALSE)
    SampleGroups$BloodNormal <- SampleNames[Matches == 1]

    # 11: Solid tissue derived normal
    Matches <- regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]11[.|-]*", SampleNames,
        perl = FALSE, useBytes = FALSE)
    SampleGroups$SolidNormal <- SampleNames[Matches == 1]

    # 20 Cellines
    Matches <- regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]20[.|-]*", SampleNames,
        perl = FALSE, useBytes = FALSE)
    SampleGroups$CellLines <- SampleNames[Matches == 1]

    # 06 Cellines
    Matches <- regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]06[.|-]*", SampleNames,
        perl = FALSE, useBytes = FALSE)
    SampleGroups$Metastatic <- SampleNames[Matches == 1]

    return(SampleGroups)
}

#' The TCGA_GENERIC_CleanUpSampleNames function
#'
#' Internal. Cleans the samples IDs into the 12 digit format and removes doubles.
#' @param GEN_Data data matrix.
#' @param IDlength length of samples ID.
#' @return data matrix with cleaned sample names.
#' @keywords internal
#'
TCGA_GENERIC_CleanUpSampleNames <- function(GEN_Data, IDlength = 12) {
    SampleNames <- colnames(GEN_Data)
    SampleNamesShort <- as.character(apply(as.matrix(SampleNames), 2, substr, 1,
        IDlength))
    if (length(SampleNamesShort) != length(unique(SampleNamesShort))) {
        # remove the doubles
        Counts <- table(SampleNamesShort)
        Doubles <- rownames(Counts)[which(Counts > 1)]

        cat("Removing doubles for", length(Doubles), "samples.\n")
        for (i in 1:length(Doubles)) {
            CurrentDouble <- Doubles[i]
            pos <- grep(CurrentDouble, SampleNames)
            # GEN_Data[1:10,pos] cor(GEN_Data[,pos])
            GEN_Data <- GEN_Data[, -pos[2:length(pos)]]
            SampleNames <- colnames(GEN_Data)  # need to update samplenames because pos is relative to this
        }
        SampleNames <- colnames(GEN_Data)
        SampleNamesShort <- as.character(apply(as.matrix(SampleNames), 2, substr,
            1, IDlength))

        # now set the samplenames
        colnames(GEN_Data) <- SampleNamesShort
    } else {
        colnames(GEN_Data) <- SampleNamesShort
    }
    return(GEN_Data)
}

#' The TCGA_Process_EstimateMissingValues function
#'
#' Internal. Removes patients and genes with more missing values than the MissingValueThreshold, and imputes remaining missing values using Tibshirani's KNN method.
#' @param MET_Data data matrix.
#' @param MissingValueThreshold threshold for removing samples and genes with too many missing values.
#' @return the data set with imputed values and possibly some genes or samples deleted.
#' @keywords internal
#'
TCGA_Process_EstimateMissingValues <- function(MET_Data, MissingValueThreshold = 0.2) {

    # FIRST REMOVING BAD PATIENTS removing patients with too many missings
    # values
    NrMissingsPerSample <- apply(MET_Data, 2, function(x) sum(is.na(x)))/nrow(MET_Data)
    cat("Removing", sum(NrMissingsPerSample > MissingValueThreshold), "patients with more than",
        MissingValueThreshold * 100, "% missing values.\n")
    if (sum(NrMissingsPerSample > MissingValueThreshold) > 0)
        MET_Data <- MET_Data[, NrMissingsPerSample < MissingValueThreshold, drop = FALSE]

    # removing clones with too many missing values
    NrMissingsPerGene <- apply(MET_Data, 1, function(x) sum(is.na(x)))/ncol(MET_Data)
    cat("Removing", sum(NrMissingsPerGene > MissingValueThreshold), "probes with more than",
        MissingValueThreshold * 100, "% missing values.\n")
    if (sum(NrMissingsPerGene > MissingValueThreshold) > 0)
        MET_Data <- MET_Data[NrMissingsPerGene < MissingValueThreshold, , drop = FALSE]

    # knn impute using Tibshirani's method
    if (length(colnames(MET_Data)) > 1) {
        k <- 15
        KNNresults <- impute::impute.knn(as.matrix(MET_Data), k)
        MET_Data_KNN <- KNNresults$data

        # cleaning up sample names
        return(MET_Data_KNN)

    } else {
        # when only 1 sample,need to make a matrix again
        # MET_Data=as.matrix(MET_Data)
        return(MET_Data)
    }
}

#' The TCGA_GENERIC_CheckBatchEffect function
#'
#' Internal. Checks if batch correction is needed.
#' @param GEN_Data matrix with data to be corrected for batch effects.
#' @param BatchData Batch data.
#' @return the p value from ANOVA test on PCA values.
#' @keywords internal
#'
TCGA_GENERIC_CheckBatchEffect <- function(GEN_Data, BatchData) {
    # first match the samples to the batch
    Order <- match(colnames(GEN_Data), BatchData[, 1])
    BatchDataSelected <- BatchData[Order, ]
    BatchDataSelected[, 1] <- factor(BatchDataSelected[, 1])

    # PCA analysis alternatively use fast.prcomp from package gmodels, but
    # tests do not show this is faster
    PCAanalysis <- prcomp(t(GEN_Data))
    PCdata <- PCAanalysis$x
    # plot(PCdata[,1]~BatchDataSelected[,2])
    pval_init <- 1
    if (length(unique(BatchDataSelected[, 1][!is.na(BatchDataSelected[, 1])])) >
        1) {
        tmp <- aov(PCdata[, 1] ~ BatchDataSelected[, 2])
        test.results <- summary(tmp)
        pval_init <- test.results[[1]][["Pr(>F)"]][[1]]
    }
    return(pval_init)
}

#' The CorrectBatchEffect  function
#' @description  top-level wrapper function for batch correction.
#' @details
#' (1) filters the batch data and the molecular data to keep only the overlapped samples.
#' (2) removes extremely small batches.
#' (3) if the molecular data have over 50,000 features (rows), it splits the data into subsets, with 50,000 features in each subset, and perform batch correction on each subset.
#' (4) identify overlapped samples in batch corrected subsets, and merge the subsets into one matrix.
#' @param GEN_Data matrix with methylation.data or gene.expression.data with genes in rows and samples in columns
#' @param BatchData dataframe with two columns: the first column indicates the sample names, and the second column indicates the batch ids.
#' @param batch.correction.method character string. Should be either 'Seurat' or 'Combat'.
#' @param MinInBatch integer indicating the batch size threshold. Batches smaller than this threshold will be removed. Default: 5
#' @param featurePerSet integer indicating the row numbers to split the GEN_Data into small subsets. Default: 50,000
#' @return matrix with corrected data
#' @keywords internal

CorrectBatchEffect <- function(GEN_Data, BatchData, batch.correction.method, MinInBatch = 5,
    featurePerSet = 50000) {

    # Perform some input check
    if (is.null(BatchData) | BatchData == "") {
        stop("BatchData is NULL, please provide batch information.\n")
    }

    if (!batch.correction.method %in% c("Seurat", "Combat")) {
        stop("batch.correction.method must be either 'Seurat' or 'Combat'.\n")
    }

    # ---------------------------------------------------------------------------------------------
    # Step 1: Filter the batch data and the GEN_data to keep only samples
    # present in both sets
    # ---------------------------------------------------------------------------------------------
    overlapSamples <- intersect(colnames(GEN_Data), BatchData[, 1])
    cat("Found", length(overlapSamples), "samples with batch information\n")
    if (length(overlapSamples) < ncol(GEN_Data)) {
        warning("There are", ncol(GEN_Data) - length(overlapSamples), "samples that do not have batch information. These samples will be excluded in the downstream analyses !!!\n")
    }
    if (length(overlapSamples) == 0) {
        stop("No overlap samples were found between the batch data and the methylation data or gene expression data!.\nThe first column of the BatchData must overlap the column names for the methylation data or gene expression data.\n")
    }

    BatchDataSelected <- BatchData[which(BatchData[, 1] %in% overlapSamples), ]
    GEN_Data <- GEN_Data[, overlapSamples, drop = FALSE]

    # ---------------------------------------------------------------------------------------------
    # Step 2: Filter out small-sized batches (batches with size < MinInBatch)
    # ---------------------------------------------------------------------------------------------

    BatchDataSelected[, 2] <- factor(BatchDataSelected[, 2])

    NrPerBatch <- table(BatchDataSelected[, 2])
    SmallBatches <- NrPerBatch < MinInBatch
    BatchesToBeRemoved <- names(SmallBatches)[which(SmallBatches == TRUE)]
    SamplesToBeRemoved <- as.character(BatchDataSelected[which(BatchDataSelected[,
        1] %in% BatchesToBeRemoved), 1])

    if (length(SamplesToBeRemoved) > 0) {
        cat("Removing", length(which(colnames(GEN_Data) %in% SamplesToBeRemoved)),
            "samples because their batches are too small.\n")
        GEN_Data <- GEN_Data[, -which(colnames(GEN_Data) %in% SamplesToBeRemoved)]
        BatchDataSelected <- BatchDataSelected[-which(BatchDataSelected[, 1] %in%
            SamplesToBeRemoved), ]
    }

    if (ncol(GEN_Data) < MinInBatch) {
        # just checking if we have enough samples after removing the too small
        # batches
        cat("The number of samples becomes to small, no batch correction possible.\n")
        return(GEN_Data)
    }

    if (length(unique(BatchDataSelected[, 2])) == 1) {
        cat("Only one batch, no batch correction possible.\n")
        return(GEN_Data)
    }

    # Seurat expects more than 30 samples in a batch to perform batch
    # correction
    if (ncol(GEN_Data) <= 30 & batch.correction.method == "Seurat") {
        warning("The total number of samples is smaller than 30, using Combat for batch correction instead of Seurat.")
        batch.correction.method <- "Combat"
    }

    # ---------------------------------------------------------------------------------------------
    # Step 3: Split data by rows into small subsets, and perform batch
    # correction on each subset
    # ---------------------------------------------------------------------------------------------
    n <- featurePerSet
    data_subsets <- list()
    counter <- 1
    for (i in seq(from = 0, to = nrow(GEN_Data), by = n)) {
        if (i == nrow(GEN_Data))
            break
        if (i + 1 <= nrow(GEN_Data) & n + i > nrow(GEN_Data)) {
            data_subsets[[counter]] <- GEN_Data[(i + 1):nrow(GEN_Data), ]
            break
        }
        data_subsets[[counter]] <- GEN_Data[(i + 1):(n + i), ]
        counter <- counter + 1
    }
    cat("Found", nrow(GEN_Data), "features.\n")
    cat("Split the dataset into", length(data_subsets), "subsets.\n")

    rm(GEN_Data)
    gc()

    # Perform batch correction for each subset of data
    data_counter <- 1
    integrated_data <- list()
    for (i in 1:length(data_subsets)) {
        cat("Batch correction on subset", i, "\n")
        data <- data_subsets[[data_counter]]
        if (batch.correction.method == "Seurat") {
            data <- BatchCorrection_Seurat(data, BatchDataSelected)
        } else {
            data <- BatchCorrection_Combat(data, BatchDataSelected)
        }
        integrated_data[[data_counter]] <- data
        data_counter <- data_counter + 1
    }

    rm(data_subsets)
    gc()

    # ---------------------------------------------------------------------------------------------
    # Step 4: Combine batch corrected data
    # ---------------------------------------------------------------------------------------------
    BatchCorrectedData <- NULL
    if (length(integrated_data) > 1) {
        cat("Finished batch correction, merging subsets...\n")
        overlapSamples <- colnames(integrated_data[[1]])
        for (i in 2:length(integrated_data)) {
            overlapSamples <- intersect(overlapSamples, colnames(integrated_data[[i]]))
        }
        cat("Found", length(overlapSamples), "overlap samples in different subsets.\n")
        BatchCorrectedData <- integrated_data[[1]][, overlapSamples, drop = FALSE]
        for (i in 2:length(integrated_data)) {
            BatchCorrectedData <- rbind(BatchCorrectedData, integrated_data[[i]][,
                overlapSamples, drop = FALSE])
        }
    } else (BatchCorrectedData <- integrated_data[[1]])
    return(BatchCorrectedData)
}


#' The BatchCorrection_Seurat function
#' @details correct batch effects with the Seurat data integration functions.
#' @param GEN_Data matrix with methylation.data or gene.expression.data
#' @param BatchDataSelected BatchData after filtering out the small batches and selecting for overlapped samples.
#' @return corrected data matrix
#' @keywords internal

BatchCorrection_Seurat <- function(GEN_Data, BatchDataSelected) {

    if (!requireNamespace("Seurat")) {
        message("To do batch correction with Seurat, you need to install the 'Seurat' R package")
        return(invisible())
    }

    # Setting some hyper-parameters
    feature_percentage <- 0.1
    dims <- 30
    k.weight <- 30

    # Seurat requires the sample size to be greater than 30 per batch,
    # therefore we have to merge small batches into some super-clusters
    min_batch <- 30
    NrPerBatch <- table(BatchDataSelected[, 2])
    NrPerBatch <- sort(NrPerBatch, decreasing = TRUE)
    cat("Found", length(NrPerBatch), "batches.\n")

    SingleBatches <- names(NrPerBatch)[which(NrPerBatch > min_batch)]

    # Create BatchInfo dictionary for Merging Batches
    Batch_Counter <- 1
    BatchInfo <- c()
    j_max <- length(NrPerBatch)
    for (i in 1:length(NrPerBatch)) {
        # When the remaining batches do not add up to 31,incorporate all
        # remaining batches to “BATCH 1'
        if (sum(NrPerBatch[i:j_max]) < 31) {
            for (k in i:j_max) {
                BatchtoMerge <- toString(names(NrPerBatch[k]))
                BatchInfo[[BatchtoMerge]] <- 1
            }
            break
        }
        CurrentBatch <- toString(names(NrPerBatch[i]))
        print(paste0("Working on Batch ", CurrentBatch))
        if (CurrentBatch %in% SingleBatches) {
            print(paste0("Batch ", CurrentBatch, " Processed"))
            BatchInfo[[CurrentBatch]] <- Batch_Counter
            Batch_Counter <- Batch_Counter + 1
            processed <- TRUE
            (next)()
        }
        for (j in j_max:(i + 1)) {
            BatchtoMerge <- toString(names(NrPerBatch[j]))
            num <- as.integer(as.numeric((NrPerBatch[i]) + as.numeric(NrPerBatch[j])))
            print(paste0("The total number of samples in ", CurrentBatch, " and ",
                BatchtoMerge, " is ", num))
            if (num > 30) {
                print(paste0("Batch ", CurrentBatch, " Processed"))
                BatchInfo[[CurrentBatch]] <- Batch_Counter
                BatchInfo[[BatchtoMerge]] <- Batch_Counter
                Batch_Counter <- Batch_Counter + 1
                j_max <- j - 1
                processed <- TRUE
                break
            }
        }
        if (processed == FALSE) {
            BatchInfo[[CurrentBatch]] <- 0
        }
        if (j_max == i)
            break
    }

    BatchDataSelected$MergedBatch <- unlist(BatchInfo)[match(BatchDataSelected[,
        2], names(BatchInfo))]

    # Remove samples
    SamplesToBeRemoved <- as.character(BatchDataSelected[which(BatchDataSelected$MergedBatch ==
        "0"), 1])
    if (length(SamplesToBeRemoved) > 0) {
        cat("Removing", length(which(colnames(GEN_Data) %in% SamplesToBeRemoved)))
        GEN_Data <- GEN_Data[, -which(colnames(GEN_Data) %in% SamplesToBeRemoved)]
    }
    BatchDataSelected <- BatchDataSelected[which(BatchDataSelected[, "MergedBatch"] !=
        "0"), ]

    BatchDataSelected <- BatchDataSelected[order(BatchDataSelected$MergedBatch),
        ]
    rownames(BatchDataSelected) <- seq(1:nrow(BatchDataSelected))

    # Perform batch correction on merged batches
    if (length(unique(BatchDataSelected$MergedBatch)) == 1) {
        # only one batch after merging the small-sized batches, no batch
        # correction is possible with Seurat, using Combat for batch
        # correction.
        warning("There is only 1 batch after merging the small-sized batches, no batch correction with Seurat is possible. Using Combat for batch correction instead...\n")
        GEN_Data <- BatchCorrection_Combat(GEN_Data, BatchDataSelected)
        return(GEN_Data)
    }

    seurat_objects <- list()
    batch_counter <- 1
    for (batch in unique(BatchDataSelected$MergedBatch)) {
        batch_samples <- BatchDataSelected[BatchDataSelected$MergedBatch == batch,
            1]
        data <- GEN_Data[, as.character(batch_samples), drop = FALSE]
        seurat_objects[[batch_counter]] <- Seurat::CreateSeuratObject(counts = data,
            project = paste0("MergedBatch_", batch))
        batch_counter <- batch_counter + 1
    }

    feature_number <- nrow(GEN_Data)
    features <- rownames(GEN_Data)
    rm(GEN_Data)
    gc()

    cat("==========================================================\n")
    cat("Finding variable features on Seurat Objects...\n")
    start <- Sys.time()
    for (i in 1:length(seurat_objects)) {
        seurat_objects[[i]] <- Seurat::FindVariableFeatures(seurat_objects[[i]],
            selection.method = "vst", nfeatures = feature_number * feature_percentage,
            verbose = FALSE)
    }
    end <- Sys.time()
    find.variable.time <- as.numeric(end - start, units = "mins")
    cat("Find variable features in minutes: ", find.variable.time, "\n")

    cat("==========================================================\n")
    cat("Finding Integration Anchors and Integrating Data on Seurat Objects.\n")
    start <- Sys.time()
    data_anchors <- Seurat::FindIntegrationAnchors(object.list = seurat_objects,
        dims = 1:dims)
    end <- Sys.time()
    find.anchor.time <- as.numeric(end - start, units = "mins")
    cat("Find anchor time in minutes: ", find.anchor.time, "\n")

    # The k.weight parameter needs to be less than the number of samples.
    start <- Sys.time()
    data_integrated <- Seurat::IntegrateData(anchorset = data_anchors, dims = 1:dims,
        k.weight = k.weight, features.to.integrate = features)
    end <- Sys.time()
    integrate.data.time <- as.numeric(end - start, units = "mins")
    cat("Intergrate data time in minutes: ", integrate.data.time, "\n")
    data_integrated <- Seurat::GetAssayData(object = data_integrated, slot = "data")

    data_integrated <- as.matrix(data_integrated)

    # The dataset may have NA values after batch correction, perform missing
    # value estimation again.
    if (sum(is.na(data_integrated)) > 0)
        data_integrated <- TCGA_Process_EstimateMissingValues(data_integrated)
    return(data_integrated)
}

#' The BatchCorrection_Combat function
#' @details correct batch effects with Combat
#' @param GEN_Data matrix with methylation.data or gene.expression.data
#' @param BatchDataSelected BatchData after filtering out the small batches and selecting for overlapped samples
#' @return corrected data matrix
#' @keywords internal

BatchCorrection_Combat <- function(GEN_Data, BatchDataSelected) {

    BatchDataSelected[, 2] <- factor(BatchDataSelected[, 2])
    BatchDataSelected[, 1] <- factor(BatchDataSelected[, 1])

    # reordering samples (not really necessary as Combat does this too)
    order <- match(colnames(GEN_Data), BatchDataSelected[, 1])
    BatchDataSelected <- BatchDataSelected[order, ]
    BatchDataSelected[, 2] <- factor(BatchDataSelected[, 2])

    # running combat
    CombatResults <- ComBat_NoFiles(GEN_Data, BatchDataSelected)
    data_integrated <- CombatResults[, -1]
    class(data_integrated) <- "numeric"

    # The dataset may have NA values after batch correction, perform missing
    # value estimation again.
    if (sum(is.na(data_integrated)) > 0)
        data_integrated <- TCGA_Process_EstimateMissingValues(data_integrated)

    return(data_integrated)
}

#' The TCGA_Download_GeneExpression function
#'
#' @description Download gene expression data from TCGA.
#' @param CancerSite character string indicating the TCGA cancer code.
#' @param TargetDirectory character with directory where a folder for downloaded files will be created.
#' @param mode character string indicating whether we should download the gene expression data for miRNAs or lncRNAs, instead of for protein-coding genes. See details for more information.
#' @param downloadData logical indicating if the data should be downloaded (default: TRUE). If False, the url of the desired data is returned.
#' @return list with paths to downloaded files for gene expression.
#' @details
#' mode: when mode is set to 'Regular', this function downloads the level 3 RNAseq data (file tag 'mRNAseq_Preprocess.Level_3'). Since there is not enough RNAseq data for OV and GBM, the micro array data is
#' downloaded. If you plan to run the EpiMix on miRNA- or lncRNA-coding genes, please specify the 'mode' parameter to 'miRNA' or 'lncRNA'.
#' @export
#' @keywords download
#' @examples
#' \donttest{
#' # Example #1 : download regular gene expression data for ovarian cancer
#' GEdirectories <- TCGA_Download_GeneExpression(CancerSite = 'OV', TargetDirectory = tempdir())
#'
#' # Example #2 : download miRNA expression data for ovarian cancer
#' GEdirectories <- TCGA_Download_GeneExpression(CancerSite = 'OV',
#'                                               TargetDirectory = tempdir(),
#'                                               mode = 'miRNA')
#'
#' # Example #3 : download lncRNA expression data for ovarian cancer
#' GEdirectories <- TCGA_Download_GeneExpression(CancerSite = 'OV',
#'                                                TargetDirectory = tempdir(),
#'                                                mode = 'lncRNA')
#' }

#'
TCGA_Download_GeneExpression <- function(CancerSite, TargetDirectory, mode = "Regular",
    downloadData = TRUE) {

    mode <- tolower(mode)
    if (!mode %in% c("regular", "enhancer", "mirna", "lncrna")) {
        stop("'mode' must be one of the followings: 'Regular', 'Enhancer', 'miRNA', 'lncRNA'.Please specify a correct input.\n")
    }

    options(timeout = 1e+05)
    dir.create(TargetDirectory, showWarnings = FALSE)
    dataType <- dataFileTag <- NULL

    # Settings
    TCGA_acronym_uppercase <- toupper(CancerSite)

    dataType <- "stddata"

    if (mode == "regular" || mode == "enhancer") {
        cat("Searching MA data for:", CancerSite, "\n")
        # get RNA seq data (GBM does not have much RNAseq data.)
        dataFileTag <- "mRNAseq_Preprocess.Level_3"
        # special case for GBM and OV, not enough RNAseq data, so using the
        # microarray data instead
        if (CancerSite == "GBM") {
            dataFileTag <- c("Merge_transcriptome__agilentg4502a_07_1__unc_edu__Level_3__unc_lowess_normalization_gene_level__data",
                "Merge_transcriptome__agilentg4502a_07_2__unc_edu__Level_3__unc_lowess_normalization_gene_level__data")
        } else if (CancerSite == "OV") {
            dataFileTag <- "Merge_transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data"
        }
        if (length(dataFileTag) == 1) {
            MAdirectories <- get_firehoseData(downloadData, saveDir = TargetDirectory,
                TCGA_acronym_uppercase = TCGA_acronym_uppercase, dataFileTag = dataFileTag)
        } else {
            # a few data sets have multiple gene expression data sets.
            MAdirectories <- c()
            for (i in 1:length(dataFileTag)) {
                MAdirectories <- c(MAdirectories, get_firehoseData(downloadData,
                  saveDir = TargetDirectory, TCGA_acronym_uppercase = TCGA_acronym_uppercase,
                  dataFileTag = dataFileTag[i]))
            }
        }
    } else if (mode == "mirna") {
        cat("Searching miRNA data for:", CancerSite, "\n")
        dataFileTag <- "miRseq_Preprocess.Level_3"
        if (length(dataFileTag) == 1) {
            MAdirectories <- get_firehoseData(downloadData, saveDir = TargetDirectory,
                TCGA_acronym_uppercase = TCGA_acronym_uppercase, dataFileTag = dataFileTag)
        } else {
            # a few data sets have multiple gene expression data sets.
            MAdirectories <- c()
            for (i in 1:length(dataFileTag)) {
                MAdirectories <- c(MAdirectories, get_firehoseData(downloadData,
                  saveDir = TargetDirectory, TCGA_acronym_uppercase = TCGA_acronym_uppercase,
                  dataFileTag = dataFileTag[i]))
            }
        }
    } else if (mode == "lncrna") {
        cat("Downloading lncRNA data for:", CancerSite, "\n")
        MAdirectories <- getLncRNAData(CancerSite)
    }
    return(MAdirectories = MAdirectories)
}

#' The TCGA_Preprocess_GeneExpression function
#'
#' @description Pre-processes gene expression data from TCGA.
#' @param CancerSite character string indicating the TCGA cancer code.
#' @param MAdirectories character vector with directories with the downloaded data. It can be the object returned by the GEO_Download_GeneExpression function.
#' @param mode character string indicating whether the genes in the gene expression data are miRNAs or lncRNAs. Should be either 'Regular', 'Enhancer', 'miRNA' or 'lncRNA'. This value should be consistent with the same parameter in the TCGA_Download_GeneExpression function. Default: 'Regular'.
#' @param doBatchCorrection logical indicating whether to perform batch effect correction. Default: False.
#' @param batch.correction.method character string indicating the method to perform batch correction. The value should be either 'Seurat' or 'Combat'. Default: 'Seurat'. Seurat is much fatster than the Combat.
#' @param MissingValueThresholdGene threshold for missing values per gene. Genes with a percentage of NAs greater than this threshold are removed. Default is 0.3.
#' @param MissingValueThresholdSample threshold for missing values per sample. Samples with a percentage of NAs greater than this threshold are removed. Default is 0.1.
#' @param cores integer indicating the number of cores to be used for performing batch correction with Combat
#' @details
#' Pre-process includes eliminating samples and genes with too many NAs, imputing NAs, and doing Batch correction. If the rownames of the gene expression data are ensembl ENSG names or ENST names, the function will convert them to the human gene symbol (HGNC).
#' @return pre-processed gene expression data matrix.
#' @export
#' @keywords preprocess
#' @examples
#' \donttest{
#'
#' # Example #1: Preprocessing gene expression for Regular mode
#'
#'  GEdirectories <- TCGA_Download_GeneExpression(CancerSite = 'OV', TargetDirectory = tempdir())
#'  GEProcessedData <- TCGA_Preprocess_GeneExpression(CancerSite = 'OV',  MAdirectories = GEdirectories)
#'
#' # Example #2: Preprocessing gene expression for miRNA mode
#'
#'  GEdirectories <- TCGA_Download_GeneExpression(CancerSite = 'OV',
#'                                                TargetDirectory = tempdir(),
#'                                                mode = 'miRNA')
#'
#'  GEProcessedData <- TCGA_Preprocess_GeneExpression(CancerSite = 'OV',
#'                                                    MAdirectories = GEdirectories,
#'                                                    mode = 'miRNA')
#'
#' # Example #3: Preprocessing gene expression for lncRNA mode
#'
#'  GEdirectories <- TCGA_Download_GeneExpression(CancerSite = 'OV',
#'                                                TargetDirectory = tempdir(),
#'                                                mode = 'lncRNA')
#'
#'  GEProcessedData <- TCGA_Preprocess_GeneExpression(CancerSite = 'OV',
#'                                                    MAdirectories = GEdirectories,
#'                                                    mode = 'lncRNA')
#'
#' }
#'
TCGA_Preprocess_GeneExpression <- function(CancerSite, MAdirectories, mode = "Regular",
    doBatchCorrection = FALSE, batch.correction.method = "Seurat", MissingValueThresholdGene = 0.3,
    MissingValueThresholdSample = 0.1, cores = 1) {

    # set up some paramerters
    if (cores > 1 & batch.correction.method == "Combat") {
        # unregister()
        cat("Registering sockets on multiple CPU cores...\n")
        cl <- parallel::makeCluster(cores)
    }

    # ---------------------------------------------------------------------------------------------
    # Step 1: Load gene expression data
    # ---------------------------------------------------------------------------------------------

    BatchData <- EpiMix_GetData("TCGA_BatchData")
    MinPerBatchCancer <- 5
    MinPerBatchNormal <- 2

    mode <- tolower(mode)
    if (!mode %in% c("regular", "enhancer", "mirna", "lncrna")) {
        stop("'mode' must be one of the followings: 'Regular', 'Enhancer', 'miRNA', 'lncRNA'.Please specify a correct input.\n")
    }

    if (mode == "regular" | mode == "mirna") {
        if (mode == "mirna") {
            MAstring <- "miRseq_RPKM_log2.txt"
        } else {
            # Processing MA data, special case for OV and GBM where no RNA seq
            # data is available
            if (CancerSite == "OV" || CancerSite == "GBM") {
                MAstring <- "transcriptome__agilent"
            } else if (CancerSite == "STAD" || CancerSite == "ESCA") {
                # for these cancers RSEM data does not exist.
                MAstring <- "mRNAseq_RPKM_log2.txt"
            } else {
                MAstring <- "mRNAseq_RSEM_normalized_log2.txt"
            }
        }

        if (grepl("Windows", Sys.info()["sysname"])) {
            # If Windows I'll create a virtual drive to handle the long file
            # names issue Create a virtual drive to overcome long names issue
            # in Windows
            MAdirectoriesOrig <- MAdirectories
            virtualDir <- MAdirectories
            virtualDir <- gsub("\\\\", "/", virtualDir)
            virtualDir <- substr(virtualDir, 1, nchar(virtualDir) - 1)
            system(paste("subst x:", virtualDir))
            MAdirectories <- "x://"
        }

        # ---------------------------------------------------------------------------------------------
        # Step 2: Preprocess gene expression data
        # ---------------------------------------------------------------------------------------------

        if (length(MAdirectories) > 1) {
            cat("\tFound multiple MA data sets.\n")
            DataSetsCancer <- list()
            GeneListsCancer <- list()
            SampleListsCancer <- list()
            MetaBatchDataCancer <- data.frame()

            DataSetsNormal <- list()
            GeneListsNormal <- list()
            SampleListsNormal <- list()
            MetaBatchDataNormal <- data.frame()
            for (i in 1:length(MAdirectories)) {
                cat("\tProcessing data set", i, "\n")
                MAfiles <- dir(MAdirectories[i])
                MatchedFile <- grep(MAstring, MAfiles)
                if (length(MatchedFile) > 0) {
                  # Getting the cancer data first
                  DataSetsCancer[[i]] <- Preprocess_MAdata_Cancer(CancerSite, MAdirectories[i],
                    MAfiles[MatchedFile], MissingValueThresholdGene = MissingValueThresholdGene,
                    MissingValueThresholdSample = MissingValueThresholdSample, doBatchCorrection = doBatchCorrection,
                    batch.correction.method = batch.correction.method, BatchData = BatchData)
                  GeneListsCancer[[i]] <- rownames(DataSetsCancer[[i]])
                  SampleListsCancer[[i]] <- colnames(DataSetsCancer[[i]])
                  currentBatchCancer <- matrix(i, length(colnames(DataSetsCancer[[i]])),
                    1)  # growing a batch data object
                  currentBatchDataCancer <- data.frame(SampleName = colnames(DataSetsCancer[[i]]),
                    Batch = currentBatchCancer)
                  MetaBatchDataCancer <- rbind(MetaBatchDataCancer, currentBatchDataCancer)

                  # Getting the normal data as well.
                  DataSetsNormal[[i]] <- Preprocess_MAdata_Normal(CancerSite, MAdirectories[i],
                    MAfiles[MatchedFile], MissingValueThresholdGene = MissingValueThresholdGene,
                    MissingValueThresholdSample = MissingValueThresholdSample, doBatchCorrection = doBatchCorrection,
                    batch.correction.method = batch.correction.method, BatchData = BatchData)
                  GeneListsNormal[[i]] <- rownames(DataSetsNormal[[i]])
                  SampleListsNormal[[i]] <- colnames(DataSetsNormal[[i]])
                  currentBatchNormal <- matrix(i, length(colnames(DataSetsNormal[[i]])),
                    1)  # growing a batch data object
                  currentBatchDataNormal <- data.frame(SampleName = colnames(DataSetsNormal[[i]]),
                    Batch = currentBatchNormal)
                  MetaBatchDataNormal <- rbind(MetaBatchDataNormal, currentBatchDataNormal)

                } else {
                  cat("MA file not found for this cancer.\n")
                }
            }
            # combine data sets with Combat.
            cat("Combining data sets.\n")
            OverlapProbesCancer <- Reduce(intersect, GeneListsCancer)
            OverlapProbesNormal <- Reduce(intersect, GeneListsNormal)
            OverlapSamplesCancer <- Reduce(intersect, SampleListsCancer)
            OverlapSamplesNormal <- Reduce(intersect, SampleListsNormal)
            if (length(OverlapSamplesCancer) > 0 | length(OverlapSamplesNormal) >
                0) {
                cat("This should not happen. There is overlap between cancer or normal samples. No solution yet.\n")
            }

            for (i in 1:length(MAdirectories)) {
                DataSetsCancer[[i]] <- DataSetsCancer[[i]][OverlapProbesCancer, ]
                DataSetsNormal[[i]] <- DataSetsNormal[[i]][OverlapProbesNormal, ]
            }
            # combine cancer data sets.
            MA_TCGA_Cancer <- Reduce(cbind, DataSetsCancer)


            MA_TCGA_Cancer <- CorrectBatchEffect(GEN_Data = MA_TCGA_Cancer, BatchData = MetaBatchDataCancer,
                batch.correction.method = batch.correction.method, MinInBatch = MinPerBatchCancer)

            # combine normal data sets.
            MA_TCGA_Normal <- Reduce(cbind, DataSetsNormal)
            MA_TCGA_Normal <- CorrectBatchEffect(GEN_Data = MA_TCGA_Normal, BatchData = MetaBatchDataNormal,
                batch.correction.method = batch.correction.method, MinInBatch = MinPerBatchNormal)

        } else {
            MAfiles <- dir(MAdirectories)
            MatchedFile <- grep(MAstring, MAfiles)
            if (length(MatchedFile) > 0) {
                MA_TCGA_Cancer <- Preprocess_MAdata_Cancer(CancerSite, MAdirectories,
                  MAfiles[MatchedFile], MissingValueThresholdGene, MissingValueThresholdSample,
                  doBatchCorrection = doBatchCorrection, batch.correction.method = batch.correction.method,
                  BatchData = BatchData)

                MA_TCGA_Normal <- Preprocess_MAdata_Normal(CancerSite, MAdirectories,
                  MAfiles[MatchedFile], MissingValueThresholdGene, MissingValueThresholdSample,
                  doBatchCorrection = doBatchCorrection, batch.correction.method = batch.correction.method,
                  BatchData = BatchData)

                cat("There are", length(colnames(MA_TCGA_Cancer)), "cancer samples and",
                  length(colnames(MA_TCGA_Normal)), "normal samples in gene expression data for",
                  CancerSite, "\n")
            } else {
                stop("MA file not found for this cancer.\n")
            }
        }
        MA_TCGA_Cancer <- TCGA_GENERIC_CleanUpSampleNames(MA_TCGA_Cancer, 15)
        if (ncol(MA_TCGA_Normal) == 0) {
            MA_TCGA_Normal <- NULL
        } else {
            MA_TCGA_Normal <- TCGA_GENERIC_CleanUpSampleNames(MA_TCGA_Normal, 15)
        }

        if (grepl("Windows", Sys.info()["sysname"]))
            system("subst x: /D")  #stop virtual drive

        # ---------------------------------------------------------------------------------------------
        # Step 3: Merge normal and cancer datasets into one matrix
        # ---------------------------------------------------------------------------------------------
        # EpiMix only takes one gene expression matrix as input, so we combine
        # cancer and normal into one matrix
        MA_TCGA_Data <- NULL
        if (!is.null(MA_TCGA_Normal)) {
            overlapGenes <- intersect(rownames(MA_TCGA_Cancer), rownames(MA_TCGA_Normal))
            MA_TCGA_Cancer <- MA_TCGA_Cancer[overlapGenes, , drop = FALSE]
            MA_TCGA_Normal <- MA_TCGA_Normal[overlapGenes, , drop = FALSE]
            MA_TCGA_Data <- cbind(MA_TCGA_Cancer, MA_TCGA_Normal)
        } else {
            MA_TCGA_Data <- MA_TCGA_Cancer
        }
    } else {
        # if target.region = 'lncRNA', just perform missing value estimation
        # MA_TCGA_Cancer = unzip(paste0(MAdirectories, '/', CancerSite,
        # '_lncRNA.txt.zip'), exdir = MAdirectories)
        MA_TCGA_Cancer <- read.table(MAdirectories, check.names = FALSE)
        MA_TCGA_Cancer <- TCGA_EstimateMissingValues_MolecularData(MA_TCGA_Cancer,
            MissingValueThresholdGene = MissingValueThresholdGene, MissingValueThresholdSample = MissingValueThresholdSample)
        # The TCGA lncRNA data were preprocessed by the Kallisto-Sleuth
        # pipeline. Sleuth has already corrected technical variations in
        # different experiments, so no batch correction is needed.
        MA_TCGA_Data <- MA_TCGA_Cancer
    }
    return(MA_TCGA_Data)
}

#' The Preprocess_MAdata_Cancer function
#'
#' Internal. Pre-process gene expression data for cancer samples.
#' @param CancerSite TCGA code for the cancer site.
#' @param Directory Directory.
#' @param File File.
#' @param MissingValueThresholdGene threshold for missing values per gene. Genes with a percentage of NAs greater than this threshold are removed. Default is 0.3.
#' @param MissingValueThresholdSample threshold for missing values per sample. Samples with a percentage of NAs greater than this threshold are removed. Default is 0.1.
#' @return The data matrix.
#' @keywords internal
#'
Preprocess_MAdata_Cancer <- function(CancerSite, Directory, File, MissingValueThresholdGene = 0.3,
    MissingValueThresholdSample = 0.1, doBatchCorrection, batch.correction.method,
    BatchData) {

    MinPerBatch <- 5
    cat("Loading cancer mRNA data.\n")
    cat("\tMissing value estimation.\n")
    MA_TCGA <- TCGA_Load_MolecularData(paste(Directory, File, sep = ""))
    MA_TCGA <- TCGA_EstimateMissingValues_MolecularData(MA_TCGA, MissingValueThresholdGene,
        MissingValueThresholdSample)
    Samplegroups <- TCGA_GENERIC_GetSampleGroups(colnames(MA_TCGA))
    if (CancerSite == "LAML") {
        MA_TCGA <- MA_TCGA[, Samplegroups$PeripheralBloodCancer, drop = FALSE]
    } else {
        MA_TCGA <- MA_TCGA[, Samplegroups$Primary, drop = FALSE]
    }

    if (doBatchCorrection) {
        cat("\tBatch correction.\n")
        # Remove samples with batch number 0
        if (length(-which(BatchData[, 2] == 0)) > 0) {
            BatchData <- BatchData[-which(BatchData[, 2] == 0), ]
        }
        MA_TCGA <- CorrectBatchEffect(GEN_Data = MA_TCGA, BatchData = BatchData,
            batch.correction.method = batch.correction.method, MinInBatch = MinPerBatch)
    }

    cat("\tProcessing gene ids and merging.\n")
    Genes <- rownames(MA_TCGA)
    SplitGenes <- limma::strsplit2(Genes, "\\|")
    rownames(MA_TCGA) <- SplitGenes[, 1]
    MA_TCGA <- MA_TCGA[!rownames(MA_TCGA) %in% "?", , drop = FALSE]
    MA_TCGA <- TCGA_GENERIC_MergeData(unique(rownames(MA_TCGA)), MA_TCGA)

    return(MA_TCGA = MA_TCGA)
}

#' The Preprocess_MAdata_Normal function
#'
#' Internal. Pre-process gene expression data for normal samples.
#' @param CancerSite TCGA code for the cancer site.
#' @param Directory Directory.
#' @param File File.
#' @param MissingValueThresholdGene threshold for missing values per gene. Genes with a percentage of NAs greater than this threshold are removed. Default is 0.3.
#' @param MissingValueThresholdSample threshold for missing values per sample. Samples with a percentage of NAs greater than this threshold are removed. Default is 0.1.
#' @return The data matrix.
#' @keywords internal
#'
Preprocess_MAdata_Normal <- function(CancerSite, Directory, File, MissingValueThresholdGene,
    MissingValueThresholdSample, doBatchCorrection, batch.correction.method, BatchData) {
    MinPerBatch <- 2  # less samples in one batch when dealing with normals.

    cat("Loading normal mRNA data.\n")
    cat("\tMissing value estimation.\n")
    MA_TCGA <- TCGA_Load_MolecularData(paste(Directory, File, sep = ""))
    MA_TCGA <- TCGA_EstimateMissingValues_MolecularData(MA_TCGA, MissingValueThresholdGene,
        MissingValueThresholdSample)
    Samplegroups <- TCGA_GENERIC_GetSampleGroups(colnames(MA_TCGA))
    if (CancerSite == "LAML") {
        MA_TCGA <- MA_TCGA[, Samplegroups$BloodNormal, drop = FALSE]
    } else {
        # MA_TCGA=MA_TCGA[,Samplegroups$SolidNormal]
        MA_TCGA <- MA_TCGA[, Samplegroups$SolidNormal, drop = FALSE]
    }

    if (doBatchCorrection) {
        cat("\tBatch correction.\n")
        # Remove samples with batch number 0
        if (length(-which(BatchData[, 2] == 0)) > 0) {
            BatchData <- BatchData[-which(BatchData[, 2] == 0), ]
        }
        MA_TCGA <- CorrectBatchEffect(GEN_Data = MA_TCGA, BatchData = BatchData,
            batch.correction.method = batch.correction.method, MinInBatch = MinPerBatch)
    }

    cat("\tProcessing gene ids and merging.\n")
    Genes <- rownames(MA_TCGA)
    SplitGenes <- limma::strsplit2(Genes, "\\|")
    rownames(MA_TCGA) <- SplitGenes[, 1]
    # MA_TCGA=MA_TCGA[!rownames(MA_TCGA) %in% '?',]
    MA_TCGA <- MA_TCGA[!rownames(MA_TCGA) %in% "?", , drop = FALSE]
    MA_TCGA <- TCGA_GENERIC_MergeData(unique(rownames(MA_TCGA)), MA_TCGA)

    return(MA_TCGA = MA_TCGA)
}

#' The TCGA_Load_MolecularData function
#'
#' Internal. Reads in gene expression data. Deletes samples and genes with more NAs than the respective thresholds. Imputes other NAs values.
#' @param Filename name of the file with the data.
#' @param MissingValueThresholdGene threshold for missing values per gene. Genes with a percentage of NAs greater than this threshold are removed. Default is 0.3.
#' @param MissingValueThresholdSample threshold for missing values per sample. Samples with a percentage of NAs greater than this threshold are removed. Default is 0.1.
#' @return gene expression data.
#' @keywords internal
#'
TCGA_Load_MolecularData <- function(Filename) {

    # loading the data in R
    Matching <- regexpr(".mat$", Filename, perl = TRUE)
    if (Matching[[1]] > 0) {
        # MET_Data=TCGA_GENERIC_ReadDataMatrixMatFile(Filename)
        DataList <- R.matlab::readMat(Filename)
        MATdata <- as.matrix(DataList$RawData)
        rownames(MATdata) <- DataList[[2]]
        colnames(MATdata) <- DataList[[3]]

    } else {
        MET_Data <- read.csv(Filename, sep = "\t", row.names = 1, header = TRUE,
            na.strings = c("NA", "null"))
    }
    if (rownames(MET_Data)[1] == "Composite Element REF") {
        cat("Removing first row with text stuff.\n")
        MET_Data <- MET_Data[-1, ]
        Genes <- rownames(MET_Data)
        MET_Data <- apply(MET_Data, 2, as.numeric)
        rownames(MET_Data) <- Genes
    }

    SampleNames <- colnames(MET_Data)
    SampleNames <- gsub("\\.", "-", SampleNames)
    colnames(MET_Data) <- SampleNames
    return(MET_Data)
}

#' The TCGA_EstimateMissingValues_MolecularData function
#'
#' Internal.Deletes samples and genes with more NAs than the respective thresholds. Imputes other NAs values.
#' @param MET_Data matrix of gene expression data
#' @param MissingValueThresholdGene threshold for missing values per gene. Genes with a percentage of NAs greater than this threshold are removed. Default is 0.3.
#' @param MissingValueThresholdSample threshold for missing values per sample. Samples with a percentage of NAs greater than this threshold are removed. Default is 0.1.
#' @return gene expression data with no missing values.
#' @keywords internal
#'
TCGA_EstimateMissingValues_MolecularData <- function(MET_Data, MissingValueThresholdGene = 0.3,
    MissingValueThresholdSample = 0.1) {
    # removing clones with too many missing values
    NrMissingsPerGene <- apply(MET_Data, 1, function(x) sum(is.na(x)))/ncol(MET_Data)
    cat("Removing", sum(NrMissingsPerGene > MissingValueThresholdGene), "genes with more than",
        MissingValueThresholdGene * 100, "% missing values.\n")
    if (sum(NrMissingsPerGene > MissingValueThresholdGene) > 0)
        MET_Data <- MET_Data[NrMissingsPerGene < MissingValueThresholdGene, ]

    # removing patients with too many missings values
    NrMissingsPerSample <- apply(MET_Data, 2, function(x) sum(is.na(x)))/nrow(MET_Data)
    cat("Removing", sum(NrMissingsPerSample > MissingValueThresholdSample), "patients with more than",
        MissingValueThresholdSample * 100, "% missing values.\n")
    if (sum(NrMissingsPerSample > MissingValueThresholdSample) > 0)
        MET_Data <- MET_Data[, NrMissingsPerSample < MissingValueThresholdSample]

    # knn impute using Tibshirani's method
    if (length(colnames(MET_Data)) > 1) {
        k <- 15
        KNNresults <- impute::impute.knn(as.matrix(MET_Data), k)
        MET_Data_KNN <- KNNresults$data

        # cleaning up sample names
        MET_Data_KNN_Clean <- TCGA_GENERIC_CleanUpSampleNames(MET_Data_KNN, 15)
        return(MET_Data_KNN_Clean)

    } else {
        # when only 1 sample,need to make a matrix again
        # MET_Data=as.matrix(MET_Data)
        MET_Data_Clean <- TCGA_GENERIC_CleanUpSampleNames(MET_Data, 15)
        return(MET_Data_Clean)
    }
}

#' The TCGA_GENERIC_MergeData function
#'
#' Internal.
#' @param NewIDListUnique unique rownames of data.
#' @param DataMatrix data matrix.
#' @return data matrix.
#' @keywords internal
#'
TCGA_GENERIC_MergeData <- function(NewIDListUnique, DataMatrix) {

    NrUniqueGenes <- length(NewIDListUnique)
    MergedData <- matrix(0, NrUniqueGenes, length(colnames(DataMatrix)))
    for (i in 1:NrUniqueGenes) {
        currentID <- NewIDListUnique[i]
        # tmpData=DataMatrix[which(rownames(DataMatrix) %in% currentID),]
        tmpData <- DataMatrix[which(rownames(DataMatrix) %in% currentID), , drop = FALSE]
        if (length(rownames(tmpData)) > 1) {
            MergedData[i, ] <- colMeans(tmpData)
        } else {
            MergedData[i, ] <- tmpData
        }
    }
    rownames(MergedData) <- NewIDListUnique
    colnames(MergedData) <- colnames(DataMatrix)

    return(MergedData)
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
    ProbeAnnotation, MET_Cancer, MET_Normal = NULL, CorThreshold = 0.4) {

    # first lookup the probes matching a single gene. DO NOT USE grep, it does
    # not do exact matching, but looks for the pattern anywhere !!!
    Probes <- ProbeAnnotation[which(ProbeAnnotation[, 2] == Gene), 1]
    Probes <- Probes[which(Probes %in% rownames(MET_Cancer))]

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

#' The TCGA_GetSampleInfo function
#' @details Generate the 'sample.info' dataframe for TCGA data.
#' @param METProcessedData Matrix of preprocessed methylation data.
#' @param CancerSite Character string of TCGA study abbreviation.
#' @param TargetDirectory Path to save the sample.info. Default: ''.
#'
#' @return A dataframe for the sample groups. Contains two columns: the first column (named: 'primary') indicating the sample names, and the second column (named: 'sample.type') indicating whether each sample is a Cancer or Normal tissue.
#' @export
#' @examples
#' {
#' data(MET.data)
#' sample.info <- TCGA_GetSampleInfo(MET.data, CancerSite = 'LUAD')
#' }

TCGA_GetSampleInfo <- function(METProcessedData, CancerSite = "LUAD", TargetDirectory = "") {
    # Split up normal and cancer data
    Samplegroups <- cancer_samples <- normal_samples <- NULL
    Samplegroups <- TCGA_GENERIC_GetSampleGroups(colnames(METProcessedData))
    if (CancerSite == "LAML") {
        cancer_samples <- Samplegroups$PeripheralBloodCancer
    } else {
        cancer_samples <- Samplegroups$Primary
    }
    if (CancerSite == "LAML") {
        normal_samples <- Samplegroups$BloodNormal
    } else {
        normal_samples <- Samplegroups$SolidNormal
    }
    df.cancer <- data.frame(primary = cancer_samples, sample.type = rep("Cancer",
        length(cancer_samples)))
    df.normal <- data.frame(primary = normal_samples, sample.type = rep("Normal",
        length(normal_samples)))
    sample.info <- rbind(df.cancer, df.normal)
    if (TargetDirectory != "") {
        utils::write.csv(sample.info, paste0(TargetDirectory, "/", "sample.info.csv"),
            row.names = FALSE)
    }
    cat("There are a total of", length(cancer_samples), "cancer samples and", length(normal_samples),
        "normal samples with sample information for", CancerSite, "\n")
    return(sample.info)
}


#' The Preprocess_CancerSite_Methylation27k function
#'
#' Internal. Pre-processes DNA methylation data from TCGA from Illymina 27k arrays.
#' @param CancerSite character of length 1 with TCGA cancer code.
#' @param METdirectory character with directory where a folder for downloaded files will be created. Can be the object returned by the Download_DNAmethylation function.
#' @param MissingValueThreshold threshold for removing samples or genes with missing values.
#' @return List with pre processed methylation data for cancer and normal samples.
#' @keywords internal
#'
#'
Preprocess_CancerSite_Methylation27k <- function(CancerSite, METdirectory, doBatchCorrection,
    batch.correction.method, MissingValueThreshold) {

    if (grepl("Windows", Sys.info()["sysname"])) {
        # If Windows I'll create a virtual drive to handle the long file names
        # issue Create a virtual drive to overcome long names issue in Windows
        virtualDir <- METdirectory
        virtualDir <- gsub("\\\\", "/", virtualDir)
        virtualDir <- substr(virtualDir, 1, nchar(virtualDir) - 1)
        system(paste("subst x:", virtualDir))

        # Load data
        METfiles <- dir("x:")
        MatchedFilePosition <- grep("methylation__humanmethylation27", METfiles)
        Filename <- paste0("x://", METfiles[MatchedFilePosition])
        MET_Data <- TCGA_GENERIC_LoadIlluminaMethylationData(Filename)

        system("subst x: /D")  #stop virtual drive
    } else {
        # Not windows Load data
        METfiles <- dir(METdirectory)
        MatchedFilePosition <- grep("methylation__humanmethylation27", METfiles)
        Filename <- paste0(METdirectory, METfiles[MatchedFilePosition])
        MET_Data <- TCGA_GENERIC_LoadIlluminaMethylationData(Filename)
    }

    # Split up normal and cancer data
    Samplegroups <- TCGA_GENERIC_GetSampleGroups(colnames(MET_Data))
    if (CancerSite == "LAML") {
        MET_Data_Cancer <- MET_Data[, Samplegroups$PeripheralBloodCancer, drop = FALSE]
    } else {
        MET_Data_Cancer <- MET_Data[, Samplegroups$Primary, drop = FALSE]
    }
    MET_Data_Cancer <- TCGA_GENERIC_CleanUpSampleNames(MET_Data_Cancer, 15)
    if (CancerSite == "LAML") {
        MET_Data_Normal <- MET_Data[, Samplegroups$BloodNormal, drop = FALSE]
    } else {
        MET_Data_Normal <- MET_Data[, Samplegroups$SolidNormal, drop = FALSE]
    }
    if (length(MET_Data_Normal) > 0) {
        MET_Data_Normal <- TCGA_GENERIC_CleanUpSampleNames(MET_Data_Normal, 15)
    }

    # Clear space
    rm(MET_Data)
    gc()

    # Missing value estimation
    cat("\tMissing value estimation for the cancer samples.\n")
    MET_Data_Cancer <- TCGA_Process_EstimateMissingValues(MET_Data_Cancer, MissingValueThreshold)
    if (length(MET_Data_Normal) > 0) {
        cat("\tMissing value estimation for the normal samples.\n")
        MET_Data_Normal <- TCGA_Process_EstimateMissingValues(MET_Data_Normal, MissingValueThreshold)
    }

    if (doBatchCorrection) {
        # Batch correction for cancer and normal.
        BatchData <- EpiMix_GetData("TCGA_BatchData")
        MinPerBatch <- 5
        cat("\tBatch correction for the cancer samples.\n")
        MET_Data_Cancer <- CorrectBatchEffect(GEN_Data = MET_Data_Cancer, BatchData = BatchData,
            batch.correction.method = batch.correction.method, MinInBatch = MinPerBatch,
            featurePerSet = 50000)


        if (length(MET_Data_Normal) > 0) {
            cat("\tBatch correction for the normal samples.\n")
            MET_Data_Normal <- CorrectBatchEffect(GEN_Data = MET_Data_Normal, BatchData = BatchData,
                batch.correction.method = batch.correction.method, MinInBatch = MinPerBatch,
                featurePerSet = 50000)
        } else {
            MET_Data_Normal <- c()
        }
        # Set values <0 to 0 and >1 to 1, because of batch correction
        MET_Data_Cancer[MET_Data_Cancer < 0] <- 0
        MET_Data_Cancer[MET_Data_Cancer > 1] <- 1
        if (length(MET_Data_Normal) > 0) {
            MET_Data_Normal[MET_Data_Normal < 0] <- 0
            MET_Data_Normal[MET_Data_Normal > 1] <- 1
        }
    }

    return(list(MET_Data_Cancer = MET_Data_Cancer, MET_Data_Normal = MET_Data_Normal))
}


#' The getLncRNAData function
#' @description Helper function to retrieve the lncRNA expression data from Experiment Hub
#' @param CancerSite TCGA cancer code
#'
#' @return local file path where the lncRNA expression data are saved
#' @keywords internal

getLncRNAData <- function(CancerSite) {
    eh <- ExperimentHub::ExperimentHub()
    data <- AnnotationHub::query(eh, "EpiMix.data")
    hub_id <- data$ah_id[which(data$title == paste0(CancerSite, "_lncRNA.txt"))]
    path <- eh[[hub_id]]
    return(path)
}
