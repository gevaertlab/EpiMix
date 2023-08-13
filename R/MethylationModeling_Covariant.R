#' The SingleGene_CovariantModeling function
#'
#'Internal. For a given gene, this function fits the mixture model and defines the respective methylation states.
#' @param dependent string holding the name of the probe for that run
#' @param CancerSamp vector containing the Cancer M-values for the probe of study
#' @param Covariant vector containing Covariant names
#' @param NormalSamp vector of values containing the Normal data points for the CpG site
#' @param Cov vector containing Covariant values to be adjusted for
#' @param PvalueThreshold threshold to consider results significant.
#' @param MeanDifferenceTreshold threshold in beta value scale from which two methylation means are considered different.
#' @param minSamplePerGroup minimum number of samples required to belong to a new mixture component in order to accept it.
#' @param NoNormalMode logical, if TRUE no comparison to normal samples is performed. Defaults to FALSE.
#' @return NrComponents number of components identified.
#' @return Models an object with the parameters of the model fitted.
#' @return MethylationStates vector with DM values for each sample.
#' @return MixtureStates vector with DMvalues for each component.
#' @return Classifications a vector indicating to which component each sample was assigned.
#' @return FlipOverState FlipOverState
#' @details maxComp, PvalueThreshold, METDiffThreshold, minSamplesPerGroup are arguments for this function but are fixed in their default values for the user
#' because they are not available in the main MethylMix function, to keep it simple. It would be easy to make them available to the user if we want to.
#' @keywords internal
#'

SingleGene_CovariantModeling <- function(dependent, CancerSamp, Covariant, NormalSamp, Cov, 
                                         PvalueThreshold = 0.01, MeanDifferenceTreshold = 0.10, 
                                         minSamplesPerGroup = 1, NoNormalMode = FALSE, maxComp = 3){
  
  METdataVector <- CancerSamp
  
  # 1 component model
  w0.m = matrix(1, length(METdataVector), 1)
  mods = vector("list", maxComp + 1)
  mods[[1]] = blc_2(matrix(METdataVector, ncol = 1), w = w0.m, maxiter = 100, tol = 1e-06, verbose = FALSE)
  bic = numeric(maxComp + 1)
  bic[1] = -2 * mods[[1]]$llike + 2 * log(length(METdataVector))
  
  # 2 to maxComp components model
  for (comp in 2:maxComp) {
    
    # Divide initial groups using quantiles
    prob = seq(1/comp, (comp-1)/comp, 1/comp)
    tmpQuantiles = quantile(METdataVector, prob)
    w0.m = matrix(0, nrow = length(METdataVector), ncol = comp)
    tmpQuantiles = c(tmpQuantiles, Inf)
    w0.m[METdataVector < tmpQuantiles[1], 1] <- 1
    for (i in 2:comp) {
      w0.m[METdataVector >= tmpQuantiles[i - 1] & METdataVector < tmpQuantiles[i], i] <- 1
    }
    
    # Fit beta mixture model
    mods[[comp]] <- blc_2(matrix(METdataVector, ncol=1), w = w0.m, maxiter = 100, tol = 1e-06, verbose = FALSE)
    if (sum(is.na(mods[[comp]]$mu)) > 0) mods[[comp]]$llike=0
    df = comp * 3 - 1
    bic[[comp]] = -2 * mods[[comp]]$llike + df * log(length(METdataVector))
    
    # See differences between model's means to compare them to threshold
    model.means = sort(mods[[comp]]$mu)
    different.means = ifelse(all(abs(diff(model.means)) > MeanDifferenceTreshold), TRUE, FALSE)
    
    # Check if smallest group has at least minSamplesPerGroup observations:
    if (minSamplesPerGroup < 0) {
      minSamplesPerGroup = max(5, 0.05 * length(METdataVector))
    }
    minOK = ifelse(min(table(apply(mods[[comp]]$w, 1, which.max))) >= minSamplesPerGroup, TRUE, FALSE)
    
    # We try adding another component if the following 2 conditions are satisfied:
    #   A: Adding one component reduces BIC
    #   B: All absolute differences between methylation means in model with one extra component are above the MeanDifferenceThreshold
    # But I also check C = smallest group has at least minSamplesPerGroup observations:
    # So, continue with another component if A & B & C, else not continue, which is the same as saying not continue if !A OR !B OR !C
    # If the model was improved try one more, else break
    if (bic[[comp]] >= bic[[comp - 1]] || !different.means || !minOK) {
      NrComponents = comp - 1
      break
    } else {
      # It improved, try one more (unless comp already is maxComp, in that case we end here)
      NrComponents = comp
    }
  }  
  rm(bic)
  
  ## Building "data" matrix, as bringing whole METcancer_df matrix in is incredibly RAM inefficient. If dealing with
  ## smaller datasets, can remove the following and parameterize METcancer_df
  
  CancerSamp <- t(CancerSamp)
  METdata <- rbind(CancerSamp, Cov)
  METdata <- t(METdata)
  colnames(METdata)[1] <- dependent
  colnames(METdata)[2] <- Covariant
  METdata <- as.data.frame(METdata)
  
  ## Adjusting data with EM Algorithm to regress out the Covariant
  
  CAMAN_glm_obj <- Cov_Adjustment(dep = dependent, fixed =  Covariant, data = METdata,
                                  k = NrComponents, pop.at.risk = NULL, family = "gaussian")
  
  rm(METdata)
  
  #Initializing data for Modeling from the CAMAN_obj and the list of Normal Sample names
  
  if (NrComponents == 1){
    fitted_obs <- CAMAN_glm_obj[[1]][["fitted.values"]]
  }
  else{
    fitted_obs <- CAMAN_glm_obj@fittedObs
  }
  
  
  METdataVector <- fitted_obs
  METdataNormalVector <- NormalSamp
  
  ##removing RAM heavy objects
  rm(CAMAN_glm_obj)
  gc()
  
  
  Model = mods[[NrComponents]]
  MethylationState = matrix(0, 1, length(METdataVector))
  FlipOverState = 0
  res = list(p.value = 1)
  MixtureStates = matrix(0, NrComponents, 1)
  classification = apply(Model$w, 1, which.max)
  
  if (NrComponents == 1) {
    if (!is.null(METdataNormalVector)) {
      res = wilcox.test(METdataVector, METdataNormalVector)
      Difference = mean(METdataVector) - mean(METdataNormalVector)
    } else {
      res = list(p.value = 1)
      Difference = mean(METdataVector)
    }
    if ((res$p.value < PvalueThreshold & abs(Difference) > MeanDifferenceTreshold) | NoNormalMode) {
      MethylationState[1, ] = Difference
      MixtureStates[1, 1] = Difference
    }
    #cat(c(GeneName,": 1 component is best.\n"))
  } else {
    for (comp in seq_len(NrComponents)) {
      METdataVector_comp = METdataVector[classification == comp]
      if (!is.null(METdataNormalVector)) {
        if (length(METdataVector_comp) > 0) res = wilcox.test(METdataVector_comp, METdataNormalVector) else res$p.value = 1
        Difference = mean(METdataVector_comp) - mean(METdataNormalVector)
      } else {
        res = list(p.value = 1)
        Difference = mean(METdataVector_comp)
      }
      if ((res$p.value < PvalueThreshold & abs(Difference) > MeanDifferenceTreshold) | NoNormalMode) {
        MethylationState[1, classification == comp] = Difference
        MixtureStates[comp, 1] = Difference
      }
    }
    # Flipover correction (in first MethylMix package was done only for 2 mixture states)
    if (NrComponents == 2 || NrComponents == 3) {
      OrigOrder = order(METdataVector)
      FlipOverResults = MethylMix_RemoveFlipOver(OrigOrder, MethylationState, classification, METdataVector, NrComponents)
      MethylationState = FlipOverResults$MethylationState
      classification = FlipOverResults$classification
      FlipOverState = FlipOverResults$LearnedState
    }
    #cat(c(GeneName,": ", NrComponents, " components are best.\n"))
  }
  print("Done")
  
  return(list(MethylationStates = MethylationState,
              NrComponents = NrComponents,
              Models = list(Model),
              MixtureStates = list(MixtureStates),
              Classifications = classification,
              FlipOverStates = FlipOverState))
}



#' The MethylMix_MixtureModel_Covariant function
#'
#' Internal. Prepares all the structures to store the results and calls in a foreach loop a function that fits the mixture model in each gene
#' for Covariant-adjusted data.
#' @param GeneNames vector containing the names of all Genes of study for each patient
#' @param METcancer matrix with methylation data for cancer samples (genes in rows, samples in columns).
#' @param METnormal matrix with methylation data for normal samples (genes in rows, samples in columns). If NULL no comparison to normal samples will be done.
#' @param FunctionalProbes vector with genes names to be considered for the mixture models.
#' @param Covariant vector containing the name(s) of the Covariants in the data
#' @param NoNormalMode logical, if TRUE no comparison to normal samples is performed. Defaults to FALSE.
#' @return MethylationStates matrix of DM values, with driver genes in the rows and samples in the columns.
#' @return NrComponents matrix with the number of components identified for each driver gene.
#' @return Plot_data an object with the data necessary for plotting a histogram
#' @return MethylationDrivers character vector with the genes found by MethylMix as differentially methylated and transcriptionally predictive (driver genes).
#' @return MixtureStates a list with a matrix for each driver gene containing the DM values.
#' @return Classifications a vector indicating to which component each sample was assigned.
#' @import doSNOW
#' @keywords internal
#'

MethylMix_MixtureModel_Covariant <- function(GeneNames, METcancer = NULL, METnormal = NULL,
                                             FunctionalProbes, NoNormalMode = FALSE, Covariant) {
  
  # Reformatting Cancer Sample and Normal Sample Matrices
  METcancer = as.matrix(METcancer)
  Cov <- METcancer[nrow(METcancer),]
  METcancer <- METcancer[-nrow(METcancer),]
  rownames(METcancer) <- GeneNames[1:length(GeneNames)-1]
  if (nrow(METnormal) < nrow(METcancer)){
    METnormal = as.matrix(METnormal[-nrow(METnormal),])
  }
  
  
  
  # Removing all Beta Values of 0, as converting to M-Value will yield to "-Inf".
  METcancer[METcancer==0.00000000] <- 0.00000001
  METnormal[METnormal==0.00000000] <- 0.00000001
  METcancer[METcancer == 1.0000000000] <- 0.9999999999
  METnormal[METnormal == 1.0000000000] <- 0.9999999999
  
  
  FunctionalGenes <- FunctionalProbes
  
  
  # overlap of samples
  if (!is.null(METnormal)) {
    GeneOverlap = intersect(intersect(rownames(METcancer), rownames(METnormal)), FunctionalGenes)
    METnormal = METnormal[GeneOverlap, , drop=FALSE]
    METcancer = METcancer[GeneOverlap, , drop = FALSE]
  } else {
    GeneOverlap = intersect(rownames(METcancer), FunctionalGenes)
    METcancer = METcancer[GeneOverlap, ,drop=FALSE]
  }
  
  
  cat("\nStarting M-Value mixture modeling.\n")
  MethylationStates = matrix(0, length(rownames(METcancer)), length(colnames(METcancer)))
  rownames(MethylationStates) = rownames(METcancer)
  colnames(MethylationStates) = colnames(METcancer)
  
  
  cat("Running M-value mixture model on",length(rownames(METcancer)),"probes and on",length(colnames(METcancer)),"samples.\n")
  
  # set a progress bar
  
  pb = NULL
  iterations = length(rownames(METcancer))
  pb <- utils :: txtProgressBar(max = iterations, style = 3)
  progress <- function(n) utils :: setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  i <- NULL # to avoid "no visible binding for global variable" in R CMD check
  
  # Reformatting Cancer Matrix and Normal Matrix for input in "Cov_Adjustment"  and 
  # "SingleGene_CovariantModeling" functions
  
  METcancer_df <- as.data.frame(t(beta2m(METcancer)))
  METnormal <- as.data.frame(t(beta2m(METnormal)))
  rm(METcancer)
  
  #
  
  if (!is.null(METnormal)) {
    res <- foreach::foreach(a = iter(METcancer_df, by='col'), i = 1:ncol(METcancer_df), .combine = "combineForEachOutput",
                            .export = c("SingleGene_CovariantModeling", "Cov_Adjustment", "mix.perform_glm", "computeDensities",
                                        "MethylMix_RemoveFlipOver", "blc_2","betaEst_2", "mix.densDistr_Epi"),
                            .packages = "CAMAN", .options.snow=opts) %dopar% 
      SingleGene_CovariantModeling(dependent = colnames(METcancer_df)[i], Covariant = Covariant, Cov = Cov,
                                   CancerSamp = a, NormalSamp = METnormal[,i])
    
  } else {
    res <- foreach::foreach(a = iter(METcancer_df, by='col'), i = 1:ncol(METcancer_df), .combine = "combineForEachOutput",
                            .export = c("SingleGene_CovariantModeling2", "Cov_Adjustment", "mix.perform_glm", "computeDensities",
                                        "MethylMix_RemoveFlipOver", "blc_2","betaEst_2", "mix.densDistr_Epi"),
                            .packages = "CAMAN", .options.snow=opts) %dopar% 
      SingleGene_CovariantModeling2(dependent = colnames(METcancer_df)[i], Covariant = Covariant, Cov = Cov,
                                    CancerSamp = a, NormalSamp = METnormal[,i])
    
  }
  
  close(pb)
  
  METcancer <- as.matrix(t(METcancer_df))
  rm(METcancer_df)
  
  if (is.null(dim(res$Classifications))) res$Classifications <- matrix(res$Classifications, nrow = 1)
  if (is.null(dim(res$MethylationStates))) res$Classifications <- matrix(res$Classifications, nrow = 1)
  
  rownames(res$MethylationStates) <- rownames(METcancer)
  if(nrow(res$Classifications) == nrow(METcancer)){
    rownames(res$Classifications) <- rownames(METcancer)
  }
  colnames(res$MethylationStates) <- colnames(res$Classifications) <- colnames(METcancer)
  
  # Removing the genes without any differential methylation.
  #cat("Removing the genes without any differential methylation...\n")
  if (!NoNormalMode) {
    # If NoNormalMode == T, no comparison to normal is made, and we don't remove genes with methylation states equal to 0
    # (for example, in pancancer analysis, running MethylMix only with normal samples, we use NoNormalMode = TRUE, and genes with only one state equal to 0 are kept.)
    NonZeroPositions = rowSums(res$MethylationStates) != 0
    res$NrComponents = res$NrComponents[NonZeroPositions, drop=FALSE]
    res$MixtureStates = res$MixtureStates[NonZeroPositions, drop=FALSE]
    res$Models = res$Models[NonZeroPositions, drop=FALSE]
    res$MethylationStates = res$MethylationStates[NonZeroPositions, , drop=FALSE]
    if(nrow(res$Classifications) == nrow(METcancer)){
      res$Classifications = res$Classifications[NonZeroPositions, , drop=FALSE]
    }
  }
  
  # Adding names and removing things that we don't want to return
  #cat("Adding names and removing things that we don't want to return...\n")
  res$MethylationDrivers = rownames(res$MethylationStates)
  res$FlipOverStates <- res$FunctionalGenes <- NULL
  names(res$NrComponents) = rownames(res$MethylationStates)
  names(res$MixtureStates) = rownames(res$MethylationStates)
  names(res$Models) = rownames(res$MethylationStates)
  
  
  return(res)
}