test_TCGA_getData <- function(){
  EpiMixResults <- TCGA_GetData(CancerSite = "OV",
                                mode = "Regular",
                                roadmap.epigenome.ids =  "E097", # only required for running the enhancer mode
                                outputDirectory = tempdir(),
                                cores = 10)

  RUnit::checkEquals(nrow(EpiMixResults$FunctionalPairs), 789)
}








