test_miRNA_mode <- function(){

  data(MET.data)
  data(microRNA.data)
  data(LUAD.sample.annotation)

  EpiMixResults <- EpiMix(methylation.data = MET.data,
                         gene.expression.data = microRNA.data,
                         mode = "miRNA",
                         sample.info = LUAD.sample.annotation,
                         group.1 = "Cancer",
                         group.2 = "Normal",
                         met.platform = "HM450",
                         OutputRoot =  tempdir())

  RUnit::checkEquals(nrow(EpiMixResults$FunctionalPairs), 6)

}








