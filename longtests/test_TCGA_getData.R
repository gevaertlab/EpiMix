test_TCGA_getData <- function(){

  data(MET.data)
  data(mRNA.data)
  data(LUAD.sample.annotation)

  EpiMixResults_Regular <- EpiMix(methylation.data = MET.data,
                                  gene.expression.data = mRNA.data,
                                  sample.info = LUAD.sample.annotation,
                                  group.1 = "Cancer",
                                  group.2 = "Normal",
                                  met.platform = "HM450",
                                  OutputRoot = tempdir())

  RUnit::checkEquals(nrow(EpiMixResults$FunctionalPairs), 789)
}








