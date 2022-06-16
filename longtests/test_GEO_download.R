test_GEO_download <- function(){
  METdirectories <- GEO_Download_DNAMethylation(AccessionID = 'GSE114134', targetDirectory = tempdir())
  RUnit::checkTrue(length(METdirectories)>0)
}








