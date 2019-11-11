library("BSgenome.Brapa.3.0")
library("IsoformSwitchAnalyzeR")
library("sleuth")
library("rhdf5")
options(stringsAsFactors = FALSE)

eq <- importIsoformExpression("/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results")

myDesign <- data.frame(sampleID = c("B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11","B12","B13","B14","B15"),condition = c("PMC","PMC","PMC","Tet","Tet","Tet","Mic","Mic","Mic","Bin","Bin","Bin","Pol","Pol","Pol"))

comps <- data.frame(condition_1 = c("PMC", "Tet", "Mic", "Bin"), condition_2 = c("Tet", "Mic", "Bin", "Pol"))

aSwitchList <- importRdata(isoformCountMatrix = eq$counts, isoformRepExpression = eq$abundance, designMatrix = myDesign, comparisonsToMake = comps, isoformExonAnnoation = "/scratch/punim0687/agolicz/rapa.isoforms.08032019/7.IsoformSwitchAnalyzeR/analysis/Rapa.03042019.merged.renamed.nm.r.gtf", showProgress = TRUE)

sub.t <-scan("Rapa.03042019.merged.renamed.nm.r.ni.nh.expressed.gtf.transcripts", what="")
aSwitchListS<-subsetSwitchAnalyzeRlist(aSwitchList,aSwitchList$isoformFeatures$isoform_id %in% sub.t)

aSwitchListS <- analyzeORF(aSwitchListS, genomeObject = BSgenome.Brapa.3.0,  minORFlength=300, orfMethod = "longest")

write.table(aSwitchListS$orfAnalysis, file="orf.analysis.tsv", sep="\t", quote=F, row.names=T)

