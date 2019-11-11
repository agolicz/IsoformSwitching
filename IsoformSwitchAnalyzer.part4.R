library("BSgenome.Brapa.3.0")
library("IsoformSwitchAnalyzeR")
load("aSwitchListFilteredA2.obj")

write.table(subset(aSwitchListFilteredA2$isoformFeatures, isoform_switch_q_value < 0.05 & abs(dIF) > 0.2), file="isoform.features.significant.tsv", sep="\t", quote=F, row.names=F)
write.table(aSwitchListFilteredA2$AlternativeSplicingAnalysis, file="alternative.splicing.tsv", sep="\t", quote=F, row.names=F)
write.table(aSwitchListFilteredA2$switchConsequence, file="switch.consequnce.tsv", sep="\t", quote=F, row.names=F)
write.table(aSwitchListFilteredA2$isoformSwitchAnalysis, file="switch.analsis.tsv", sep="\t", quote=F, row.names=F)

pdf("splicing.summary.pdf")
ss<-extractSplicingSummary(aSwitchListFilteredA2, splicingToAnalyze = c('A3', 'A5', 'ES', 'IR', 'MES'), alpha=0.05, dIFcutoff = 0.2, onlySigIsoforms = TRUE, plot = TRUE, returnResult = TRUE, asFractionTotal = TRUE, localTheme = theme_minimal())
dev.off()
write.table(ss, file="splicing.summary.tsv", sep="\t", quote=F, row.names=F)

ass<-aggregate(nrIsoWithConsequences ~ AStype, ss, sum)
write.table(ass, file="splicing.summary.sum.tsv", sep="\t", quote=F, row.names=F)

pdf("splicing.enrichment.pdf")
se<-extractSplicingEnrichment(aSwitchListFilteredA2, splicingToAnalyze=c('A3', 'A5', 'ES', 'IR', 'MES'), alpha=0.05, dIFcutoff = 0.2, onlySigIsoforms = TRUE, returnResult=TRUE, returnSummary=FALSE, localTheme = theme_minimal())
dev.off()
write.table(se, file="splicing.enrichment.tsv", sep="\t", quote=F, row.names=F)

pdf("consequnce.enrichment.pdf")
cr<-extractConsequenceEnrichment(aSwitchListFilteredA2, consequencesToAnalyze = 'all' , alpha=0.05, dIFcutoff = 0.2, returnResult=TRUE,returnSummary=FALSE, localTheme = theme_minimal())
dev.off()
write.table(cr, file="consequnce.enrichment.tsv", sep="\t", quote=F, row.names=F)
