### Add annotation
library("BSgenome.Brapa.3.0")
library("IsoformSwitchAnalyzeR")
library("sleuth")
library("rhdf5")

load("aSwitchListFilteredA1.obj")

aSwitchListFilteredA1 <- analyzeCPC2(aSwitchListFilteredA1, "/scratch/punim0687/agolicz/rapa.isoforms.08032019/7.IsoformSwitchAnalyzeR/analysis/prep/cpc2.results", removeNoncodinORFs = TRUE)
aSwitchListFilteredA1 <- analyzePFAM(aSwitchListFilteredA1, "/scratch/punim0687/agolicz/rapa.isoforms.08032019/7.IsoformSwitchAnalyzeR/analysis/prep/pfam.results" )
aSwitchListFilteredA1 <- analyzeSignalP(aSwitchListFilteredA1, "/scratch/punim0687/agolicz/rapa.isoforms.08032019/7.IsoformSwitchAnalyzeR/analysis/prep/signalp.results")

aSwitchListFilteredA1<-analyzeAlternativeSplicing(aSwitchListFilteredA1, onlySwitchingGenes=TRUE, alpha=0.05, dIFcutoff = 0.2, showProgress=TRUE, quiet=FALSE)

aSwitchListFilteredA2<-analyzeSwitchConsequences(aSwitchListFilteredA1,consequencesToAnalyze=c('intron_retention','coding_potential','ORF_seq_similarity','domain_length','domains_identified','signal_peptide_identified', 'tss', 'tts', 'intron_structure'), alpha=0.05, dIFcutoff = 0.2, onlySigIsoforms=TRUE, removeNonConseqSwitches=FALSE,showProgress=TRUE,quiet=FALSE)

aSwitchListFilteredA3<-analyzeSwitchConsequences(aSwitchListFilteredA1,consequencesToAnalyze=c('intron_retention','coding_potential','ORF_seq_similarity','domain_length','domains_identified','signal_peptide_identified', 'tss', 'tts', 'intron_structure'), alpha=0.05, dIFcutoff = 0.2, onlySigIsoforms=FALSE, removeNonConseqSwitches=FALSE, showProgress=TRUE, quiet=FALSE)

write.table(aSwitchListFilteredA2$switchConsequence,file="switch.consequences.tsv", sep="\t", row.names=F, quote=F)

write.table(aSwitchListFilteredA3$switchConsequence,file="switch.consequences.withnoSig.tsv", sep="\t", row.names=F, quote=F)

save(aSwitchListFilteredA2, file="aSwitchListFilteredA2.obj")
save(aSwitchListFilteredA3, file="aSwitchListFilteredA3.obj")

