library("BSgenome.Brapa.3.0")
library("IsoformSwitchAnalyzeR")
library("sleuth")
library("rhdf5")
options(stringsAsFactors = FALSE)

eq <- importIsoformExpression("/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results", interLibNormTxPM=FALSE)

myDesign <- data.frame(sampleID = c("B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11","B12","B13","B14","B15"),condition = c("PMC","PMC","PMC","Tet","Tet","Tet","Mic","Mic","Mic","Bin","Bin","Bin","Pol","Pol","Pol"))

comps <- data.frame(condition_1 = c("PMC", "Tet", "Mic", "Bin"), condition_2 = c("Tet", "Mic", "Bin", "Pol"))

aSwitchList <- importRdata(isoformCountMatrix = eq$counts, isoformRepExpression = eq$abundance, designMatrix = myDesign, comparisonsToMake = comps, isoformExonAnnoation = "/scratch/punim0687/agolicz/rapa.isoforms.08032019/7.IsoformSwitchAnalyzeR/analysis/Rapa.03042019.merged.renamed.nm.r.gtf", showProgress = TRUE)

#Sleuth
options(stringsAsFactors = FALSE)
c1<-data.frame(sample = c("B1", "B2", "B3", "B4", "B5", "B6"), condition = c("PMC", "PMC", "PMC", "Tet", "Tet", "Tet") , path = c("/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B1", "/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B2", "/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B3", "/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B4", "/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B5", "/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B6"))
c2<-data.frame(sample = c("B4", "B5", "B6", "B7", "B8", "B9"), condition = c("Tet", "Tet", "Tet", "Mic", "Mic", "Mic") , path = c("/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B4", "/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B5", "/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B6", "/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B7", "/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B8", "/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B9"))
c3<-data.frame(sample = c("B7", "B8", "B9", "B10", "B11", "B12"), condition = c("Tet", "Tet", "Tet", "Bin", "Bin", "Bin") , path = c("/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B7", "/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B8", "/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B9", "/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B10", "/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B11", "/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B12"))
c4<-data.frame(sample = c("B10", "B11", "B12", "B13", "B14", "B15"), condition = c("Tet", "Tet", "Tet", "Pol", "Pol", "Pol") , path = c("/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B10", "/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B11", "/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B12", "/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B13", "/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B14", "/scratch/punim0687/agolicz/rapa.isoforms.08032019/6.kallisto/results/B15"))

so.c1 <- sleuth_prep(c1, ~condition)
so.c1 <- sleuth_fit(so.c1)
so.c1 <- sleuth_fit(so.c1, ~1, 'reduced')
so.c1 <- sleuth_lrt(so.c1, 'reduced', 'full')
sleuth_table.c1 <- sleuth_results(so.c1, 'reduced:full', 'lrt', show_all = TRUE)
sleuth_table.c1$condition_1<-"PMC"
sleuth_table.c1$condition_2<-"Tet"

so.c2 <- sleuth_prep(c2, ~condition)
so.c2 <- sleuth_fit(so.c2)
so.c2 <- sleuth_fit(so.c2, ~1, 'reduced')
so.c2 <- sleuth_lrt(so.c2, 'reduced', 'full')
sleuth_table.c2 <- sleuth_results(so.c2, 'reduced:full', 'lrt', show_all = TRUE)
sleuth_table.c2$condition_1<-"Tet"
sleuth_table.c2$condition_2<-"Mic"

so.c3 <- sleuth_prep(c3, ~condition)
so.c3 <- sleuth_fit(so.c3)
so.c3 <- sleuth_fit(so.c3, ~1, 'reduced')
so.c3 <- sleuth_lrt(so.c3, 'reduced', 'full')
sleuth_table.c3 <- sleuth_results(so.c3, 'reduced:full', 'lrt', show_all = TRUE)
sleuth_table.c3$condition_1<-"Mic"
sleuth_table.c3$condition_2<-"Bin"

so.c4 <- sleuth_prep(c4, ~condition)
so.c4 <- sleuth_fit(so.c4)
so.c4 <- sleuth_fit(so.c4, ~1, 'reduced')
so.c4 <- sleuth_lrt(so.c4, 'reduced', 'full')
sleuth_table.c4 <- sleuth_results(so.c4, 'reduced:full', 'lrt', show_all = TRUE)
sleuth_table.c4$condition_1<-"Bin"
sleuth_table.c4$condition_2<-"Pol"

dfst <- do.call("rbind", list(sleuth_table.c1, sleuth_table.c2, sleuth_table.c3, sleuth_table.c4))

t2g<-read.csv("transcript2gene.tsv",sep="\t",header=T, stringsAsFactors = FALSE)
sog.c1 <- sleuth_prep(c1,  ~condition, target_mapping = t2g, aggregation_column = 'gene_id', num_cores = 1)
sog.c1 <- sleuth_fit(sog.c1, ~condition, 'full')
sog.c1 <- sleuth_fit(sog.c1, ~1, 'reduced')
sog.c1 <- sleuth_lrt(sog.c1, 'reduced', 'full')
sleuth_table.g.c1 <- sleuth_results(sog.c1, 'reduced:full', 'lrt', show_all = TRUE)
sleuth_table.g.c1$condition_1<-"PMC"
sleuth_table.g.c1$condition_2<-"Tet"

sog.c2 <- sleuth_prep(c2,  ~condition, target_mapping = t2g, aggregation_column = 'gene_id', num_cores = 1)
sog.c2 <- sleuth_fit(sog.c2, ~condition, 'full')
sog.c2 <- sleuth_fit(sog.c2, ~1, 'reduced')
sog.c2 <- sleuth_lrt(sog.c2, 'reduced', 'full')
sleuth_table.g.c2 <- sleuth_results(sog.c2, 'reduced:full', 'lrt', show_all = TRUE)
sleuth_table.g.c2$condition_1<-"Tet"
sleuth_table.g.c2$condition_2<-"Mic"

sog.c3 <- sleuth_prep(c3,  ~condition, target_mapping = t2g, aggregation_column = 'gene_id', num_cores = 1)
sog.c3 <- sleuth_fit(sog.c3, ~condition, 'full')
sog.c3 <- sleuth_fit(sog.c3, ~1, 'reduced')
sog.c3 <- sleuth_lrt(sog.c3, 'reduced', 'full')
sleuth_table.g.c3 <- sleuth_results(sog.c3, 'reduced:full', 'lrt', show_all = TRUE)
sleuth_table.g.c3$condition_1<-"Mic"
sleuth_table.g.c3$condition_2<-"Bin"

sog.c4 <- sleuth_prep(c4,  ~condition, target_mapping = t2g, aggregation_column = 'gene_id', num_cores = 1)
sog.c4 <- sleuth_fit(sog.c4, ~condition, 'full')
sog.c4 <- sleuth_fit(sog.c4, ~1, 'reduced')
sog.c4 <- sleuth_lrt(sog.c4, 'reduced', 'full')
sleuth_table.g.c4 <- sleuth_results(sog.c4, 'reduced:full', 'lrt', show_all = TRUE)
sleuth_table.g.c4$condition_1<-"Bin"
sleuth_table.g.c4$condition_2<-"Pol"

dfsg <- do.call("rbind", list(sleuth_table.g.c1, sleuth_table.g.c2, sleuth_table.g.c3, sleuth_table.g.c4))

### Isoform level
aSwitchList$isoformFeatures$iso_q_value <- dfst$qval[match(
    paste0(
        aSwitchList$isoformFeatures$isoform_id,
        aSwitchList$isoformFeatures$condition_1,
        aSwitchList$isoformFeatures$condition_2
    ),
    paste0(
        dfst$target_id,
        dfst$condition_1,
        dfst$condition_2
    )
)]

### Gene level
aSwitchList$isoformFeatures$gene_q_value <- dfsg$qval[match(
    paste0(
        aSwitchList$isoformFeatures$gene_id,
        aSwitchList$isoformFeatures$condition_1,
        aSwitchList$isoformFeatures$condition_2
    ),
    paste0(
        dfsg$target_id,
        dfsg$condition_1,
        dfsg$condition_2
    )
)]

sub.t <-scan("/scratch/punim0687/agolicz/rapa.isoforms.08032019/7.IsoformSwitchAnalyzeR/analysis/Rapa.03042019.merged.renamed.nm.r.ni.nh.expressed.coding.300.gtf.transcripts", what="")
aSwitchListS<-subsetSwitchAnalyzeRlist(aSwitchList,aSwitchList$isoformFeatures$isoform_id %in% sub.t)

aSwitchListFiltered <- preFilter(aSwitchListS, geneExpressionCutoff = 5, isoformExpressionCutoff = 2, removeSingleIsoformGenes = TRUE)

aSwitchListFilteredA1 <- isoformSwitchTestDEXSeq(aSwitchListFiltered, alpha = 0.05, dIFcutoff = 0.2, correctForConfoundingFactors=FALSE, overwriteIFvalues = FALSE, reduceToSwitchingGenes = TRUE)

aSwitchListFilteredA1 <- preFilter(aSwitchListFilteredA1, reduceToSwitchingGenes=TRUE, keepIsoformInAllConditions = FALSE, alpha=0.05, dIFcutoff = 0.2)

aSwitchListFilteredA1 <- analyzeORF(aSwitchListFilteredA1, BSgenome.Brapa.3.0)

aSwitchListFilteredA1 <- extractSequence(aSwitchListFilteredA1, BSgenome.Brapa.3.0, onlySwitchingGenes = TRUE, alpha = 0.05, dIFcutoff = 0.2,  pathToOutput = "/scratch/punim0687/agolicz/rapa.isoforms.08032019/7.IsoformSwitchAnalyzeR/analysis")

save(aSwitchListFilteredA1, file="aSwitchListFilteredA1.obj")

