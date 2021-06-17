library("topGO")
library("plyr")

writeTopGO <- function(GOterm, top, txt, mg, uni){
    # give properly formatted background in format: gene	go;go;go
    annAT <- readMappings(uni, sep="\t", IDsep=";")
    allgenes <- names(annAT)
    # give file with your genes of interest, one gene_id per line
    mygenes <-scan(mg,what="")
    geneList <- factor(as.integer(allgenes %in% mygenes))
    names(geneList) <- allgenes
    GOdata <-new ("topGOdata", ontology = GOterm, allGenes = geneList, nodeSize = top, annot=annFUN.gene2GO, gene2GO=annAT)

    resultFisher <- runTest(GOdata, algorithm = "weight", statistic ="fisher")
    allRes <-GenTable(GOdata, resFisher = resultFisher, topNodes = 100)
    names(allRes)[length(allRes)] <- "p.value"
    write.table(allRes, file=txt, sep="\t",quote=FALSE,row.names=FALSE)
    allRes$nc<-rep(txt,dim(allRes)[1])
    allGO = genesInTerm(GOdata)
    print(allGO["GO:0008380"])
    sampleGO = lapply(allGO,function(x) x[x %in% mygenes] )
    df.g <- ldply (sampleGO, data.frame)
    names(df.g)<-c("GO.ID","Genes")
    df.res<-merge(allRes,df.g,by="GO.ID")
    write.table(df.res,file=paste(txt,".go.table",sep=""),sep="\t",row.names=F,quote=F)
}

#different node sizes give slightly different results, 5 or 10 are recommended as most meaningful by the author

#three command line arguments
#1. list with interesting genes
#2. annotated gene universe
#3. name of the file

args <- commandArgs(trailingOnly = TRUE)
mG <- args[1]
universe <- args[2]
name <- args[3]
bp <- paste("10topBP.", name,".txt", sep="")
mf <- paste("10topMF.", name,".txt", sep="")
cc <- paste("10topCC.", name,".txt", sep="")
writeTopGO("BP", 10, bp, mG, universe)
writeTopGO("MF", 10, mf, mG, universe)
writeTopGO("CC", 10, cc, mG, universe)

