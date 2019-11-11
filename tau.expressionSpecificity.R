tau<-function(x){
    if(any(is.na(x))) stop('NA\'s need to be 0.')
    if(any(x<0)) stop('Negative input values not permitted. Maybe data is log transformed?')
    t<-sum(1-x/max(x))/(length(x)-1)
}
e<-read.csv("mean.expression.subset.tsv", sep="\t", row.names=1)
tu <- apply(e, 1, tau)
sp<-data.frame(Transcript=row.names(e),Tau=tu)
write.table(sp,file="tau.tsv",sep="\t",row.names=F,quote=F)
