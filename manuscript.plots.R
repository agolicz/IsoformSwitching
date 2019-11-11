library(pheatmap)
library(RColorBrewer)
library(viridis)

#heatmap Fig3A

t<-read.csv("IRG.PMC.Tet.t1.tpm",sep="\t",header=T)
t2<-read.csv("IRG.PMC.Tet.t2.tpm",sep="\t",header=T)
t$Gene<-NULL
t2$Gene<-NULL
h1<- pheatmap(log1p(t), cluster_cols=F, scale="row", show_rownames=F, color=magma(100))
row.order = h1$tree_row$order
values.to.show = t2[row.order,]
pheatmap(values.to.show, cluster_cols=F, cluster_rows=F, scale="row", show_rownames=F, color=magma(100))


# Position plot Fig4A

t<-read.csv("introns.position",sep="\t", header=F)
ggplot(t, aes(t$V7))+geom_density(fill="lightgreen", alpha=0.5)+theme_minimal()+xlab("Transcript position 5' -> 3'")+ylab("Density")+theme(text = element_text(size=14))
ggsave("Density.plot.pdf", height=4, width=4)

# GC plot Fig4B

t<-read.csv("all.gc",sep="\t")
t$Type<-factor(t$Type, levels=c("IR", "non-IR AS", "CS"))
ggplot(t, aes(x=Type, y=GC, fill=Type))+geom_violin(alpha=0.5)+theme_minimal()+xlab("")+ylab("GC content")+ theme(legend.position = "none")+stat_summary(fun.y=median, geom="point", size=2, color="black")+theme(text = element_text(size=14))
ggsave("GC.pdf", height=4, width=4)

# Length plot Fig4C

t<-read.csv("all.len",sep="\t")
t$Type<-factor(t$Type, levels=c("IR", "non-IR AS", "CS"))
ggplot(t, aes(x=Type, y=log(Length), fill=Type))+geom_violin(alpha=0.5)+theme_minimal()+xlab("")+ylab("Length [log(bp)]")+ theme(legend.position = "none") +stat_summary(fun.y=median, geom="point", size=2, color="black")+theme(text = element_text(size=14))

ggsave("Length.pdf", height=4, width=4)

# Expression plots Fig5A

t<-read.csv("expression.transcripts.all",sep="\t")
t$Type<-factor(t$Type, levels=c("IR transcripts", "non-IR transcripts"))
ggplot(t, aes(log1p(Expression), fill=Type))+geom_density(aes(linetype=Type), alpha=0.25)+theme_minimal()+xlab("log1p(Expression)")+ylab("Density")+ theme(legend.title = element_blank())+theme(text = element_text(size=14))

ggsave("Expression.transcripts.pdf", height=4, width=5)

t<-read.csv("expression.transcripts.tet",sep="\t")
t$Type<-factor(t$Type, levels=c("IR transcripts", "non-IR transcripts"))
ggplot(t, aes(log1p(Expression), fill=Type))+geom_density(aes(linetype=Type), alpha=0.25)+theme_minimal()+xlab("log1p(Expression)")+ylab("Density")+ theme(legend.title = element_blank())+theme(text = element_text(size=14))

ggsave("Expression.transcripts.tet.pdf", height=4, width=5)

t<-read.csv("expression.genes.all",sep="\t")
t$Type<-factor(t$Type, levels=c("IR genes", "non-IR genes"))
ggplot(t, aes(log1p(Expression), fill=Type))+geom_density(aes(linetype=Type), alpha=0.25)+theme_minimal()+xlab("log1p(Expression)")+ylab("Density")+ theme(legend.title = element_blank())+theme(text = element_text(size=14))

ggsave("Expression.genes.pdf", height=4, width=5)

t<-read.csv("expression.genes.tet",sep="\t")
t$Type<-factor(t$Type, levels=c("IR genes", "non-IR genes"))
ggplot(t, aes(log1p(Expression), fill=Type))+geom_density(aes(linetype=Type), alpha=0.25)+theme_minimal()+xlab("log1p(Expression)")+ylab("Density")+ theme(legend.title = element_blank())+theme(text = element_text(size=14))

ggsave("Expression.genes.tet.pdf", height=4, width=5)

#Tau plot Fig5B

t<-read.csv("all.tau",sep="\t")
t$Type<-factor(t$Type, levels=c("IR transcripts", "non-IR transcripts tetrad", "non-IR transcripts"))
ggplot(t, aes(x=Type, y=Tau, fill=Type))+geom_violin(alpha=0.5)+theme_minimal()+xlab("")+ylab("Tau")+ theme(legend.position = "none") +stat_summary(fun.y=median, geom="point", size=2, color="black")+ theme(axis.text.x = element_text(angle = 20))+theme(text = element_text(size=14))

ggsave("Tau.pdf", height=4, width=4)


# Splice site plots Fig5C

t<-read.csv("as.file1.file2.s.s.ss.s1.seq.2.score",sep="\t")
t$Type<-factor(t$Type, levels=c("IR", "non-IR AS"))
ggplot(t, aes(x=Type, y=Score, fill=Type))+geom_violin(alpha=0.5)+theme_minimal()+xlab("")+ylab("Score")+ theme(legend.position = "none") +stat_summary(fun.y=median, geom="point", size=2, color="black")+theme(text = element_text(size=14))

ggsave("AS.3.pdf", height=4, width=4)

t<-read.csv("as.file1.file2.s.s.ss.s2.seq.2.score",sep="\t")
t$Type<-factor(t$Type, levels=c("IR", "non-IR AS"))
ggplot(t, aes(x=Type, y=Score, fill=Type))+geom_violin(alpha=0.5)+theme_minimal()+xlab("")+ylab("Score")+ theme(legend.position = "none") +stat_summary(fun.y=median, geom="point", size=2, color="black")+theme(text = element_text(size=14))

ggsave("AS.5.pdf", height=4, width=4)

t<-read.csv("cs.file1.file2.s.s.ss.s1.seq.2.score",sep="\t")
t$Type<-factor(t$Type, levels=c("IR", "CS"))
ggplot(t, aes(x=Type, y=Score, fill=Type))+geom_violin(alpha=0.5)+theme_minimal()+xlab("")+ylab("Score")+ theme(legend.position = "none") +stat_summary(fun.y=median, geom="point", size=2, color="black")+theme(text = element_text(size=14))

ggsave("CS.3.pdf", height=4, width=4)

t<-read.csv("cs.file1.file2.s.s.ss.s2.seq.2.score",sep="\t")
t$Type<-factor(t$Type, levels=c("IR", "CS"))
ggplot(t, aes(x=Type, y=Score, fill=Type))+geom_violin(alpha=0.5)+theme_minimal()+xlab("")+ylab("Score")+ theme(legend.position = "none") +stat_summary(fun.y=median, geom="point", size=2, color="black")+theme(text = element_text(size=14))

ggsave("CS.5.pdf", height=4, width=4)

#SDG8 FigS5

t<-read.csv("SDG8.expression",sep="\t", row.names=1)
pheatmap(log1p(t), cluster_cols=F, cluster_rows=F, scale="row", show_rownames=T, color=magma(100))

# Splicing summary Fig1C

t<-read.csv("pairs.fraction.txt", sep="\t")
t2 <- t[order(t$Fraction),] 
t2$Type <- factor(t2$Type, levels = t2$Type)
ggplot(data=t2, aes(x=Type, y=Fraction, fill=Type)) + geom_bar(stat="identity", color="black")+coord_flip()+theme_minimal()+xlab("")
ggsave("splicing.summary.global.pdf", width=5, height=4)

#Expression trajectory lineplot Fig3C

t<-read.csv("all.tpm", sep="\t")
t$Gene<-NULL
mdata <- melt(t, id=c("Type"))
mdata$value<-log1p(mdata$value)
aggdata <-aggregate(value~Type+variable, data=mdata, FUN=median, na.rm=TRUE)
ggplot(aggdata, aes(x=variable, y=value, group=Type)) +  geom_line(aes(linetype=Type, color=Type), size=1)+ geom_point(aes(color=Type), size=2)+theme_minimal()+xlab("")+ylab("Median of expression [log1p(TPM)]")
ggsave("expression.trajectories.pdf", width=5, height=4)

t<-read.csv("all.m.tpm", sep="\t")
t$Gene<-NULL
mdata <- melt(t, id=c("Type"))
mdata$value<-log1p(mdata$value)
aggdata <-aggregate(value~Type+variable, data=mdata, FUN=mean, na.rm=TRUE)
ggplot(aggdata, aes(x=variable, y=value, group=Type)) +  geom_line(aes(linetype=Type, color=Type), size=1)+ geom_point(aes(color=Type), size=2)+theme_minimal()+xlab("")+ylab("Median of expression [log1p(TPM)]")+theme(text = element_text(size=14), axis.text.x = element_text(angle=45, hjust=1))

ggsave("expression.m.trajectories.pdf", width=6, height=4)

#Expression trajectory boxplot FigS3

t<-read.csv("all.m.tpm", sep="\t")
t$Gene<-NULL
mdata <- melt(t, id=c("Type"))
mdata$value<-log1p(mdata$value)
ggplot(mdata, aes(x=variable, y=value, fill=Type)) +  geom_boxplot() +theme_minimal()+xlab("")+ylab("Expression [log1p(TPM)]")+theme(text = element_text(size=14), axis.text.x = element_text(angle=45, hjust=1))
ggsave("FigS3.pdf", width=6, height=4)

