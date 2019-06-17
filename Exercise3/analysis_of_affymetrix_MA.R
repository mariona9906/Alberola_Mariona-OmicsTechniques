## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ------------------------------------------------------------------------
setwd("~/Alberola_Mariona-OmicsTechniques/Exercise3")
targets <- read.csv("~/Alberola_Mariona-OmicsTechniques/Exercise3/dades/targets.txt", sep=";")
targets

workingDir <-getwd()
dataDir <- file.path(workingDir, "dades")
resultsDir <- file.path(workingDir, "results")
setwd(resultsDir)


## ------------------------------------------------------------------------
if (!require(BiocManager)) install.packages("BiocManager")

installifnot <- function (pkg){
  if (!require(pkg, character.only=T)){
    BiocManager::install(pkg)
}else{
  require(pkg, character.only=T)
  }
}

installifnot("pd.mogene.1.0.st.v1")
installifnot("mogene10sttranscriptcluster.db")
installifnot("oligo")
installifnot("limma")
installifnot("Biobase")
installifnot("arrayQualityMetrics")
installifnot("genefilter")
installifnot("multtest")
installifnot("annotate")
installifnot("xtable")
installifnot("gplots")
installifnot("scatterplot3d")


## ------------------------------------------------------------------------
#TARGETS
workingDir <-getwd()
dataDir <- file.path(workingDir, "dades")
resultsDir <- file.path(workingDir, "results")
setwd(resultsDir)
class(targets)
expressions <- read.delim("./dades/GSE106402_series_matrix.txt", header = TRUE, row.names = 1)

rawData <-new("ExpressionSet", exprs=as.matrix(expressions))

#CELFILES
#CELfiles <- list.celfiles(file.path(dataDir))
#CELfiles
#f <- read.celfiles(file.path(dataDir,CELfiles))


#DEFINE SOME VARIABLES FOR PLOTS
sampleNames <- as.character(targets$Sample_Name)
sampleColor <- as.character(targets$Colors)



## ------------------------------------------------------------------------
boxplot(rawData, which="all",las=2, main="Intensity distribution of RAW data", 
        cex.axis=0.6, col=sampleColor, names=sampleNames)

#HIERARQUICAL CLUSTERING
clust.euclid.average <- hclust(dist(t(exprs(rawData))),method="average")
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of RawData", 
     cex=0.7,  hang=-1)

#PRINCIPAL COMPONENT ANALYSIS
plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE, formapunts=NULL, myCex=0.8,...)
{
  pcX<-prcomp(t(X), scale=scale) # o prcomp(t(X))
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colors, pch=formapunts, 
       xlim=c(min(pcX$x[,1])-100000, max(pcX$x[,1])+100000),ylim=c(min(pcX$x[,2])-100000, max(pcX$x[,2])+100000))
  text(pcX$x[,1],pcX$x[,2], labels, pos=3, cex=myCex)
  title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep=" "), cex=0.8)
}

plotPCA(exprs(rawData), labels=sampleNames, dataDesc="raw data", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)

#SAVE TO A FILE
pdf(file.path(resultsDir, "QCPlots_Raw.pdf"))
boxplot(rawData, which="all",las=2, main="Intensity distribution of RAW data", 
        cex.axis=0.6, col=sampleColor, names=sampleNames)
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of samples of RawData", 
     cex=0.7,  hang=-1)
plotPCA(exprs(rawData), labels=sampleNames, dataDesc="raw data", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)
dev.off()


## ------------------------------------------------------------------------

#eset<-rma(rawData)

#write.exprs(eset, file.path(resultsDir, "NormData.txt"))

installifnot("affyPLM")
eset<-normalize.ExpressionSet.invariantset(rawData)
?rma


write.exprs(eset, file.path(resultsDir, "NormData.txt"))




## ------------------------------------------------------------------------
#BOXPLOT
boxplot(eset, las=2, main="Intensity distribution of Normalized data", cex.axis=0.6, 
        col=sampleColor, names=sampleNames)

#HIERARQUICAL CLUSTERING
clust.euclid.average <- hclust(dist(t(exprs(eset))),method="average")
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of NormData", 
     cex=0.7,  hang=-1)

#PRINCIPAL COMPONENT ANALYSIS
plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE, formapunts=NULL, myCex=0.8,...)
{
  pcX<-prcomp(t(X), scale=scale) # o prcomp(t(X))
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colors, pch=formapunts, 
       xlim=c(min(pcX$x[,1])-10, max(pcX$x[,1])+10),ylim=c(min(pcX$x[,2])-10, max(pcX$x[,2])+10))
  text(pcX$x[,1],pcX$x[,2], labels, pos=3, cex=myCex)
  title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep=" "), cex=0.8)
}

plotPCA(exprs(eset), labels=sampleNames, dataDesc="NormData", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)

#SAVE TO A FILE
pdf(file.path(resultsDir, "QCPlots_Norm.pdf"))
boxplot(eset, las=2, main="Intensity distribution of Normalized data", cex.axis=0.6, 
        col=sampleColor, names=sampleNames)
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of NormData", 
     cex=0.7,  hang=-1)
plotPCA(exprs(eset), labels=sampleNames, dataDesc="selected samples", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)
dev.off()


#ARRAY QUALITY METRICS

arrayQualityMetrics(expressionset = eset, reporttitle="QualityControl", force=TRUE)



## ------------------------------------------------------------------------
annotation(eset) <- "org.Mm.eg.db"
eset_filtered <- nsFilter(eset, var.func=IQR,
         var.cutoff=0.75, var.filter=TRUE,
         filterByQuantile=TRUE)
#NUMBER OF GENES OUT
print(eset_filtered$filter.log$numLowVar)

#NUMBER OF GENES IN
print(eset_filtered$eset)


## ------------------------------------------------------------------------
#CONTRAST MATRIX.lINEAR MODEL

design <-model.matrix(~ 0+targets$Experiment_Types)
design
colnames(design)<-c("cn", "un", "cs", "us")
rownames(design)<- targets$Sample_Name 
print(design)

#COMPARISON
cont.matrix1 <- makeContrasts (
  cnVSun = cn-un,
  cnVScs = cn-cs,
  cnVSus = cn-us,
  unVScs = un-cs,
  unVSus = un-us,
  csVSus = cs-us,
  levels=design)
## Where: 
## cn = siCONT_noTSA, 
## un = siUHRF1_noTSA, 
## cs = siCONT_TSA, 
## us = siUHRF1_TSA
comparison1 <- "siCONT_noTSA_vs_siUHRF1_noTSA"
comparison2 <- "siCONT_noTSA_vs_siCONT_TSA"
comparison3 <- "siCONT_noTSA_vs_siUHRF1_TSA"
comparison4 <- "siUHRF1_noTSA_vs_siCONT_TSA"
comparison5 <- "siUHRF1_noTSA_vs_siUHRF1_TSA"
comparison6 <- "siCONT_TSA_vs_siUHRF1_TSA"

cont.matrix1 

#MODEL FIT

fit1 <- lmFit(eset, design)
fit.main1 <- contrasts.fit(fit1, cont.matrix1)
fit.main1 <- eBayes(fit.main1)

#######




## ------------------------------------------------------------------------

#FILTER BY FALSE DISCOVERY RATE AND FOLD CHANGE
#topTab <-  topTable (fit.main1, number=nrow(fit.main1), coef="Induced.vs.WT", adjust="fdr",lfc=abs(3))
topTab_cnVSun <- topTable (fit.main1, number=nrow(fit.main1), coef="cnVSun", adjust="fdr"); head(topTab_cnVSun)
topTab_cnVScs <- topTable (fit.main1, number=nrow(fit.main1), coef="cnVScs", adjust="fdr"); head(topTab_cnVScs)
topTab_cnVSus <- topTable (fit.main1, number=nrow(fit.main1), coef="cnVSus", adjust="fdr"); head(topTab_cnVSus)
topTab_unVScs <- topTable (fit.main1, number=nrow(fit.main1), coef="unVScs", adjust="fdr"); head(topTab_unVScs)
topTab_unVSus <- topTable (fit.main1, number=nrow(fit.main1), coef="unVSus", adjust="fdr"); head(topTab_unVSus)
topTab_csVSus <- topTable (fit.main1, number=nrow(fit.main1), coef="csVSus", adjust="fdr"); head(topTab_csVSus)


## ------------------------------------------------------------------------
#EXPORTED TO CSV AND HTML FILE
write.csv2(topTab_cnVSun, file= file.path(resultsDir,paste("Selected.Genes.in.comparison.",
                                                    comparison1, ".csv", sep = "")))
write.csv2(topTab_cnVScs, file= file.path(resultsDir,paste("Selected.Genes.in.comparison.",
                                                    comparison2, ".csv", sep = "")))
write.csv2(topTab_cnVSus, file= file.path(resultsDir,paste("Selected.Genes.in.comparison.",
                                                    comparison3, ".csv", sep = "")))
write.csv2(topTab_unVScs, file= file.path(resultsDir,paste("Selected.Genes.in.comparison.",
                                                    comparison4, ".csv", sep = "")))
write.csv2(topTab_unVSus, file= file.path(resultsDir,paste("Selected.Genes.in.comparison.",
                                                    comparison5, ".csv", sep = "")))
write.csv2(topTab_csVSus, file= file.path(resultsDir,paste("Selected.Genes.in.comparison.",
                                                    comparison6, ".csv", sep = "")))




## ------------------------------------------------------------------------
volcanoplot(fit.main1, highlight=10, names=fit.main1$ID, 
            main = paste("Differentially expressed genes", colnames(cont.matrix1), sep="\n"))
abline(v = c(-3, 3))


pdf(file.path(resultsDir,"Volcanos.pdf"))
volcanoplot(fit.main1, highlight = 10, names = fit.main1$ID, 
            main = paste("Differentially expressed genes", colnames(cont.matrix1), sep = "\n"))
abline(v = c(-3, 3))
dev.off()
class(topTab_cnVScs)


## ------------------------------------------------------------------------
#PREPARE THE DATA
my_frame <- data.frame(exprs(eset))
head(my_frame)
HMdata <- merge(my_frame, topTab_cnVSun, by.x = 0, by.y = 0)
rownames(HMdata) <- HMdata$Row.names
HMdata <- HMdata[, -c(1,10:15)]
head(HMdata)
HMdata2 <- data.matrix(HMdata, rownames.force=TRUE)
head(HMdata2)
write.csv2(HMdata2, file = file.path(resultsDir,"Data2HM.csv"))

#HEATMAP PLOT
my_palette <- colorRampPalette(c("blue", "red"))(n = 299)

heatmap.2(HMdata2,
          Rowv=TRUE,
          Colv=TRUE,
          main="HeatMap siCONT_noTSA_vs_siUHRF1_noTSA",
          scale="row",
          col=my_palette,
          sepcolor="white",
          sepwidth=c(0.05,0.05),
          cexRow=0.5,
          cexCol=0.9,
          key=TRUE,
          keysize=1.5,
          density.info="histogram",
          ColSideColors=sampleColor,
          tracecol=NULL,
          srtCol=30)

#EXPORT TO PDF FILE
pdf(file.path(resultsDir,"HeatMap InducedvsWT.pdf"))
heatmap.2(HMdata2,
          Rowv=TRUE,
          Colv=TRUE,
          main="HeatMap Induced.vs.WT FC>=3",
          scale="row",
          col=my_palette,
          sepcolor="white",
          sepwidth=c(0.05,0.05),
          cexRow=0.5,
          cexCol=0.9,
          key=TRUE,
          keysize=1.5,
          density.info="histogram",
          ColSideColors=sampleColor,
          tracecol=NULL,
          srtCol=30)
dev.off()


## ------------------------------------------------------------------------
BiocManager::install("hgug4110b.db")



require(hgug4110b.db)
require(hugene20sttranscriptcluster.db)

class(hgug4110b.db)
probes_itvsct<-rownames(head(topTab_itvsct))
antot_itvsct<-select(hugene20sttranscriptcluster.db,probes_itvsct,
                     columns = c("ENTREZID","SYMBOL","GENENAME"))
anot_itvsct


## ------------------------------------------------------------------------
require(hgu133a.db)
columns(hgu133a.db)


