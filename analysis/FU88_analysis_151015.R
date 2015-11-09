## A First look at the oxBS/BS pilot. 

## PL 151014
## 

library(minfi)
setwd("~/projects/FU88/analysis")

#Set main dir for gt data
baseDir <- '~/projects/FU88/FU88_150916_ResultReport/FU88_150916_IDAT'


#--- Preprocessing ----------------------------------------

#Read sample info
targets <- read.450k.sheet(baseDir)

#Add BS/oxBS status
targets$bs <- c('BS','oxBS','BS','oxBS','BS','oxBS','BS','oxBS','BS','oxBS','BS','oxBS','CTL','BS','oxBS','BS','oxBS','BS','oxBS','BS','oxBS','BS','oxBS','BS','oxBS')

#Read experiment data
RGset <- read.450k.exp(targets = targets)
#initial qc report
qcReport(RGset, sampNames = targets$Sample_Name, sampGroups= targets$bs, pdf='qcReport.rgset.bs_grp.pdf')

#filter on detP 0.01
detP <- detectionP(RGset)
failed <- detP > 0.01
maxFail <- 0
mset <- preprocessRaw(RGset)
#drop all rows with a failed probe in any sample
mset <- mset[rowSums(failed) <= maxFail, ]

#Run SWAN norm. and write pdf/save Rdata
mset.swan <-preprocessSWAN(RGset, mSet = mset)
pdf('densityPlot.msetSwan.pdf')
densityPlot(mset.swan, sampGroups=targets$bs, main= sprintf('Beta values for p<0.01 probes (n=%s)', nrow(mset.swan)))
dev.off()

#save.image('mset.swan.20151019.Rdata')


#--- PCA --------------------------------------------------
betaVals <- getBeta(mset.swan)
pcaResults <- prcomp(t(betaVals))
pdf('pca.plot.pdf', w=14/2.54, h=14/2.54)
par(las=1, mgp=c(2, 0.7, 0))
plot(pcaResults$x, main='PCA on beta values', 
    xlab=sprintf('PC1 (sd: %s%%)', round (100 * (pcaResults$sdev[1] / sum(pcaResults$sdev)))),
    ylab=sprintf('PC2 (sd: %s%%)', round (100 * (pcaResults$sdev[2] / sum(pcaResults$sdev)))),
    pch=19, col=ifelse(targets$bs=="BS", "blue", "red"), xlim=c(-100, 150))
grid(col='grey30')
text(x=pcaResults$x[, "PC1"] + 35, y=pcaResults$x[, "PC2"], labels=targets$Sample_Name, cex=0.5)
dev.off()

#--- Hierarchical clustering
pdf('hierarchical_clustering.plot.pdf')
plot(hclust(dist(t(betaVals))), labels=targets$Sample_Name)
dev.off()

#save.image('post.PCA.HClust.Rdata')

#--- Correlations in Beta values --------------------------
library(corrplot)
correlation.percent <- round(cor(getBeta(mset.swan)) * 100, 1)
colnames(correlation.percent) <- targets$Sample_Name
rownames(correlation.percent) <- targets$Sample_Name

#remove diag of 1 from max calc for plot
diag(correlation.percent) <- min(correlation.percent)

#set colormap
col <- colorRampPalette(c("yellow", "yellow", "orange", "red"))

pdf('corrplot.beta.pdf', pointsize=12)
par(xpd=NA, cex=0.35)
corrplot(correlation.percent, type="upper", 
         order="alphabet", cl.lim= c(min(correlation.percent), 
                                     max(correlation.percent)), diag=FALSE, 
         col=col(20), addCoefasPercent = FALSE,
         addCoef.col="black", tl.col="black", is.corr=FALSE,
         method='circle', mar=c(0,0,4,0), main='Pairwise correlation of beta values', tl.cex=1, cl.cex=1, pch.cex=0.7)
dev.off()

###--- Mapping --------------------------------------------
#Currently not implemented


###--- 5hmC analysis
# Using M-values for better statistics
M <- getM(mset.swan, type="beta", betaThreshold = 0.001)

#split analysis on sample groups
ALL.set <- which(targets$Sample_Group == 'ALL')
CEP.set <- which(targets$Sample_Group == 'CEP')
BG.set <- which(targets$Sample_Group == 'BG')
CEGX.set <- which(targets$Sample_Group == 'CEGX')

#Test for difference in methylation between BS and oxBS = hydroxymethylation
hydroxy.probes.ALL <- dmpFinder(M[,ALL.set], pheno=targets$bs[ALL.set], type="categorical", shrinkVar=TRUE)
hydroxy.probes.CEP <- dmpFinder(M[,CEP.set], pheno=targets$bs[CEP.set], type="categorical", shrinkVar=TRUE)
hydroxy.probes.BG <- dmpFinder(M[,BG.set], pheno=targets$bs[BG.set], type="categorical", shrinkVar=TRUE)
hydroxy.probes.CEGX <- dmpFinder(M[,CEGX.set], pheno=targets$bs[CEGX.set], type="categorical", shrinkVar=TRUE)

#Function to return a matrix of means per samplegroup
aggSamp<-function(X, sampGroup){
  cnames<-sort(unique(sampGroup))
  aggMat<-matrix(data= NA, nrow= nrow(X), ncol= length(cnames))
  colnames(aggMat)<- cnames
  rownames(aggMat)<- rownames(X)
  for(i in 1:length(cnames)){
    xcol<-which(sampGroup == cnames[i])
    aggMat[, i] <- rowMeans(X[, xcol])
  }
  return(aggMat)
}

#Make tables with mean beta per group and probe
mean.beta.CEP <- aggSamp(X= getBeta(mset.swan)[,CEP.set], sampGroup= targets$bs[CEP.set])
mean.beta.CEGX <- aggSamp(X= getBeta(mset.swan)[,CEGX.set], sampGroup= targets$bs[CEGX.set])

hydroxy.probes.CEP <- cbind(hydroxy.probes.CEP, mean.beta.CEP[rownames(hydroxy.probes.CEP),])


