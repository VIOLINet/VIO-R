#!/usr/bin/env Rscript
args = commandArgs( trailingOnly=TRUE )

#Parse 
if ( length(args)<2 ) {
  stop( "At least two groups must be provided.\nExample: Rscript limma_one_gse.R GSM733843,GSM733844 GSM733852,GSM733853", call.=FALSE )
}

library(Biobase)
library(GEOquery)
library(GEOmetadb)
library(limma)

groups = c()
samples = c()
count = 0
for ( i in 1:length(args) ) {
  tokens = strsplit(args[i],',')[[1]]
  for ( j in 1:length(tokens) ) {
    count = count + 1
    groups[count] = i
    samples[count] = tokens[j]
  }
}

if ( !file.exists( "GEOmetadb.sqlite" ) ) {
  getSQLiteFile()
}

geoMeta = geoConvert(samples, out_type = c("gpl","gse"))

gplList = unique(geoMeta$gpl[,2])
if ( length(gplList)>1 ) {
  stop( "Multiple platform detected. Please check your input.", call.=FALSE )
}

gseList = unique(geoMeta$gse[,2])
if ( length(gseList)>1 ) {
  for ( i in 1:length(gseList) ) {
    gseTmp = getGEO( gseList[1] )
    if ( length(gseTmp) > 1 ) {
      index = grep( gplList[1], attr( gseTmp, "names" ) )
    } else {
      index = 1
    }
    gse = gseTmp[[index]]
  }
}

fvarLabels(gse) = make.names( fvarLabels(gse) )

gsmList = strsplit( gse@experimentData@other$sample_id, ' ' )[[1]]
if ( all( samples %in% gsmList ) ) {
  gsmMatch = which(gsmList %in% samples)
  groupList = rep('X', length(gsmList))
  for ( i in 1:count ) {
    groupList[gsmMatch[i]] = groups[i]
  }
  gsms = paste( groupList, collapse="" )
} else {
  stop( "Samples(GSM) are not in the same series(GSE). Please check your input.", call.=FALSE )
}

groups = groupList[gsmMatch]
gse = gse[,gsmMatch]

ex = exprs(gse)
qx = as.numeric( quantile( ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=TRUE ) )
LogC = (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if ( LogC ) {
  ex[which(ex <= 0)] <- NaN
  exprs(gse) <- log2(ex) 
}

groups = paste( "G", groups, sep="" )
fl = as.factor( groups )

gse$description = fl
design = model.matrix( ~ description + 0, gse )
colnames( design ) = levels( fl )
fit = lmFit( gse, design )

if ( length( unique(groups) ) == 2 ) {
  contrast = c(paste( unique(groups), collapse="-" ) )
  
  cont.matrix = makeContrasts( contrasts=contrast, levels=design )
  fit2 = contrasts.fit( fit, cont.matrix )
  fit2 = eBayes( fit2, 0.01 )
  tT = topTable( fit2, adjust="fdr", sort.by="B", number=Inf )
  tT = tT[tT[,"P.Value"]<0.05 & abs(tT[,"logFC"])>1,]
  tT = subset( tT, select=c( "ID", "Gene.Symbol", "Gene.Title", "logFC", "P.Value", "adj.P.Val" ) )
  print( tT )
} else {
  contrast = c()
  contrastNum = 0
  groupNum = length(unique(groups))
  for ( i in 1:(groupNum-1) ) {
    for ( j in (i+1):groupNum ) {
      contrastNum = contrastNum + 1
      contrast[contrastNum] = paste( unique(groups)[i], unique(groups)[j], sep="-" )
    }
  }
  cont.matrix = makeContrasts( contrasts = contrast, levels=design )
  fit2 = contrasts.fit( fit, cont.matrix )
  fit2 = eBayes(fit2, 0.01)
  tT = topTable( fit2, adjust="fdr", sort.by="B", number=Inf )
  tT = tT[tT[,"P.Value"]<0.05,]
  tT <- subset(tT, select=c( "ID", "Gene.Symbol", "Gene.Title", "F", "P.Value", "adj.P.Val" ) )
  print( tT )
}