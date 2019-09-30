# File: 14_gseaSpatialAllContrasts.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: gene set enrichment analysis for the datasets
# Date: 30/09/2019


## set variables and source libraries
## libraries to load
library(gage)

lFiles = list.files('results/', pattern='DEAnalysis*', full.names = T, ignore.case = T)

ldfData = lapply(lFiles, function(x) as.data.frame(read.csv(x, header=T, row.names=1, stringsAsFactors = F)))
names(ldfData) = lFiles
sapply(ldfData, nrow)

# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

sapply(ldfData, function(df) identical(rownames(df), rn))


cvTitle = gsub('results//DEAnalysis', '', names(ldfData))
cvTitle = gsub('.xls', '', cvTitle)

## load spatial coordinates data
oMsigGS.c2 = readList('~/Data/MetaData/mm10GeneChromosome.gmt')

## choose a contrast to work with loop through
for (i in 1:length(ldfData)){
  dfContrast = ldfData[[i]]
  # for a contrats of choice create the list
  iContFc = dfContrast$logFC
  ## add enterez ids
  names(iContFc) = as.character(dfContrast$ind)
  head(iContFc)
  head(dfContrast)
  oGage = gage(iContFc, oMsigGS.c2)
  
  dfGreater = data.frame(oGage$greater)
  #str(dfGreater)
  #i = which(dfGreater$p.val < 0.01)
  #rownames(dfGreater[i,])
  
  dfLess = data.frame(oGage$less)
  #str(dfLess)
  #i = which(dfLess$p.val < 0.01)
  #rownames(dfLess[i,])
  
  write.csv(dfGreater[,c('p.val', 'q.val', 'set.size')], file=paste('results/', cvTitle[i], '_upregulated_spatial.xls', sep=''))
  write.csv(dfLess[,c('p.val', 'q.val', 'set.size')], file=paste('results/', cvTitle[i], '_downregulated_spatial.xls', sep=''))
}
