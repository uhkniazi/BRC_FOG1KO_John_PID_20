# File: bioType.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: assign DE genes to their biotypes
# Date: 04/11/2019

source('header.R')


## load the common binary matrix of DE genes created in earlier results
dfCommonGenes = read.csv('results/commonDEGenes.xls', header=T, row.names=1)
head(dfCommonGenes)

library(org.Mm.eg.db)

dfEns = select(org.Mm.eg.db, keys=rownames(dfCommonGenes), keytype = 'ENTREZID', columns = 'ENSEMBL')
dfEns = na.omit(dfEns)
rn = rownames(dfCommonGenes)
rn = rn[rn %in% dfEns$ENTREZID]

dfEns = select(org.Mm.eg.db, keys=rn, keytype = 'ENTREZID', columns = 'ENSEMBL')
dfEns = na.omit(dfEns)
#dfEns = dfEns[!duplicated(dfEns$ENTREZID), ]
i = match(rn, dfEns$ENTREZID)
dfEns = dfEns[i,]
identical(dfEns$ENTREZID, rn)

library(EnsDb.Mmusculus.v79)

columns(EnsDb.Mmusculus.v79)

## total types of biotypes
listTxbiotypes(EnsDb.Mmusculus.v79)
df = select(EnsDb.Mmusculus.v79, keys = rn, keytype = 'ENTREZID', columns = 'TXBIOTYPE')
as.data.frame(table(df$TXBIOTYPE))
barplot(table(df$TXBIOTYPE), main='de expressed genes biotype', cex.names=0.7, las=1)

### all count matrix
library(RMySQL)

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# get the query
g_did
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 41) AND (MetaFile.comment like "%count%")')
dfSample = dbGetQuery(db, q)
dfSample
n = paste0(dfSample$location, dfSample$name)
load(n[2])

## load the metadata i.e. covariates
q = paste0('select Sample.* from Sample where Sample.idData = 41')
dfSample = dbGetQuery(db, q)
dim(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)

## make count matrix
names(lCounts)
mCounts = do.call(cbind, lCounts)
colnames(mCounts) = names(lCounts)

# sanity check
identical(dfSample$id, as.integer(colnames(mCounts)))

mData = mCounts
dim(mData)

f = ifelse(dfSample$group1 == 'Differentiated', 'D', 'ND')
fReplicates = paste0(f, '-', dfSample$group3)
table(fReplicates)
dfSample$fReplicates = factor(fReplicates)
# combine the technical replicates
i = seq_along(1:ncol(mData))
m = tapply(i, dfSample$fReplicates, function(x) {
  return(x)
})

mData = sapply(m, function(x){
  return(rowSums(mCounts[,x]))
})

# get a shorter version of dfSample after adding technical replicates
dfSample.2 = dfSample[sapply(m, function(x) return(x[1])), ]
identical(colnames(mData), as.character(dfSample.2$fReplicates))

## first normalise the data
# drop the samples where average across rows is less than 3
i = rowMeans(mData)
table( i < 3)
# FALSE  TRUE 
# 13658 10809 
mData = mData[!(i< 3),]
dim(mData)

df = select(EnsDb.Mmusculus.v79, keys = rownames(mData), keytype = 'ENTREZID', columns = 'TXBIOTYPE')
barplot(table(df$TXBIOTYPE), main='all expressed genes biotype', cex.names=0.7, las=1)
