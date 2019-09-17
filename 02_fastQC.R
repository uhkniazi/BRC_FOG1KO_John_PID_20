# File: 02_fastQC.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: quality checks on the fastq files before trimming
# Date: 16/09/2019


## set variables and source libraries
source('header.R')
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CFastqQuality/experimental/CFastqQuality.R'
download(url, 'CFastqQuality.R')

# load the required packages
source('CFastqQuality.R')
# delete the file after source
unlink('CFastqQuality.R')

## connect to mysql database to get sample information
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'Sample')
# another way to get the query, preferred
g_did
dfSample = dbGetQuery(db, "select * from Sample where idData=41;")
# remove any whitespace from the names
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)
head(dfSample)
## get the file names from the files table
dbListFields(db, 'File')

## write query to get file names
q = paste0('select File.*, Sample.id as sid, Sample.* from File, Sample where Sample.idData=41 AND Sample.id=File.idSample AND File.type="fastq"', ';')

# df = lapply(q, function(x) dbGetQuery(db, x))
# dfFiles = do.call(rbind, df)
dfFiles = dbGetQuery(db, q)
# close connection after getting data
dbDisconnect(db)
str(dfFiles)
# remove redundant columns
dfFiles = dfFiles[,-c(5,7)]
#### get the names of the fastq files for first sequencing run
setwd('dataExternal/data/raw/')
csFiles = list.files('.', pattern = '*.gz', recursive = T)

# sanity check if all files present in directory
table(dfFiles$name %in% csFiles)

## create factors and variables for fastqqualitybatch object
i = grep('_1.fastq.gz', dfFiles$name)
fReadDirection = rep('2', times=nrow(dfFiles))
fReadDirection[i] = '1'
fReadDirection = factor(fReadDirection)
i = dfFiles$group1
i[i == 'Non-differentiated'] = 'ND'
i[i == "Differentiated"] = 'D'
cNames = paste(i, dfFiles$group3, 'Rep', dfFiles$group2, as.character(fReadDirection),  sep='_')
lMetaData = list(files=dfFiles, samples=dfSample)

ob = CFastqQualityBatch(csFiles, cNames, fReadDirection, lMetaData)


# # split the file names into batches by list, to run analysis on each sample
# lFilesIndex = split(dfFiles$name, dfFiles$idSample)
# names(lFilesIndex) = dfSample$title





## perform the analysis one sample at a time
## function to write the qa files
write.qa = function(fls, indir, title){
  wd = getwd()
  setwd(indir)
  ob = CFastqQuality(fls, title)
  setwd(wd)
  cat(paste('done', title, '\n'))
  return(ob)
}

ivFilesIndex = seq_along(csFiles)

lOb = lapply(ivFilesIndex, function(x){
  tryCatch(write.qa(dfFiles$name[x], getwd(), as.character(dfFiles$id[x])), error=function(e) NULL)
})

names(lOb) = as.character(dfFiles$id)
setwd(gcswd)
n = make.names(paste('CFastqQuality data id 39 gui'))
lOb$meta.1 = dfSample
lOb$meta.2 = dfFiles
lOb$desc = paste('CFastqQuality data id 39 gui', date())
n2 = paste0('~/Data/MetaData/', n)
#save(lOb, file=n2)

## note: comment out as this entry has been made in db
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/', comment='pre trim FASTQ quality checks on Gui mouse data')
dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)

### create the plots of interest
getwd()
lOb$desc = NULL
lOb$meta.1 = NULL
lOb$meta.2 = NULL

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
g_did
q = paste0('select File.id as fid, File.name, File.type, Sample.* from File, Sample where Sample.idData = 39 AND (File.idSample = Sample.id AND File.type like "%fastq%")')
dfBatches = dbGetQuery(db, q)
dbDisconnect(db)
identical(dfFiles$id, dfBatches$fid)
identical(names(lOb), as.character(dfBatches$fid))

# assign sample id (which is unique key) to sample name
cvSampleNames = as.character(dfBatches$id)
cvSampleNames = paste(cvSampleNames, gsub('.+(_[1|2])\\.fastq.gz$', '\\1', dfBatches$name), sep='')
names(lOb) = cvSampleNames

pdf(file='results/qa.fastq.pdf')

iReadCount = sapply(lOb, CFastqQuality.getReadCount)
iReadCount = iReadCount/1e+6

barplot(iReadCount, las=2, main='Pre-Trim Read Count', ylab = 'No. of Reads in Millions', cex.names =0.8, col=grey.colors(2))

mQuality = sapply(lOb, function(x){
  m = mGetReadQualityByCycle(x)
  m = colMeans(m, na.rm = T)
  return(m)
})

matplot(mQuality, type='l', main='Pre-trim base quality', ylab = 'Mean Score', xlab='Position in Read')

lReadWidth = lapply(lOb, iGetReadWidth)
boxplot(lReadWidth, las=2, main='Pre-trim Read Width', ylab = 'Read Width', col=grey.colors(2), outline=F, xaxt='n')
axis(1, at=1:length(lReadWidth), labels = names(lReadWidth), cex.axis=0.8, las=2)

## plot all the alphabets by cycle for each forward and reverse reads
i = grep('_1', cvSampleNames)

lAlphabets = lapply(i, function(x){
  m = t(mGetAlphabetByCycle(lOb[[x]]))
  m = m[,c('A', 'T', 'G', 'C')]
  r = rowSums(m)
  m = sweep(m, 1, r, '/')
  return(m)
})

matplot(lAlphabets[[1]], type='l', main='Pre-trim Sequence Content - Forward Strands', ylab = 'Proportion of Base count', xlab='Position in Read')
temp = lapply(lAlphabets[-1], function(x)
  matlines(x, type='l'))
legend('topleft', legend = colnames(lAlphabets[[1]]), lty=1:4, col=1:4, ncol=2, lwd=2)

# reverse strands
i = grep('_2', cvSampleNames)

lAlphabets = lapply(i, function(x){
  m = t(mGetAlphabetByCycle(lOb[[x]]))
  m = m[,c('A', 'T', 'G', 'C')]
  r = rowSums(m)
  m = sweep(m, 1, r, '/')
  return(m)
})

matplot(lAlphabets[[1]], type='l', main='Pre-trim Sequence Content - Reverse Strands', ylab = 'Proportion of Base count', xlab='Position in Read')
temp = lapply(lAlphabets[-1], function(x)
  matlines(x, type='l'))
legend('topleft', legend = colnames(lAlphabets[[1]]), lty=1:4, col=1:4, ncol=2, lwd=2)

### some diagnostic plots
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

fReadDirection = rep(NA, times=nrow(dfBatches))
i = grepl('_1', cvSampleNames)
fReadDirection[i] = '_1'
fReadDirection[!i] = '_2'
fReadDirection = factor(fReadDirection)

str(dfBatches)
mBatch = mQuality
#colnames(mBatch) = paste0(dfBatches$fid, '-', dfBatches$title, '-', as.character(fReadDirection))

oDiag = CDiagnosticPlots(mBatch, 'Pre-trim Base Quality')
l = CDiagnosticPlotsGetParameters(oDiag)
l$PCA.jitter = F; l$HC.jitter = F;
oDiag = CDiagnosticPlotsSetParameters(oDiag, l)

plot.mean.summary(oDiag, fReadDirection)
plot.sigma.summary(oDiag, fReadDirection)
boxplot.median.summary(oDiag, fReadDirection)
plot.PCA(oDiag, fReadDirection)
plot.dendogram(oDiag, fReadDirection, labels_cex = 0.8)

## try a new batch
fBatch = factor(dfBatches$group3)
levels(fBatch)
plot.mean.summary(oDiag, fBatch)
plot.sigma.summary(oDiag, fBatch)
boxplot.median.summary(oDiag, fBatch)
plot.PCA(oDiag, fBatch)
plot.dendogram(oDiag, fBatch, labels_cex = 0.8)

dev.off(dev.cur())

