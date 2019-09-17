# File: 01_createDbEntry.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: list the samples and create appropriate db entries
# Date: 16/9/2019


## set variables and source libraries
source('header.R')

## connect to mysql database 
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# sample and file table
dbGetQuery(db, paste('describe Sample;'))
cSampleCol = dbGetQuery(db, paste('describe Sample;'))$Field[-1]

dbGetQuery(db, paste('describe File;'))
cFileCol = dbGetQuery(db, paste('describe File;'))$Field[-1]

# setwd(gcRemoteDir)
setwd('dataExternal/')
setwd('data/raw/')

# list the files
cvFiles = list.files(pattern = 'fastq.gz', recursive = T)

# each sample has 2 files 
fSplit = gsub('_[1|2].fastq.gz', '', cvFiles)

lFiles = split(cvFiles, fSplit)

setwd(gcswd)
## load the metadata file
dfMeta = read.csv('dataExternal/metaData.csv', header=T, stringsAsFactors = F)
str(dfMeta)
dfMeta$FileName = gsub('./RNASeq_\\d+/', '', dfMeta$FileName)

# sanity check
table(as.character(dfMeta$FileName) %in% unique(fSplit))

## order the table in the same sequence as file names
i = match(names(lFiles), as.character(dfMeta$FileName))
dfMeta = dfMeta[i,]
identical(as.character(dfMeta$FileName), names(lFiles))
identical(as.character(dfMeta$FileName), unique(fSplit))

## create the entry for samples
cSampleCol

dfSamples = data.frame(idProject=g_pid, idData=g_did, title=dfMeta$SampleLabel, 
                       description= paste('sample_name', as.character(dfMeta$FileName),
                                          'genotype', as.character(dfMeta$Genotype),
                                          'group1 is Treatment',
                                          'group2 is technical replicate',
                                          'group3 is clone', sep=';'),
                       group1 = dfMeta$Differentiated, group2= dfMeta$TechnicalReplicate, group3=dfMeta$Clone)
# write this data to the database
rownames(dfSamples) = NULL

### NOTE: Do not execute this anymore as entry created
# # write this table to database
#dbWriteTable(db, name='Sample', value=dfSamples, append=T, row.names=F)
# get this table again from database with ids added
g_did
dfSamples = dbGetQuery(db, paste('select * from Sample where Sample.idData = 41;'))

# create entries for these files in the database
dbListTables(db)
cn = dbListFields(db, 'File')[-1]
cn
## get the file names from the description column
temp = strsplit(dfSamples$description, ';')
cvNames = sapply(temp, function(x) x[2])
identical(cvNames, names(lFiles))
dfSamples$title = cvNames
# get the names of the samples
temp = lapply(dfSamples$title, function(x){
  # get the file names
  df = data.frame(name=lFiles[[x]], type='fastq', idSample=dfSamples[dfSamples$title == x, 'id'])
  return(df)
})

dfFiles = do.call(rbind, temp)
rownames(dfFiles) = NULL

# write this table to database
## note: do not execute as it is already done
#dbWriteTable(db, name='File', value=dfFiles, append=T, row.names=F)

dbDisconnect(db)
