# File: 18_bam2ExpressedSequence.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: use the bam file to load a gene of choice and extract reads in that region
#       and align those reads with the reference sequence for the gene
# Date: 20/12/2019

source('header.R')
library('RMySQL')
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)
library(GenomicAlignments)

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'File')
dbListFields(db, 'Sample')
# get the query
g_did
q = paste0('select Sample.id as sid, Sample.group1, Sample.group2, Sample.title, File.* from Sample, File
           where (Sample.idData = 41) AND (File.idSample = Sample.id AND File.type like "mcgill quality 10 sorted bam")')
dfSample = dbGetQuery(db, q)
nrow(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)
# remove any whitespace from the names
dfSample$name = gsub(" ", "", dfSample$name, fixed = T)
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)
i = grep('q10_sort.bam', dfSample$name)
dfSample = dfSample[i,]
#### set working directory to appropriate location with bam files
setwd('dataExternal/bams/')
csFiles = list.files('.', pattern = '*.bam$', recursive = T)
# check if these files match the file names in database
table(dfSample$name %in% csFiles)

dfSample = dfSample[dfSample$name %in% csFiles, ]
dim(dfSample)
csFiles = dfSample$name

# enterez id for the gene of interest
# reset format to NCBI from UCSC
tx = TxDb.Mmusculus.UCSC.mm10.knownGene
seqlevelsStyle(tx)
seqlevelsStyle(tx) = 'NCBI'
oGRgene = genes(tx)
oGRgene = oGRgene['22761']
oGRgene = GRanges(as.character(seqnames(oGRgene)), ranges(oGRgene), as.character(strand(oGRgene)))

bs = BSgenome.Mmusculus.UCSC.mm10
seqlevelsStyle(bs) = 'NCBI'
oDSRef = getSeq(bs, oGRgene)

f_getSeq = function(pile){
  # split the pile up data frame on positions
  df = data.frame(P=pile$pos, N=pile$nucleotide, C=pile$count)
  l = split(df, df$P)
  seq = array(0, dim=c(length(l), 4), dimnames = list(names(l), c('A', 'T', 'G', 'C')))
  for (i in 1:nrow(seq)){
    s2 = seq[i,]
    # sequence at current nucleotide position
    x = l[[i]]
    s = x$C
    names(s) = x$N
    # match the positions of names
    p = match(names(s), names(s2))
    # if a non standard residue is present then the position will not match
    # and will have an NA there
    p2 = !is.na(p)
    p = p[p2]
    s2[p] = s[p2]
    seq[i,] = s2
  }
  return(seq)
}

# NAME: f_oDNAStringSetConvertPWAMatrix
# ARGS: pwa.matrix = a pairwise alignment object in matrix form, generated via a call like
#       as.matrix(pairwiseAlignment(patten, subject))
# DESC: converts the sequences in the matrix to a DNAStringSet object which can be exported to a FASTA
#       file to be viewed in a viewer of choice e.g. jalview
# RETS: a DNAStringSet object with all sequences from the pwa matrix
f_oDNAStringSetConvertPWAMatrix= function(pwa.matrix){
  require(Biostrings)
  # sequence to return
  seq.export = DNAStringSet()
  # extract each sequence from the matrix
  for (i in 1:nrow(pwa.matrix)){
    s = DNAStringSet(pwa.matrix[i,])
    s = DNAStringSet(unlist(s))
    seq.export = append(seq.export, s)
  }
  # set names for sequences
  names(seq.export) = rownames(pwa.matrix)
  return(seq.export)
}

oPile = lapply(csFiles, function(x) pileup(x, scanBamParam = ScanBamParam(which = oGRgene), 
               pileupParam = PileupParam(distinguish_strands = F, max_depth = 1000, min_base_quality = 30, 
                                         min_mapq = 30)))

lParam = lapply(oPile, f_getSeq)

lConsensus = lapply(lParam, function(param){
  cvConsensus = apply(param[,c('A', 'T', 'G', 'C')], 1, function(x){
    return(c('A', 'T', 'G', 'C')[which.max(x)])
  })
  names(cvConsensus) = rownames(param)
  return(cvConsensus)
})

names(lConsensus) = paste(dfSample$idSample, dfSample$title, sep = '_')

#param = data.frame(param, consensus=cvConsensus, stringsAsFactors = F)
## reformat the data matrix 
dfResult = data.frame(reference=unlist(strsplit(toString(oDSRef[[1]]), split = '')),
                      row.names=c(start(oGRgene):end(oGRgene)), stringsAsFactors = F)

for (i in 1:length(lConsensus)){
  m = match(names(lConsensus[[i]]), rownames(dfResult))
  dfResult = cbind(dfResult, '-', stringsAsFactors=F)
  dfResult[m, i+1] = lConsensus[[i]]
  colnames(dfResult)[i+1] = names(lConsensus)[i]
}

## create appropriate file name
cvExport = paste0('../../results/alignments/', 'all.fasta')
oDSalignment = f_oDNAStringSetConvertPWAMatrix(t(dfResult))
export(oDSalignment, cvExport)


