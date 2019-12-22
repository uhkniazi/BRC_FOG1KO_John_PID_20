# File: 18_bam2ExpressedSequence.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: use the bam file to load a gene of choice and extract reads in that region
#       and align those reads with the reference sequence for the gene
# Date: 20/12/2019

library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# enterez id for the gene of interest
oGRgene = genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
oGRgene = oGRgene['22761']

library(GenomicAlignments)
bam = readGAlignments(file.choose(), param=ScanBamParam(which=oGRgene))
sequenceLayer