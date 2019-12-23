# File: 18_bam2ExpressedSequence.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: use the bam file to load a gene of choice and extract reads in that region
#       and align those reads with the reference sequence for the gene
# Date: 20/12/2019

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)

# enterez id for the gene of interest
oGRgene = genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
oGRgene = oGRgene['22761']

library(GenomicAlignments)
bam = readGAlignments(file.choose(), param=ScanBamParam(which=oGRgene, what='seq'))

# sequenceLayer
# s = mcols(bam[1])$seq
# sequenceLayer(s, cigar(bam)[1])
s2 = stackStringsFromBam(file.choose(), param=ScanBamParam(which=oGRgene, what='seq'))

oPile = pileup(file.choose(), scanBamParam = ScanBamParam(which = oGRgene), 
               pileupParam = PileupParam(distinguish_strands = F, max_depth = 1000, min_base_quality = 30, 
                                         min_mapq = 30))

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

getSequenceParameters = function(ivSeq, cRefBase, prior=c(A=1/2, T=1/2, G=1/2, C=1/2), iSize=1000){
  if(!require(LearnBayes)) stop('R Package LearnBayes required')
  ## internal functions
  # get alpha values for dirichlet posterior
  getAlpha = function(seq, prior=c(A=1/2, T=1/2, G=1/2, C=1/2)){
    #a = letterFrequency(seq, letters='ATGC', OR=0, as.prob=F)
    alpha = seq + prior #a + prior
    return(alpha)
  }
  # calculate dirichlet variance 
  getDirichletVariance = function(alpha){
    al0 = sum(alpha)
    denom = al0^2 * (al0+1)
    ret = sapply(alpha, function(x){
      return(x * (al0 - x) / denom)
    })
    return(ret)
  }
  # get posterior theta from posterior dirichlet
  getPosterior = function(alpha, n=1000){
    p = rdirichlet(n, alpha)
    colnames(p) = names(alpha)
    #m = colMeans(p)
    return(p)
  }
  # posterior preditive distribution based on the theta from dirichlet posterior
  # set n to adjust size of sum alpha
  getPosteriorPredict = function(theta, n=1){
    ret = t(apply(theta, 1, function(x) rmultinom(1, n, x)))
    colnames(ret) = colnames(theta)
    return(ret)
  }
  # posterior gamma, a component of the dirichlet distribution
  ## see gelman P 583 and bayesian computations with R page 66
  getPosteriorGamma = function(alpha.scale, base, n=1000, prior){
    # adjust alpha by removing dirichlet prior
    alpha.scale = alpha.scale - prior
    i = which(names(alpha.scale) == base)
    alpha.new = c(alpha.scale[i], sum(alpha.scale[-i]))
    names(alpha.new) = c('Base', 'Other')
    ## set a non-informative jeffery's prior for gamma distributed rate
    jef.prior = c(alpha=0.5, beta=1)
    # adjust the new alpha for gamma posterior
    alpha.new = alpha.new + jef.prior['alpha']
    # convert the alpha to rate per 1000
    alpha.new = (alpha.new/sum(alpha.new)) * 1000
    rg = sapply(seq_along(alpha.new), function(x) {
      return(rgamma(n, alpha.new[x], jef.prior['beta']))
    })
    colnames(rg) = names(alpha.new)
    return(rg)
  }
  ###### processing steps
  ## get posterior values
  a = getAlpha(ivSeq, prior)
  # get posterior dirichlet variance
  var = getDirichletVariance(a)
  # get posterior via simulation
  p = getPosterior(a, iSize)
  #   r = getPosteriorPredict(colMeans(p), 1000)
  #   # get variance by adding variance of each binomial component of the posterior predictive data
  #   (sum(apply(r, 2, var)))
  var = sum(var)
  theta = colMeans(p)[cRefBase]
  if (is.na(theta)) {rate=c('Base'= NA, 'Other'=NA) } else {
    rate = round(colMeans(getPosteriorGamma(a, cRefBase, iSize, prior)), 2)}
  return(c(theta=theta, var=var, lambda.base=rate['Base'], lambda.other=rate['Other'], ivSeq))
}

s = f_getSeq(oPile)

oDSRef = getSeq(BSgenome.Mmusculus.UCSC.mm10, oGRgene)
param = sapply(rownames(s), function(x){
  gr = GRanges(seqnames(oGRgene), ranges = IRanges(start = as.numeric(x), width=1), strand=strand(oGRgene))
  ref = getSeq(BSgenome.Mmusculus.UCSC.mm10, gr)
  getSequenceParameters(s[x,], as.character(ref[[1]]))
})

colnames(param) = rownames(s)
rownames(param) = c('theta', 'var', 'lambda.base', 'lambda.other', 'A', 'T', 'G', 'C')
param = t(param)
dim(param)
# create factors for quantiles
param = na.omit(param)
dim(param)
# var.q = cut(param[,'var'],breaks = quantile(param[,'var'], c(0, 0.05, 0.95, 1)), include.lowest = T, 
#             labels = c('q.low', 'q,m', 'q.high'))
theta.q = cut(param[,'theta'],breaks = quantile(param[,'theta'], c(0, 0.01, 0.95, 1)), include.lowest = T, 
              labels = c('q.low', 'q.m', 'q.high'))
param = cbind(param, theta.q)

## reformat the data matrix 
mFile = matrix(NA, nrow=width(oDSRef), ncol=ncol(param), dimnames=list(start(oGRgene):end(oGRgene), colnames(param)))
m = match(rownames(param), rownames(mFile))
mFile[m,] = param
cvConsensus = apply(mFile[,c('A', 'T', 'G', 'C')], 1, function(x){
  return(c('A', 'T', 'G', 'C')[which.max(x)])
})
cvConsensus = unlist(cvConsensus)
dfResult = data.frame(mFile)
dfResult$consensus = NA
m = match(names(cvConsensus), rownames(mFile))
dfResult[m,'consensus'] = cvConsensus

pat = paste0(cvConsensus, collapse = '')
pwa = pairwiseAlignment(pat, oDSRef, type='global')

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

c = unlist(strsplit(toString(oDSRef[[1]]), split = ''))
m = rbind(c, as.matrix(pwa))
rownames(m) = c('Reference', 'Consensus')
oDSexport = f_oDNAStringSetConvertPWAMatrix(m)
export(oDSexport, 'temp/alignment.fasta')

m = rbind(c, dfResult$consensus)
