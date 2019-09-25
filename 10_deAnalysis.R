# File: 10_deAnalysis.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: Modelling and selecting DE genes
# Date: 24/09/2019

## load the data
source('header.R')
library(RMySQL)

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# select the right table using data and project id
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

ivProb = apply(mData, 1, function(inData) {
  inData[is.na(inData) | !is.finite(inData)] = 0
  inData = as.logical(inData)
  lData = list('success'=sum(inData), fail=sum(!inData))
  return(mean(rbeta(1000, lData$success + 0.5, lData$fail + 0.5)))
})

hist(ivProb)

library(DESeq2)
sf = estimateSizeFactorsForMatrix(mData)
mData.norm = sweep(mData, 2, sf, '/')

identical(colnames(mData.norm), as.character(dfSample.2$fReplicates))

## delete sample section after testing
mData.norm = round(mData.norm, 0)

set.seed(123)
i = sample(1:nrow(mData.norm), 10, replace = F)
dfData = data.frame(t(mData.norm[i,]))

#dfData = data.frame(t(mData.norm))
dfData = stack(dfData)
dfSample.2$group1 = ifelse(dfSample.2$group1 == 'Differentiated', 'D', 'ND')
str(dfSample.2)

dfData$fTreatment = factor(dfSample.2$group1, levels = c('ND', 'D'))
dfData$fGenotype = factor(dfSample.2$group3)
dfData$fBatch = dfData$fTreatment:dfData$fGenotype

dfData$Coef = factor(dfData$fBatch:dfData$ind)

dfData = droplevels.data.frame(dfData)
dfData = dfData[order(dfData$Coef), ]
str(dfData)

# # setup the model
library(lme4)
fit.lme1 = glmer.nb(values ~ 1 + (1 | Coef), data=dfData)
summary(fit.lme1)
ran = ranef(fit.lme1, condVar=F)

plot(log(fitted(fit.lme1)), resid(fit.lme1), pch=20, cex=0.7)
lines(lowess(log(fitted(fit.lme1)), resid(fit.lme1)), col=2)

## setup the stan model
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='nbinomResp1RandomEffectsMultipleScales.stan')

## calculate hyperparameters for variance of coefficients
l = gammaShRaFromModeSD(sd(log(dfData$values+0.5)), 2*sd(log(dfData$values+0.5)))
# # ## set initial values
# ran = ranef(fit.lme1)
# r1 = ran$Coef
# r2 = ran$Coef.adj1
# r3 = ran$Coef.adj2
# 
# initf = function(chain_id = 1) {
#   list(sigmaRan1 = 1, sigmaRan2=1)
# }

## subset the data to get the second level of nested parameters
## this is done to avoid loops in the stan script to map the scale parameters
## of each ind/gene to the respective set of coefficients for jitters
d = dfData[!duplicated(dfData$Coef), ]

lStanData = list(Ntotal=nrow(dfData), 
                 Nclusters1=nlevels(dfData$Coef),
                 NScaleBatches1 = nlevels(dfData$ind), # to add a separate scale term for each gene
                 NgroupMap1=as.numeric(dfData$Coef),
                 NBatchMap1=as.numeric(d$ind), # this is where we use the second level mapping
                 Nphi=nlevels(dfData$ind),
                 NphiMap=as.numeric(dfData$ind),
                 y=dfData$values, 
                 gammaShape=l$shape, gammaRate=l$rate,
                 intercept = mean(log(dfData$values+0.5)), intercept_sd= sd(log(dfData$values+0.5))*3)

ptm = proc.time()

fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=4,
                    pars=c('sigmaRan1',
                           'phi',
                           #'mu',
                           'rGroupsJitter1'
                           #'betas',
                           #'phi_scaled'
                           ),
                    cores=4, control=list(adapt_delta=0.99, max_treedepth = 11))#, init=initf)
# save(fit.stan, file='temp/fit.stan.nb_13Mar.rds')
ptm.end = proc.time()
print(fit.stan, c('sigmaRan1[1]', 'phi'), digits=3)
print(fit.stan, c('rGroupsJitter1'))
#traceplot(fit.stan, 'betas')
traceplot(fit.stan, c('sigmaRan1[1]'))
traceplot(fit.stan, c('sigmaRan1[2]'))
traceplot(fit.stan, c('rGroupsJitter1[1]', 'sigmaRan1[1]'))

## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$rGroupsJitter1
dim(mCoef)
# # ## get the intercept at population level
# iIntercept = as.numeric(extract(fit.stan)$betas)
# ## add the intercept to each random effect variable, to get the full coefficient
# mCoef = sweep(mCoef, 1, iIntercept, '+')

## function to calculate statistics for differences between coefficients
getDifference = function(ivData, ivBaseline){
  stopifnot(length(ivData) == length(ivBaseline))
  # get the difference vector
  d = ivData - ivBaseline
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(list(z=z, p=p))
}

## split the data into the comparisons required
d = data.frame(cols=1:ncol(mCoef), mods=levels(dfData$Coef))
# the split is done below on : symbol, but factor name has a : symbol due
# to creation of interaction earlier, do some acrobatics to sort that issue
## split this factor into sub factors
f = strsplit(as.character(d$mods), ':')
d = cbind(d, do.call(rbind, f))
head(d)
d$`1` = d$`1`:d$`2`
d = d[,-4]
colnames(d) = c(colnames(d)[1:2], c('fBatch', 'ind'))
str(d)
d$split = factor(d$ind)

levels(d$fBatch)
## repeat this for each comparison

## get a p-value for each comparison
l = tapply(d$cols, d$split, FUN = function(x, base='ND:WT', deflection='D:WT') {
  c = x
  names(c) = as.character(d$fBatch[c])
  dif = getDifference(ivData = mCoef[,c[deflection]], ivBaseline = mCoef[,c[base]])
  r = data.frame(ind= as.character(d$ind[c[base]]), coef.base=mean(mCoef[,c[base]]), 
                 coef.deflection=mean(mCoef[,c[deflection]]), zscore=dif$z, pvalue=dif$p)
  r$difference = r$coef.deflection - r$coef.base
  #return(format(r, digi=3))
  return(r)
})

dfResults = do.call(rbind, l)
dfResults$adj.P.Val = p.adjust(dfResults$pvalue, method='BH')

### plot the results
dfResults$logFC = dfResults$difference
dfResults$P.Value = dfResults$pvalue
library(org.Mm.eg.db)
## remove X from annotation names
dfResults$ind = gsub('X', '', as.character(dfResults$ind))
df = AnnotationDbi::select(org.Mm.eg.db, keys = as.character(dfResults$ind), columns = 'SYMBOL', keytype = 'ENTREZID')
i = match(dfResults$ind, df$ENTREZID)
df = df[i,]
dfResults$SYMBOL = df$SYMBOL
identical(dfResults$ind, df$ENTREZID)
## produce the plots 
f_plotVolcano(dfResults, 'WT:15.5 vs WT:13.5')#, fc.lim=c(-2.5, 2.5))
f_plotVolcano(dfResults, 'WT:15.5 vs WT:13.5', fc.lim=range(dfResults$logFC))

m = tapply(dfData$values, dfData$ind, mean)
i = match(rownames(dfResults), names(m))
m = m[i]
identical(names(m), rownames(dfResults))
plotMeanFC(log(m), dfResults, 0.01, 'WT:14.5 vs WT:13.5')
table(dfResults$adj.P.Val < 0.01)
## save the results 
write.csv(dfResults, file='results/DEAnalysisMut:15.5VsMut:14.5.xls')

######### do a comparison with deseq2
dfDesign = data.frame(Treatment = factor(dfSample.2$group1, levels=c('WT', 'Mut')):factor(dfSample.2$group3)
                      , row.names=colnames(mData))

oDseq = DESeqDataSetFromMatrix(mData, dfDesign, design = ~ Treatment)
oDseq = DESeq(oDseq)
plotDispEsts(oDseq)
oRes = results(oDseq, contrast=c('Treatment', 'Mut:15.5', 'Mut:14.5'))
plotMA(oRes)
temp = as.data.frame(oRes)
i = match((dfResults$ind), rownames(temp))
temp = temp[i,]
identical((dfResults$ind), rownames(temp))
plot(dfResults$logFC, temp$log2FoldChange, pch=20)
table(oRes$padj < 0.01)