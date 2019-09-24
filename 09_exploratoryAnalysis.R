# File: 09_exploratoryAnalysis.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: quality checks on the count matrix with covariates
# Date: 24/09/2019

source('header.R')

## load the data
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
## do both data sets from mcgill and bam2fq route, choose the right index
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

## load CDiagnostics and test
## compare the normalised and raw data
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

# perform both data sets together, mcgill and bam2fq route
oDiag.1 = CDiagnosticPlots(log(mData.fq+0.5), 'bam2fq')
oDiag.2 = CDiagnosticPlots(log(mData.mc+0.5), 'mcgill')
# the batch variable we wish to colour by, 
# this can be any grouping/clustering in the data capture process
# e.g. in this case it is different lanes/machines
dfSample.2$group1 = ifelse(dfSample.2$group1 == 'Differentiated', 'D', 'ND')
str(dfSample.2)
fBatch = factor(dfSample.2$group1)
fBatch = factor(dfSample.2$group3)
fBatch = factor(dfSample.2$group3):factor(dfSample.2$group1)

## compare the 2 methods using various plots
par(mfrow=c(1,2))
boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.5)
boxplot.median.summary(oDiag.2, fBatch, legend.pos = 'topright', axis.label.cex = 0.5)

plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.5)
plot.mean.summary(oDiag.2, fBatch, axis.label.cex = 0.5)

plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.5)
plot.sigma.summary(oDiag.2, fBatch, axis.label.cex = 0.5)

plot.missing.summary(oDiag.1, fBatch, axis.label.cex = 0.5, cex.main=1)
plot.missing.summary(oDiag.2, fBatch, axis.label.cex = 0.5, cex.main=1)

# plot.PCA(oDiag.1, fBatch, cex.main=1)
# plot.PCA(oDiag.2, fBatch, cex.main=1)
# # 
# plot.dendogram(oDiag.1, fBatch, labels_cex = 0.8, cex.main=0.7)
# plot.dendogram(oDiag.2, fBatch, labels_cex = 0.8, cex.main=0.7)

## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.2)
l$PCA.jitter = F
l$HC.jitter = F

oDiag.1 = CDiagnosticPlotsSetParameters(oDiag.1, l)
oDiag.2 = CDiagnosticPlotsSetParameters(oDiag.2, l)
plot.PCA(oDiag.1, fBatch, legend.pos = 'topleft')
plot.PCA(oDiag.2, fBatch, legend.pos='topleft')
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.7)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.7)

# ## extreme values
# oDiag.1 = calculateExtremeValues(oDiag.1)
# oDiag.2 = calculateExtremeValues(oDiag.2)
# m1 = mGetExtremeValues(oDiag.1)
# m2 = mGetExtremeValues(oDiag.2)

# ## samples with most extreme values
# apply(m1, 2, function(x) sum(x > 0))
# apply(m2, 2, function(x) sum(x > 0))
# 
# ## variables that are contributing to this
# v1 = apply(m1, 1, function(x) sum(x > 0))
# v2 = apply(m2, 1, function(x) sum(x > 0))
# 
# which(v1 > 0)
# which(v2 > 0)
par(mfrow=c(1,1))
plot(oDiag.2@lData$PCA$sdev)
plot.PCA(oDiag.2, fBatch)
mPC = oDiag.2@lData$PCA$x[,1:4]
mPC = scale(mPC)
## try a linear mixed effect model to account for varince
library(lme4)
dfData = data.frame(mPC)
dfData = stack(dfData)
str(dfData)
library(lattice)
densityplot(~ values, data=dfData)
densityplot(~ values | ind, data=dfData)

str(dfSample.2)
dfData$fTreatment = factor(dfSample.2$group1)
dfData$fGroup = factor(dfSample.2$group3)
dfData$fBatch = dfData$fTreatment:dfData$fGroup

densityplot(~ values | ind, groups=fBatch, data=dfData, auto.key = list(columns=3))

# format data for modelling
dfData$Coef.1 = factor(dfData$fTreatment:dfData$ind)
dfData$Coef.2 = factor(dfData$fGroup:dfData$ind)
dfData$Coef.int = factor(dfData$Coef.1:dfData$Coef.2)

str(dfData)

fit.lme1 = lmer(values ~ 1  + (1 | Coef.1), data=dfData)
fit.lme2 = lmer(values ~ 1  + (1 | Coef.1) + (1 | Coef.2), data=dfData)
#fit.lme3 = lmer(values ~ 1  + (1 | Coef.1) + (1 | Coef.2) + (1 | Coef.int), data=dfData)
#fit.lme4 = lmer(values ~ 1  + (1 | Coef.int), data=dfData)

anova(fit.lme1, fit.lme2)

summary(fit.lme1)
summary(fit.lme2)

plot((fitted(fit.lme1)), resid(fit.lme1), pch=20, cex=0.7)
lines(lowess((fitted(fit.lme1)), resid(fit.lme1)), col=2)
hist(dfData$values, prob=T)
lines(density(fitted(fit.lme1)))

## fit model with stan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='tResponsePartialPooling.stan')

str(dfData)
m1 = model.matrix(values ~ Coef.1 - 1, data=dfData)
m2 = model.matrix(values ~ Coef.2 - 1, data=dfData)
m3 = model.matrix(values ~ Coef.int - 1, data=dfData)
m = cbind(m1, m2)#, m3)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
                 NscaleBatches=2, NBatchMap=c(rep(1, times=nlevels(dfData$Coef.1)),
                                              rep(2, times=nlevels(dfData$Coef.2))),
                 #rep(3, times=nlevels(dfData$Coef.int))),
                 y=dfData$values)

fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=2, pars=c('betas', 'populationMean', 'sigmaPop', 'sigmaRan',
                                                                         'nu', 'mu'),
                    cores=2)#, control=list(adapt_delta=0.99, max_treedepth = 12))
print(fit.stan, c('betas', 'populationMean', 'sigmaPop', 'sigmaRan', 'nu'), digits=3)

traceplot(fit.stan, 'populationMean')
traceplot(fit.stan, 'sigmaPop')
traceplot(fit.stan, 'sigmaRan')

# ## just using the one coefficient
# m = model.matrix(values ~ Coef.int - 1, data=dfData)
# 
# lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
#                  NscaleBatches=1, NBatchMap=c(rep(1, times=nlevels(dfData$Coef.int))),
#                  y=dfData$values)
# 
# fit.stan.2 = sampling(stanDso, data=lStanData, iter=5000, chains=2, pars=c('betas', 'populationMean', 'sigmaPop', 'sigmaRan',
#                                                                            'nu', 'mu'),
#                       cores=2)
# print(fit.stan.2, c('betas', 'populationMean', 'sigmaPop', 'sigmaRan', 'nu'), digits=3)
# 
# traceplot(fit.stan.2, 'populationMean')
# traceplot(fit.stan.2, 'sigmaPop')
# traceplot(fit.stan.2, 'sigmaRan')

m = cbind(extract(fit.stan)$sigmaRan, extract(fit.stan)$sigmaPop) 
dim(m)
m = log(m)
colnames(m) = c('Treatment', 'Genotype', 'Residual')
pairs(m, pch=20, cex=0.5, col='grey')

df = stack(data.frame(m[,-3]))
histogram(~ values | ind, data=df, xlab='Log SD')


hist(dfData$values, prob=T)
plot(density(dfData$values))
## use appropriate model name e.g. fit.stan or fit.stan.2 to do model checks, in this case
# fit.stan.2 is a better model.
mFitted = extract(fit.stan)$mu

apply(mFitted[sample(1:nrow(mFitted), size = 100), ], 1, function(x) lines(density(x)))

# quick residual check
m = colMeans(mFitted)
r = dfData$values - m
plot(m, r, pch=20)
lines(lowess(m, r), col=2)

### plot the posterior predictive values
m = extract(fit.stan, c('mu', 'nu', 'sigmaPop'))
dim(m$mu)
i = sample(1:5000, 1000)
muSample = m$mu[i,]
nuSample = m$nu[i]
sigSample = m$sigmaPop[i]

## t sampling functions
dt_ls = function(x, df, mu, a) 1/a * dt((x - mu)/a, df)
rt_ls <- function(n, df, mu, a) rt(n,df)*a + mu

## use point-wise predictive approach to sample a new value from this data
ivResp = dfData$values
mDraws = matrix(NA, nrow = length(ivResp), ncol=1000)

for (i in 1:ncol(mDraws)){
  mDraws[,i] = rt_ls(length(ivResp), nuSample[i], muSample[i,], sigSample[i])
}

yresp = density(ivResp)
plot(yresp, xlab='', main='Fitted distribution', ylab='density', lwd=2)#, ylim=c(0, 1))
temp = apply(mDraws[,1:30], 2, function(x) {x = density(x)
#x$y = x$y/max(x$y)
lines(x, col='red', lwd=0.2)
})

## calculate bayesian p-value for this test statistic
getPValue = function(Trep, Tobs){
  left = sum(Trep <= Tobs)/length(Trep)
  right = sum(Trep >= Tobs)/length(Trep)
  return(min(left, right))
}
## define some test quantities to measure the lack of fit
## define a test quantity T(y, theta)
## variance
T1_var = function(Y) return(var(Y))

## min quantity
T1_min = function(Y){
  return(min(Y))
} 

## max quantity
T1_max = function(Y){
  return(max(Y))
} 

## mean quantity
T1_mean = function(Y){
  return(mean(Y))
} 

## mChecks
mChecks = matrix(NA, nrow=4, ncol=1)
rownames(mChecks) = c('Variance', 'Max', 'Min', 'Mean')
colnames(mChecks) = c('model 1')

t1 = apply(mDraws, 2, T1_var)
mChecks['Variance', 1] = getPValue(t1, var(ivResp))

## testing for outlier detection i.e. the minimum value show in the histograms earlier
t1 = apply(mDraws, 2, T1_min)
t2 = T1_min(ivResp)
mChecks['Min',1] = getPValue(t1, t2)

## maximum value
t1 = apply(mDraws, 2, T1_max)
t2 = T1_max(ivResp)
mChecks['Max', 1] = getPValue(t1, t2)

## mean value
t1 = apply(mDraws, 2, T1_mean)
t2 = T1_mean(ivResp)
mChecks['Mean', 1] = getPValue(t1, t2)

mChecks
