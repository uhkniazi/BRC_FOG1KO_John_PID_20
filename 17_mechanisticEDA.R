# File: 17_mechanisticEDA.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: use sets of interacting genes to explore mechanistic models
# Date: 22/11/2019

source('header.R')

###########################################################
############ load the count matrix and normalise
###########################################################
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

###########################################################

###########################################################
########### scale and format data for modelling
###########################################################
dim(mData.norm)
mData.norm = log(mData.norm+0.5)
mData.scaled = t(scale(t(mData.norm)))

cvScan = scan(what=character())

## convert entrez ids to symbols
library(org.Mm.eg.db)
cvGenes.keyword = unique(cvScan)
length(cvGenes.keyword)
df = AnnotationDbi::select(org.Mm.eg.db, keys = cvGenes.keyword, columns = 'SYMBOL', keytype = 'ENTREZID')
df = na.omit(df)
i = match(df$ENTREZID, rownames(mData.scaled))
mData.scaled = mData.scaled[i,]
identical(rownames(mData.scaled), df$ENTREZID)
rownames(mData.scaled) = df$SYMBOL

## create data.frame
dfData = data.frame(t(mData.scaled))
###########################################################

###########################################################
########## Band3/Slc4a1 expression models of various sizes
###########################################################
library(rethinking)
## choose the variable to model
colnames(dfData)

fit.1 <- quap(
  alist(
    Slc4a1 ~ dnorm(mu, sigmaPop),
    mu <- b0 + Ga*Gata1,
    b0 ~ dnorm(0, 2),
    c(Ga) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=dfData,
  start=list(b0=0)
)
summary(fit.1)
precis(fit.1)
post = extract.samples(fit.1, n=1000)
dim(post)
head(post)
precis(post)
plot(precis(post))
plot(coeftab(fit.1), pars=c('b0', 'Ga'))
pairs(post)

###### add another covariate, fog1 22761 zfpm1
fit.2 <- quap(
  alist(
    Slc4a1 ~ dnorm(mu, sigmaPop),
    mu <- b0 + Ga*Gata1 + Fog*Zfpm1,
    b0 ~ dnorm(0, 2),
    c(Ga, Fog) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=dfData,
  start=list(b0=0)
)
summary(fit.2)
plot(coeftab(fit.1, fit.2), pars=c('b0', 'Ga', 'Fog'))
plot(compare(fit.1, fit.2))
pairs(fit.2)

################## simulations from DAG
nsim = 1000
df.Sim = data.frame(1:nsim)
df.Sim$Gata1 = rnorm(nsim)
df.Sim$Zfpm1 = rnorm(nsim)
m = cbind(1, as.matrix(df.Sim[,-1]))
head(m)
coef(fit.2)
m = m %*% coef(fit.1)[-4]
df.Sim$Slc4a1 = rnorm(nsim, m, coef(fit.2)['sigmaPop'])


## recover coefficients using simulated data fit
fit.2.sim <- quap(
  alist(
    Slc4a1 ~ dnorm(mu, sigmaPop),
    mu <- b0 + Ga*Gata1 + Fog*Zfpm1,
    b0 ~ dnorm(0, 2),
    c(Ga, Fog) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=df.Sim,
  start=list(b0=0)
)

plot(coeftab(fit.2, fit.2.sim))
## simulate new response variable
str(df.Sim)
coef(fit.1.cd)
m = cbind(1, df.Sim$f13.Peanut, df.Sim$total.IgE)
head(m)
m = m %*% coef(fit.1.cd)[-4]
dim(m)
df.Sim$CD63.Act = rnorm(nsim, m[,1])

## now fit the simulated data in the same model
fit.1.cd.sim <- quap(
  alist(
    CD63.Act ~ dnorm(mu, sigmaPop),
    mu <- b0 + b1*f13.Peanut + b2*total.IgE,
    b0 ~ dcauchy(0, 2),
    c(b1, b2) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=df.Sim,
  start=list(b0=0)
)
summary(fit.1.cd.sim)
plot(coeftab(fit.1.cd, fit.1.cd.sim), pars=c('b0', 'b1', 'b2'))

fit.2.cd.sim <- quap(
  alist(
    CD63.Act ~ dnorm(mu, sigmaPop),
    mu <- b0 + b1*f13.Peanut,
    b0 ~ dcauchy(0, 2),
    c(b1) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=df.Sim,
  start=list(b0=0)
)
summary(fit.2.cd.sim)
plot(coeftab(fit.1.cd, fit.1.cd.sim, fit.2.cd, fit.2.cd.sim), pars=c('b0', 'b1', 'b2'))

plot(density(dfData.pa$CD63.Act))
plot(density(df.Sim$CD63.Act))

###########################################################
