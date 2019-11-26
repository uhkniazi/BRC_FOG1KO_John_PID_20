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
dim(mData.scaled)
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

################## simulations from DAG1
nsim = 1000
df.Sim = data.frame(1:nsim)
df.Sim$Gata1 = rnorm(nsim)
df.Sim$Zfpm1 = rnorm(nsim)
m = cbind(1, as.matrix(df.Sim[,-1]))
head(m)
coef(fit.2)
m = m %*% coef(fit.1)[-4]
df.Sim$Slc4a1 = rnorm(nsim, m)


## recover coefficients using simulated data fit
fit.1.sim <- quap(
  alist(
    Slc4a1 ~ dnorm(mu, sigmaPop),
    mu <- b0 + Ga*Gata1,
    b0 ~ dnorm(0, 2),
    c(Ga) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=df.Sim,
  start=list(b0=0)
)

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

plot(coeftab(fit.2, fit.1.sim, fit.2.sim), pars=c('b0', 'Ga', 'Fog'))

################## simulations from DAG2
nsim = 1000
df.Sim2 = data.frame(1:nsim)
df.Sim2$Gata1 = rnorm(nsim)
df.Sim2$Lmo2 = rnorm(nsim)
df.Sim2$Ldb1 = rnorm(nsim)
df.Sim2$Tcf3 = rnorm(nsim)
df.Sim2$Stil = rnorm(nsim)

fit.temp <- quap(
  alist(
    Zfpm1 ~ dnorm(mu, sigmaPop),
    mu <- b0 + Ga*Gata1 + Lm*Lmo2 + Ld*Ldb1 + Tc*Tcf3 + St*Stil,
    b0 ~ dnorm(0, 2),
    c(Ga, Lm, Ld, Tc, St) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=dfData,
  start=list(b0=0)
)
summary(fit.temp)
colnames(df.Sim2)
m = cbind(1, as.matrix(df.Sim2[,-1]))
head(m)
coef(fit.temp)
m = m %*% coef(fit.temp)[-7]
df.Sim2$Zfpm1 = rnorm(nsim, m)

fit.temp.sim <- quap(
  alist(
    Zfpm1 ~ dnorm(mu, sigmaPop),
    mu <- b0 + Ga*Gata1 + Lm*Lmo2 + Ld*Ldb1 + Tc*Tcf3 + St*Stil,
    b0 ~ dnorm(0, 2),
    c(Ga, Lm, Ld, Tc, St) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=df.Sim2,
  start=list(b0=0)
)

plot(coeftab(fit.temp, fit.temp.sim))

### generate the Slc data now
colnames(df.Sim2)
m = cbind(1, as.matrix(df.Sim2[,c('Gata1', 'Zfpm1')]))
head(m)
coef(fit.2)
m = m %*% coef(fit.2)[-4]
df.Sim2$Slc4a1 = rnorm(nsim, m)

## recover coefficients using simulated data fit to Dag2
fit.2.sim2 <- quap(
  alist(
    Slc4a1 ~ dnorm(mu, sigmaPop),
    mu <- b0 + Ga*Gata1 + Fog*Zfpm1,
    b0 ~ dnorm(0, 2),
    c(Ga, Fog) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=df.Sim2,
  start=list(b0=0)
)

plot(coeftab(fit.2, fit.2.sim, fit.2.sim2), pars=c('b0', 'Ga', 'Fog'))

#### fog1 alone
fit.1.f <- quap(
  alist(
    Slc4a1 ~ dnorm(mu, sigmaPop),
    mu <- b0 + Fog*Zfpm1,
    b0 ~ dnorm(0, 2),
    c(Fog) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=dfData,
  start=list(b0=0)
)

plot(coeftab(fit.1, fit.1.f, fit.2, fit.2.sim, fit.2.sim2), pars=c('b0', 'Ga', 'Fog'))
plot(coeftab(fit.1, fit.1.f, fit.2, fit.2.sim2))


#### add a section on fake-data simulation 
#### compare with the original data
## simulate data
fake.fit.1 = sim(fit.1, n = 300)
fake.fit.1.sim = sim(fit.1.sim, n=300)

fake.fit.2 = sim(fit.2, n = 300)
fake.fit.2.sim = sim(fit.2.sim, n=300)
fake.fit.2.sim2 = sim(fit.2.sim2, n=300)
hist(dfData$Slc4a1, prob=T)
lines(density(colMeans(fake.fit.2)), col=1)
lines(density(colMeans(fake.fit.2.sim)), col=2)
lines(density(colMeans(fake.fit.2.sim2)), col=3)

lines(density(df.Sim$Slc4a1), col=2)
lines(density(df.Sim2$Slc4a1), col=3)
lines(density(dfData$Slc4a1), col=1)
###########################################################

###########################################################
########## Hbb.b1 expression models of various sizes
###########################################################
fit.1a <- quap(
  alist(
    Hbb.b1 ~ dnorm(mu, sigmaPop),
    mu <- b0 + Ga*Gata1,
    b0 ~ dnorm(0, 2),
    c(Ga) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=dfData,
  start=list(b0=0)
)

fit.1b <- quap(
  alist(
    Hbb.b1 ~ dnorm(mu, sigmaPop),
    mu <- b0 + Fog*Zfpm1,
    b0 ~ dnorm(0, 2),
    c(Fog) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=dfData,
  start=list(b0=0)
)

###### both covariates Gata1 and fog1 22761 zfpm1
fit.2 <- quap(
  alist(
    Hbb.b1 ~ dnorm(mu, sigmaPop),
    mu <- b0 + Ga*Gata1 + Fog*Zfpm1,
    b0 ~ dnorm(0, 2),
    c(Ga, Fog) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=dfData,
  start=list(b0=0)
)
summary(fit.2)
plot(coeftab(fit.1a, fit.1b, fit.2), pars=c('b0', 'Ga', 'Fog'))
plot(compare(fit.1a, fit.1b, fit.2))
pairs(fit.2)

### nurd complex and its affect on fog
mNurd = as.matrix(dfData[,c('Hdac1', 'Hdac2', 'Mta1', 'Mta2', 'Mta3', 'Mbd2', 'Mbd3', 'Chd3', 'Chd4', 'Rbbp4', 'Rbbp7')])
mPc = prcomp(t(mNurd), scale=F)
mPc = mPc$rotation[,1:4]

dfNurd = data.frame(Zfpm1=dfData$Zfpm1, mPc)

# fit.nurd <- quap(
#   alist(
#     Zfpm1 ~ dnorm(mu, sigmaPop),
#     mu <- b0 + hd1*Hdac1 + hd2*Hdac2 + mt1*Mta1 + mt2*Mta2 + mt3*Mta3 + mb2*Mbd2 + mb3*Mbd3 + ch3*Chd3 + ch4*Chd4 + rb4*Rbbp4 + rb7*Rbbp7,
#     b0 ~ dnorm(0, 1),
#     c(hd1, hd2, mt1, mt2, mt3, mb2, mb3, ch3, ch4, rb4, rb7) ~ dnorm(0, 0.5),
#     sigmaPop ~ dexp(1)
#   ), data=dfData,
#   start=(list(b0=0, hd1=0, hd2=0, mt1=0, mt2=0, mt3=0, mb2=0, mb3=0, ch3=0, ch4=0, rb4=0, rb7=0))
# )

fit.nurd <- quap(
  alist(
    Zfpm1 ~ dnorm(mu, sigmaPop),
    mu <- b0 + p1*PC1 + p2*PC2 + p3*PC3 + p4*PC4, 
    b0 ~ dnorm(0, 1),
    c(p1, p2, p3, p4) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=dfNurd,
  start=(list(b0=0))
)

summary(fit.nurd)
plot(coeftab(fit.nurd))

fake.fit.nurd = sim(fit.nurd, n = 300)
hist(dfData$Zfpm1, prob=T)
lines(density(colMeans(fake.fit.nurd)), col=1)

#### simulate new hbb using original model and 2 step dag
#### sim 1 = using original model
nsim = 1000
df.Sim = data.frame(1:nsim)
df.Sim$Gata1 = rnorm(nsim)
df.Sim$Zfpm1 = rnorm(nsim)
m = cbind(1, as.matrix(df.Sim[,-1]))
head(m)
coef(fit.2)
m = m %*% coef(fit.2)[-4]
df.Sim$Hbb.b1 = rnorm(nsim, m)

### sim 2 = first step to simulate zfpm1
nsim = 1000
df.Sim2 = data.frame(1:nsim)
df.Sim2$PC1 = rnorm(nsim)
df.Sim2$PC2 = rnorm(nsim)
df.Sim2$PC3 = rnorm(nsim)
df.Sim2$PC4 = rnorm(nsim)
m = cbind(1, as.matrix(df.Sim2[,-1]))
head(m)
coef(fit.nurd)
m = m %*% coef(fit.nurd)[-6]
df.Sim2$Zfpm1 = rnorm(nsim, m)

## now simulate the second step of the model
df.Sim2$Gata1 = rnorm(nsim)
colnames(df.Sim2)
m = cbind(1, as.matrix(df.Sim2[,c('Gata1', 'Zfpm1')]))
head(m)
coef(fit.2)
m = m %*% coef(fit.2)[-4]
df.Sim2$Hbb.b1 = rnorm(nsim, m)

#### sim 3 = Gata1 complex, affecting Fog1
nsim = 1000
df.Sim3 = data.frame(1:nsim)
df.Sim3$Gata1 = rnorm(nsim)
df.Sim3$Lmo2 = rnorm(nsim)
df.Sim3$Ldb1 = rnorm(nsim)
df.Sim3$Tcf3 = rnorm(nsim)
df.Sim3$Stil = rnorm(nsim)

## estimate the coefficients to generate fog from original data
fit.gataPlex <- quap(
  alist(
    Zfpm1 ~ dnorm(mu, sigmaPop),
    mu <- b0 + Ga*Gata1 + Lm*Lmo2 + Ld*Ldb1 + Tc*Tcf3 + St*Stil,
    b0 ~ dnorm(0, 2),
    c(Ga, Lm, Ld, Tc, St) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=dfData,
  start=list(b0=0)
)
summary(fit.gataPlex)
colnames(df.Sim3)
m = cbind(1, as.matrix(df.Sim3[,-1]))
head(m)
coef(fit.gataPlex)
m = m %*% coef(fit.gataPlex)[-7]
df.Sim3$Zfpm1 = rnorm(nsim, m)

## now simulate the second step of the model
colnames(df.Sim3)
m = cbind(1, as.matrix(df.Sim3[,c('Gata1', 'Zfpm1')]))
head(m)
coef(fit.2)
m = m %*% coef(fit.2)[-4]
df.Sim3$Hbb.b1 = rnorm(nsim, m)

### now check which of the 2 simulated data sets recover the 
### coefficients for the original model
fit.2.sim1 <- quap(
  alist(
    Hbb.b1 ~ dnorm(mu, sigmaPop),
    mu <- b0 + Ga*Gata1 + Fog*Zfpm1,
    b0 ~ dnorm(0, 2),
    c(Ga, Fog) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=df.Sim,
  start=list(b0=0)
)

fit.2.sim2 <- quap(
  alist(
    Hbb.b1 ~ dnorm(mu, sigmaPop),
    mu <- b0 + Ga*Gata1 + Fog*Zfpm1,
    b0 ~ dnorm(0, 2),
    c(Ga, Fog) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=df.Sim2,
  start=list(b0=0)
)

fit.2.sim3 <- quap(
  alist(
    Hbb.b1 ~ dnorm(mu, sigmaPop),
    mu <- b0 + Ga*Gata1 + Fog*Zfpm1,
    b0 ~ dnorm(0, 2),
    c(Ga, Fog) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=df.Sim3,
  start=list(b0=0)
)

plot(coeftab(fit.2, fit.2.sim1, fit.2.sim2, fit.2.sim3), pars=c('Ga', 'Fog'))
plot(coeftab(fit.1a, fit.1b, fit.2, fit.2.sim1, fit.2.sim2, fit.2.sim3), pars=c('b0', 'Ga', 'Fog'))

hist(dfData$Hbb.b1, prob=T)
lines(density(colMeans(sim(fit.2.sim1, n=100))))
plot(dfData$Hbb.b1, dfData$Zfpm1, pch=20, col=c(1,2, 3, 4, 5)[as.numeric(factor(dfSample.2$group3))])
points(colMeans(sim(fit.2, n=100)), dfData$Zfpm1, col=c(1,2, 3, 4, 5)[as.numeric(factor(dfSample.2$group3))], pch='*')
legend('bottomright', legend = levels(factor(dfSample.2$group3)), fill=c(1,2,3,4, 5))

plot(colMeans(sim(fit.2.sim1, n=100)), df.Sim$Zfpm1, pch=20)
plot(df.Sim$Hbb.b1, df.Sim$Zfpm1, pch=20)

plot(colMeans(sim(fit.2.sim2, n=100)), df.Sim2$Zfpm1, pch=20)
plot(df.Sim2$Hbb.b1, df.Sim2$Zfpm1, pch=20)

plot(colMeans(sim(fit.2.sim1, n=100)), df.Sim$Gata1, pch=20)
plot(df.Sim$Hbb.b1, df.Sim$Gata1, pch=20)

##########################################################

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
