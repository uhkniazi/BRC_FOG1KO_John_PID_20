# File: 11_vennDiagram.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: using the results for each contrast, create venn diagrams and some plots
# Date: 30/9/2019

source('header.R')

lFiles = list.files('results/', pattern='DEAnalysis*', full.names = T, ignore.case = T)

ldfData = lapply(lFiles, function(x) as.data.frame(read.csv(x, header=T, row.names=1, stringsAsFactors = F)))
names(ldfData) = lFiles
sapply(ldfData, nrow)

# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

sapply(ldfData, function(df) identical(rownames(df), rn))

cvTitle = gsub('results//DEAnalysis', '', names(ldfData))
cvTitle = gsub('.xls', '', cvTitle)

for (i in seq_along(cvTitle)){
  df = ldfData[[i]]
  hist(df$logFC, xlab='Log Fold Change', ylab='', main=paste('Fold Changes ', cvTitle[i], sep=''))
  f_plotVolcano(df, cvTitle[i], 0.01, fc.lim=range(df$logFC))
}

names(ldfData)
## select significant genes
dfContrast1.sub = ldfData[[1]][ldfData[[1]]$adj.P.Val < 0.1,]
dfContrast2.sub = ldfData[[2]][ldfData[[2]]$adj.P.Val < 0.1,]
dfContrast3.sub = ldfData[[3]][ldfData[[3]]$adj.P.Val < 0.1,]
dfContrast4.sub = ldfData[[4]][ldfData[[4]]$adj.P.Val < 0.1,]
dfContrast5.sub = ldfData[[5]][ldfData[[5]]$adj.P.Val < 0.1,]
dfContrast6.sub = ldfData[[6]][ldfData[[6]]$adj.P.Val < 0.1,]
dfContrast7.sub = ldfData[[7]][ldfData[[7]]$adj.P.Val < 0.1,]
dfContrast8.sub = ldfData[[8]][ldfData[[8]]$adj.P.Val < 0.1,]
dfContrast9.sub = ldfData[[9]][ldfData[[9]]$adj.P.Val < 0.1,]
dfContrast10.sub = ldfData[[10]][ldfData[[10]]$adj.P.Val < 0.1,]


library(VennDiagram)

# create a list for overlaps
lVenn = list(rownames(dfContrast1.sub), rownames(dfContrast2.sub), rownames(dfContrast3.sub),
             rownames(dfContrast4.sub), rownames(dfContrast5.sub), rownames(dfContrast6.sub),
             rownames(dfContrast7.sub), rownames(dfContrast8.sub), rownames(dfContrast9.sub),
             rownames(dfContrast10.sub)
)
names(ldfData)
names(lVenn) = cvTitle
# calculate overlaps
#lVenn.overlap = calculate.overlap(lVenn)
venn.diagram(lVenn[c(1, 3, 5)], filename = 'results/venn_all_contrasts_DiffMutVsDiffWT.tif', margin=0.1)
venn.diagram(lVenn[c(2, 4, 6, 7)], filename = 'results/venn_all_contrasts_DiffVsUndiff.tif', margin=0.1)
venn.diagram(lVenn[c(8, 9, 10)], filename = 'results/venn_all_contrasts_UndiffMutVsUndiffWT.tif', margin=0.1)

## skip this section of separating them into up and down as this makes too many
## combinations for this data set
### repeat the analysis but separate up and down regulated genes
## select significant genes
# dfContrast1.up = dfContrast1.sub[dfContrast1.sub$logFC > 0, ]
# dfContrast1.down = dfContrast1.sub[dfContrast1.sub$logFC < 0, ]
# 
# dfContrast2.up = dfContrast2.sub[dfContrast2.sub$logFC > 0, ]
# dfContrast2.down = dfContrast2.sub[dfContrast2.sub$logFC < 0, ]
# 
# dfContrast3.up = dfContrast3.sub[dfContrast3.sub$logFC > 0, ]
# dfContrast3.down = dfContrast3.sub[dfContrast3.sub$logFC < 0, ]
# 
# dfContrast4.up = dfContrast4.sub[dfContrast4.sub$logFC > 0, ]
# dfContrast4.down = dfContrast4.sub[dfContrast4.sub$logFC < 0, ]
# 
# dfContrast5.up = dfContrast5.sub[dfContrast5.sub$logFC > 0, ]
# dfContrast5.down = dfContrast5.sub[dfContrast5.sub$logFC < 0, ]
# 
# dfContrast6.up = dfContrast6.sub[dfContrast6.sub$logFC > 0, ]
# dfContrast6.down = dfContrast6.sub[dfContrast6.sub$logFC < 0, ]
# 
# dfContrast7.up = dfContrast7.sub[dfContrast7.sub$logFC > 0, ]
# dfContrast7.down = dfContrast7.sub[dfContrast7.sub$logFC < 0, ]
# 
# dfContrast8.up = dfContrast8.sub[dfContrast8.sub$logFC > 0, ]
# dfContrast8.down = dfContrast8.sub[dfContrast8.sub$logFC < 0, ]
# 
# dfContrast9.up = dfContrast9.sub[dfContrast9.sub$logFC > 0, ]
# dfContrast9.down = dfContrast9.sub[dfContrast9.sub$logFC < 0, ]
# 
# dfContrast10.up = dfContrast10.sub[dfContrast10.sub$logFC > 0, ]
# dfContrast10.down = dfContrast10.sub[dfContrast10.sub$logFC < 0, ]
# 
# 
# # create a list for overlaps
# lVenn = list(rownames(dfContrast1.up), rownames(dfContrast1.down),
#              rownames(dfContrast2.up), rownames(dfContrast2.down), 
#              rownames(dfContrast3.up), rownames(dfContrast3.down),
#              rownames(dfContrast4.up), rownames(dfContrast4.down),
#              rownames(dfContrast5.up), rownames(dfContrast5.down),
#              rownames(dfContrast6.up), rownames(dfContrast6.down),
#              rownames(dfContrast7.up), rownames(dfContrast7.down),
#              rownames(dfContrast8.up), rownames(dfContrast8.down),
#              rownames(dfContrast9.up), rownames(dfContrast9.down),
#              rownames(dfContrast10.up), rownames(dfContrast10.down)
# )
# 
# 
# #cvTitle = gsub('results//DEAnalysis(\\w+)VsControl.xls', '\\1', names(ldfData))
# cvTitle.up = paste(cvTitle, 'up', sep='-')
# cvTitle.down = paste(cvTitle, 'down', sep='-')
# o = c(1, 11, 2, 12, 3, 13, 4, 14, 5, 15, 6, 16, 7, 17, 8, 18, 9, 19, 10, 20)
# ## sanity check
# matrix(c(cvTitle.up, cvTitle.down)[o], ncol = 2, byrow = T)
# cvTitle = c(cvTitle.up, cvTitle.down)[o]
# names(lVenn) = cvTitle

## create a binary matrix
cvCommonGenes = unique(do.call(c, lVenn))
mCommonGenes = matrix(NA, nrow=length(cvCommonGenes), ncol=length(lVenn))
for (i in 1:ncol(mCommonGenes)){
  mCommonGenes[,i] = cvCommonGenes %in% lVenn[[i]]
}
rownames(mCommonGenes) = cvCommonGenes
colnames(mCommonGenes) = names(lVenn)

# create groups in the data based on permutation with repetition: 2 ^ ncol 
# https://www.mathsisfun.com/combinatorics/combinations-permutations.html
mCommonGenes.grp = mCommonGenes
set.seed(123)
dm = dist(mCommonGenes.grp, method='binary')
hc = hclust(dm)
plot(hc)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.1)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mCommonGenes.grp = cbind(mCommonGenes.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mCommonGenes.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## put these results together
dfCommonGenes = data.frame(mCommonGenes, sig.pvals=rowSums(mCommonGenes), groups=cp, Symbol=ldfData[[1]][cvCommonGenes, 'SYMBOL'])
## gene names have an X before them remove those if present
rownames(dfCommonGenes) = (gsub(pattern = '^X', replacement = '', x = rownames(dfCommonGenes)))
head(dfCommonGenes)

write.csv(dfCommonGenes, file='results/commonDEGenes.xls')
