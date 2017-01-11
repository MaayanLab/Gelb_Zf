library(limma)
library(edgeR)
library(statmod)

setwd('~/Documents/Zichen_Projects/Flavia2/')
## help functions
splitTopTable <- function(table, filename) {
  # to split a topTable into up and down genes
  # and write them into csv files
  maskUp <- table$logFC > 0
  maskDn <- table$logFC < 0
  tableUp <- table[maskUp, ]
  tableDn <- table[maskDn, ]
  write.csv(tableUp, file=paste0(filename, '-up.csv'), row.names=F)
  write.csv(tableDn, file=paste0(filename, '-dn.csv'), row.names=F)
}

write.limma.topTables <- function(file, tables) {
  # write multiple toptables to a excel file
  # tables is a list of topTables
  require(xlsx, quietly = TRUE)
  splited.topTables <- list()
  for (table.name in names(tables)) {
    table <- tables[[table.name]]
    maskUp <- table$logFC > 0
    maskDn <- table$logFC < 0
    tableUp <- table[maskUp, ]
    tableDn <- table[maskDn, ]
    table.up.name <- paste0(table.name, '-up')
    table.dn.name <- paste0(table.name, '-dn')
    splited.topTables[[table.dn.name]] <- tableDn
    splited.topTables[[table.up.name]] <- tableUp
  }
  
  splited.topTables.names <- names(splited.topTables)
  for (i in 1:length(splited.topTables.names)) {
    table.name <- splited.topTables.names[i]
    table <- splited.topTables[[table.name]]
    if (i == 1) {
      write.xlsx(table, file, sheetName = table.name, row.names=TRUE)
    } else {
      write.xlsx(table, file, sheetName = table.name, row.names=TRUE, append=TRUE)
    }
  }
}

## prepare count.mat
df <- read.csv('featureCounts_matrix.csv',check.names=F)
count.mat <- as.matrix(df[2:length(df)])
genes <- as.vector(df[[1]])
lengths.df <- read.csv('featureCounts_gene_lengths.csv', check.names=F)
gene.lengths <- lengths.df[[2]]
rownames(count.mat) <- genes
# colnames(count.mat) <- c('dKO_1', 'dKO_2', 'WT_1', 'WT_2')
meta_df <- read.csv('meta_df.csv')

## prepare for the design matrix
genotypes <- factor(meta_df$genotype)
tissues <- factor(meta_df$tissue)
conditions <- factor(paste(meta_df$genotype, meta_df$tissue, sep="."))
design <- model.matrix(~0+conditions)
colnames(design) <- levels(conditions)


### filter count matrix
filter <- apply(count.mat, 1, function(x) length(x[x>20])>=2)
count.mat.filtered <- count.mat[filter,]
genes <- genes[filter]

## voom-norm
dge <- DGEList(counts=count.mat.filtered)
dge$genes <- genes
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot=TRUE)

## estimate the correlation between measurements made on the same subject
corfit <- duplicateCorrelation(v,design,block=meta_df$animal_id)
corfit$consensus

## run limma to call DEGs with design matrix
fit <- lmFit(v,design, block=meta_df$animal_id, correlation=corfit$consensus)
# fit <- lmFit(v,design)
cont.matrix <- makeContrasts(KOvsWT.BASAL = dko.BASAL-wt.BASAL,
                             KOvsWT.ER = dko.ER-wt.ER,
                             KOvsWT.LP = dko.LP-wt.LP,
                             levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

## examine number of DEGs
qval.cutoff <- 0.05
foldChange.cutoff <- 1
table.DEGs <- topTable(fit2,coef=1, number=dim(fit2)[1], p.value=qval.cutoff, lfc=foldChange.cutoff)
cat('DEGs:', dim(table.DEGs)[1])

## write excel file
tables <- list()
table.DEGs <- topTable(fit2,coef=1, number=dim(fit2)[1])
tables[['dKOvsWT']] <- table.DEGs
write.limma.topTables('DEGs_voom-limma.xlsx', tables)
