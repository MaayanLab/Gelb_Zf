library(edgeR)
library(limma)

RowVar <- function(x) {
	## to filter out zero variance genes
	rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}


runLimma <- function(data, rids, sampleclass) {
	## run limma to identify DEGs.
	## data: matrix of expression data
	## rids: vector of genes
	## sampleclass: factor with 1 being ctrls and 2 being perts, 0 being samples to exclude
	mask <- sampleclass != 0
	sampleclass <- sampleclass[mask]
	data <- data[, mask]

	design <- cbind(Grp1=1, Grp2vs1=as.numeric(sampleclass) - 1)
	rownames(data) <- rids
	
	fit <- lmFit(data, design)
	fit <- eBayes(fit)
	p.values <- as.numeric(fit$p.value[,2])
	q.values <- p.adjust(p.values, method="BH")
	logFC <- as.numeric(fit$coefficients[,2]) 
	results <- list(pvalues=p.values, qvalues=q.values, logFC=logFC)
	results
}


runVoomLimma <- function(data, rids, sampleclass) {
	## run voom-limma to identify DEGs for RNA-seq read counts.
	## data: read count matrix of expression data
	## rids: vector of genes
	## sampleclass: factor with 1 being ctrls and 2 being perts	
	mask <- sampleclass != 0
	sampleclass <- sampleclass[mask]
	data <- data[, mask]
	
	design <- cbind(Grp1=1, Grp2vs1=as.numeric(sampleclass) - 1)
	rownames(data) <- rids

	dge <- DGEList(counts=data)
	dge$genes <- rids
	dge <- calcNormFactors(dge)
	v <- voom(dge, design, plot=FALSE)

	fit <- lmFit(v, design)
	fit <- eBayes(fit)
	p.values <- as.numeric(fit$p.value[,2])
	q.values <- p.adjust(p.values, method="BH")
	logFC <- as.numeric(fit$coefficients[,2]) 
	results <- list(pvalues=p.values, qvalues=q.values, logFC=logFC)
	results
}
