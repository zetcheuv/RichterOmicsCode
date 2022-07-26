#!/usr/bin/Rscript --default-packages=stats,utils,methods,grDevices,graphics
#@seb2016 (sebastien.hergalant@inserm.fr)
#@seb2018 - revised for linear models, moderated t-tests and covariates, added comBat analysis for batch effects

### transcriptomic analyses at gene level
### matrix prep, filtering, annotation, normalization, comBat, DE (moderated t-test)

### 5 arguments - gene_count_matrix.csv|transcript_count_matrix.csv TR2IDs.txt|ENST2IDs.txt sample.table readcutoff[~20se|~40pe] fdrcutoff[0.05]
# 1- gene_count_matrix.csv = raw read count compil output from pipeline, gene- or transcript-wise, as output by stringtie
# 2- TR2IDs.txt = ID annotations with corresponding Ensembl and Gene symbols, as output by stringtie
# 3- sample.table = phenodata with sample names, groups and transciptome names
# 4- READcutoff = by gene mean raw-read-depth threshold to apply for filtering low expression values across all (or control) samples

### bioc libraries
suppressPackageStartupMessages(library(limma)) # for linear models, bayes, voom and so much else
suppressPackageStartupMessages(library(edgeR)) # for TMM normalization

### args
args <- commandArgs(trailingOnly=TRUE)
genecountFile <- ifelse(is.na(args[1]),"./3-transcripts/all-samples_counts/gene_count_matrix.csv",args[1])
stopifnot(file.exists(genecountFile))
tridFile <- ifelse(is.na(args[2]),"./3-transcripts/all-samples_counts/TR2IDs.txt",args[2])
stopifnot(file.exists(tridFile))
phenoFile <- ifelse(is.na(args[3]),"sample.table",args[3])
stopifnot(file.exists(phenoFile))
readcut <- ifelse(is.na(args[4]), 20, as.numeric(args[4])) # 20 for single-read, 40 for paired-end and 80 for paired-end transcripts???
stopifnot(is.numeric(readcut))
fdrcut <- ifelse(is.na(args[5]), 0.05, as.numeric(args[5]))
stopifnot(is.numeric(fdrcut))


######################################
### 
whole.table <- read.csv(genecountFile, header=T, stringsAsFactor=F, check.names=F)
genes <- read.delim(tridFile, header=F, stringsAsFactors=F, sep="\t", check.names=F)
phenodata <- read.csv(phenoFile, header=TRUE, sep="\t")
# add an automatic prefix to output file names (usually "genes_" or "transcripts_")
filePrefix <- paste0(sub("_.*", "", basename(genecountFile)), "s_")

## annotating genes and samples
gs <- genes[match(whole.table[,1], genes[,1]),2]
gs[is.na(gs)] <- as.character(whole.table[is.na(gs),1])
rownames(whole.table) <- gs
whole.table <- whole.table[,-1]
colnames(whole.table) <- sub("_stats", "", colnames(whole.table))
#now reorder samples according to the phenodata file, and make a quick check by the way
stopifnot(all(sort(colnames(whole.table)) == sort(as.character(phenodata$Sample_Name))))
whole.table <- whole.table[, as.character(phenodata$Sample_Name)]
#change colnames by TransName column if present
if("TransName" %in% colnames(phenodata)) colnames(whole.table) <- as.character(phenodata$TransName)

## filtering
whole.table <- whole.table[!apply(whole.table,1,function(x) length(which(x==0)) > length(x)-2),]
whole.table <- whole.table[,!apply(whole.table,2,function(x) length(which(x==0)) > length(x)-2)]
#allow only genes with sufficient mean expression level across samples
whole.table <- whole.table[which(apply(whole.table, 1, mean) >= readcut),] # ne pas remplacer les valeurs par 0

# set the default colors for plotting by generating random colors not taken from the grey shades, as they are too similar
options(device="png")
palette(sample(colors()[grep('gr(a|e)y', colors(), invert = T)], ncol(whole.table)))
x <- dev.off()
## normalization with TMM and VOOM, do some plotting
dge <- DGEList(whole.table)
dge <- calcNormFactors(dge, method="TMM")
# voom can optionally normalize with quantiles if the data are too noisy
x <- dev.next()
png(filename=paste0(filePrefix,"meanvar.png"))
voom <- voom(dge, plot=TRUE, normalize.method="none")
x <- dev.off()
logcpmadj <- as.matrix(as.data.frame(voom$E))
# densities
x <- dev.next()
png(filename=paste0(filePrefix,"densities.png"))
plotDensities(logcpmadj, legend=FALSE)
x <- dev.off()
# sample correlations
median_vector_norm <- apply(logcpmadj, 1, median)
correl_norm <- apply(logcpmadj, 2, cor, median_vector_norm)
x <- dev.next()
png(filename=paste0(filePrefix,"correlations.png"))
plot(correl_norm,type="l",ylim=c(0,1),ylab="correlation coeff.",xlab="samples",lwd=5,col="blue")
x <- dev.off()

# save these normalized and filtered data into a table that can be used in statistics, clusterings, etc...
write.table(logcpmadj, paste0(filePrefix,"norm-log2cpm.txt"), quote=F, sep="\t", row.names=T, col.names=NA)


