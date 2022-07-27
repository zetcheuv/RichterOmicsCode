#@seb2022 sebastien.hergalant@inserm.fr
#workflow for computing LCS (linear classifier scores) from expression data
#applicable to public datasets (like the ones found on Gene Expression Omnibus) from microarray or RNAseq data

#requirement1: normalized/filtered/prepared transcriptome table with genes in rows and samples in columns
#requirement2: the 215 gene signature file with associated up and down tags for weighting the scores

# TRANS: a matrix containing transcriptomic data, with gene names as rownames (can be duplicated)
stopifnot(is.matrix(TRANS))
# scoringSig215: a table with 215 gene names in column 1 and weights column 2
scoringSig215 <- c(rep(-1,93), rep(1,122))
names(scoringSig215) <- read.delim("LCS_215genesig.txt", stringsAsFactors=F)[,2]

#reduce the transcriptome to the 215 genes 
TRANS_scoring <- TRANS[rownames(TRANS) %in% names(scoringSig215),]
scoringSig <- scoringSig215[names(scoringSig215) %in% rownames(TRANS_scoring),]

#aggregate rows matching the same gene in one row by summing its expression
TRANS_scoring <- aggregate(TRANS_scoring, by=list(rownames(TRANS_scoring)), sum)
rownames(TRANS_scoring) <- TRANS_scoring[,1]
TRANS_scoring <- TRANS_scoring[,-1]

#finally ensure that both object are in the same gene order
TRANS_scoring <- TRANS_scoring[names(scoringSig),]

#now everything is set, compute the scores
# 1) standardize the data
TRANS_scoring <- t(scale(t(TRANS_scoring)))
# 2) remove the last permille (upper outlier effect dampening, the lowest outliers have no other effect than diminishing the scores)
TRANS_scoring[TRANS_scoring > quantile(,0.999)] <- quantile(TRANS_scoring,0.999)
# 3) apply LCS over samples. LCS is dependant on the dataset (constitution, size, quality, ...)
scores <- apply(TRANS_scoring, 2, function(x) mean(x * scoringSig) )
# 4) compute Zscores from the scores, these will be more comparable across datasets
zscores <- sapply(scores, function(x) (x-mean(scores))/sd(scores))
# 5) associate a pvalue to the Zscores, and so to the scores
scores_pvalues <- pnorm(-zscores)

#output these results to a table
write.table(cbind(sort(scores, decreasing=T), sort(zscores, decreasing=T), sort(scores_pvalues)), "LCS.txt", sep="\t", row.names=T, col.names=F, quote=F)





