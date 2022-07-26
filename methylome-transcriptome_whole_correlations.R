#!/usr/bin/Rscript
# Script for methylation-expression correlation
# @seb2017 (sebastien.hergalant@inserm.fr)
# Adapted for RNASeq data (transcript level), added some features, some major revisions, rewrote chunks of code
# @seb2018 rewrote the promoter region part, some annotation parts for expanding results
# @seb2020 rewrote the split transcript/gene annots and split CpG annots for a better genome-wide correlation, works with EPIC, works with genes and not just at transcript level, check params and set default values, commented code with loads of stuff where relevant, added integrations for positive correlations within gene bodies, optimized some code for time consuming tasks (no parallel computing required)

#args[1] = transcriptome prepared/normalized table (transcripts with ENST, or genes with ENSG)
#args[2] = methylome beta-values/prepared table (all probes can be kept for computation, cleaning SNP sites and else can occur later)
#args[3] = transcripts or genes ? (transcriptome)
#args[4] = EPIC or 450k annotation file (methylome)

# Load biomaRt and start it
suppressPackageStartupMessages(library(biomaRt))
ensembl <- useMart("ensembl")
# if some orthology with Homo sapiens needs to be done with transcriptome data, do it beforehand !
ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl) # should stay 'as is' for now because methylation platforms are annotated for Hs

# Load tables, ordering expression/methylation table and sanity checks
args <- commandArgs(trailingOnly=TRUE)
trans <- ifelse(is.na(args[1]), "correl_trans.txt.gz", args[1]) 
stopifnot(file.exists(trans))
meth <- ifelse(is.na(args[2]), "correl_meth.txt.gz", args[2])
stopifnot(file.exists(meth))
transLevel <- tolower(ifelse(is.na(args[3]), "transcript", args[3]))
stopifnot(transLevel %in% c("transcript", "gene"))
ucscAnnot <- ifelse(transLevel == "transcript", "Accession", "Name")
methAnnotFile <- ifelse(is.na(args[4]), "/home/chac/methyl/Illumina_EPIC850k_simplifiedAnnots.txt.gz", args[4])
stopifnot(file.exists(methAnnotFile))

# let's go!
log2exprs <- read.table(trans, row.names=1, header=TRUE, sep="\t")
Bvalues <- read.table(meth, row.names=1, header=TRUE, sep="\t") # M instead of B here? not really important for correlations...

log2exprs <- log2exprs[,colnames(log2exprs) %in% colnames(Bvalues)]
Bvalues <- Bvalues[,colnames(Bvalues) %in% colnames(log2exprs)]

log2exprs <- log2exprs[,order(colnames(log2exprs))]
Bvalues <- Bvalues[,order(colnames(Bvalues))]

cat("rnaseq data curated and reordered samples:", colnames(log2exprs), "\n")
cat("methyl data curated and reordered samples:", colnames(Bvalues), "\n")
stopifnot(colnames(log2exprs) == colnames(Bvalues))

# Annotation for transcritome with stable RefSeq transcripts (NM_XXX)
cat("annotating the transcriptome...\n")
# delete trID only transcripts/genes (with no ensembl_id linked)
log2exprs <- log2exprs[grep("ENS[GT][0-9]+", rownames(log2exprs)),] # cannot go with -grep("^trID", ...) as it could return integer(0) and lead to error
# split all transcript/gene ENST/ENSG annotations when multiple ENST/ENSG by row. Only one ENST per row, ENSG can be many, this complicates the approach. Cannot go by the easier grep "ENST/G" as well because with ENST they can be mixed (with their ENSG), which will lead to the wrong output
col_split <- lapply(strsplit(rownames(log2exprs), "\\|"), function(x) unlist(lapply(strsplit(x,"//"), "[", 2))) # the 'good' ensID is always located second in split
col_split <- lapply(col_split, function(x) if(length(x) == 0){x <- " "} else x) # in case some elements are empty. This should not happen here !
# and now finally duplicate rows with multiple annots
log2exprs <- data.frame(ensembl_id=unlist(col_split), apply(log2exprs, 2, function(x) rep(x, sapply(col_split, FUN=length))))

annotations_expression <- getBM(attributes=c(paste("ensembl", transLevel, "id", sep="_"), ifelse(transLevel == "transcript", "refseq_mrna", "external_gene_name"), "hgnc_symbol"), filters=paste("ensembl", transLevel, "id", sep="_"), values = unique(log2exprs$ensembl_id), mart=ensembl)
colnames(annotations_expression) <- c("ensembl_id", "refseq_id", "HGNC_symbol") # create the refseq_id column that will be used to merge later with methylome
log2exprs_annotated <- merge(annotations_expression, log2exprs, by="ensembl_id")
rm(log2exprs, annotations_expression)

# Annotation for methylation microarray. BEWARE: UCSC_RefGene_Group levels are 5'UTR, 3'UTR, TSS1500, TSS200, 1st Exon and Body. Simple quotes in 5'UTR and 3'UTR are interpreted with read.table default options as quotes separators so quote is set as "\"" in order to avoid bad importation
cat("annotating the methylome...\n")
annotations_methylation <- read.table(file=methAnnotFile, header=TRUE, stringsAsFactors=FALSE, comment.char="", quote="\"", sep="\t") # do not put cg as rownames !
colnames(annotations_methylation)[1] <- "methylProbe_id"
# split "UCSC_RefGene_Name" or "UCSC_RefGene_Accession" annotations and duplicate rows when multiple annots
col_split <- lapply(strsplit(as.character(annotations_methylation[,paste("UCSC_RefGene", ucscAnnot, sep="_")]),";"), unique) # way too much duplicated annots by row, make them unique
col_split <- lapply(col_split, function(x) if(length(x) == 0){x <- " "} else x) # for empty annots, yes here it happens !
annotations_methylation <- data.frame(refseq_id=unlist(col_split), apply(annotations_methylation, 2, function(x) rep(x, sapply(col_split, FUN=length)))) # refseq_id column created here
Bvalues_annotated <- merge(annotations_methylation, Bvalues, by.x="methylProbe_id", by.y=0) # merge by cg, for annotations it shall not be 0! (as they demultiply when splitting)
rm(Bvalues, annotations_methylation)

# Merge tables, by anything you want! By refseq_id column :)
cat("merging transcriptome and methylome big tables...\n")
merged_table <- merge(Bvalues_annotated, log2exprs_annotated, by="refseq_id", suffixes=c(".meth", ".tran"))
rm(Bvalues_annotated, log2exprs_annotated, col_split)

# Spearman's correlation between Bvalues and log2exprs from big_merged_annotated_table (wholeTable) and returning pvalues and Spearman's rho coefficients
corrMethTrans <- function(wholeTable){
  Bvalues_correlation <- as.numeric(wholeTable[grep("\\.meth$", names(wholeTable))])
  log2exprs_correlation <- as.numeric(wholeTable[grep("\\.tran$", names(wholeTable))])
  cor_spearman <- cor.test(Bvalues_correlation, log2exprs_correlation, method="spearman")
  results <- c(cor_spearman$p.value, cor_spearman$estimate)
  names(results) <- c("pvalue", "rho")
  return(results)
}

# Apply conditional rho aggregate on unique gene names for defined subsets of probes (promoter regions, gene bodies, ...) as found in subTable
# minRho: if > 0 the condition is TRUE when rho >= minRho ; if < 0 the condition is TRUE when rho <= minRho
# n : above condition must be TRUE for at least n CpG / gene
# computes about ~5 min with a 80k row subTable and ~8k unique genes
extractSigCorrs <- function(subTable, minRho = 0.7, n = 1){
  uniqueGenes <- unique(subTable$refseq_id)
  numSigPerGene <- if(minRho < 0){
     sapply(uniqueGenes, function(x){ nrow(subTable[subTable$refseq_id == x & subTable$rho <= minRho, ]) }) }
	 else { sapply(uniqueGenes, function(x){ nrow(subTable[subTable$refseq_id == x & subTable$rho >= minRho, ]) }) }
  names(numSigPerGene) <- uniqueGenes
  # return all CpGs with the significant region, they can easily be sorted after
  subTable[subTable$refseq_id %in% names(numSigPerGene[numSigPerGene >= n]),]
}

# big_merged_table: the whole annotated table with all samples [.meth for methylation status] [.tran for expression levels]
# rm_samp_pat: group(s) to exclude. Can represent a unique group (ie "C") or multiple groups "(WT|C)" or no group ""
correlate_with_promoters <- function(big_merged_table, rm_samp_pat){
  # Removing unwanted samples from big_merged_table
  cat(paste0("removing samples corresponding to general pattern -- ", rm_samp_pat))
  whole_table <- big_merged_table[, colnames(big_merged_table) %in% grep(paste0("^", rm_samp_pat, "_[A-Z]+[0-9]+\\.(meth|tran)$"), colnames(big_merged_table), value = TRUE, invert = TRUE)]
  
  # Fetching the remaining group names once the unwanted groups are removed
  groups_remaining <- unique(sub("_[A-Z]+[0-9]+\\.(meth|tran)$", "", grep("^[A-Za-z]+_[A-Z]+[0-9]+\\.(meth|tran)$", colnames(whole_table), value = TRUE)))
  
  # Compute means and stuff for these groups
  cat(paste0(" -- and compute some stats for the remaining groups ", paste(groups_remaining, collapse="-"), "... "))
  for (sample_group in groups_remaining){
    mean.meth <- apply(whole_table[, grep(paste0("^", sample_group, "_[A-Z]+[0-9]+\\.meth$"), colnames(whole_table))], 1, mean)
    mean.tran <- apply(whole_table[, grep(paste0("^", sample_group, "_[A-Z]+[0-9]+\\.tran$"), colnames(whole_table))], 1, function(x) mean(2^x))
    whole_table <- cbind(mean.meth, mean.tran, whole_table)
    colnames(whole_table)[1:2] <- c(paste0("mean_meth_",sample_group), paste0("mean_tran_",sample_group))
  }

  cat("done\ncomputing correlations (may take a little while)... ")
  # parallelize? maybe not needed, computation time is ~5 min for 500k rows
  whole_table <- cbind(t(as.data.frame(apply(whole_table, 1, corrMethTrans))), whole_table)
  # Conditional negative correlation aggregates for promoter related probes ~another 5 min
  cat("done\nextracting significant results (rho <= -0.33 for 3 CpGs or more) within promoter/TSS ranges... ")
  table_TSS <- subset(whole_table, UCSC_RefGene_Group == "TSS1500" | UCSC_RefGene_Group == "TSS200" | UCSC_RefGene_Group == "5'UTR" | UCSC_RefGene_Group == "1stExon")
  table_TSSneg <- extractSigCorrs(table_TSS, minRho = -1/3, n = 3)
  cat("done\nextracting significant results (rho <= -0.33 for 3 CpGs or more) within promoter/TSS ranges... ")
  table_TSSpos <- extractSigCorrs(table_TSS, minRho = 1/3, n = 3)
  cat("done\nextracting significant results (rho >= 0.33 for 3 CpGs or more) within promoter/TSS ranges... ")
  
  # Writing some tables of results
  cat("done\ndumping integration results in files... ")
  gz <- gzfile(paste("integrome", transLevel, paste(groups_remaining, collapse = "-"), "wholeResults.txt.gz", sep="_"), "w")
  write.table(whole_table, gz, quote=F, row.names=FALSE, sep="\t")
  close(gz)
  gz <- gzfile(paste("integrome", transLevel, paste(groups_remaining, collapse = "-"), "sigTSSpromoters3NEG.txt.gz", sep="_"), "w")
  write.table(table_TSSneg, gz, quote=F, row.names=FALSE, sep="\t")
  close(gz)
  gz <- gzfile(paste("integrome", transLevel, paste(groups_remaining, collapse = "-"), "sigTSSpromoters3POS.txt.gz", sep="_"), "w")
  write.table(table_TSSpos, gz, quote=F, row.names=FALSE, sep="\t")
  close(gz)
  cat(" done\n")
  rm(table_TSS, table_Body, whole_table)
}

## for testing only: Dlike with inter, then Clike with inter, then Clike with Dlike, and finally with everything (=without excluding anything)
### GIVE SAMPLE NAMES LIKE: (Clike|Dlike|inter)_[A-Z]+[0-9]+ (here especially) or anything like [A-Za-z]+_[A-Z]+[0-9]+ to be more general
for (pattern in c("Clike", "Dlike", "inter", "")) correlate_with_promoters(merged_table, pattern)

# Finished!
cat("all done\n")
