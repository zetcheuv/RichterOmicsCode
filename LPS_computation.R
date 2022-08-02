#@seb2022 sebastien.hergalant@inserm.fr
#workflow for computing LPS (linear predictor scores) from methylation data

#requirement1: normalized/filtered/prepared methylome beta-values table with CpG probes in rows and samples in columns
#requirement2: the 4863 CpG signature file with associated t-statistics for weighting the scores

# FULL_final: a matrix or data frame (see RS_methylome_datasets_and_stats.R to obtain it)
FULL_final <- as.data.frame(FULL_final)
# scoringSig4863: a table with 4863 CpG probe names in first column and weighted-t statistics in second 
scoringSig4863 <- read.delim("LPS_4863CpGsig.txt", stringsAsFactors=F)

# reduce to shared CpGs between objects
FULL_scoring <- FULL_final[rownames(FULL_final) %in% scoringSig4863[,1],]
scoringSig <- scoringSig4863[scoringSig4863[,1] %in% rownames(FULL_scoring),]
# ensure that all CpGs are ordered the same in both objects
FULL_scoring <- FULL_scoring[scoringSig[,1],]

# now compute LPS, in aggregated (for complete 4863 CpGs when present), then average mode (if CpGs are missing in some datasets)
scoresSUM <- apply(FULL_scoring, 2, function(x) sum(x*scoringSig[,2]))
scoresMEAN <- apply(FULL_scoring, 2, function(x) mean(x*scoringSig[,2]))

# mean and standart deviation scores for DLBCL and CLL spaces, in MEAN mode rather than SUM as it is fully comparable with other datasets if all 4863 CpGs are not all present
# these statistics come from LPS that were pre-computed on a 68 DLBCL and 215 CLL cohort
DLBCL_score_mean <- -8.085893
DLBCL_score_sd <- 1.426302
CLL_score_mean <- -2.074467
CLL_score_sd <- 0.5114767

# for each RS sample, compute the P-value to classify in DLBCL space p(S in DLBCL), ie to be "DLBCL-like RS" 
DLBCLlike_pvalue <- sapply(sort(scoresMEAN), function(x) dnorm(x, DLBCL_score_mean, DLBCL_score_sd) / (dnorm(x, DLBCL_score_mean, DLBCL_score_sd) + dnorm(x, CLL_score_mean, CLL_score_sd)))
# alternatively, one can compute the 1 - p(S in DLBCL), ie for each RS sample, compute the associated P-value to be "highCLL-derived RS" 
CLLderived_pvalue <- sapply(sort(scoresMEAN), function(x) dnorm(x, CLL_score_mean, CLL_score_sd) / (dnorm(x, DLBCL_score_mean, DLBCL_score_sd) + dnorm(x, CLL_score_mean, CLL_score_sd)))
