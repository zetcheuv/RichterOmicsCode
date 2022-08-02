#@seb2022 - sebastien.hergalant@inserm.fr
# 1) steps to follow to generate both FULL and EPIC final datasets from raw idat files
# QC - merging - cleaning - filtering - normalization - imputation - batch effects - evaluation - visualization ...
# 2) Some statistics, linear modelling, DMR search and functional normalization - as an example

###################################################
###################################################

# package requirements
library(minfi) # main dataset manipulations
require(minfiData)
library(DMRcate) # useful functions / DMR analyses 
library(RColorBrewer) # for MDS plots and PCA
library(sva) # batch effect assessment and memoval
library(impute) # missing values imputation
library(FactoMineR) # PCA
library(factoextra)
library(corrplot) #correlation plots and heat maps

# requirement: phenodata files need to be complete with all biological, clinical, statistical, omics... variables
# requirement: phenodata samples need to be presented in the same order as in the datasets
# FULL: complete dataset with merged EPIC and 450k samples 
# EPIC dataset: 850k (EPIC chip) - idat files in idats_850k/ directory
# BC450k (temporary) dataset: 450k chip - idat files in idats_450k/ directory 

options(width=250)
#load phenodata files
finalPHENO450k <- read.delim("phenodataWith450samples.txt", stringsAsFactors=F)
EPICpheno <- read.delim("phenodataWithEPICsamples.txt", stringsAsFactors=F)
FULLpheno <- rbind(finalPHENO450k, EPICpheno)

#load RGS datasets
BC450k_rgs <- read.metharray(paste0("idats_450k/", finalPHENO450k$Basename))
EPIC_rgs <- read.metharray(paste0("idats_850k/", finalPHENO850k$Basename), force=T) #force=T in case of multiple EPIC chip versions
FULL_rgs <- combineArrays(BC450k_rgs, EPIC_rgs, outType="IlluminaHumanMethylation450k")

# after this step the FULL and EPIC datasets basically undergo the same manipulations (with more things going on for FULL), so the following code will only focus on the FULL dataset
#generate full QC report
qcReport(FULL_rgs, pdf= "FULL_QC.pdf")
qcReport(EPIC_rgs, pdf="EPIC_QC.pdf")
#if samples or CpGs fail, have a bad profile, ... -> remove

#compute failed positions
detP <- detectionP(FULL_rgs)
failed <- detP>0.01
failedSamples <- which(colMeans(failed)>0.1) # if more than 10% CpGs fail, the sample fails
failedCpGs <- which(rowMeans(failed)>0.1) #if for than 10% samples fail, the CpG fails
failed <- failed[-failedCpGs, -failedSamples]

#compute methylSets and remove failed samples and CpGs
#RAW data
FULL_raw <- preprocessRaw(FULL_rgs)[-failedCpGs,-failedSamples]
#SWAN normalization
FULL_swan <- preprocessSWAN(FULL_rgs)[-failedCpGs, -failedSamples]

#for what follows, repeat the same for RAW
#beta-value matrices 
FULL_swanBeta <- getBeta(FULLplacenta_swan)
#replace remaining failed positions with missing values
FULL_swanBeta[failed] <- NA
#impute missing values with k nearest neighbors
FULL_swanBeta <- impute.knn(data=FULL_swanBeta, k=20, maxp=10000, rng.seed=123456789)$data # seed for reproducibility
#remove SNP, cross-hybridization and gonosome probes
FULL_swanBeta <- rmSNPandCH(FULLplacenta_swanBeta, dist=2, mafcut=0.05, rmcrosshyb=TRUE, rmXY=TRUE)
#step for removing batch effects between arrays
FULL_final <- ilogit2(ComBat(dat=logit2(FULL_swanBeta), batch=FULLpheno$ArrayBatch, ref.batch="450k"))

#correlation structure for duplicate samples after batch effect removal
FULL_corplot <- cor(FULL_final[,FULLpheno[FULLpheno$dupliPlot != "", "Sample_Name"]])
corrplot(FULL_corplot, order="hclust", method="square", type="lower", is.corr=F)

#MDS plots for SWAN and RAW data (with 1000, 10k, 100k and all CpGs, with any variable to estimate)
#below with 10k and variable ArrayType to explore batch effects on replicate samples for example
mdsPlot(dat=FULL_swanBeta, numPositions=10000, sampGroups=FULLpheno$ArrayType, pal=c("lightgray", colorRampPalette(brewer.pal(8, "Paired"))(40)), legendPos="bottomleft", legendNCol=10, pch=19)
mdsPlot(dat=FULL_final, numPositions=10000, sampGroups=FULLpheno$ArrayType, pal=c("lightgray", colorRampPalette(brewer.pal(8, "Paired"))(40)), legendPos="bottomleft", legendNCol=10, pch=19)
mdsPlot(dat=FULL_rawBeta, numPositions=10000, sampGroups=FULLpheno$ArrayBatch, pal=c("lightgray", colorRampPalette(brewer.pal(8, "Paired"))(40)), legendPos="bottomleft", legendNCol=10, pch=19)

#blood cell deconvolution
FULL_bloodCellCounts <- estimateCellCounts(rgSet=FULL_rgs, verbose=T)
FULLpheno <- cbind(FULLpheno, FULL_bloodCellCounts)

#gender prediction (for further gender check)
FULL_predictedGender <- as.data.frame(getSex(object=mapToGenome(FULL_raw)))
FULLpheno <- cbind(FULLpheno, FULL_predictedGender)

# average over replicates once replicate structure is estimated - remove bad replicates if needed
# need replicate sample name to be equal for this to work
FULLnorepl_final <- sapply(unique(colnames(FULL_final)), function(x) rowMeans(FULL_final[,grepl(x, colnames(FULL_final)), drop=F]))

# outliers and final explorations - given as an example with replicates and sample names
FULL_PCA <- PCA(t(FULL_final), ncp=10, graph=F, scale.unit=F)
fviz_pca_ind(FULL_PCA, geom.ind="text", pointsize=2.5, pointshape=21, addEllipses=F, axes=c(1,2), repel=T)

# final dataset & phenodata here
write.table(FULLpheno, "FULLphenodata_final.txt", sep="\t", quote=F, row.names=T, col.names=NA)
write.table(FULLnorepl_final, "FULLmethylation_final.txt", sep="\t", quote=F, row.names=T, col.names=NA)

###########################################
###########################################
###########################################

### And now for some statistics / linear modelling - DMRs - Some functional annotations

# packages - bioc libraries
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(liftOver))
suppressPackageStartupMessages(library(missMethyl))
suppressPackageStartupMessages(library(ExperimentHub))

# some functions
### basic statistics (mean/median/sd/diffmethyl): a bit long and forceful is computed on the whole chip (>450k)
# test.matrix must have probes as rownames and ordered samples: group1 (length=n1) then group2 (length=n2) 
add_stats <- function(methyl.bvalues, n1, n2, name1="1", name2="2"){
   stopifnot(ncol(methyl.bvalues) == n1+n2)
   mymat <- methyl.bvalues

   mean1 <- apply(mymat,1,function(x) mean(x[1:n1])) # if one wanted to apply mean to log2 values, this is also the way
   sd1 <- apply(mymat,1,function(x) sd(x[1:n1])) # for sd on log2 values I don't really know if this is more adapted
   mean2 <- apply(mymat,1,function(x) mean(x[(n1+1):(n1+n2)]))
   sd2 <- apply(mymat,1,function(x) sd(x[(n1+1):(n1+n2)]))
   foldchange <- mean1-mean2 # it is not really a fc here, more a diff in methylation status
   
   myres <- data.frame(row.names=rownames(mymat), round(foldchange,3), round(mean1,3), round(sd1,3), round(mean2,3), round(sd2,3))
   colnames(myres) <- c("diffMethyl",paste(name1,"_mean",sep=""),paste(name1,"_stdev",sep=""),paste(name2,"_mean",sep=""),paste(name2,"_stdev",sep=""))
   return(myres)
}
### this function extracts all useful stats and write the results in a file named [test_name]x"_stats.txt"
# [what.matrix] = Bvalues or Mvalues matrix, sample order matter!!! (coherence between that and fit$design)
# [cutoff] = max qvalue allowed - cut the crap we rarely want the whole 450k table anyway
# [adjust] = inherited from limma, topTable and the likes
# [contrast_coef] = column number in contrast.matrix to identify the test and extract the right samples from design and ultimately from [what]
# [ig1 and ig2] = indexes representing the sample positions in fit$design rows of interest or in what.matrix columns
# [groupName1 and 2]  = strings representing the group names as found in fit$contrasts and fit$design
# [annots.matrix] = annotation file, mandatory: probenames in first column
# [mod] = fit bayesian model - for convenience
generate_stats <- function(test_name="TESTvsCTRL", contrast_coef=1, adjust="BH", cutoff=0.2, what.matrix=Bvalues, annots.matrix=annotation, mod=fit, ig1=index1, ig2=index2, groupName1="1", groupName2="2"){
    stopifnot(contrast_coef > 0 && contrast_coef <= ncol(mod$contrasts))
    stopifnot(adjust %in% c("none", "BH", "BY", "holm", NULL))
    stopifnot(cutoff > 0 && cutoff <= 1)
    stopifnot(is.matrix(what.matrix))
    stopifnot(ncol(what.matrix) == nrow(mod$design))
    stopifnot(length(ig1)+length(ig2) <= ncol(what.matrix))    

    results.matrix <- topTable(mod, coef=contrast_coef, adjust=adjust, number=nrow(what.matrix), sort.by="p")
    write.table(results.matrix,paste(test_name,"_stats.whole-raw.txt",sep=""),quote=F,sep="\t",row.names=T,col.names=NA)
    results.matrix <- results.matrix[results.matrix$adj.P.Val <= cutoff,]
    write.table(results.matrix,paste(test_name,"_stats-fdrcut-raw.txt",sep=""),quote=F,sep="\t",row.names=T,col.names=NA)
    what.matrix <- what.matrix[rownames(results.matrix),, drop=FALSE] # will be cut accordingly, drop as to keep matrix str if it has only one row left

    results.matrix <- merge(results.matrix, add_stats(what.matrix[,c(ig1,ig2), drop=FALSE],n1=length(ig1),n2=length(ig2),name1=groupName1,name2=groupName2), by.x=0, by.y=0)
    results.matrix <- merge(results.matrix, annots.matrix, by.x=1, by.y=1)[,-c(2,3,4,7)] # remove logFC, aveExp, t-stat and B values

    results.matrix <- results.matrix[order(results.matrix$P.Value),]
    write.table(results.matrix,paste(test_name,"_stats.fdrcut-annotated-dup-splitonly.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)
}
#rewrite of this function from DMRcate, for hg38 better handle - updated in its new version, env handling
extractRanges <- function(dmrcoutput, genome = c("hg19", "hg38", "mm10")){
    genome <- match.arg(genome)
    if (!is(dmrcoutput, "DMResults")) { stop("Error: dmrcoutput is not a DMResults object. Please create one with dmrcate().") }
    coords <- extractCoords(dmrcoutput@coord)
    # before we had: coords <- cbind(coords, dmrcoutput$results[, c("no.cpgs", "minfdr", "Stouffer", "maxbetafc", "meanbetafc")])
    coords <- cbind(coords, dmrcoutput@no.cpgs, dmrcoutput@min_smoothed_fdr, dmrcoutput@Stouffer, dmrcoutput@HMFDR, dmrcoutput@Fisher, dmrcoutput@maxdiff, dmrcoutput@meandiff)
    coords$chromStart <- as.integer(as.character(coords$chromStart))
    coords$chromEnd <- as.integer(as.character(coords$chromEnd))
    ranges <- makeGRangesFromDataFrame(coords, keep.extra.columns = TRUE)

    # some work has been done here, first with xLiftOver, and now much much better with liftOver -- in its new version, DMRcate handles things differently from now on, this part has been adapted to fit my old requirements
    # chain infos are to be downloaded from UCSC -> genome -> download -> hg19 -> liftOver -> gunzip the file somewhere and import it to lift the ranges
    # CompressedGRangesList is created so there is a need to unlist that, but it may create some duplicated/split coordinates so remove the duplicates - no need anymore with new version, apparently
    eh = ExperimentHub()
    switch(genome, hg19 = { grt = eh[["EH3132"]] }, hg38 = { ranges = unlist(liftOver(ranges, import.chain("/home/chac/methyl/hg19ToHg38.over.chain")), use.names=F) ; grt = eh[["EH3134"]] }, mm10 = { grt = eh[["EH3136"]] })
    genesidx <- as.data.frame(findOverlaps(ranges, grt))
    genesover <- tapply(genesidx$subjectHits, genesidx$queryHits,
        function(x) grt$symbol[x])
    op.A <- sapply(genesover, function(l) paste(l, collapse = " ")) # changes here
    name.A <- names(genesover)
    m.A <- as.numeric(name.A)
    M <- length(ranges)
    overlapping.genes <- rep(NA_character_, M)
    overlapping.genes[m.A] <- op.A
    ranges$overlapping.genes <- overlapping.genes
    colnames(values(ranges)) <- sub("dmrcoutput@", "", colnames(values(ranges)))
    ranges
}

############ let's go!
#simplify var names a bit
Mvalues <- as.matrix(logit2(FULLnorepl_final))
Bvalues <- FULLnorepl_final
phenodata <- FULLpheno

# designs and contrasts over 3 columns of variables  provided as an example here
design.model <- model.matrix(~0 + phenodata$Sample_Group)
colnames(design.model) <- sub("^phenodata\\$Sample_Group", "", colnames(design.model))
design.model2 <- model.matrix(~0 + phenodata$Sample_Group2)
colnames(design.model2) <- sub("^phenodata\\$Sample_Group2", "", colnames(design.model2))
design.model3 <- model.matrix(~0 + phenodata$Sample_Group3)
colnames(design.model3) <- sub("^phenodata\\$Sample_Group3", "", colnames(design.model3))
contrast.matrix <- makeContrasts(CLL-Richter, DLBCL-Richter, CLL-DLBCL, levels=design.model)
contrast.matrix2 <- makeContrasts(U_CLL-M_CLL, U_Richter-M_Richter, U_CLL-U_Richter, M_CLL-M_Richter, levels=design.model2)
contrast.matrix3 <- makeContrasts(Richter_Clike-Richter_Dlike, DLBCL-Richter_Dlike, CLL-Richter_Clike, DLBCL-Richter_Inter, CLL-Richter_Inter, levels=design.model3)

#linear modelling with empirical Bayes
fitMod <- lmFit(Mvalues, design.model)
fitMod2 <- lmFit(Mvalues, design.model2)
fitMod3 <- lmFit(Mvalues, design.model3)
fitMod <- contrasts.fit(fitMod, contrast.matrix)
fitMod2 <- contrasts.fit(fitMod2, contrast.matrix2)
fitMod3 <- contrasts.fit(fitMod3, contrast.matrix3)
fitMod <- eBayes(fitMod)
fitMod2 <- eBayes(fitMod2)
fitMod3 <- eBayes(fitMod3)

# pvalues adjustments with Benjamini-Hochberg on the M-values and some added infos on the B-values + DMR search and annotations
lapply(list(fitMod, fitMod2, fitMod3), function(fm){ for(i in 1:ncol(fm$contrasts)){
   group1 <- rownames(fm$contrasts)[which(fm$contrasts[,i] == 1)]
   group2 <- rownames(fm$contrasts)[which(fm$contrasts[,i] == -1)]
   indgr1 <- which(fm$design[,group1] == 1)
   indgr2 <- which(fm$design[,group2] == 1)
   testgroupname <- paste0(length(indgr1), group1, "_vs_", length(indgr2), group2)

   # basic stats and stuff / cutoff is the fdr threshold for annotating the results 
   generate_stats(test_name=testgroupname, contrast_coef=i, adjust="BH", cutoff=annotcut, annots.matrix=annots, mod=fm, ig1=indgr1, ig2=indgr2, groupName1=group1, groupName2=group2)

   # DMR search / if fdrcut==0 then the whole thing will abort
   cpg_annots <- cpg.annotate("array", Mvalues, what="M", arraytype=chipType, analysis.type="differential", design=fm$design, contrasts=T, cont.matrix=fm$contrasts, coef=colnames(fm$contrasts)[i], fdr=fdrcut)

   if(length(which(cpg_annots@ranges@elementMetadata@listData$ind.fdr <= fdrcut)) < 5) next

   # demarcate regions the extract the coordinates, overlap with TSS, promoters and transcripts, and work a bit on the final tables
   DMRs <- dmrcate(cpg_annots, lambda=1000, C=2, min.cpgs=3)
   if(length(DMRs@coord) < 3) next

   # hg19 must be last because the hg19 related objects will be used again below the loop
   for(hg in c("hg38","hg19")){
      dmrRanges <- extractRanges(DMRs, genome = hg)

      # GO enrichment analysis via missMethyl - does not always work, problems with BiasedUrn when too few candidates
      tryCatch({
         for(db in c("GO", "KEGG")){
            #hypomethylated AND hypermethylated together
            FAranges <- goregion(dmrRanges[abs(dmrRanges$maxdiff) > 0.1,], all.cpg = rownames(Mvalues), collection = db, array.type = chipType)
            FAranges <- FAranges[order(FAranges$P.DE),]
            write.table(FAranges, paste(testgroupname,"_DMRs.",hg,"_both.FA-",db,".txt",sep=""), sep="\t", quote = F, row.names = T, col.names = NA)
            #hypomethylated alone
            FAranges <- goregion(dmrRanges[abs(dmrRanges$maxdiff) > 0.1 & dmrRanges$meandiff < 0,], all.cpg = rownames(Mvalues), collection = db, array.type = chipType)
            FAranges <- FAranges[order(FAranges$P.DE),]
            write.table(FAranges, paste(testgroupname,"_DMRs.",hg,"_hypo.FA-",db,".txt",sep=""), sep="\t", quote = F, row.names = T, col.names = NA)
            #hypermethylated alone
            FAranges <- goregion(dmrRanges[abs(dmrRanges$maxdiff) > 0.1 & dmrRanges$meandiff > 0,], all.cpg = rownames(Mvalues), collection = db, array.type = chipType)
            FAranges <- FAranges[order(FAranges$P.DE),]
            write.table(FAranges, paste(testgroupname,"_DMRs.",hg,"_hyper.FA-",db,".txt",sep=""), sep="\t", quote = F, row.names = T, col.names = NA)
         }
      }, error = function(e) { print("Error in FA for DMRs with missMethyl::goregion()") ; message(e) ; cat("\n") })

      dmrRanges <- dmrRanges[order(abs(dmrRanges$meandiff),decreasing=T),]
      dmrTable <- as.data.frame(dmrRanges)
      dmrTable$label <- paste0("DMR",1:nrow(dmrTable))
      write.table(dmrTable[,-5], paste(testgroupname,"_DMRs.",hg,".txt",sep=""), quote=F, sep="\t", row.names=F, col.names=T)
   }

} })







