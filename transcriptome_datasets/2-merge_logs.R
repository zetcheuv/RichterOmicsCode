#!/usr/bin/Rscript
# @seb2016 (sebastien.hergalant@inserm.fr)
# usage: merge_logs.R "path/files.pattern"
# !!! protect path/files.pattern inside "" to avoid pattern from being interpreted by the shell !!!
# files.pattern is a true REGEXP!! think of replacing * by .* and regular dots (.) by \\.  etc...
# numbers in pointed files must be formed with dots (.), not commas (,)
# median between numbers is computed

args <- commandArgs(trailingOnly = TRUE)
if(is.na(args[1])){
   filepattern <- "*.log"
   filedir <- "."
} else{
   filedir <- dirname(args[1])
   filepattern <- basename(args[1])
}

logfiles <- list.files(path=filedir,pattern=filepattern,full.names=TRUE)
loglist <-lapply(logfiles,readLines,n=-1)
splitlist <- lapply(loglist,function(x) strsplit(gsub("([0-9]+(.[0-9]+)?)","~\\1~",x),"~"))

resultlog <- splitlist[[1]]

for(i in 1:length(resultlog)){
   for(j in 1:length(resultlog[[i]])){
      suppressWarnings(res <- round(mean(as.numeric(lapply(splitlist,function(x) x[[i]][j]))),2))
      if(!is.na(res)) resultlog[[i]][j] <- res
   }
}

writeLines(unlist(lapply(resultlog,paste,collapse="")))

