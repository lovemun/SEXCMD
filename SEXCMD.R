#!/bin/env Rscript

## Total read count (EXOME based)
TOTAL_READ_COUNT <- 1E6

## Coefficients by Sequencing type
WES_COEF <- c(1, 1)
WGS_COEF <- c(1, 150)
RNA_COEF <- c(1, 5)

## Argument handle
args <- commandArgs(trailingOnly = TRUE)
usage_message <- "USAGE: Rscript SEXCMD.R SEX_Marker.fasta total_read_counts Seq_type[1/2/3] [XY/ZW] input.fastq.gz
    Sequencing Type 1:Whole Exome Seq, 2:RNA-Seq, 3:Whole Genome Seq, Sex type : XY or ZW"

if (length(args)!= 3) {
  write(usage_message, stderr())
  quit(status = 1)
}
marker_file <- as.character(args[1])
type_flag <- as.numeric(args[2])
sextype <- as.character(args[3])
input_fastq <- as.character(args[4])

## check input files
check_input <- function(input_file){
  if (!file.exists(input_file)) {
    write(usage_message, stderr())
    write(paste("ERROR: Can not find", input_file), stderr())
    quit(status = 1)
  }
}
check_input(marker_file)
check_input(input_fastq)

## check sequencing type flag
if (type_flag==1){
  COEF <- WES_COEF
} else if (type_flag==2) {
  COEF <- RNA_COEF
} else if (type_flag==3) {
  COEF <- WGS_COEF
} else {
  write(usage_message, stderr())
  quit(status = 1)
}
TOTAL_READ_COUNT <- as.integer(TOTAL_READ_COUNT*COEF[2])

## check BWA, samtools and index marker sequence
indexing_command <- paste("bwa index -a is", marker_file)
bwt_file <- paste(marker_file,".bwt", sep="")
if (Sys.which('bwa')=="" | Sys.which('samtools')==""){
  write("ERROR: Can not find bwa or samtools in $PATH", stderr())
  quit(status = 1)
} else if (!file.exists(bwt_file) | file.info(bwt_file)$size==0){
  if (system(indexing_command)==0){
    write(paste("LOG!!:", indexing_command, "has been completed."), stderr())
  } else {
    write(paste("ERROR!!:", indexing_command, "has been failed."), stderr())
    quit(status = 1)
  }
}

## check bwa mem
#mem_flag <- system(paste("bash -c 'bwa mem", marker_file,
#				"<(gzip -dc",input_fastq,"|head -n 4) &> /dev/null'"))
#if (mem_flag!=0) {
#	write("ERROR: unrecognized command 'mem'. Update bwa over v0.7.", stderr())
#	quit(status = 1)
#}

## required R packages
if (!require("e1071", character.only = TRUE)){
  install.packages("e1071", dep=TRUE, repos='http://cran.rstudio.com/')
  library("e1071")
}
if (!require("seqinr", character.only = TRUE)){
  install.packages("seqinr", dep=TRUE, repos='http://cran.rstudio.com/')
  library("seqinr")
}
library(parallel)
cpu.thres <- detectCores()/2

## calculate stat and count
stat_file <- paste(input_fastq, ".stat", sep='')
qt_stat_file <- paste('"',stat_file,'"', sep='')
count_file <- paste(input_fastq, ".count", sep='')
calc_cmd <- paste("gzip -dc",input_fastq,"|head -n",format(TOTAL_READ_COUNT*4, scientific=F),
                  "|awk '{print substr($0,1,101)}'|awk '{if (NR%4==2 && $0!~/N/) {sum+=length;i+=1};print $0}",
                  "END{print sum/i,i >",qt_stat_file,"}'|bwa mem -t 1",marker_file,"-",
                  "|samtools view -S -F 4 -q 30 -|cut -f3|sort|uniq -c|awk '{print $2\"\t\"$1}' >",count_file)
system(calc_cmd)
write(paste("LOG!!: Read mapping and count has been completed."), stderr())

stat <- read.table(stat_file, header=F)
avg_read_len <- stat[[1]]
total_read_count <- stat[[2]]

## check minimum total count
err_msg <- "ERROR!!: Total mapped read count is below 3. 
Increase input data size or TOTAL_READ_COUNT value."
tryCatch({
  count <- read.table(count_file, header=F, as.is=T)
  if (sum(count$V2) <= 3){
    write(err_msg, stderr())
    quit(status = 1)
  }
}, error = function(e) {
  write(err_msg, stderr())
  quit(status = 1)		
}) 

## build count table
marker_seq_list <- read.fasta(marker_file, set.attributes=F)
marker_length_list <- lapply(marker_seq_list, length)
marker_seq_list[] <- 0
marker_seq_list[count$V1] <- count$V2
ct <- data.frame(Length=do.call(rbind, marker_length_list[startsWith(x = names(marker_seq_list), prefix = "chrX")]), Count=do.call(rbind, marker_seq_list[startsWith(x = names(marker_seq_list), prefix = "chrX")]))
ct <- cbind(ct, do.call(rbind, marker_seq_list[startsWith(x = names(marker_seq_list), prefix = "chrY")]))
colnames(ct) <- c("Length", "Xcount", "Ycount")

## gender prediction
rs <- sum(ct$Ycount)/sum(ct$Xcount)
if (sextype = "XY"){
  if (rs < 0.2){
    sex <- "F"
  } else {
    sex <- "M"
  }
  ct <- rbind(ct, c(paste0("Ratio of X and Y counts is ", rs), paste0("Sex of this sample ", input_fastq, " is ", sex), ""))
}
if (sextype = "ZW"){
  if (rs < 0.2){
    sex <- "M"
  } else {
    sex <- "F"
  }
  ct <- rbind(ct, c(paste0("Ratio of Z and W counts is ", rs), paste0("Sex of this sample ", input_fastq, " is ", sex), ""))
}
rownames(ct)[nrow(ct)] <- "Sex_Determination"

# Save the results
result_file <- paste(input_fastq, ".OUTPUT", sep='')
write.table(ct, result_file, quote=F, col.names=NA, sep="\t")

## make clean
if (all(file.remove(stat_file, count_file))) {
  write(paste("LOG!!: SEXCMD calculation result has been saved on ",result_file,".",sep=''),stderr())	
}



