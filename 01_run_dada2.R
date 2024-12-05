#!/usr/bin/env Rscript
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

## run dada2 for 16S amplicon sequences
#Supply the arguments from the command line

args <- commandArgs(TRUE)
if (length(args) < 8){
    stop("Usage: run_dada2.R inputDir workdir QCLIB SCRIPTS THREADS DECONTAM_THRESHOLDS DECONTAM_COLUMN DECONTAM_NEGATIVE METADATA(optional)");
}

path <- normalizePath(args[1]) #inputDir
wd <- normalizePath(args[2]) #workdir
QCLIB <- normalizePath(args[3]) #the path to qc_analysis.r (metagenomics-pipeline)
SCRIPTS <- normalizePath(args[4])
DECONTAM_THRESHOLDS <- args[5]
DECONTAM_COLUMN  <- args[6]
DECONTAM_NEGATIVE <- args[7]
THREADS <- as.numeric(args[8]) #threads to use in the pipeline
METADATA <- args[9]

if (DECONTAM_COLUMN != "False") {
    if (!is.na(METADATA)) {
        metadata <- read.table(normalizePath(args[9]), sep = "\t", header = TRUE, row.names = 1)
    } else {
        stop("You must provide a valid metadata file.")
    }
} else if (DECONTAM_COLUMN == "False" & DECONTAM_THRESHOLDS != "False") {
    stop("You must provide a valid column name for decontamination.")
}

suppressMessages(library(dada2))
suppressMessages(source(QCLIB))
suppressMessages(source(paste0(SCRIPTS, "/decontam.R")))

fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_L001"), `[`, 1)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# create samplesheet
sample_sheet <- data.frame(
    sampleID = sample.names,
    forwardReads = fnFs,
    reverseReads = fnRs
)

# qc and trunc estimation
raw_fqc <- fastq_quality_control(sample_sheet, wd, raw = TRUE)
trunc_parameters <- lapply(raw_fqc, trunc_estimation, min_qual = 30)
message(paste0("Trunc parameters:", "\nFW: ", trunc_parameters$forwardReads, "\nRV: ", trunc_parameters$reverseReads))

# launch filterAndTrim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(trunc_parameters$forwardReads, trunc_parameters$reverseReads), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=THREADS)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#Sample inference with standard all sample pooling
dadaFs <- dada(filtFs, err=errF, multithread=THREADS, pool=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=THREADS, pool=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, justConcatenate=TRUE)

#get asv sequence and count table and write to a tab-delimited file
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=THREADS, verbose=TRUE)
#Write seqtab.nochim
write.table(seqtab.nochim, file=paste(wd, 'seqtab_nochim.tsv', sep='/'), sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

# Decontamination
output.dir <- file.path(wd, "decontam")
if (DECONTAM_THRESHOLDS != "False" & DECONTAM_THRESHOLDS != "default") {

    dir.create(output.dir, showWarnings = FALSE)

    threshold.vector <- as.numeric(unlist(strsplit(DECONTAM_THRESHOLDS, ",")))

    seqtab.nochim.nocontam <- runDecontamBatch(seqtab.nochim, metadata, column.name=DECONTAM_COLUMN, neg.key=DECONTAM_NEGATIVE, threshold.vector=threshold.vector, output.dir=output.dir, threads=THREADS)

    decontamStats(seqtab.nochim, seqtab.nochim.nocontam)

} else if (DECONTAM_THRESHOLDS == "default") {

    dir.create(output.dir, showWarnings = FALSE)

    seqtab.nochim.nocontam <- runDecontamBatch(seqtab.nochim, metadata, column.name=DECONTAM_COLUMN, neg.key=DECONTAM_NEGATIVE, output.dir=output.dir,threads=THREADS)

    decontamStats(seqtab.nochim, seqtab.nochim.nocontam)

} else {

    message("No decontamination will be performed.")
    seqtab.nochim.nocontam <- seqtab.nochim

    if (DECONTAM_COLUMN != "False") {
        dir.create(output.dir, showWarnings = FALSE)
        runDecontamBatch(seqtab.nochim, metadata, column.name=DECONTAM_COLUMN, neg.key=DECONTAM_NEGATIVE, output.dir=output.dir, threads=THREADS, del.contaminants=FALSE)
    }
}

#Write asv sequence and count table
write.csv(t(seqtab.nochim.nocontam), file=paste(wd, 'asv_seqNcount.csv', sep='/'), quote = FALSE)

#Get process stats and write to a tab-delimited file
#getN <- function(x) sum(getUniques(x))
#track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
#colnames(track) <- c("Primer_trimmed", "Filtered", "Denoised_FWD", "Denoised_REV", "Merged", "Non_chimera")
#rownames(track) <- sample.names
#write.table(track, file=paste(wd, 'DADA2_summary.csv', sep="/"), sep=" ") #save processing stats to a file