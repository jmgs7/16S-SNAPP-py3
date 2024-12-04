#!/usr/bin/env Rscript
## Swift Biosciences 16S snapp workflow - Adapted by Bionformatics Service @ Genyo (Granada, Spain)
## Author (GitHub: jmgs7) 24102024
splitDataFrame <- function(main.df, metadata.df, column.name, output.dir=".", threads=0) {

    suppressMessages(library(parallel))
    
    if (threads == 0) {
        threads <- detectCores()
    }

    # Check if column.name exists in metadata.df
    if (!column.name %in% colnames(metadata.df)) {
        stop("The specified column.name does not exist in the metadata dataframe.")
        }
    # Ensure rownames of main.df match rownames in metadata.df
    if (!all(rownames(main.df) %in% rownames(metadata.df))) {
        stop("Not all rownames in the main dataframe are present in the metadata dataframe.")
        }

    # Ensure main.df is a dataframe and not a matrix.
    if (is.matrix(main.df)) {
        main.df <- as.data.frame(main.df)
        message("Warning: main.df was input as a matrix. Converting to data.frame...")
    }

    # Extract the metadata column corresponding to rownames
    metadata_values <- metadata.df[rownames(main.df), column.name]

    # Split the main dataframe by metadata column
    split.list <- split(main.df, metadata_values)

    if (output.dir != FALSE) {
        # Write each dataframe to a TSV file
        mclapply(names(split.list), function(name) {
            filename <- paste0(output.dir, "/", "decontam_batch_", name, ".tsv")
            write.table(split.list[[name]], file = filename, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
            }, mc.cores = threads)
    }

    return(split.list)
}

runDecontam <- function(seqtab.nochim, neg.key = "Knegativo", threshold=0.4, output.stats=FALSE) {
    suppressMessages(library(decontam))

    #You need to adjust the number of FALSES and TRUES and their order according to you sample distribution.
    vector.decontam <- grepl(neg.key, rownames(seqtab.nochim), ignore.case = TRUE) # TRUE is the negative control.
    contam.df <- isContaminant(as.matrix(seqtab.nochim), neg = vector.decontam, threshold=threshold) # Set a stricter threshold.
    contam_asvs <- row.names(contam.df[contam.df$contaminant == TRUE, ])
    seqtab.nochim.nocontam <- seqtab.nochim[,!colnames(seqtab.nochim) %in% contam_asvs]

    if (output.stats != FALSE) {
        write.table(contam.df, file=output.stats, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
    }

    return(seqtab.nochim.nocontam)
}

runDecontamBatch <- function(seqtab.nochim, metadata, column.name="extraction_batch", neg.key="Knegativo", threshold.vector=NULL, output.dir=".", threads=0) {

    batch.list <- splitDataFrame(seqtab.nochim, metadata, column.name, output.dir, threads)

    if (is.null(threshold.vector)) {
        threshold.vector <- rep(0.4, length(batch.list))
    } else if (length(batch.list) != length(threshold.vector)) {
        stop("The length of the thresholds vector must match the number of dataframes in the list.")
    }

    suppressMessages(library(parallel))
    suppressMessages(library(data.table))

    if (threads == 0) {
        threads <- detectCores()
    }

    temp.df <- as.data.frame(rbindlist(mclapply(seq_along(batch.list), function(i) {
        seqtab.nochim <- batch.list[[i]]
        threshold <- threshold.vector[i]
        temp.df <- runDecontam(seqtab.nochim, neg.key, threshold, paste0(output.dir, "/", "decontam_stats_", names(batch.list)[i], ".tsv"))
        temp.df$id <- rownames(temp.df)
        return(temp.df)
    }, mc.cores = threads))) # It outputs a data table which does not accept rownames. We have to store the rownames in a columna and transform to dataframe.
    
    rownames(temp.df) <- temp.df$id
    seqtab.nochim.nocontam <- temp.df[, !colnames(temp.df) %in% c("id")]

    return(seqtab.nochim.nocontam)
}

decontamStats <- function(seqtab.nochim, seqtab.nochim.nocontam) {
    message(
        paste0(
            "Total inferred ASVs: ", ncol(seqtab.nochim), 
            "\nTotal contam ASVs: ", (ncol(seqtab.nochim) - ncol(seqtab.nochim.nocontam)), 
            "\nTotal non-contam ASVs: ", ncol(seqtab.nochim.nocontam), 
            " (", round(ncol(seqtab.nochim.nocontam)/ncol(seqtab.nochim)*100, 2), "%)"
            )
    ) 
}