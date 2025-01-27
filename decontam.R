#!/usr/bin/env Rscript
## Swift Biosciences 16S snapp workflow - Adapted by Bionformatics Service @ Genyo (Granada, Spain)
## Author (GitHub: jmgs7) 24102024

loadLibrary <- function(lib) {
    if (!require(lib, character.only = TRUE)) {
        message(paste0("Package ", lib, " not found, trying to install..."))
        install.packages(lib)
    if (!require(lib, character.only = TRUE)) {
        stop(paste0("Package ", lib, " not found."))
        }
    } 
}

splitDataFrame <- function(main.df, metadata.df, batch.column, output.dir=".", threads=0) {
    loadLibrary("parallel")
    
    if (threads == 0) {
        threads <- detectCores()
    }

    # Check if batch.column exists in metadata.df
    if (!batch.column %in% colnames(metadata.df)) {
        stop("The specified batch.column does not exist in the metadata dataframe.")
        }
    # Ensure rownames of main.df match rownames in metadata.df
    if (!all(rownames(main.df) %in% rownames(metadata.df))) {
        stop("Not all rownames in the main dataframe are present in the metadata dataframe.")
        }

    # Ensure main.df is a dataframe and not a matrix.
    if (is.matrix(main.df)) {
        main.df <- as.data.frame(main.df)
        message("Warning: main.df was inputted as a matrix. Converting to data.frame...")
    }

    # Extract the metadata column corresponding to rownames
    metadata_values <- metadata.df[rownames(main.df), batch.column]

    # Split the main dataframe by metadata column
    split.list <- split(main.df, metadata_values)

    if (output.dir != FALSE) {
        # Write each dataframe to a TSV file
        mclapply(names(split.list), function(name) {
            filename <- paste0(output.dir, "/", "decontamBatch_", name, ".tsv")
            write.table(split.list[[name]], file = filename, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
            }, mc.cores = threads)
    }

    return(split.list)
}

plotDecontamHist <- function(contam.df, output.file = "./decontam_hist.pdf", title = "Histogram of p values") {
    loadLibrary("ggplot2")

    ggplot(contam.df, aes(x = p)) +
        geom_histogram(bins = 20, fill = "steelblue", color = "black") +
        scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
        labs(title = title,
            x = "p statistic",
            y = "Frequency (# of ASVs)") +
        theme_classic()

    ggsave(output.file, width = 8, height = 6, device = "pdf")

}

runDecontam <- function(seqtab.nochim, neg.key = "Knegativo", metadata = NULL, neg.column = NULL, threshold=0.4, output.stats=FALSE, del.contaminants = TRUE) {
    loadLibrary("decontam")

    #You need to adjust the number of FALSES and TRUES and their order according to you sample distribution.
    if (!is.null(metadata)) {
        if (is.null(neg.column)) {
            message("WARNING: If metadata is provided, neg.column must be specified.")
            message(paste0("Trying to create decontam vector by detecting the neg.key (", neg.key, ") on sample names..."))
            vector.decontam <- grepl(neg.key, rownames(seqtab.nochim), ignore.case = TRUE) # TRUE is the negative control.
        } else {
            vector.decontam <- metadata[[neg.column]] == neg.key
        }
    } else {
        vector.decontam <- grepl(neg.key, rownames(seqtab.nochim), ignore.case = TRUE) # TRUE is the negative control.
    }

    contam.df <- isContaminant(as.matrix(seqtab.nochim), neg = vector.decontam, threshold=threshold) # Set a stricter threshold.
    contam_asvs <- row.names(contam.df[contam.df$contaminant == TRUE, ])

    if (output.stats != FALSE) {
        write.table(contam.df, file = output.stats, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
        plotName = gsub("\\.\\w+$", "", output.stats)
        plotDecontamHist(contam.df, paste0(plotName, "_plot.pdf"), title = basename(plotName))
    }

    if (del.contaminants) {
        seqtab.nochim.nocontam <- seqtab.nochim[,!colnames(seqtab.nochim) %in% contam_asvs]
        return(seqtab.nochim.nocontam)
    }
}

runDecontamBatch <- function(seqtab.nochim, metadata, batch.column="extraction_batch", neg.key="Knegativo", neg.column = NULL, threshold.vector=NULL, output.dir=".", del.contaminants=TRUE, threads=0) {

    batch.list <- splitDataFrame(seqtab.nochim, metadata, batch.column, output.dir, threads)
    metadata.list <- split(metadata, metadata[batch.column])
    if (is.null(threshold.vector)) {
        threshold.vector <- rep(0.4, length(batch.list))
    } else if (length(batch.list) != length(threshold.vector)) {
        stop("The length of the thresholds vector must match the number of decontam batches.")
    }

    loadLibrary("parallel")
    loadLibrary("data.table")

    if (threads == 0) {
        threads <- detectCores()
    }


    # Run decontam and deletes contam asvs or not in function of the del.contaminants argument.
    if (del.contaminants) {

        temp.df <- as.data.frame(
            rbindlist(
                mclapply(
                    seq_along(batch.list), 
                        function(i) {
                            seqtab.nochim <- batch.list[[i]]
                            metadata <- metadata.list[[i]]
                            threshold <- threshold.vector[i]
                            temp.df <- runDecontam(seqtab.nochim, neg.key, metadata, neg.column, threshold, paste0(output.dir, "/", "decontamBatch_", names(batch.list)[i], "_stats.tsv"), del.contaminants)
                            temp.df$id <- rownames(temp.df)
                            return(temp.df)
                            }, 
                        mc.cores = threads), # Number of cores for mcapply.
                fill = TRUE) # rbindlist: Fill missing columns after decontamination with NAs to later substitute with 0.
            ) # It outputs a data table which does not accept rownames. We have to store the rownames in a column and transform to dataframe.

        # Reset the rownames.
        rownames(temp.df) <- temp.df$id
        temp.df <- temp.df[, !colnames(temp.df) %in% c("id")]

        # Handle NAs and all-0 columns.
        temp.df[is.na(temp.df)] <- 0 # Replace NAs with 0
        seqtab.nochim.nocontam <- temp.df[, !(colSums(temp.df) == 0)] # Remove columns with all zeros

        return(seqtab.nochim.nocontam)

    } else {

        mclapply(
                seq_along(batch.list), 
                    function(i) {
                        seqtab.nochim <- batch.list[[i]]
                        threshold <- threshold.vector[i]
                        runDecontam(seqtab.nochim = seqtab.nochim, neg.key = neg.key, threshold = threshold, output.stats = paste0(output.dir, "/", "decontamBatch_", names(batch.list)[i], "_stats.tsv"), del.contaminants = del.contaminants)
                        }, 
            mc.cores = threads) # Number of cores for mcapply.

    }
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