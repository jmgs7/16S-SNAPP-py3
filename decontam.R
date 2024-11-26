#!/usr/bin/env Rscript
## Swift Biosciences 16S snapp workflow - Adapted by Bionformatics Service @ Genyo (Granada, Spain)
## Author (GitHub: jmgs7) 24102024

split_dataframe_by_metadata <- function(main_df, metadata_df, column_name, ouput_dir=".", threads=0) {

    suppressMessages(library(parallel))
    
    if (threads == 0) {
        threads <- detectCores()
    }

    # Check if column_name exists in metadata_df
    if (!column_name %in% colnames(metadata_df)) {
        stop("The specified column_name does not exist in the metadata dataframe.")
        }
    # Ensure rownames of main_df match rownames in metadata_df
    if (!all(rownames(main_df) %in% rownames(metadata_df))) {
        stop("Not all rownames in the main dataframe are present in the metadata dataframe.")
        }

    # Ensure main_df is a dataframe and not a matrix.
    if (is.matrix(main_df)) {
        main_df <- as.data.frame(main_df)
        message("Warning: main_df was input as a matrix. Converting to data.frame...")
    }

    # Extract the metadata column corresponding to rownames
    metadata_values <- metadata_df[rownames(main_df), column_name]

    # Split the main dataframe by metadata column
    split_list <- split(main_df, metadata_values)

    if (ouput_dir != FALSE) {
        # Write each dataframe to a TSV file
        mclapply(names(split_list), function(name) {
            filename <- paste0(ouput_dir, "/", "decontam_batch_", name, ".tsv")
            write.table(split_list[[name]], file = filename, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
            }, mc.cores = threads)
    }

    return(split_list)
}

run_decontam <- function(seqtab.nochim, neg_key = "Knegativo", threshold=0.4) {
    suppressMessages(library(decontam))

    #You need to adjust the number of FALSES and TRUES and their order according to you sample distribution.
    vector_for_decontam <-  grepl(neg_key, rownames(seqtab.nochim), ignore.case = TRUE) # TRUE is the negative control.
    contam_df <- isContaminant(as.matrix(seqtab.nochim), neg = vector_for_decontam, threshold=threshold) # Set a stricter threshold.
    contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])
    seqtab.nochim.nocontam <- seqtab.nochim[,!colnames(seqtab.nochim) %in% contam_asvs]
    return(seqtab.nochim.nocontam)
}

run_decontam_batches <- function(seqtab.nochim, metadata, column_name="extraction_batch", neg_key="Knegativo", threshold_vector=NULL, ouput_dir=".", threads=0) {

    batch_list <- split_dataframe_by_metadata(seqtab.nochim, metadata, column_name, ouput_dir, threads)

    if (is.null(threshold_vector)) {
        threshold_vector <- rep(0.4, length(batch_list))
    } else if (length(batch_list) != length(threshold_vector)) {
        stop("The length of the thresholds vector must match the number of dataframes in the list.")
    }

    suppressMessages(library(parallel))
    suppressMessages(library(data.table))

    if (threads == 0) {
        threads <- detectCores()
    }

    seqtab.nochim.nocontam <- rbindlist(mclapply(seq_along(batch_list), function(i) {
        seqtab.nochim <- batch_list[[i]]
        threshold <- threshold_vector[i]
        run_decontam(seqtab.nochim, neg_key, threshold)
    }, mc.cores = threads))

    return(seqtab.nochim.nocontam)
}

decontam_stats <- function(seqtab.nochim, seqtab.nochim.nocontam) {
    message(
        paste0(
            "Total inferred ASVs: ", ncol(seqtab.nochim), 
            "\nTotal contam ASVs: ", (ncol(seqtab.nochim) - ncol(seqtab.nochim.nocontam)), 
            "\nTotal non-contam ASVs: ", ncol(seqtab.nochim.nocontam), 
            " (", round(ncol(seqtab.nochim.nocontam)/ncol(seqtab.nochim)*100, 2), "%)"
            )
    ) 
}