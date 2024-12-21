#!/bin/bash
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

#Run blastn three times with different DBs and queries to make increasingly
#more complete sets of matches between reads and reference.
# The first and second blasn allow to filter the RDP DB to include only sequences that are likely to be found in the samples.
# The third one is the one that will be used to determine which of the reference sequences could be considered to be present in the sample and be used as a bona fide full 16S sequence.

##Format a blastDB of asv for second blast in order to capture a more complete
#set of matches

if ! compgen -G "asv.*" > /dev/null; then

    echo -e "\nMaking ASV blast DB...\n    Starts: $(date)">>$log
    start=$(date +%s.%N)
    $MAKEBLASTDB \
        -dbtype nucl \
        -in asv_uniq.fasta \
        -out asv

    echo "    Ends: $(date)">>$log
    end=$(date +%s.%N)
    runtime=$(python -c "print(${end} - ${start})")
    echo "    Making ASV blastDB Runtime: $runtime sec" >> $log

fi


##Blast the pre-formatted reference DB to find the hits above the cutoff
# In this first blast step, we use the full reference database from RDP and filter the templates that have a high quality hit.
if [ ! -f reflist.txt ]; then

    echo -e "\nRunning first blast...\n    Starts: $(date)">>$log
    start=$(date +%s.%N)
    # Output format: 6 (BLAST standard output format), std (Standard output), include: qlen (query length), slen (subject length).
    $BLASTN \
        -num_threads $THREADS \
        -db ${RDP_FULL_BLAST} \
        -query asv_uniq.fasta \
        -out blastn_1.txt \
        -outfmt "6 std qlen slen" \
        -perc_identity 97 \
        -qcov_hsp_perc 97 # % identity and query coverage cutoffs

    echo "    Ends: $(date)">>$log
    end=$(date +%s.%N)
    runtime=$(python -c "print(${end} - ${start})")
    echo "    First Runtime: $runtime sec" >> $log

    #obtain top hits of each query sequence from all samples
    echo -e "\nFiltering the first blast...\n    Starts: $(date)">>$log
    start=$(date +%s.%N)
    # Filters by perc_identity and qcov_hsp_perc
    cat blastn_1.txt \
        | awk -F '\t' '{if($3 >= 99 && ($4/$13) > 0.99) print $2}' \
        | sort -u \
        > reflist.txt # Reflist contains the IDs from the FILTERED matched reference sequences.
    template_count=$(< reflist.txt wc -l)

fi


# Filters the fastq of the refset to only include the filtered matched reference sequences.
# This file is used in the next step. The objective of the first two blastn is to "enrich" the reference set to include only sequences that are likely to be found in the samples.
if [ ! -f exp_refset.fasta ]; then

    cut -f2 blastn_1.txt | sort -u > exp_reflist.txt # exp_reflist contains the IDs from ALL the matched reference sequences..
    ${RDPHOME}/ReadSeq \
        select-seqs exp_reflist.txt  \
        exp_refset.fasta fasta Y ${RDP_FULL_SEQ}

fi


# Filters the fastq of the refset to only include the matched reference sequences (no filtering).
if [ ! -f seqmatch/train_*.trainee ]; then

    ${RDPHOME}/ReadSeq \
        select-seqs reflist.txt \
        refset.fasta fasta Y exp_refset.fasta
    echo "    Ends: $(date)">>$log
    end=$(date +%s.%N)
    runtime=$(python -c "print(${end} - ${start})")
    echo -e "\nExtracting ${template_count} preliminary candidate template sequences\n    Starts: $(date)">>$log
    echo "    Completed parsing first blast. Runtime: $runtime sec" >> $log

fi 


# Prepare SeqMatch DBecho -e "\nFormatting SeqMatch DB ...\n    Starts: $(date)">>$log
# This is a clustering step.
if ! compgen -G "blastn_2_derep.*" > /dev/null; then

    start=$(date +%s.%N)
    [  -d seqmatch ] && echo "Directory seqmatch does exist, please delete..." && exit
    mkdir seqmatch
    ${RDPHOME}/SequenceMatch train exp_refset.fasta seqmatch/train # Trains with all the matched reference sequences.
    echo "    Ends: $(date)">>$log
    end=$(date +%s.%N)
    runtime=$(python -c "print(${end} - ${start})")
    echo "    Completed training SeqMatch. Runtime: $runtime sec" >> $log

    # Run blast using preliminary candidate templetes against asv DB in order to obtain a
    #more complete set of template-asv matches
    echo -e "\nRunning second blast...\n    Starts: $(date)">>$log
    start=$(date +%s.%N)
    $BLASTN \
        -num_threads $THREADS \
        -db asv \
        -max_target_seqs 2000 \
        -query refset.fasta \
        -out blastn_2.txt \
        -outfmt "6 std qlen slen" \
        -perc_identity 99

    echo "    Ends: $(date)">>$log
    end=$(date +%s.%N)
    runtime=$(python -c "print(${end} - ${start})")
    echo "    Second blast Runtime: $runtime sec" >> $log

    cat blastn_2.txt \
        | awk -F '\t' '{if($3 >= 97 && ($4/$14) > 0.97) print $0}' \
        > blastn_2_filtered.txt

    echo -e "\nRunning dereplicating blast_2...\n    Starts: $(date)">>$log
    start=$(date +%s.%N)
    ${SCRIPTS}/derep_hitsets.py blastn_2_filtered.txt \
        > blastn_2_derep_IDs.txt
    echo "    Ends: $(date)">>$log
    end=$(date +%s.%N)
    runtime=$(python -c "print(${end} - ${start})")
    echo "    Dereplicate Second blast Runtime: $runtime sec" >> $log

    echo -e "\nRunning extracting dereplicated blastn_2...\n    Starts: $(date)">>$log
    start=$(date +%s.%N)
    ${RDPHOME}/ReadSeq \
        select-seqs blastn_2_derep_IDs.txt \
        blastn_2_derep.fasta fasta Y refset.fasta
    echo "    Ends: $(date)">>$log
    end=$(date +%s.%N)
    runtime=$(python -c "print(${end} - ${start})")
    echo "    Extracted dereplicated blast_2 hits Runtime: $runtime sec" >> $log

    echo -e "\nRunning making blastDB with dereplicated blastn_2 hits ...\n    Starts: $(date)">>$log
    start=$(date +%s.%N)
    $MAKEBLASTDB \
        -dbtype nucl \
        -in blastn_2_derep.fasta \
        -out blastn_2_derep
    echo "    Ends: $(date)">>$log
    end=$(date +%s.%N)
    runtime=$(python -c "print(${end} - ${start})")
    echo "    Making blastDB with dereplicated blastn_2 hits Runtime: $runtime sec" >> $log

fi


# Third blastn
if [ ! -d asv_tmp ]; then

    echo -e "\nRunning blastn with all asv seqs against dereplicated blastn_2 hits ...\n    Starts: $(date)">>$log
    start=$(date +%s.%N)
    $BLASTN \
        -num_threads $THREADS \
        -db blastn_2_derep \
        -max_target_seqs 2000 \
        -query asv_uniq.fasta \
        -out blastn_3.txt \
        -outfmt "6 std qlen slen" \
        -perc_identity 97 \
        -qcov_hsp_perc 97
    echo "    Ends: $(date)">>$log
    end=$(date +%s.%N)
    runtime=$(python -c "print(${end} - ${start})")
    echo "    Done blastn with all asv seqs against dereplicated blastn_2 hits Runtime: $runtime sec" >> $log

    cat blastn_3.txt \
        | awk -F '\t' '{if($3 >= 97 && ($4/$13) > 0.97) print $0}' \
        > blastn_3_filtered.txt

    #yank blastn results and create pickle dictionary file for each of all samples
    ${SCRIPTS}/pickle_blastn_by_sample.py asv_count.csv blastn_3_filtered.txt asv.uc


    #Split ASV file into multiple ones containing 500 sequences in each
    mkdir asv_tmp
    cat asv_PE.fasta | \
        (cd asv_tmp; split -a 8 -l 1000)
        #(cd asv_tmp; split -a 8 --additional-suffix=.fasta -l 1000)

fi