#!/bin/bash
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502
## Swift Biosciences 16S snapp.sh workflow

## The main wrapper to run the wordflow
## Remember to edit/set the parameters in config.txt file
## Run as snapp.sh config.txt inputDir workDir
## make sure work dir exists before running the pipeline

START=$(date +%s.%N)

if [ $# -ne 3 ]
    then
         echo 'snapp.sh config.txt inputdir workdir'
         exit
fi
set -e
set -x

export SCRIPTS=$(dirname $(readlink -f $0)) #the path to scripts
source $1 #load config.txt
echo `cat ${1}`
INPUTDIR=$(readlink -f $2) # the absolute dir where the fastq.gz file are located
export WD=$(readlink -f $3) #work folder
runlog='log'
current_time=$(date "+%Y.%m.%d-%H.%M.%S")
echo $current_time
export log=${runlog}.${current_time}

[ ! -d "$WD" ] && echo "Directory $WD does not exist, please create one..." && exit

##Make a directory for primer-trimmed sequence files
mkdir ${WD}/trimmed
mkdir RESDIR
export RESDIR=$(readlink -f $PWD/RESDIR)
cd ${WD}

##Match and trim primers from PE reads
echo -e "Matching/trimming primers...\n    Starts: $(date)" >> $log
start=$(date +%s.%N)
trimStats='trimStats.txt'
printf "SampleID\tStarting read count\tPrimer trimmed\n" >> ${trimStats}
for R1 in ${INPUTDIR}/*_R1_*fastq.gz; do
    echo $PRIMERS
    R2=${R1/_R1_/_R2_} #the path to the read 2 file
    basenameR1=${R1##*/}
    basenameR2=${R2##*/}
    prefixR1=${basenameR1%".fastq.gz"}
    prefixR2=${basenameR2%".fastq.gz"}
    prefix=${prefixR1%_L001_R1_*}
    totalCount=$(bc <<<  $(zcat ${R1}|wc -l)/4)

    #PE option to trim primers off the reads
    $CUTADAPT -e 0.10 -g file:$PRIMERS -G file:$PRIMERS \
              -o trimmed/${basenameR1} \
              -p trimmed/${basenameR2} \
              --untrimmed-output ${prefixR1}_NotTrimmed.fastq \
              --untrimmed-paired-output ${prefixR2}_NotTrimmed.fastq \
              $R1 $R2 \
              --max-n 0 \
              --minimum-length ${READLEN}

    [[ ! -s trimmed/${basenameR1} ]] \
        && rm trimmed/${basenameR1} \
        && rm trimmed/${basenameR2} \
        && echo "${prefix} zero match primers" \
        && rm ${prefix}* \
        && printf "${prefix}\t${totalCount}\t0}\n" >> ${trimStats}\
        && continue

    trimmedCount=$(bc <<<  $(zcat trimmed/${basenameR1}|wc -l)/4)
    trimPCT=$(bc <<< "scale=4 ; (${trimmedCount}/${totalCount})*100" )
    printf "${prefix}\t${totalCount}\t${trimmedCount}\n" >> ${trimStats}
    echo "$prefix trimmed: ${trimPCT}%" >> $log

done
rm *NotTrimmed.fastq
echo "    Ends: $(date)">>$log
end=$(date +%s.%N)
runtime=$(python -c "print(${end} - ${start})")
echo "    Primer trimming Runtime: $runtime sec" >> $log

##Run Dada2 to obtain ASVs and remove chimeric reads from primer-trimmed
##fastq PE files of all samples
echo -e "\nRunning DADA2...\n    Starts: $(date)">>$log
start=$(date +%s.%N)
printf '\n' >> $log
echo -e "\nDADA2 processing stats:" >> $log
#Prepare asv count tables and sequences in PEs, single formats
${SCRIPTS}/run_dada2.R trimmed/ . ${READLEN} >> $log #run dada2 with R
# cat "DADA2_summary.csv" | sed -e 's/"//g' >> $log
${SCRIPTS}/get_asv_files.py asv_seqNcount.csv asv

echo -e "\n    Unique ASV pairs from DADA2: $(bc <<< "$(grep -c ">" \
    asv_seq.fasta)/2")"  >>$log
echo "    Ends: $(date)">>$log
end=$(date +%s.%N)
runtime=$(python -c "print(${end} - ${start})")
echo "    DADA2 Runtime: $runtime sec" >> $log

##Classify asv PEs with RDP classifier
echo -e "\nRunning RDP classifier...\n    Starts: $(date)">>$log
start=$(date +%s.%N)
if [[ -z $RDP_CLASSIFIER ]]
then
    ${RDPHOME}/classifier \
        -f fixrank \
        -o asv_PE.cls \
        asv_PE.fasta
else
    ${RDPHOME}/classifier \
        -t $RDP_CLASSIFIER \
        -o asv_PE.cls \
        asv_PE.fasta
fi

end=$(date +%s.%N)
runtime=$(python -c "print(${end} - ${start})")
echo "    Read classfication Runtime: $runtime sec" >> $log

##Data reduction, i.e. reduce the number of reads for downstream processing
#Retain only unique sequences comparing both strands
echo -e "\nDereplicating ASVs...\n    Starts: $(date)" >>$log
start=$(date +%s.%N)
$VSEARCH --cluster_size asv_seq.fasta \
         --strand both \
         --iddef 1 \
         --id 1.00 \
         --uc asv.uc \
         --centroid asv_uniq.fasta
end=$(date +%s.%N)
runtime=$(python -c "print(${end} - ${start})")
echo "    Dereplication Runtime: $runtime sec" >> $log

${SCRIPTS}/blastn_for_templates.sh

##Determine the candiate reference sequenes, associate reads, allocate read counts,
##and classify each consensus sequence representing each associated amplicon set.
echo -e "\nConverging candidate template sequences...\n    Starts: $(date)">>$log
start=$(date +%s.%N)
${SCRIPTS}/converge.py \
        asv_count.csv \
        asv_PE.fasta \
        asv_PE.cls \
        $log

end=$(date +%s.%N)
runtime=$(python -c "print(${end} - ${start})")
echo "    Converging Runtime: $runtime sec" >> $log
echo "\n" >> $log

##Cluster all sequences in consensus files using kmer in R
start=$(date +%s.%N)
cat *consensus.fasta > all_cons.fasta
${SCRIPTS}/sum_count.py . all

$VSEARCH --sortbysize all_cons.fasta --output all_cons_sorted.fasta

mv all_cons_sorted.fasta all_cons.fasta

[ ! -d "tmp" ] && echo "create tmp diretory" && mkdir tmp

#split the sorted fasta file into smaller ones for running minimap2 alignment
${SCRIPTS}/splitFasta.py all_cons.fasta 10000 tmp
for fas in tmp/*.fas; do
    prefix=${fas%.fas}
    prefix=${prefix##*/all_cons_}
    minimap2 -a --eqx -o tmp/${prefix}.sam all_cons.fasta ${fas}
done

cat tmp/*.sam | grep -v "^@" > all_cons.sam
cutoff=0.03
${SCRIPTS}/parse_minimap2_clr.py all_cons.sam all_cons.clr $cutoff
${SCRIPTS}/get_OTU_table.py . all_cons.fasta all_cons.clr OTUs

END=$(date +%s.%N)
runtime=$(python -c "print(${END} - ${start})")
echo "    Clustering time: $runtime sec" >> $log
echo "\n" >> $log
echo "\n" >> $log
runtime=$(python -c "print(${END} - ${START})")
echo -e "\nWhole process completed in: $runtime sec">>$log

#clean up
rm -rf tmp
