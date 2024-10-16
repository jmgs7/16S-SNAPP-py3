#!/usr/bin/env python
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

# to obtain a ref-read dataframe with count value for each sample
## Ported to Python 3 on 20210106


def get_ref_read_df(
    refset, count_table
):  # to obtain a ref-read dataframe with count value for each sample
    import pandas as pd

    id_dict = {ref.ID: {} for ref in refset}
    for refseq in refset:
        ref_id = refseq.ID
        read_ids = refseq.getReadIDs()
        id_dict[ref_id] = {read_id: count_table[read_id] for read_id in read_ids}
    return pd.DataFrame.from_dict(id_dict)  # dictioary to a DF


# add reference`sequence strings to refseq objects
def fetch_refseq(RDPHOME, id_file_name, outfile_name, reffile_name):
    import subprocess
    import os

    subprocess.check_call(
        [
            os.path.join(RDPHOME, "ReadSeq"),
            "select-seqs",
            id_file_name,
            outfile_name,
            "fasta",
            "Y",
            reffile_name,
        ],
        stdin=None,
        stdout=None,
        stderr=None,
        shell=False,
    )
    return 1


def update_refseq(
    DF, reffile_name, refset
):  # add additional attributes to the refseq objects
    recs = open(reffile_name, "r").read().strip(">").split("\n>")
    # Iterate all refseq objects for multiple tasks
    for rec in recs:
        lines = rec.split("\n")
        ID = lines[0].split()[0]
        seq = "".join(lines[1:]).replace(" ", "")
        refset[ID].addSeq(seq)
        refset[ID].addReadCounts(DF)
        refset[ID].addRegs()
    return refset


## make a unique id for each consensus by adding a serial number to its template ID
class Name_proxy:
    counter = {}

    def get_assumed_id(self, ID):
        if not ID in self.__class__.counter:
            self.__class__.counter[ID] = 1
        else:
            self.__class__.counter[ID] += 1
        number = self.__class__.counter[ID]
        name = ID + "_" + str("0" * (3 - len(str(number))) + str(number))
        return name


## classify consensus sequences in batch, fetch, and add assignments to Refseq objects
def classify_proxy(sample_id, RDPHOME, WD):
    import subprocess
    import os

    # classfy the consensus sequences
    try:
        train_set = os.environ["RDP_CLASSIFIER"]  # use the specified training set
        subprocess.check_call(
            [
                os.path.join(RDPHOME, "classifier"),
                "-t",
                train_set,
                "-o",
                os.path.join(WD, sample_id + ".cls"),
                os.path.join(WD, sample_id + "_consensus.fasta"),
            ]
        )
    except KeyError:  # use default training set
        subprocess.check_call(
            [
                os.path.join(RDPHOME, "classifier"),
                "-f",
                "fixrank",  # Change -f fixrank to -f allrank when uses a custon classifier with species-level resolution.
                "-o",
                os.path.join(WD, sample_id + ".cls"),
                os.path.join(WD, sample_id + "_consensus.fasta"),
            ]
        )
    return 1


def build_tree(
    seqfile_name, WD, RESDIR
):  # make a phylogenetic tree from template sequences
    import os

    aligned = os.path.join(WD, "templates_mafft.fasta")
    tree = os.path.join(RESDIR, "templates_mafft.tree")
    os.system(
        "mafft --quiet --thread 4 FASTA > ALIGNED".replace(
            "FASTA", seqfile_name
        ).replace("ALIGNED", aligned)
    )

    os.system(
        "fasttree -quiet -nopr -nt -gtr < ALIGN > \
               TREE".replace(
            "ALIGN", aligned
        ).replace(
            "TREE", tree
        )
    )
    return 1


def run_seqmatch(folder_name, WD):  # run seqmatch of PE in parallel processing mode
    import os
    import concurrent.futures

    # the split files names of pe_seq
    fnames = [
        os.path.join(folder_name, f)
        for f in os.listdir(folder_name)
        if os.path.isfile(os.path.join(folder_name, f))
    ]
    # all_seq_matches = []
    rs_dict = {}  # all SequenceMatch results
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        results = [executor.submit(seq_match, WD, fname) for fname in fnames]
        for f in concurrent.futures.as_completed(results):
            # all_seq_matches.append(f.result())
            rs_dict.update(f.result())
    return rs_dict


def seq_match(WD, QUERY):  # function to run seqmatch
    import os

    DB = os.path.join(WD, "seqmatch")
    cmd = "${RDPHOME}/SequenceMatch \
    seqmatch -k 1 -s 0.4 DB QUERY".replace(
        "DB", DB
    ).replace(
        "QUERY", QUERY
    )
    rs = os.popen(cmd).readlines()
    rs_dict = {}
    for line in rs[1:]:
        asv_id, hit_id = line.strip().split("\t")[0:2]
        rs_dict[asv_id] = hit_id
    return rs_dict


def read_seq(seqfile_name):
    seq_dict = {}
    recs = open(seqfile_name, "r").read().strip(">").split("\n>")
    for rec in recs:
        lines = rec.split("\n")
        ID = lines[0].split()[0]
        seq = "".join(lines[1:]).replace(" ", "")
        seq_dict[ID] = seq
    return seq_dict


def rev_complement(seq):
    anticodon = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "U": "A",
        "M": "K",
        "K": "M",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "V": "B",
        "B": "V",
        "H": "D",
        "D": "H",
        "X": "X",
        "N": "N",
        "-": "-",
        ".": ".",
    }
    revseq = ""
    while seq:  # iterate over every character in the string
        base = seq[-1].upper()
        try:
            revseq = revseq + anticodon[base]
        except KeyError:
            revseq = revseq + base
        seq = seq[
            :-1
        ]  # delete the character positioned at[-1], which has been processed
    return revseq
