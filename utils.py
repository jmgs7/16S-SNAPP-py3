#!/usr/bin/env python
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

# to obtain a ref-read dataframe with count value for each sample
## Ported to Python 3 on 20210106


def get_ref_read_df(
    refset, count_table
):  # to obtain a ref-read dataframe with count value for each sample
    """
    Construct a ref-read dataframe with count value for each sample

    Parameters
    ----------
    refset : list of Refseq objects
        List of Refseq objects
    count_table : dict
        Dictionary of read IDs to their counts

    Returns
    -------
    pandas.DataFrame
        Dataframe with columns as reference IDs and index as read IDs, and
        values as the count of each read mapped to the reference
    """

    import pandas as pd

    id_dict = {ref.ID: {} for ref in refset}
    for refseq in refset:
        ref_id = refseq.ID
        read_ids = refseq.getReadIDs()
        id_dict[ref_id] = {read_id: count_table[read_id] for read_id in read_ids}
    return pd.DataFrame.from_dict(id_dict)  # dictionary to a DF


# add reference`sequence strings to refseq objects
def fetch_refseq(RDPHOME, id_file_name, outfile_name, reffile_name):
    """
    Fetch the reference sequences from the reference database based on the IDs in the
    file, and write the sequences to a new file.

    Parameters
    ----------
    RDPHOME : str
        The path to the RDPTools directory.
    id_file_name : str
        The name of the file containing the reference IDs.
    outfile_name : str
        The name of the output file to write the sequences to.
    reffile_name : str
        The name of the reference file to fetch the sequences from.

    Returns
    -------
    int
        1 if the command was successfully called.

    """

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
    """
    Add additional attributes to the Refseq objects.

    This function reads a reference file to extract sequences and updates
    the Refseq objects in the provided refset with the full-length sequence,
    read counts, and mapped regions.

    Parameters
    ----------
    DF : pandas.DataFrame
        The dataframe containing read counts for each reference.
    reffile_name : str
        The name of the reference file containing sequences.
    refset : dict
        A dictionary of Refseq objects keyed by reference ID.

    Returns
    -------
    dict
        The updated dictionary of Refseq objects with added attributes.
    """

    recs = open(reffile_name, "r").read().strip(">").split("\n>")
    # Iterate all refseq objects for multiple tasks
    for rec in recs:
        lines = rec.split("\n")
        if lines:
            ID = lines[0].split()[0]
            seq = "".join(lines[1:]).replace(" ", "")
            refset[ID].addSeq(seq)
            refset[ID].addReadCounts(DF)
            refset[ID].addRegs()
        else:
            continue
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
def classify_proxy(sample_id, RDPHOME, RDPHOME_CUSTOM, WD):
    """
    Classify consensus sequences using the RDP classifier.

    This function classifies consensus sequences for a given sample using the
    RDP classifier. It attempts to use a custom training set specified in the
    environment variable `RDP_CLASSIFIER`. If not available, it falls back to
    using the default training set.

    Parameters
    ----------
    sample_id : str
        The sample ID for which consensus sequences are classified.
    RDPHOME : str
        The path to the RDPTools directory containing the default classifier.
    RDPHOME_CUSTOM : str
        The path to the custom RDPTools directory containing a custom classifier.
    WD : str
        The working directory where input consensus sequences are stored and
        output classifications will be saved.

    Returns
    -------
    int
        Returns 1 to indicate the classification was executed successfully.
    """

    import subprocess
    import os

    # classfy the consensus sequences
    try:
        train_set = os.environ["RDP_CLASSIFIER"]  # use the specified training set
        subprocess.check_call(
            # Use local git instalation to avoid memory overload of conda's RDP classifier binary.
            [
                "java",
                "-jar",
                "-Xmx8g",
                os.path.join(RDPHOME_CUSTOM, "classifier.jar"),
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
    """
    Make a phylogenetic tree from template sequences.

    Parameters
    ----------
    seqfile_name : str
        File name of the fasta file containing the template sequences.
    WD : str
        Working directory where the aligned file will be saved.
    RESDIR : str
        Directory where the tree will be saved.

    Returns
    -------
    int
        Returns 1 to indicate the tree building was successful.
    """

    import os

    aligned = os.path.join(WD, "templates_mafft.fasta")
    tree = os.path.join(RESDIR, "templates_mafft.tree")
    os.system(
        "mafft --quiet --thread $THREADS FASTA > ALIGNED".replace(
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
    """
    Run seqmatch in parallel processing mode for paired end sequences.

    Parameters
    ----------
    folder_name : str
        The folder name containing the split files of paired end sequences.
    WD : str
        The working directory where the seqmatch DB is located.

    Returns
    -------
    dict
        A dictionary of SequenceMatch results.

    Notes
    -----
    This function will run seqmatch for all the split files of paired end
    sequences in the specified folder in parallel mode. The results are stored
    in a dictionary where the keys are the file names and the values are the
    SequenceMatch results.
    """

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
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=int(os.environ.get("THREADS"))
    ) as executor:
        results = [executor.submit(seq_match, WD, fname) for fname in fnames]
        for f in concurrent.futures.as_completed(results):
            # all_seq_matches.append(f.result())
            rs_dict.update(f.result())
    return rs_dict


def seq_match(WD, QUERY):  # function to run seqmatch
    """
    Run seqmatch for a query fasta file.

    Parameters
    ----------
    WD : str
        The working directory where the seqmatch DB is located.
    QUERY : str
        The path to the query fasta file.

    Returns
    -------
    dict
        A dictionary of SequenceMatch results.
    """

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
    """
    Read a fasta file and return a dictionary of sequences where the keys are
    the sequence IDs and the values are the sequences.

    Parameters
    ----------
    seqfile_name : str
        The path to the fasta file.

    Returns
    -------
    dict
        A dictionary of sequences.
    """

    seq_dict = {}
    recs = open(seqfile_name, "r").read().strip(">").split("\n>")
    for rec in recs:
        lines = rec.split("\n")
        ID = lines[0].split()[0]
        seq = "".join(lines[1:]).replace(" ", "")
        seq_dict[ID] = seq
    return seq_dict


def rev_complement(seq):
    """
    Compute the reverse complement of a given DNA sequence.

    Parameters
    ----------
    seq : str
        The DNA sequence.

    Returns
    -------
    str
        The reverse complement of the given DNA sequence.
    """

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
