#!/usr/bin/env python
# used to aggregage sequences and counts from all processed samples

import os
import sys
import math
import json


def get_lineages(tax_filename, CONF):
    """
    Parse taxonomic assignments from a file and return a mapping of sequence IDs to their corresponding lineages.

    This function reads a tab-separated file containing taxonomic information for sequences.
    It processes each line to extract the sequence ID and its associated taxonomic lineage,
    filtering based on a confidence threshold. The lineage is constructed as a semicolon-separated
    string of taxonomic ranks, prefixed by their rank abbreviation, and is stored in a dictionary
    keyed by the sequence ID.

    Parameters:
    tax_filename (str): The path to the file containing taxonomic assignments.
    CONF (float): The confidence threshold for accepting taxonomic assignments.

    Returns:
    dict: A dictionary mapping sequence IDs to their taxonomic lineages.
          Sequence IDs that cannot be classified to the domain level are marked as "Unclassified".
    """
    Hash = {}
    lines = open(tax_filename.strip(), "r").readlines()
    for row in lines:
        cols = row.strip().split("\t")
        ID = cols[0]
        if ID.find("asv") == -1:
            ID = integer_size(ID)[0]
        levels = cols[2:]
        i = 0
        lineage = []
        while i < len(levels):
            level = levels[i : i + 3]
            name, rank, conf = level
            name = rank[0] + "__" + name.replace('"', "")
            if float(conf) < CONF:
                break
            else:
                lineage.append(name.replace('"', ""))
            i += 3
        lineage = ";".join(lineage).strip()
        if lineage == "":  # rare cases the sequence can't be classified to domain level
            lineage = "Unclassified"
        Hash[ID] = lineage
    return Hash


def combine_consensus(f, sample_name):
    """
    Combine consensus sequences from all samples and write them to a single file.

    Parameters:
    f (str): The path to the consensus file for this sample.
    sample_name (str): The name of the sample being processed.

    Returns:
    int: A count of the number of sequences processed.
    """
    global counter
    count_table[sample_name] = {}
    cons = open(f, "r").read().strip(">").split("\n>")
    for con in cons:
        lines = con.split("\n")
        ID = lines[0].strip()
        seq = lines[1].strip()
        SID, sample_name, size = ID.split(";")
        sample_name = sample_name.split("=")[1]
        newid, size = integer_size(ID)
        # if ID.find('S00') != -1: #consensuse sequences
        if not ID.find("asv") != -1:  # consensuse sequences
            consensus_out.write(
                ">" + newid + "\n" + seq + "\n"
            )  # output the consensus sequences with integer size
        count_table[sample_name][newid] = int(size)
    return 1


def integer_size(name):
    """
    Take a string of format "ID;size=123" and return a tuple where the first element is the new string
    with the size rounded up to the nearest integer and the second element is the integer size.

    Parameters:
    name (str): The string of format "ID;size=123" to be processed.

    Returns:
    tuple: A tuple containing the new string with the size rounded up to the nearest integer and
           the integer size.
    """
    ID, size = name.split(";size=")
    size = int(math.ceil(float(size)))
    newid = ID + ";size=" + str(size)
    return newid, size


count_table = {}
idmap = {}
consensus_out = open("all_cons.fasta", "w")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("sum_count.py wd outname")
        sys.exit()

wd = sys.argv[1]  # work directory
outname = sys.argv[2]  # name prefix
outname = os.path.join(wd, outname)

counter = 0  # for renaming consensus sequences
for f in os.listdir(wd):
    if not f.endswith("_consensus.fasta"):
        continue
    sample_name = f.split("_consensus.fasta")[0]
    path = os.path.join(wd, f)
    combine_consensus(path, sample_name)
consensus_out.close()

with open(outname + "_counts.json", "w") as outfile:
    json.dump(count_table, outfile)

# load taxonomic assignment of all consensus and asv sequences and save it to json
taxonomy_table = {}
for f in os.listdir(wd):
    if not f.endswith(".cls"):
        continue
    sample_name = f.split(".cls")[0]
    path = os.path.join(wd, f)
    taxonomy_table.update(get_lineages(path, 0.7))
with open(outname + "_assignments.json", "w") as outfile:
    json.dump(taxonomy_table, outfile)
