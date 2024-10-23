#!/usr/bin/env python
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

## to parse the classifier results into a hash of counts keyed by lineage name

## Ported to Python 3 on 20210106

import os


def get_lineages(tax_filename, CONF):
    """
    Parse the classifier results into a hash of counts keyed by lineage name.

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
        ID = ID.split(";")[0]
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


def get_best_tax(ref_tax, read_taxs):  # choose the better between reftax and readtax
    """
    Choose the better between reftax and readtax based on the resolution of taxonomic lineage.

    Parameters:
    ref_tax (str): The taxonomic lineage of the reference sequence.
    read_taxs (list): A list of taxonomic lineages of the associated reads.

    Returns:
    str: The most resolved lineage between the reference and associated reads.
    """
    taxs = {}
    for tax in read_taxs:
        level = tax.count(";")
        if not level in taxs:
            taxs[level] = []
        taxs[level].append(tax)
    resolution = taxs.keys()
    # resolution.sort() #changed to next line for Python 3
    resolution = sorted(resolution)
    resolution.reverse()
    top_lineage = taxs[resolution[0]][0]
    if ref_tax.count(";") >= resolution[0]:
        top_lineage = ref_tax
    return top_lineage  # the most resolved lineage


def add_lineages(
    wd, sample_id, tax_dict, refset
):  # refset pass through to add lineage information
    """
    Add the taxonomic lineage of the reference sequences to the Refseq objects.

    Parameters:
    wd (str): The working directory containing the taxonomic assignment files.
    sample_id (str): The sample ID associated with the taxonomic assignments.
    tax_dict (dict): A dictionary mapping sequence IDs to their taxonomic lineages.
    refset (dict): The dictionary of Refseq objects.

    Returns:
    dict: The updated dictionary of Refseq objects with added taxonomic lineages.
    """
    ref_tax_file = os.path.join(wd, sample_id + ".cls")
    ref_tax_hash = get_lineages(ref_tax_file, 0.7)
    for ref_id in refset.keys():
        read_ids = list(refset[ref_id].getReadIDs())
        ref_lineage = ref_tax_hash[ref_id]
        read_lineage_list = list(set([tax_dict[ID][0] for ID in read_ids]))
        best_lineage = get_best_tax(ref_lineage, read_lineage_list)
        refset[ref_id].addAssign(best_lineage)
    return refset
