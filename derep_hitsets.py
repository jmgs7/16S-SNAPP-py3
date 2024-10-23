#!/usr/bin/env python
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

# the mini workflow to dereplidate blastn file prior to converge step
## simply dereplicate read sets by IDs to make a list of unique sets
## Ported to Python 3 on 20210106

import sys


# To make reference-to-read set hash keyed by reference seq IDs
def getSets(blastn):
    """
    Generate a hash of reference-to-read set from a blastn output file.

    Parameters
    ----------
    blastn : str
        Path to the blastn output file.

    Returns
    -------
    dict
        A dictionary where keys are reference sequence IDs and values are
        sets of read IDs that hit the reference sequence.
    """
    f = open(blastn, "r")
    Hash = {}
    while 1:
        try:
            line = next(f)
            sID, qID = line.strip().split("\t")[0:2]
            if not sID in Hash:
                Hash[sID] = set([])
            Hash[sID].add(qID)
        except StopIteration:
            break
    return Hash


# Simply dereplicate the list of sets of read IDs
def getUniqSets(Hash):
    """
    Dereplicate the list of sets of read IDs, ensuring uniqueness.

    Parameters
    ----------
    Hash : dict
        A dictionary where keys are reference sequence IDs and values are
        sets of read IDs that hit the reference sequence.

    Returns
    -------
    dict
        A dictionary where keys are unique reference sequence IDs and values
        are unique sets of read IDs.
    """
    uniqSets = []
    dereped = {}
    for sID, qSet in Hash.items():
        if not qSet in uniqSets:
            uniqSets.append(qSet)
            dereped[sID] = qSet
    return dereped


def getDerepHitSets(blastn):
    """
    Dereplicate a blastn output file by reference-to-read set to obtain
    a list of unique reference sequence IDs.

    Parameters
    ----------
    blastn : str
        Path to the blastn output file.

    Returns
    -------
    list
        A list of unique reference sequence IDs.
    """
    sets = getSets(blastn)
    uniqSets = getUniqSets(sets)
    return uniqSets.keys()


if __name__ == "__main__":
    derepIDs = getDerepHitSets(sys.argv[1])
    for ID in derepIDs:
        print(ID)
