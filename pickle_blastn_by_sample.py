#!/usr/bin/env python
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

# extract blastn hits as the dictionary and save as pickle for later use
## Ported to Python 3 on 20210106

import sys
import os
import pandas as pd
import blastn_parser
import timeit
import pickle
import derep_hitsets


def sliceSampleHits(Dict, IDs):  # extract blast info for this sample
    """
    Slice a subset of the blastn hits dictionary Dict for a given sample by
    extracting only the blast info associated with the given list of IDs.

    Parameters
    ----------
    Dict : dict
        A dictionary of blastn hits where each key is a reference ID and
        the value is a dictionary mapping ASV IDs to lists of two lists
        representing the read 1 and read 2 coordinates.
    IDs : list
        A list of ASV IDs to extract from Dict.

    Returns
    -------
    Hash : dict
        A dictionary of blastn hits where each key is a reference ID and
        the value is a dictionary mapping ASV IDs to lists of two lists
        representing the read 1 and read 2 coordinates. The dictionary will
        only contain the references and their associated ASVs found in the
        given list of IDs.
    """
    Hash = {}
    IDs = set(IDs)
    for refid in Dict.keys():
        if not refid in Hash:
            Hash[refid] = {}
        asvids = set(Dict[refid].keys())
        count = 0  # count the qualified asv matches
        for asvid in IDs.intersection(asvids):
            R1_info = Dict[refid][asvid][0]
            R2_info = Dict[refid][asvid][1]
            if R1_info == [] or R2_info == []:
                continue
            else:
                Hash[refid][asvid] = [R1_info, R2_info]
                count += 1
        if count == 0:  # no qualified asv PE matches for this refseq
            del Hash[refid]
    return Hash


def pickle_blastn(countTable, blastn, ucFile):
    """
    Process BLASTN results and ASV count data to generate pickle files for each sample.

    This function creates a directory named 'pickle' and processes BLASTN results
    along with an ASV count table. It extracts and dereplicates hit information for
    each sample and saves it as a pickle file. Additionally, it writes reverse
    complement ASV IDs to a text file.

    Parameters
    ----------
    countTable : str
        Path to the ASV count CSV file.
    blastn : str
        Path to the BLASTN output file.
    ucFile : str
        Path to the UC file containing sequence clustering information.

    Outputs
    -------
    - Pickle files for each sample containing dereplicated hit information.
    - A text file 'reverse_complement_IDs.txt' containing reverse complement ASV IDs.
    """
    os.system("mkdir -p pickle")
    start_time = timeit.default_timer()
    rDF = pd.read_csv(countTable, sep=",", header=0, index_col=0)
    blastn_info, rc_list = blastn_parser.get_blastn_hits(blastn, ucFile)
    print(len(blastn_info))
    blastn_info_time = timeit.default_timer()
    print(blastn_info_time - start_time, "Dictionary time")
    for sampleID in rDF.columns:
        sample_time = timeit.default_timer()
        sample_count = rDF.loc[:, sampleID]
        sample_count = sample_count[
            sample_count > 0
        ]  # drop 0's count series of this sample_PE_IDs
        sample_PE_IDs = list(sample_count.index)
        sampleInfo = sliceSampleHits(blastn_info, sample_PE_IDs)
        print(len(sampleInfo))
        # print sampleInfo, 'sampleInfo'
        slice_time = timeit.default_timer()
        print(slice_time - sample_time, sampleID, "slice_time")
        derepInfo = derep_hitsets.getUniqSets(sampleInfo)
        f = open("pickle/" + sampleID + ".pkl", "wb")
        pickle.dump(derepInfo, f)
        derep_time = timeit.default_timer()
        print(len(derepInfo), "derepInfo")
        print(derep_time - slice_time, sampleID, "derep_time")
    with open("reverse_complement_IDs.txt", "w") as out:
        out.write("\n".join(rc_list))
    print(derep_time - start_time, "total_time")


if __name__ == "__main__":
    if not len(sys.argv) == 4:
        print("pickle_blastn_by_sample.py asv_count.csv blastn asv.uc")
        sys.exit()

    table = sys.argv[1]
    blast = sys.argv[2]
    uc = sys.argv[3]
    pickle_blastn(table, blast, uc)
