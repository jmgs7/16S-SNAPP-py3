#!/usr/bin/env python
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

# Class of templates to be used for consensus sequences
## Ported to Python 3 on 20210106

from operator import itemgetter
from itertools import groupby
import numpy as np


class Refseq:
    def __init__(self, ID):
        self.ID = ID
        self.seq = ""  # full length reference sequence
        self.reads = {}
        self.positions = []
        self.consensus = ""

    def addRead(self, read_info):  # add read info: mapped positions
        """
        Add a read to the reference sequence.

        Parameters
        ----------
        read_info : dict
            A dictionary of read information with two keys: 'ID' and 'pos'.
            'ID' is the read ID, and 'pos' is a tuple of two integers, the start
            and end positions of the read mapped to the reference sequence.

        Returns
        -------
        None
        """
        ID = read_info["ID"]
        start, end = read_info["pos"]
        positions = range(start, end)  # mapped positions
        self.reads[ID] = positions  # mapped positions
        for position in positions:
            if not position in self.positions:
                self.positions.append(position)
        self.positions.sort()

    def addPEread(
        self, ID, read_info_R1, read_info_R2
    ):  # add PE read info: two sets of positions
        """
        Add a paired-end read to the reference sequence.

        Parameters
        ----------
        ID : str
            The read ID.
        read_info_R1 : tuple
            A tuple of two integers, the start and end positions of the read
            mapped to the reference sequence for the first read in the pair.
        read_info_R2 : tuple
            A tuple of two integers, the start and end positions of the read
            mapped to the reference sequence for the second read in the pair.

        Returns
        -------
        None
        """
        s1, e1 = read_info_R1
        s2, e2 = read_info_R2
        positions = set(range(s1, e1)).union(set(range(s2, e2)))
        positions = list(positions)
        positions.sort()
        self.reads[ID] = positions  # mapped positions
        for position in positions:
            if not position in self.positions:
                self.positions.append(position)
        self.positions.sort()

    def getReadIDs(self):
        """
        Return a set of all read IDs that have been mapped to this reference
        sequence.

        Returns
        -------
        set
            A set of all read IDs that have been mapped to this reference
            sequence.
        """
        return set(self.reads.keys())

    ## Execute the following operations when read mapping is complete and read count normalization is done
    def addSeq(self, seq):  # the full length sequence of reference
        """
        Set the full length sequence of the reference.

        Parameters
        ----------
        seq : str
            The full length sequence of the reference.

        Returns
        -------
        None
        """
        self.seq = seq.upper()

    def addRegs(
        self,
    ):  # this Python 3 version is slightly different from Python 2 version but has the same return, i.e. consecutively covered  regions once the read mapping is complete
        """
        Divide the reference sequence into regions of consecutive coverage.

        This function takes the list of positions where reads have been mapped to
        the reference sequence and divides it into regions of consecutive
        coverage. The output is a list of lists of integers, where each inner list
        represents a region of consecutive coverage and each integer is the
        position of a base in the reference sequence.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.baseRegs = []
        for k, g in groupby(enumerate(self.positions), lambda x: x[0] - x[1]):
            group = map(itemgetter(1), g)
            group = list(map(int, group))
            self.baseRegs.append(group)

    def addReadCounts(
        self, DF
    ):  # add minimized read counts attributable to this reference of this sample
        """
        Add minimized read counts attributable to this reference of this sample.

        Parameters
        ----------
        DF : pandas.DataFrame
            Dataframe with columns as reference IDs and index as read IDs, and
            values as the count of each read mapped to the reference

        Returns
        -------
        None
        """
        readIDs = self.reads.keys()
        self.count = DF.loc[readIDs, self.ID].to_dict()  # count of reads by readID

        self.baseFreq = {pos: 0 for pos in self.positions}  # all mapped positions
        for readID, positions in self.reads.items():
            for pos in positions:
                self.baseFreq[pos] += self.count[readID]

    def addAssign(
        self, tax
    ):  # replace tax assignment with better resolved reads and eventually the emsemble's assignment
        """
        Assign a taxonomic classification to the reference sequence.

        Parameters
        ----------
        tax : str
            The taxonomic classification to be assigned to the reference sequence.

        Returns
        -------
        None
        """
        self.tax = tax

    def getCountSum(self):  # return the count sum of reads mapped to this reference
        """
        Calculate the total count of reads mapped to this reference.

        Returns
        -------
        int
            The sum of read counts mapped to the reference.
        """
        return sum(self.count.values())

    def getMeanBaseCov(
        self,
    ):  # return the average time a mapped base is covered by reads
        """
        Calculate the average time a mapped base is covered by reads.

        Returns
        -------
        float
            The average coverage of mapped bases.
        """
        return np.array(self.baseFreq.values()).mean()

    def getAlignLen(self):
        """
        Calculate the total length of aligned regions.

        Returns
        -------
        int
            The total length of aligned regions.
        """
        length = 0
        for reg in self.baseRegs:
            length += len(reg)
        return length

    def getAlignPct(self):  # fraction of all aligned positions
        """
        Calculate the fraction of all aligned positions.

        Returns
        -------
        float
            The fraction of all aligned positions.
        """
        return round(len(self.positions) / float(len(self.seq)), 2)

    def getProxy(
        self,
    ):  # return the proxy sequence formed by concatenating all mapped regions
        """
        Return the proxy sequence formed by concatenating all mapped regions.

        This function creates a proxy sequence by extracting and concatenating
        segments of the reference sequence that correspond to mapped regions.
        The segments are joined with "NNNNNNN" representing gaps between them.

        Returns
        -------
        str
            The proxy sequence with concatenated mapped regions, separated by "NNNNNNN".
        """
        proxy = [self.seq[reg[0] : reg[-1]] for reg in self.baseRegs]
        return "NNNNNNN".join(proxy)


if __name__ == "__main__":
    myref = Refseq(ID="r1", seq="GCAACTGGACTGGAA")
    myref.addReads(["asv1", "asv2", "asv3"])
    myref.addReads(["asv8"])
    myref.addBases(set([4, 5, 9, 3, 10]))
    myref.addBases(set([4, 5, 40, 59, 10]))
    print(myref.ID)
    print(myref.seq)
    print(myref.readIDs)
    print(myref.getCount())
    print(myref.bases)
    print(myref.getBaseCov())
    print(myref.getRegs())
