#!/usr/bin/env python
## Swift Biosciences 16S snapp workflow - Adapted by Bionformatics Service @ Genyo (Granada, Spain)
## Author (GitHub: jmgs7) 24102024

import pandas as pd
import sys

count_table = sys.argv[0]
taxonomy_table = sys.argv[1]

count_table = pd.read_csv(count_table, sep="\t", header=0, index_col=0)
taxonomy_table = pd.read_csv(taxonomy_table, sep="\t", header=0, index_col=0)

def convert_index_to_seq_index(df):
    max_rows = len(df)
    max_digits = len(str(max_rows))
    old_to_new_index_map = {}
    for i, old_index in enumerate(df.index):
        new_index = f"seq_{str(i).zfill(max_digits)}"
        old_to_new_index_map[old_index] = new_index
        df.index[i] = new_index
    return df, old_to_new_index_map

count_table, id_map = convert_index_to_seq_index(count_table)
taxonomy_table.index = [id_map[index] for index in taxonomy_table.index]