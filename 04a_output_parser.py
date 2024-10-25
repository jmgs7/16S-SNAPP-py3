#!/usr/bin/env python
## Swift Biosciences 16S snapp workflow - Adapted by Bionformatics Service @ Genyo (Granada, Spain)
## Author (GitHub: jmgs7) 24102024

import pandas as pd
import sys
import os
import timeit
from utils import build_tree


def create_index_dict(df):
    """
    Creates a dictionary mapping old index values to new index values of a pandas DataFrame.


    The new index values are generated as "seq_" followed by a zero-padded
    integer corresponding to the row position (1-indexed).

    Parameters
    ----------
        df (pandas.DataFrame): Input DataFrame.

    Returns
    -------
        dict: Dictionary mapping old index values to new index values.

    Example:
        >>> df = pd.DataFrame({'A': [1, 2, 3]})
        >>> index_map = create_index_dict(df)
        >>> print(index_map)
        {0: 'seq_001', 1: 'seq_002', 2: 'seq_003'}
    """

    import pandas as pd

    max_rows = len(df)
    max_digits = len(str(max_rows))
    old_to_new_index_map = {}
    for i, old_index in enumerate(df.index):
        new_index = f"seq_{str(i+1).zfill(max_digits)}"
        old_to_new_index_map[old_index] = new_index
    return old_to_new_index_map


def parse_fasta(fasta_in, fasta_out, id_map):
    """
    Parse a FASTA file, change its heeaders based on a dictionary
    mapping old index values to new index values,
    and sort the sequences alphabetically based on their headers.

    Parameters
    ----------
        fasta_in (str): The path to the input FASTA file.
        fasta_out (str): The path to the output FASTA file.
        id_map (dict): A dictionary mapping old header to new headers.

    Returns
    -------
        None
    """

    with open(fasta_in, "r") as fasta:
        # The recurrent use of .strip() avoids trailing newlines and spaces problems.
        parsed_fasta = []
        for line in fasta:
            if line.startswith(">"):
                seq_id = line.replace(">", "").strip()
                parsed_fasta.append(id_map[seq_id])
            else:
                parsed_fasta.append(
                    line.upper().strip()
                )  # Make sure the sequence is uppercase

    def sort_fasta(fasta_list):
        """
        This is a helper function that sorts the fasta file based on the headers.

        Parameters
        ----------
        fasta_list : list
            A list of strings representing the fasta file.

        Returns
        -------
        sorted_fasta : Pandas DataFrame
            A pandas DataFrame mapping headers to sequences.
        """

        import pandas as pd

        sorted_fasta = pd.DataFrame.from_dict(
            dict(zip(fasta_list[::2], fasta_list[1::2])), "index"
        )
        sorted_fasta.sort_index(inplace=True)

        return sorted_fasta

    # Sort the parsed fasta
    sorted_fasta = sort_fasta(parsed_fasta)

    # Transform to a dict and write the fasta file
    output = sorted_fasta.to_dict()[0]
    with open(fasta_out, "w") as f:
        for seq_id, seq in output.items():
            f.write(f">{seq_id}\n{seq}\n")

    # Return the sorted fasta DataFrame
    return sorted_fasta


def build_seq_table(fasta_df, id_map):
    """
    Constructs a sequence table from a DataFrame of FASTA sequences and an ID mapping.

    This function takes a DataFrame containing FASTA sequences and a dictionary
    mapping new sequence IDs to their original IDs. It creates a new DataFrame
    with two columns: 'Origin', which contains the original sequence IDs, and
    'Sequence', which contains the corresponding sequences.

    Parameters
    ----------
    fasta_df : pandas.DataFrame
        A DataFrame where each row represents a sequence with its ID as the index.
    id_map : dict
        A dictionary mapping new sequence IDs to original sequence IDs.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with 'Origin' and 'Sequence' columns, mapping original IDs
        to their respective sequences.
    """

    import pandas as pd

    # Make a deep copy to avoid modifying the original DataFrame
    sequence_table = fasta_df.copy(deep=True)
    # Rename index column and sequence column
    sequence_table.index.names = ["ID"]
    sequence_table.columns = ["Sequence"]

    # Invert the id_mapping so it accepts the simplified ID as key
    inv_map = {v: k for k, v in id_map.items()}

    # Add the original ID to the sequence table
    sequence_table["Origin"] = [inv_map[header] for header in fasta_df.index]
    # Reorder the columns
    sequence_table = sequence_table[["Origin", "Sequence"]]

    return sequence_table


##main##
if __name__ == "__main__":
    if not len(sys.argv) == 4:
        print(len(sys.argv))
        print("output_parser.py feature-table.tsv taxonomy-table.tsv templates.fasta")
        sys.exit()


start_time = timeit.default_timer()

# environmental variables
WD = os.environ["WD"]  # work directory
RESDIR = os.environ["RESDIR"]  # results directory

# Arguments
count_table = sys.argv[1]
taxonomy_table = sys.argv[2]
templates_fasta = sys.argv[3]

# Creates specific directory for results
results_dir = os.path.join(RESDIR, "results")
os.makedirs(results_dir, exist_ok=True)

# Import the count and taxonomy tables
count_table = pd.read_csv(count_table, sep="\t", header=0, index_col=0)
taxonomy_table = pd.read_csv(taxonomy_table, sep="\t", header=0, index_col=0)

# Map new and old index values
id_map = create_index_dict(count_table)

# Rename the index of the count and taxonomy tables
count_table = count_table.rename(index=id_map)
taxonomy_table = taxonomy_table.rename(index=id_map)

# Export the count table
count_table.to_csv(
    os.path.join(results_dir, "sequences_counts.tsv"),
    sep="\t",
    header=True,
    index=True,
)

# Format the taxonomy string in the taxonomy table
taxonomy_table["lineage"] = taxonomy_table["lineage"].str.replace(
    "\\w__", "", regex=True
)
taxonomy_table = taxonomy_table["lineage"].str.split(";", expand=True)
taxonomy_table.columns = [
    "Root",
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
]

# Export the taxonomy table
taxonomy_table.to_csv(
    os.path.join(results_dir, "sequences_taxonomy.tsv"),
    sep="\t",
    header=True,
    index=True,
)

# Parse the templates fasta and outputs the new fasta file
parsed_fasta_df = parse_fasta(
    templates_fasta, os.path.join(results_dir, "sequences.fasta"), id_map
)

# Render the sequence table
sequence_table = build_seq_table(parsed_fasta_df, id_map)
sequence_table.to_csv(
    os.path.join(results_dir, "sequences_table.tsv"), sep="\t", header=True, index=True
)


# Build the phylogenetic tree with the parsed sequences
build_tree(os.path.join(results_dir, "sequences.fasta"), results_dir, results_dir)
# Rename the output files
os.rename(
    os.path.join(results_dir, "templates_mafft.tree"),
    os.path.join(results_dir, "phylogenetic_tree.tree"),
)
os.rename(
    os.path.join(results_dir, "templates_mafft.fasta"),
    os.path.join(results_dir, "sequences_alinged.fasta"),
)


elapsed = timeit.default_timer() - start_time
print("Total output_parser.py time: ", elapsed)
