#!/usr/bin/env python
## Swift Biosciences 16S snapp workflow - Adapted by Bionformatics Service @ Genyo (Granada, Spain)
## Author (GitHub: jmgs7) 12112024

import sys
import os
import timeit


def run_fastqc(fastq_dir, output_dir=None):
    """
    Runs FastQC on fastq files in a given directory.

    Parameters
    ----------
    fastq_dir : str
        Directory containing the fastq files to be processed.
    output_dir : str, optional
        Directory to write the output reports. If not specified, a directory named
        "QC" will be created under fastq_dir.

    Returns
    -------
    output_dir : str
        The directory where the output files were written.
    """

    if not output_dir:
        output_dir = os.path.join(fastq_dir, "QC")
        os.makedirs(output_dir, exist_ok=True)
    else:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

    os.system(
        f"fastqc -q -t {os.environ.get("THREADS")} -o {output_dir} {fastq_dir}/*.fastq.gz"
    )  # run fastqc on fastq files in fastq_dir
    return output_dir


def run_multiqc(fastqc_dir, output_dir=None):
    """
    Runs MultiQC on FastQC output files in a given directory.

    Parameters
    ----------
    fastqc_dir : str
        Directory containing the output files from FastQC.
    output_dir : str, optional
        Directory to write the output report. If not specified, a directory named
        "multiqc" will be created under fastqc_dir.

    Returns
    -------
    output_dir : str
        The directory where the output file was written.
    """

    if not output_dir:
        output_dir = os.path.join(fastqc_dir, "multiqc")
        os.makedirs(output_dir, exist_ok=True)
    else:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

    os.system(
        f"multiqc -f -q -o {output_dir} {fastqc_dir}"
    )  # run multiqc on fastqc files in fastqc_dir


##main##
if __name__ == "__main__":
    if not len(sys.argv) == 2:
        print(len(sys.argv))
        print("QC.py fastq_directory output_directory(optional)")
        sys.exit()

start_time = timeit.default_timer()

fastq_dir = sys.argv[1]
try:
    output_dir = sys.argv[2]
except IndexError:
    output_dir = None

fastqc_dir = run_fastqc(fastq_dir, output_dir)
run_multiqc(fastqc_dir)

elapsed = timeit.default_timer() - start_time
print("Total QC.py time: ", elapsed)
