#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import os
import sys
import gzip
from mimetypes import guess_type
import argparse
import subprocess
from functools import partial

import pandas as pd
from Bio import SeqIO


def run_minimap2(reference, query, output_file, query_extension,
                 max_loaded_bases):
    """Execute minimap2 to map short sequences to a reference genome.

    Parameters
    ----------
    reference : str
        Path to the FASTA file with the reference.
    map_fasta : str
        Path to FASTA file with the short sequences
        to map against the reference.
    output_file : str
        Path to the output file with mapping results.

    Returns
    -------
    List with following elements:
        stdout : list
            List with stdout from minimap2 in bytes.
        stderr : list
            List with the stderr from minimpa2 in bytes.
    """
    # -I parameter to control number of target bases loaded into memory
    # --secondary=yes to output secondary alignments that might have high score
    minimap_args = ['minimap2 -I {0} --cs -cx sr --secondary=yes {1} {2} > '
                    '{3}'.format(max_loaded_bases, reference, query, output_file)]

    minimap_proc = subprocess.Popen(minimap_args,
                                    shell=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)

    stderr = minimap_proc.stderr.readlines()
    stdout = minimap_proc.stdout.readlines()

    return [stdout, stderr]


reference = '/home/rmamede/Desktop/pneumo_reads/T2T-CHM13v2/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna'
query = ['/home/rmamede/Desktop/pneumo_reads/test_data/PCR351_1.fastq.gz']
output_directory = '/home/rmamede/Desktop/pneumo_reads/test_results'
max_loaded_bases = '200M'
def main(reference, query, output_directory, max_loaded_bases, no_cleanup):

    # Create output directory
    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)
    else:
        sys.exit('Output directory already exists.')

    # Query might include multiple files
    # Get basename for all queries
    # Need to detect file format
    query_basenames = {}
    query_extensions = {}
    for file in query:
        basename = os.path.basename(file)
        # Remove file extension
        basename, extension = os.path.splitext(basename)
        compressed = False
        if extension == '.gz':
            basename, extension = os.path.splitext(basename)
            compressed = True
        query_basenames[basename] = file
        # Remove dot from extension
        query_extensions[basename] = [extension[1:], compressed]

    print('Mapping...', end='')
    paf_files = {}
    for k, v in query_basenames.items():
        # Map reads against reference
        output_paf = os.path.join(output_directory, f'{k}.paf')
        stdout, stderr = run_minimap2(reference, v, output_paf,
                                      query_extensions[k], max_loaded_bases)
        paf_files[k] = output_paf
    print('done.')

    # Get list of identifiers from mapping results
    mapped_ids = {}
    for k, v in paf_files.items():
        current_ids = pd.read_csv(v, delimiter='\t', usecols=[0], names=['seqid'])
        # Set to get distinct identifiers
        current_ids = set(current_ids['seqid'].tolist())
        print(f'{len(current_ids)} sequences mapped against reference.')
        # Save identifiers that mapped against reference
        mapped_ids_outfile = os.path.join(output_directory,
                                          f'{k}_mapped_ids.txt')
        with open(mapped_ids_outfile, 'w') as outfile:
            outfile.write('\n'.join(current_ids)+'\n')
        mapped_ids[k] = mapped_ids_outfile

    print('Writing filtered seqs...', end='')
    for k, v in query_basenames.items():
        # Read list of ids that mapped against reference
        with open(mapped_ids[k], 'r') as infile:
            query_mapped = infile.readlines()

        # Read and filter records in query file
        # Check if file is compressed
        _open = partial(gzip.open, mode='rt') if query_extensions[k][1] is True else open
        with _open(v) as infile:
            records_iterator = SeqIO.parse(infile, query_extensions[k][0])
            filtered_records = (record
                                for record in records_iterator
                                if record.id not in query_mapped)

            filtered_records_outfile = os.path.join(output_directory,
                                                    f'{k}_filtered.{query_extensions[k][0]}')
            SeqIO.write(filtered_records, filtered_records_outfile,
                        f'{query_extensions[k][0]}')
    print('done.')


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-r', '--reference', type=str,
                        required=True, dest='reference',
                        help='Path to the FASTA file with the '
                             'reference/contaminant sequences to map against.')

    parser.add_argument('-q', '--query', type=str, nargs='+',
                        required=True, dest='query',
                        help='Path to FASTA or/and FASTQ files with the '
                             'sequences that will be mapped against the '
                             'reference.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to the output directory where the output '
                             'files will be created.')

    parser.add_argument('-I', '--max-loaded-bases', type=str,
                        required=False, default='200M',
                        dest='max_loaded_bases',
                        help='Maximum number of bases loaded by minimap2.')

    parser.add_argument('--no-cleanup', required=False, action='store_true',
                        dest='no_cleanup',
                        help='Keep PAF files with mapping results.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
