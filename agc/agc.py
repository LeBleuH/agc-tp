#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
import textwrap
from pathlib import Path
from collections import Counter
from typing import Iterator, Dict, List
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Mingjie HUANG"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Mingjie HUANG"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Mingjie HUANG"
__email__ = "mingjie.huang@etu.u-paris.fr"
__status__ = "Developpement"



def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for OTU clustering (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="Kmer size for OTU clustering (default 8)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=Path,
                        default=Path("OTU.fasta"), help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    """Read a compressed fasta and extract all fasta sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :return: A generator object that provides the Fasta sequences (str).
    """
    with gzip.open(amplicon_file, "rt") as monfich:
        sequence_parts = []

        for line in monfich:
            line = line.strip()

            if line.startswith(">"):
                if sequence_parts:
                    sequence = "".join(sequence_parts)
                    if len(sequence) >= minseqlen:
                        yield sequence

                sequence_parts = []

            elif line:
                sequence_parts.append(line)
        
        # Process the last sequence
        if sequence_parts:
            sequence = "".join(sequence_parts)
            if len(sequence) >= minseqlen:
                yield sequence


def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count] of sequence with a count >= mincount and a length >= minseqlen.
    """
    sequences = Counter(read_fasta(amplicon_file, minseqlen))

    for sequence, count in sorted(
        sequences.items(),
        key=lambda item: item[1],
        reverse=True
    ):
        if count >= mincount:
            yield [sequence, count]

def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list:  (list) A list of aligned sequences in the format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """
    sequence1, sequence2 = alignment_list

    if len(sequence1) != len(sequence2):
        raise ValueError("Aligned sequences must have the same length.")
    
    if len(sequence1) == 0:
        return 0.0

    identical = sum(
        nt1 == nt2
        for nt1, nt2 in zip(sequence1,sequence2)
    )

    return identical / len(sequence1) * 100.0

def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int, mincount: int, chunk_size: int, kmer_size: int) -> List:
    """Compute an abundance greedy clustering regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """
    otu_list = []

    # Dereplicated sequences are sorted by decreasing abundance
    sequences = dereplication_fulllength(
        amplicon_file, 
        minseqlen, 
        mincount
        )
    
    for sequence, count in sequences:
        is_similar = False

        for otu_sequence, otu_count in otu_list:
            alignment = nw.global_align(
                sequence,
                otu_sequence,
                gap_open=-1,
                gap_extend=-1,
                matrix=str(Path(__file__).parent / "MATCH")
            )

            identity = get_identity(alignment)

            if identity > 97:
                is_similar = True
                break
        
        # Create a new OTU if no representative is sufficiently similar
        if not is_similar:
            otu_list.append([sequence, count])
    
    return otu_list


def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    with open(output_file, "w", encoding="utf-8") as output:
        for index, (sequence, count) in enumerate(OTU_list, start=1):
            output.write(
                f">OTU_{index} occurrence:{count}\n"
            )
            output.write(textwrap.fill(sequence, width=80) + "\n")


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Perform abundance greedy clustering
    otu_list = abundance_greedy_clustering(
        args.amplicon_file,
        args.minseqlen,
        args.mincount,
        args.chunk_size,
        args.kmer_size
    )

    # Write OTU sequences to file
    write_OTU(
        otu_list, 
        args.output_file
    )



if __name__ == '__main__':
    main()
