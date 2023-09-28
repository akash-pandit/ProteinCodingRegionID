#!/usr/bin/env python3

import sys
import os

# passed input: filepaths

sequence, protein = sys.stdin.read().split()

if sequence.split('.')[1] not in ['fasta', 'fna', 'ffn', 'faa', 'frn', 'fa']:
    print("Invalid nucleotide sequence file, please use a FASTA file.")
    exit(1)
if protein.split('.')[1] not in ['fasta', 'fna', 'ffn', 'faa', 'frn', 'fa']:
    print("Invalid protein sequence file, please use a FASTA file.")
    exit(1)