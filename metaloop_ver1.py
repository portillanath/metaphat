import sys, os
from Bio import SeqIO
import numpy as np
import pandas as pd
from pathlib import Path
import networkx as nx
from networkx import DiGraph
import random
import matplotlib.pyplot as plt


path = Path("./raw_data/")
glob_path = path.glob("metagenoma_12000_1000.fastq")
reads = []


# This function reads the fastq file and creates an array with the sequences
for file_path in glob_path:
    with open(file_path, "r") as f:
        for line in f:
            if line[0] == '@':
                id=line.strip()
            if line[0]== "A" or line[0] == "T" or line[0]=="C" or line[0]=="G":
                seq=line.strip()
                reads.append(seq)
kmers = []
kmer_coverage = {}

def get_random_substrings(reads, a, b, n):
    """
    Genera n substrings de cada string en la lista reads que tengan un tamaño aleatorio, uniformemente distribuido, entre a y b.
    Calcula la cobertura de cada kmer en los substrings generados y devuelve un diccionario con la cobertura de cada kmer.
    
    Args:
    reads: La lista de strings de la que se van a generar los substrings.
    a: El tamaño mínimo de los substrings.
    b: El tamaño máximo de los substrings.
    n: El número de substrings que se van a generar por string.
    
    Returns:
    Un diccionario con la cobertura de cada kmer en los substrings generados.
    """
    for read in reads:
        read_length = len(read)
        for i in range(n):
            ksize = random.randint(a, b)
            initial = random.randint(0, read_length - ksize)
            final = initial + ksize
            kmer = read[initial:final]
            kmers.append(kmer)
            if kmer not in kmer_coverage:
                kmer_coverage[kmer] = 1
            else:
                kmer_coverage[kmer] += 1
    total_kmers = len(kmers)
    for kmer in kmer_coverage:
        kmer_coverage[kmer] /= total_kmers
    return kmer_coverage

get_random_substrings(reads, 15, 17, 25)

kmers = list(kmer_coverage.keys())
kmer_coverage = list(kmer_coverage.values())

with open("kmers.txt", "w") as f:
    # Write each kmer to a new line in the file
    for kmer in kmers:
        f.write(kmer + "\n")

