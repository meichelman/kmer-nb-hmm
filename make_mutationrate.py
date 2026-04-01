import numpy as np
from collections import defaultdict

from helper_functions import Make_folder_if_not_exists

def make_mutation_rate(obs_file, out_file, window_size):
    kmer_dict = defaultdict(lambda: defaultdict(int))
    contig_lengths = defaultdict(int)
    with open(obs_file) as data:
        for line in data:
            contig, start, end, count = line.split('\t')
            window = int(start) - int(start) % window_size
            kmer_dict[contig][window] += int(count)
            contig_lengths[contig] = int(end)
    
    kmer_arr = []
    assembly_positions = []
    for contig in kmer_dict:
        lastwindow = max(kmer_dict[contig]) + window_size
        lastwindow = int(lastwindow)

        for window in range(0, lastwindow, window_size):
            kmer_arr.append(kmer_dict[contig][window])
            assembly_positions.append([contig, window, window + window_size])

    kmer_arr = np.array(kmer_arr)
    assembly_length = sum(contig_lengths.values())
    assembly_avg = np.sum(kmer_arr) / np.sum(assembly_length)

    Make_folder_if_not_exists(out_file)
    with open(out_file, 'w') as out:
        out.write('contig\tstart\tend\tmutationrate\n')
        for pos, kmer_count in zip(assembly_positions, kmer_arr):
            contig, start, end = pos
            mutrate = round(kmer_count / window_size / assembly_avg, 5)
            out.write(f'{contig}\t{start}\t{end}\t{mutrate}\n')