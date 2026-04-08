import numpy as np
import argparse
from collections import defaultdict
from helper_functions import Make_folder_if_not_exists

        
    
def make_obs_rate(obs_file, out_file, bin_size):
    
    # Determine contig offsets and window counts
    kmer_dict = defaultdict(lambda: defaultdict(int))
    contig_lengths = defaultdict(int)
    with open(obs_file) as data:
        for line in data:
            contig, start, end, count = line.split('\t')
            window = int(start) - int(start) % bin_size
            kmer_dict[contig][window] += int(count)
            contig_lengths[contig] = max(contig_lengths[contig], int(end))

    # Calculate observation rates
    kmer_arr = []
    assembly_positions = []
    for contig in kmer_dict:
        last_window = max(kmer_dict[contig]) + bin_size
        last_window = int(last_window)

        for window in range(0, last_window, bin_size):
            kmer_arr.append(kmer_dict[contig][window])
            # Store actual window end, capped at contig length
            actual_end = min(window + bin_size, contig_lengths[contig])
            assembly_positions.append([contig, window, actual_end])

    # Calculate average k-mer count per base across the assembly
    kmer_arr = np.array(kmer_arr)
    assembly_length = sum(contig_lengths.values())
    assembly_avg = np.sum(kmer_arr) / assembly_length

    # Write output
    Make_folder_if_not_exists(out_file)
    with open(out_file, 'w') as out:
        out.write('contig\tstart\tend\tobs_rate\n')
        for position, kmer_count in zip(assembly_positions, kmer_arr):
            contig, start, end = position
            actual_window_size = end - start  # use real span, not nominal window_size
            obs_rate = round(kmer_count / actual_window_size / assembly_avg, 5)
            out.write(f'{contig}\t{start}\t{end}\t{obs_rate}\n')


def print_script_usage():
    toprint = f'''
    Hidden Markov Model for archaic human introgression inference using the ARCkmerFinder output.

    Usage:
    python observation_rate.py -obs [bed_file]
    python observation_rate.py -obs [bed_file] -out [output_file] -bin_size [bin_size]
    
    > Estimate observation rate
        -obs            Input file with observation data (required)
        -out            Output file with observation rates estimates (default: 'observation_rate.bed')
        -bin_size       Parameter defining size of bins (default: 1 Mb)
        
    '''

    return toprint


def main():
    parser = argparse.ArgumentParser(description=print_script_usage(), formatter_class=argparse.RawTextHelpFormatter)
    
    mutation_rate = parser.add_parser('Observation rate', help='Estimate the observation rate')
    mutation_rate.add_argument("-obs",help="Input file with observation data (required)", type=str, required = True)
    mutation_rate.add_argument("-out", metavar='',help="Output file with observation rates estimates (default: 'observation_rate.bed')", default = 'observation_rate.bed')
    mutation_rate.add_argument("-bin_size", metavar='',help="Parameter defining size of bins (default: 1 Mb)", type=int, default = 1_000_000)
    
    args = parser.parse_args()
    
    if hasattr(args, 'obs'):
        print('-' * 40)
        print(f'> Observations file: {args.obs}')
        print(f'> Estimated observation rates written to: {args.out}')
        print(f'> Bin size: {args.bin_size}')
        print('-' * 40)

        print('Estimating observation rates...')
        make_obs_rate(args.obs, args.out, args.bin_size)
        print('Done')
    else:
        print(print_script_usage())


if __name__ == "__main__":
    main()

