import numpy as np
import re
import argparse

parser = argparse.ArgumentParser(description="Creates a new alignment by removing all positions where the ratio of gaps to nongaps is less than the specified ratio.")

parser.add_argument('-p','--percent-nongapped',  help='Remove positions where fewer than p% of the position are not gaps')
parser.add_argument('-i','--input',  help='Alignment')
parser.add_argument('-o','--output',  help='Position-filtered alignment')

args = parser.parse_args()

p = float(args.percent_nongapped) / 100.


# Find positions to retain by summing gaps (0) and nongaps (1)
n_seqs = 0
nongap = re.compile('\w') # Regular expression for alphanumeric
count_nongaps = []

with open(args.input) as fh:
    for line in fh.readlines():
        line = line.strip()

        # For each new seq
        if '>' in line:

            n_seqs += 1

            # Store
            if len(count_nongaps) > 1:
                if n_seqs == 2: # First item in array
                    nongaps_sum = np.array(count_nongaps)
                elif n_seqs > 2:
                    nongaps_sum = nongaps_sum + np.array(count_nongaps)
                else:
                    pass            

            # Restart list
            count_nongaps = []

        else:

            # Position is amino acid (1) or not (0)
            for i in line:
                if nongap.match(i):
                    count_nongaps.append(1)
                else:
                    count_nongaps.append(0)

# Boolean for positions where nongaps are greater than or equal to the cutoff ratio
positions_to_retain = (nongaps_sum / n_seqs) >= p

# Number of expected positions in new and old alignments
m_in = positions_to_retain.shape[0]
m_out = positions_to_retain[positions_to_retain].shape[0]

# Filter alignment and write to new file
with open(args.output, 'w') as aln:

    split = 80 # Characters before creating new line

    with open(args.input) as fh1:

        for line in fh1.readlines():
            line = line.strip()

            if '>' in line:

                # New seq
                gene_id = line
                sequence = []

            else:

                # Extend sequences
                for i in line:
                    sequence.append(i)

            # If sequence is full length
            if len(sequence) == m_in:

                # Join string, split for new lines, and write to new file
                sequence = ''.join(np.array(sequence)[positions_to_retain])
                splits = []
                for i in range(0, int(m_out / split)):
                    splits.append(''.join(sequence)[int(i*split) : int((i+1)*split)])           

                aln.write('\n'.join([gene_id ]+ splits))
                aln.write('\n')

