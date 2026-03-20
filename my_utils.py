
import numpy as np
import pysam

#

def diffs(x):
  x = np.array(x)
  diff = list()
  for i in range(len(x)-1):
    diff.append(np.abs(x[i+1]-x[i]))
  diff = np.array(diff)
  return diff

def estimate_mu(signals,window_size=4,diff_thre=0.4):
  signals = signals.reshape(-1)
  sig_means = means_windowing(signals, window_size)#median
  differences = diffs(sig_means)
  duration = len(signals)/np.sum(differences > diff_thre)
  if duration > 20:
    duration = 20
  mu = 0.0927*duration+1.1454
  return mu

def means_windowing(x,size=4):
  x = np.array(x)
  x_ = np.zeros((x.size-size+1,size))
  for i in range(size-1):
    x_[:,i] = x[i:-size+i+1]
  x_[:,-1] = x[size-1:]
  means = np.median(x_,axis=1)
  return means

# reference extraction
def get_reference_alignment(bam_file, read_id):
  with pysam.AlignmentFile(bam_file, "rb") as bam:
    for read in bam.fetch():
      if read.query_name == read_id:
        return {
                    "reference_name": read.reference_name,
                    "reference_start": read.reference_start,
                    "reference_end": read.reference_end,
                    "cigar_string": read.cigarstring,
                    "mapping_quality": read.mapping_quality,
                    "is_reverse": read.is_reverse,
                        }
  return None


def read_fasta_to_dict(fa_file):
#Returns dictionary with sequence IDs as keys and sequences as values

    sequences = {}
    current_id = None
    current_sequence = []

    with open(fa_file, 'r') as f:
        for line in f:
            line = line.strip()

            if line.startswith('>'):
                # If we have a previous sequence, store it
                if current_id is not None:
                    sequences[current_id] = ''.join(current_sequence)

                current_id = line[1:]  # Remove '>'
                current_sequence = []

            elif line:  # Non-empty line (sequence data)
                current_sequence.append(line)

    # last sequence
    if current_id is not None:
        sequences[current_id] = ''.join(current_sequence)

    return sequences

def find_primer_start(bam_file, read_id, reference, primer_sequence):
  # Get read alignment details
  alignment = get_reference_alignment(bam_file, read_id)
  if not alignment:
    return None  # Read not found

  # Extract reference sequence in the read-aligned region
  ref_seq = reference.fetch(alignment["reference_name"],
                            alignment["reference_start"],
                            alignment["reference_end"])

  # Find primer in this sequence
  primer_index = ref_seq.find(primer_sequence)
  if primer_index == -1:
    return None  # Primer not found in this region

  primer_start = alignment["reference_start"] + primer_index
  return primer_start



