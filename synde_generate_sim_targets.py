import os, shutil, glob
from tqdm import tqdm
import math
import json
import pysam
import sys
import subprocess
import helper
import binascii
import struct
import numpy as np
import h5py
import warnings
import random
import argparse
from pathlib import Path
sys.path.append('/dss/dsshome1/03/ge73wor2/lokatt/seq_raw_mapping')
from local_aligner import local_aligner, align_top
from map_moves import map_moves_to_positions
from extract_segment import main_explicit, read_positions_file
from multiprocessing import Pool, Manager
from concurrent.futures import ProcessPoolExecutor, as_completed
import threading
from queue import Queue

CHOME = "/dss/dsshome1/03/ge73wor2"
PATH_TO_GUPPY = os.path.join(CHOME, "ont-guppy/bin/guppy_basecaller")
CFG_FILE='~/ont-guppy/data/dna_r9.4.1_450bps_hac.cfg'
PATH_TO_CPP_EXEC = "viterbi/viterbi_nanopore.out"
MODEL_STRIDE = 3
MAX_READS = 20000
MIN_READS = 10000
SIM_YEAR=2020
bonito_model_path = 'dna_r9.4.1'
barcode_search_extend_len = 25
barcode_extend_penalty = 0.6
model_stride = 3
# Global lock for file writing (shared across processes)
file_lock = None

def get_parser():
    parser = argparse.ArgumentParser(description='generate decoded lists from raw signal')
    parser.add_argument('--exp_num',type=int,required=True)
    parser.add_argument('--output_json_path',type=str,required=True)
    parser.add_argument('--num_threads',type=int,default=1)
    parser.add_argument('--ps_subsample',type=int,default=6)
    parser.add_argument('--ps_thresh',type=float,default=0.0)
    parser.add_argument('--model',type=str,default='dna_r9.4.1')
    args = parser.parse_args()

    return args

def get_available_gpus():
    """Return list of available GPU indices using nvidia-smi"""
    try:
        result = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=index", "--format=csv,noheader"],
            encoding="utf-8"
        )
        gpus = [int(x.strip()) for x in result.strip().split("\n") if x.strip().isdigit()]
        return gpus
    except Exception as e:
        print("WARNING: Could not detect GPUs, defaulting to GPU 0")
        return [0]


class GPUPool:
    """
    Thread/Process-safe GPU allocator.
    Each worker checks out a GPU and returns it when done.
    """
    def __init__(self, gpu_ids):
        self.queue = Queue()
        for gid in gpu_ids:
            self.queue.put(gid)

    def acquire(self):
        return self.queue.get()

    def release(self, gid):
        self.queue.put(gid)



def extract_seq(filename, lnum):
    with open(filename, 'r') as f:
        lines = f.readlines()
    if len(lines) >= lnum:
        seq = lines[lnum-1].strip()
    else:
        seq = 'C' * 110
    return seq

def read_fasta_to_dict(fa_file):
    ## Returns dictionary with sequence IDs as keys and sequences as values
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


def init_worker(lock):
    ### Initialize worker process with shared lock
    global file_lock
    file_lock = lock

def process_single_read_guppy(read_data, entry_data, BARCODES,num_reads, gpu_id):
    """Process a single read to obtain Guppy positions (only if rc=False)

    Args:
        read_data: tuple of (i, readid, raw_data)
        entry_data: dict with read metadata
        BARCODES: list of barcode sequences
        aligner_params: dict or tuple with aligner initialization parameters (instead of aligner object)
        num_reads: total number of reads (for logging)
        gpu_id: specific GPU ID assigned to this process
    """
    i = read_data[0]
    readid = read_data[1]
    raw_data = read_data[2]
    rc = entry_data.get('rc', False)

    # Set GPU for this process
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    print(f"Processing read (Guppy) {i}/{num_reads}: {readid}, GPU {gpu_id}")

    # If rc is True, set Guppy values to None
    if rc:
        return {
            "guppy_position": None,
            "guppy_aligned_primer_ref": None,
            "guppy_aligned_primer_pattern": None,
            "guppy_aligned_primer_match": None,
            "guppy_score": None,
        }

    START_BARCODE = BARCODES[0]
    END_BARCODE = BARCODES[1]
    PRIMER_LEN = len(START_BARCODE)

    rnd = str(np.random.randint(10000000))
    fast5_dir = f'tmp_input_{rnd}_{i}_{readid}'
    fast5_filename = os.path.join(fast5_dir, f'tmp.{rnd}.{i}.fast5')
    fastq_out = f'fastq_{rnd}_{i}_{readid}'

    try:
        aligner = local_aligner()

        # ================== STEP 1: create fast5 from raw signal ======================
        os.mkdir(fast5_dir)
        helper.create_fast5(raw_data, fast5_filename)

        # =========== STEP 2: get Guppy's position  =============
        # Use specific GPU (already set via CUDA_VISIBLE_DEVICES)
        guppy_cmd = (
            f"{PATH_TO_GUPPY} -i {fast5_dir} -c {CFG_FILE} "
            f"--moves_out --bam_out --save_path {fastq_out} -x cuda:0"  # Use cuda:0 since we set CUDA_VISIBLE_DEVICES
        )
        os.system(guppy_cmd)

        fastq_dir_pass = os.path.join(fastq_out, 'pass')
        fastq_filename2 = glob.glob(f"{fastq_dir_pass}/*.fastq")[0] if len(glob.glob(f"{fastq_dir_pass}/*.fastq")) > 0 else None
        bam_file = glob.glob(f"{fastq_dir_pass}/*.bam")[0] if len(glob.glob(f"{fastq_dir_pass}/*.bam")) > 0 else None

        if fastq_filename2 is None:
            fastq_dir_fail = os.path.join(fastq_out, 'fail')
            fastq_filename2 = glob.glob(f"{fastq_dir_fail}/*.fastq")[0] if len(glob.glob(f"{fastq_dir_fail}/*.fastq")) > 0 else None
            bam_file = glob.glob(f"{fastq_dir_fail}/*.bam")[0] if len(glob.glob(f"{fastq_dir_fail}/*.bam")) > 0 else None

        guppy_basecall = extract_seq(fastq_filename2, 2)

        # ------- Guppy primer position ---------
        with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam:
            for sam_record in bam:
                seq = sam_record.get_forward_sequence()
                currents_len = int(sam_record.get_tag("ns"))
                offset = int(sam_record.get_tag("ts"))
                mv = sam_record.get_tag("mv")
                stride = int(mv[0])
                moves = mv[1:]

                assert(sum(moves) == len(seq))
                read_positions = map_moves_to_positions(stride, offset, moves)
                assert(len(read_positions) == len(seq) + 1)

        alp = align_top(START_BARCODE, seq, aligner=aligner)
        score, aligned_seq1, pattern, aligned_seq2, _alignment = alp
        search_pattern = ''
        for t in aligned_seq2:
            if t != '-':
                search_pattern += t
        primer_seq_start = seq.find(search_pattern)
        primer_seq_end = primer_seq_start + PRIMER_LEN

        if primer_seq_start >= 0 and primer_seq_end < len(read_positions):
            guppy_base_positions = read_positions[primer_seq_start : primer_seq_end + 1]
        else:
            guppy_base_positions = [-1] * 20

        # Cleanup
        shutil.rmtree(fast5_dir, ignore_errors=True)
        shutil.rmtree(fastq_out, ignore_errors=True)

        return {
            "guppy_position": guppy_base_positions,
            "guppy_aligned_primer_ref": aligned_seq1,
            "guppy_aligned_primer_pattern": pattern,
            "guppy_aligned_primer_match": aligned_seq2,
            "guppy_score": score,
        }

    except Exception as e:
        print(f"Error processing read (Guppy) {i}: {str(e)}")
        import traceback
        traceback.print_exc()
        shutil.rmtree(fast5_dir, ignore_errors=True)
        shutil.rmtree(fastq_out, ignore_errors=True)
        return None  # Return None instead of error string for consistency

def cleanup_files(file_list):
    ## Helper function to clean up files safely
    for file_path in file_list:
        try:
            if os.path.exists(file_path):
                os.remove(file_path)
        except:
            pass

if __name__ == '__main__':
    args = get_parser()
    NUM_THREADS = args.num_threads
    SUBSAMPLE=args.ps_subsample
    FAST5_DIR = os.path.expanduser("~/fast5_datasets/stanford_bonito")
    READ_ID_FILE_PATH=os.path.join(FAST5_DIR, f"20200814_read_ids")

    if args.exp_num == 0:
        START_BARCODE = "GCTACATGTATACTGCGAGACAGAC"
        END_BARCODE = "CGATAGTCGCAGTCGCACATCACTC"
        PAYLOAD_LEN = 114
        MSG_LEN = 165
        MEM_CONV = 6
    elif args.exp_num == 1:
        START_BARCODE = "TCTATCTACTCGTGCTCGCTAGCTG"#"TCTATCTACTCGTGCTCGCTAGCTG"
        END_BARCODE = "TGTCTGCACTGCACTAGTCGCATGT"#"ACATGCGACTAGTGCAGTGCAGACA"
        PAYLOAD_LEN = 115
        MSG_LEN = 164
        MEM_CONV = 8
    elif args.exp_num == 2:
        START_BARCODE = "TGAGATCACAGCTACATAGTGAGAG"#"TGAGATCACAGCTACATAGTGAGAG"
        END_BARCODE = "GATAGAGCACTGTCGTGTGCAGTCA"#"TGACTGCACACGACAGTGCTCTATC"
        PAYLOAD_LEN = 117
        MSG_LEN = 164
        MEM_CONV = 11
    elif args.exp_num == 3:
        START_BARCODE="AGCGTACACGACTGAGCACACTACG"
        END_BARCODE="ATGCGCTCACACGCATCTGCTAGCA"
        PAYLOAD_LEN=112
        MSG_LEN = 180
        MEM_CONV = 6
    elif args.exp_num == 4:
        START_BARCODE = "CATCAGCAGTAGAGAGTAGCGCGAT"#"CATCAGCAGTAGAGAGTAGCGCGAT"
        END_BARCODE = "TATCATCGACGCTAGCAGTGTCTGC"#"GCAGACACTGCTAGCGTCGATGATA"
        PAYLOAD_LEN = 113
        MSG_LEN = 180
        MEM_CONV = 8
    elif args.exp_num == 5:
        START_BARCODE = "GAGTCTCTAGCGCTACGAGATATAT"#"GAGTCTCTAGCGCTACGAGATATAT"
        END_BARCODE = "AGCTCTGATGAGATCAGCAGACTGT"#"ACAGTCTGCTGATCTCATCAGAGCT"
        PAYLOAD_LEN = 115
        MSG_LEN = 180
        MEM_CONV = 11
    elif args.exp_num == 6:
        START_BARCODE = "CACGAGATCTCAGTGTCGACACGTG"
        END_BARCODE = "AGTGATAGTGATCTCGACGCAGCTA"
        PAYLOAD_LEN = 116
        MSG_LEN = 197
        MEM_CONV = 6
    else:
        print("Error! BARCODES not provided!")
        quit()

    basename = f"{SIM_YEAR}_raw_signal_{args.exp_num}"
    HDF5_FILE = os.path.join(FAST5_DIR, f"{basename}.hdf5")
    REF_FILE = os.path.join(FAST5_DIR, f"oligos_8_4_20/oligos_{args.exp_num}.fa")
    READ_ID_FILE=os.path.join(READ_ID_FILE_PATH, f"exp_{args.exp_num}/read_ids.txt")
    output_json_file = os.path.expanduser(os.path.join(args.output_json_path, f"{basename}.json"))

    assert os.path.exists(READ_ID_FILE), "READ_ID_FILE does not exist!"
    with open(READ_ID_FILE) as f:
        all_read_ids = [l.rstrip('\n') for l in f.readlines()]

    START_BARCODE_RC = helper.reverse_complement(END_BARCODE)
    END_BARCODE_RC = helper.reverse_complement(START_BARCODE)
    BARCODES = [START_BARCODE, END_BARCODE, START_BARCODE_RC, END_BARCODE_RC]

    # Check if output file exists and has all required fields
    skip_processing = False
    if os.path.exists(output_json_file):
        with open(output_json_file, 'r') as f:
            results = json.load(f)
        if len(results) >= MIN_READS:
            required_fields = [
                "guppy_position",
                "guppy_aligned_primer_ref",
                "guppy_aligned_primer_pattern",
                "guppy_aligned_primer_match",
                "guppy_score",
                "ctc_basecall_lev_primer_position"
            ]

            all_entries_complete = all(
                all(field in entry for field in required_fields)
                for entry in results
            )

            if all_entries_complete:
                skip_processing = True
                print(f"Output file has {len(results)} entries with all required fields. Skipping h5py processing.")
                quit()

    print("-------------", output_json_file, "----------------")
    print('START_BARCODE: ', START_BARCODE, ', END_BARCODE: ', END_BARCODE)
    print('START_BARCODE_RC: ', START_BARCODE_RC, ', END_BARCODE_RC: ', END_BARCODE_RC)

    # ============= PHASE 1: Generate entries with Chandak et al's primer postions =================
    print("\n" + "=" * 80)
    print("PHASE 1: Create new entries and add primer positions a la Chandak et al.")
    print("=" * 80)
    bonito_cmd = f'bonito gen_entries {args.model} {HDF5_FILE} --start_barcode {START_BARCODE} --end_barcode {END_BARCODE} --min_len {MSG_LEN+MEM_CONV} --output_json_file {output_json_file} --device cuda'
    os.system(bonito_cmd)

    # =================== PHASE 2: Add PrimerSeeker positions =========================
    print("\n" + "=" * 80)
    print("PHASE 2: add primer positions a la PrimerSeeker")
    print("params: shift=100, beam=8, subsample=6, start_offset=2, concentration threshold=0.98, examines 50% of raw read")
    print("=" * 80)
    bonito_cmd = f'bonito primer_search {args.model} {HDF5_FILE} --sim_targets {output_json_file} --results_file {output_json_file} --version opt --shift 100 --beam 8 --subsample 6 --conc_thresh 0.98 --examine_fraction 0.5 --start_offset 2'
    os.system(bonito_cmd)

    num_chunks = 200

    # ============== PHASE 3: Add Guppy positions to existing entries (Parallel with CPU Guppy) ================
    gpu_ids = get_available_gpus()
    #print(f"Detected GPUs: {gpu_ids}")
    gpu_pool = GPUPool(gpu_ids)
    with h5py.File(HDF5_FILE, 'r') as f5:
        print(f"Total reads in file: {len(all_read_ids)}")
        print("\n" + "=" * 80)
        print("PHASE 3: Adding Guppy positions")
        print("=" * 80)

        # Reload results to get all CTC entries
        with open(output_json_file, 'r') as f:
            results = json.load(f)

        # Filter entries that need Guppy processing
        entries_needing_guppy = [entry for entry in results if 'guppy_position' not in entry]

        print(f"Found {len(entries_needing_guppy)} entries needing Guppy processing")

        for chunk_start in range(0, len(entries_needing_guppy), num_chunks):
            chunk_end = min(chunk_start + num_chunks, len(entries_needing_guppy))
            chunk_entries = entries_needing_guppy[chunk_start:chunk_end]
            print(f"\nProcessing Guppy chunk {chunk_start}-{chunk_end} ({len(chunk_entries)} reads)...")

            read_data_list = []
            for entry in chunk_entries:
                readid = entry['read_id']
                i = entry['read_index']
                raw_data = f5[readid]['raw_signal'][:]
                read_data_list.append(((i, readid, raw_data), entry))

            print(f"Processing {len(read_data_list)} reads (Guppy) with {NUM_THREADS} workers...")

            with ProcessPoolExecutor(max_workers=NUM_THREADS) as executor:
                future_to_entry = {}
                for idx, read_data in enumerate(read_data_list):
                    # Assign GPU in round-robin fashion
                    gpu_id = gpu_ids[idx % len(gpu_ids)]

                    future = executor.submit(
                        process_single_read_guppy,
                        read_data[0],
                        read_data[1],
                        BARCODES,
                        len(read_data_list),
                        gpu_id  # Pass specific GPU ID
                    )
                    future_to_entry[future] = read_data[1]

                for future in as_completed(future_to_entry):
                    entry = future_to_entry[future]
                    try:
                        guppy_result = future.result(timeout=600)  # 10 min timeout
                        if guppy_result is not None and isinstance(guppy_result, dict):
                            entry.update(guppy_result)
                            print(f"Read {entry['read_id']}: Guppy Success")
                    except Exception as exc:
                        print(f"Read {entry['read_id']} generated an exception: {exc}")
                        import traceback
                        traceback.print_exc()

            # Save results after each chunk
            print(f"Saving updated results...")
            with open(output_json_file, 'w') as f:
                json.dump(results, f, indent=4)
