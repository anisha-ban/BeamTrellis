import os
import json
import subprocess
import helper
import binascii
import struct
import random as rm
import numpy as np
import h5py
import warnings
import random
import argparse
import shutil
import h5py
import pysam
from multiprocessing import Pool, Manager
from concurrent.futures import ProcessPoolExecutor, as_completed
from ont_fast5_api.fast5_interface import get_fast5_file
import threading
from queue import Queue
from my_utils import get_reference_alignment

bonito_model_path = 'dna_r9.4.1'
model_stride = 3

barcode_search_extend_len = 25
barcode_extend_penalty = 0.6

# Global lock for file writing (shared across processes)
file_lock = None

def init_worker(lock):
    """Initialize worker process with shared lock"""
    global file_lock
    file_lock = lock

def process_single_read(read_data_list, i, sim_entry, REF_FILE, BAM_FILE, PAYLOAD_LEN, num_reads):
    """Process a single read - designed to be called by worker processes"""
    readid = sim_entry["read_id"]
    print(f"Processing read {i}/{num_reads}: {readid}")
    raw_data = read_data_list[readid]
    START_BARCODE = sim_entry["fwd_primer"]
    END_BARCODE_RC = helper.reverse_complement(START_BARCODE)
    primer_length = len(START_BARCODE)
    rnd = str(np.random.randint(10000000))

    try:
        # ================== STEP 1: get raw signal ======================
        # create fast5 from raw data
        fast5_dir = f'tmp_input_{rnd}_{i}_{readid}'
        os.mkdir(fast5_dir)
        fast5_filename = os.path.join(fast5_dir, f'tmp.{rnd}.{i}.fast5')
        helper.create_fast5(raw_data, fast5_filename)

        # call bonito to generate CTC posterior table
        post_filename = f'tmp.{rnd}_{i}_{readid}.post'
        trans_filename = f'tmp.{rnd}_{i}_{readid}.trans'
        fastq_filename = f'tmp.{rnd}_{i}_{readid}.fastq'

        subprocess.run(['bonito','basecaller', bonito_model_path, fast5_dir,
                    '--post_file', post_filename, '--device','cuda:0'])

        # convert bonito post file to fastq and move (trans) files
        helper.bonito_basecall_generate_move(post_filename, fastq_filename, trans_filename)

        # get reference sequence, START_BARCODE and END_BARCODE
        alignment = get_reference_alignment(BAM_FILE, readid)
        ref_start, ref_end = alignment['reference_start'], alignment['reference_end']
        reference = pysam.Fastafile(REF_FILE)
        ref_sequence = reference.fetch(alignment["reference_name"],
                                     ref_start, ref_end)
        prim_start = ref_sequence.find(START_BARCODE)
        ref_sequence = ref_sequence[prim_start:prim_start+PAYLOAD_LEN+2*primer_length]
        END_BARCODE = ref_sequence[-primer_length:]
        START_BARCODE_RC = helper.reverse_complement(END_BARCODE)
        r"""
        if alignment is None or ((alignment and alignment['is_reverse'] == False and
                                alignment['reference_end'] - alignment['reference_start'] > PAYLOAD_LEN + 2 * primer_length + 30) == False):
            log_failed_alignment(readid)
            cleanup_files([fast5_filename, post_filename, trans_filename,
                        fastq_filename])
            #return f"Failed alignment check {i}"
            return None

        ref_start, ref_end = alignment['reference_start'], alignment['reference_end']

        if ref_start >= ref_end - PAYLOAD_LEN - 2 * primer_length:
            cleanup_files([fast5_filename, post_filename, trans_filename,
                        fastq_filename])
            #return f"Too short read"
            return None
        """

        # truncate post according to barcode
        (start_pos, end_pos, dist_start, dist_end) = helper.find_barcode_pos_in_post_ps(
            trans_filename, fastq_filename, START_BARCODE, END_BARCODE,
            barcode_search_extend_len, barcode_extend_penalty)
        (start_pos_RC, end_pos_RC, dist_start_RC, dist_end_RC) = helper.find_barcode_pos_in_post_ps(
            trans_filename, fastq_filename, START_BARCODE_RC, END_BARCODE_RC,
            barcode_search_extend_len, barcode_extend_penalty)

        cleanup_files([fast5_filename, post_filename, trans_filename,
                        fastq_filename])
        shutil.rmtree(fast5_dir, ignore_errors=True)

        if dist_start + dist_end > dist_start_RC + dist_end_RC:
            #rc = True
            return None

        if start_pos == -1 or end_pos - start_pos + 1 < PAYLOAD_LEN + 1:
            #print(f"Failure in barcode removing for read {i}")
            return None

        sim_entry["ctc_basecall_lev_primer_position"] = start_pos * model_stride
        return sim_entry

    except Exception as e:
        print(f"Error processing read {i}: {str(e)}")
        # Attempt cleanup on error
        try:
            cleanup_files([f'tmp.{rnd}.{i}', f'tmp.enc.{rnd}.{i}',
                        f'tmp.{rnd}_{i}_{readid}.post', f'tmp.{rnd}_{i}_{readid}.trans',
                        f'tmp.{rnd}_{i}_{readid}.fastq', f'tmp.{rnd}_{i}_{readid}.post.new'])
            shutil.rmtree(f'tmp_input_{rnd}_{i}_{readid}', ignore_errors=True)
        except:
            pass
        return f"Error processing {i}: {str(e)}"

def cleanup_files(file_list):
    """Helper function to clean up files safely"""
    for file_path in file_list:
        try:
            if os.path.exists(file_path):
                os.remove(file_path)
        except:
            pass

def get_read_data_list(RAW_FILE):
    print("Loading all reads...")
    read_data_list = []
    if RAW_FILE.endswith(".hdf5"):
        f5 = h5py.File(RAW_FILE, 'r')
        for i, readid in enumerate(f5.keys()):
            #ref_id = f5[readid].attrs['ref'].decode("utf-8")
            raw_data = f5[readid]['raw_signal']
            read_data_list.append((readid, raw_data))
    elif RAW_FILE.endswith(".fast5"):
        with get_fast5_file(RAW_FILE) as f5:
            for i, read in enumerate(f5.get_reads()):
                readid = read.read_id
                raw_data = read.get_raw_data()
                read_data_list.append((readid, raw_data))
    return read_data_list

def add_pos(RAW_FILE, REF_FILE, BAM_FILE, PAYLOAD_LEN, SIM_FILE, NUM_THREADS=4):

    # Load all reads first
    read_data_list = get_read_data_list(RAW_FILE)
    read_data_dict = {readid: raw_data for readid, raw_data in read_data_list}
    print(f"Found {len(read_data_list)} reads to process")
    num_reads = len(read_data_list)

    with open(SIM_FILE, 'r') as jf:
        sim_data = json.load(jf)

    results = []
    # Option 1: ProcessPoolExecutor (recommended)
    print(f"Starting parallel processing with {NUM_THREADS} workers...")
    with ProcessPoolExecutor(max_workers=NUM_THREADS) as executor:
        # Submit all jobs
        future_to_read = {executor.submit(process_single_read, read_data_dict, i, sim_entry, REF_FILE, BAM_FILE, PAYLOAD_LEN, num_reads):i  for i, sim_entry in enumerate(sim_data)}

        # Process completed jobs
        for future in as_completed(future_to_read):
            read_idx = future_to_read[future]
            try:
                result = future.result()
                if result is not None:
                    results.append(result)
                    #print(f"Read {read_idx}: Success")
                else:
                    print(f"Read {read_idx}: Skipped")
            except Exception as exc:
                print(f"Read {read_idx} generated an exception: {exc}")

        # Write results to JSON file
        print(f"\nWriting {len(results)} results to {SIM_FILE}")
        with open(SIM_FILE, 'w') as f:
            json.dump(results, f, indent=4)

r"""
def main():
    # Load all reads first
    read_data_list = get_read_data_list(RAW_FILE)
    read_data_dict = {readid: raw_data for readid, raw_data in read_data_list}
    print(f"Found {len(read_data_list)} reads to process")
    num_reads = len(read_data_list)

    with open(SIM_FILE, 'r') as jf:
        sim_data = json.load(jf)

    results = []
    # Option 1: ProcessPoolExecutor (recommended)
    print(f"Starting parallel processing with {args.process_workers} workers...")
    with ProcessPoolExecutor(max_workers=args.process_workers) as executor:
        # Submit all jobs
        future_to_read = {executor.submit(process_single_read, read_data_dict, i, sim_entry, num_reads):i  for i, entry in enumerate(sim_data)}

        # Process completed jobs
        for future in as_completed(future_to_read):
            read_idx = future_to_read[future]
            try:
                result = future.result()
                if result is not None:
                    results.append(result)
                    #print(f"Read {read_idx}: Success")
                else:
                    print(f"Read {read_idx}: Skipped")
            except Exception as exc:
                print(f"Read {read_idx} generated an exception: {exc}")

        # Write results to JSON file
        print(f"\nWriting {len(results)} results to {SIM_FILE}")
        with open(SIM_FILE, 'w') as f:
            json.dump(results, f, indent=4)

        print(f"Done! Results saved to {SIM_FILE}")

"""
