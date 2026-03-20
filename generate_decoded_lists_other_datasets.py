import os
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
from ont_fast5_api.fast5_interface import get_fast5_file
from my_utils import get_reference_alignment
from multiprocessing import Pool, Manager
from concurrent.futures import ProcessPoolExecutor, as_completed
import threading
from queue import Queue

parser = argparse.ArgumentParser(description='generate decoded lists from raw signal')
parser.add_argument('--fast5_file',type=str,required=True)
parser.add_argument('--bam_file',type=str,required=True)
parser.add_argument('--ref_file',type=str,required=True)
parser.add_argument('--out_prefix',type=str,required=True)
parser.add_argument('--mem_conv',type=int,required=True)
parser.add_argument('--msg_len',type=int,required=True)
parser.add_argument('--rate_conv',type=int,required=True)
parser.add_argument('--list_size',type=int,required=True)

parser.add_argument('--barcode_extend_penalty',type=float,default=0.6)
parser.add_argument('--num_threads',type=int,default=1)
parser.add_argument('--process_workers',type=int,default=4)  # New: number of parallel processes
parser.add_argument('--bonito_model_path',type=str)
parser.add_argument('--complexity_log_file',type=str,required=True)

args = parser.parse_args()
print(args)

# Global constants (for multiprocessing)
LIST_SIZE = args.list_size
NUM_THREADS = args.num_threads
OUT_PREFIX = args.out_prefix
FAST5_FILE = args.fast5_file
BAM_FILE = args.bam_file
REF_FILE = args.ref_file
PATH_TO_CPP_EXEC = "viterbi/viterbi_nanopore.out"
MSG_LEN = args.msg_len
MEM_CONV = args.mem_conv
RATE_CONV = args.rate_conv
CFILE = args.complexity_log_file

FAILED_BARCODE_FILE = OUT_PREFIX + "/failed_barcode_removal.txt"
FAILED_ALIGNMENT_FILE = OUT_PREFIX + "/failed_alignment.txt"

if args.bonito_model_path is None:
    bonito_model_path = 'dna_r9.4.1'
else:
    bonito_model_path = args.bonito_model_path

base_dict = {'A':0, 'C':1, 'G':2, 'T':3}
int_dict = {0:'A', 1:'C', 2:'G', 3:'T'}
primer_length = 25
barcode_search_extend_len = 25
index_len = 12
crc_len = 8

# Global lock for file writing (shared across processes)
file_lock = None

def init_worker(lock):
    """Initialize worker process with shared lock"""
    global file_lock
    file_lock = lock

def log_failed_barcode_removal(read_id):
    """Thread-safe function to log failed barcode removal"""
    global file_lock
    with file_lock:
        with open(FAILED_BARCODE_FILE, 'a') as f:
            f.write(f"{read_id}\n")

def log_failed_alignment(read_id):
    """Thread-safe function to log failed alignment"""
    global file_lock
    with file_lock:
        with open(FAILED_ALIGNMENT_FILE, 'a') as f:
            f.write(f"{read_id}\n")

def process_single_read(read_data, num_reads, cfile):
    """Process a single read - designed to be called by worker processes"""
    i, readid, raw_data = read_data
    
    print(f"Processing read {i}/{num_reads}: {readid}")
    
    decoded_filename = OUT_PREFIX + '/list_' + str(i)
    if os.path.isfile(decoded_filename):
        # Check message length
        with open(decoded_filename) as f:
            decoded_msg_list = [l.rstrip('\n') for l in f.readlines()]
            correct_msg = decoded_msg_list[-1]
        if len(correct_msg) == MSG_LEN:
            print(f'Decoded file {decoded_filename} already exists. Skipping.')
            return f"Skipped {i}"

    rnd = str(np.random.randint(10000000))

    try:
        # STEP 1: get random codeword from chosen convolutional code
        payload_bits = ''.join(rm.choice('01') for _ in range(MSG_LEN - index_len - crc_len))
        msg = helper.encode_with_crc(payload_bits)
        file_msg = f'tmp.{rnd}.{i}'
        file_seq = f'tmp.enc.{rnd}.{i}'
        
        with open(file_msg,'w') as f:
            f.write(msg+'\n')
        subprocess.run([PATH_TO_CPP_EXEC,'-m','encode','-i',file_msg,'-o',file_seq,
                       '--mem-conv',str(MEM_CONV),'-r',str(RATE_CONV),'--msg-len',str(MSG_LEN)])
        cw = helper.read_seq(file_seq)
        cw_len = len(cw)

        # STEP 2: get reference alignment
        alignment = get_reference_alignment(BAM_FILE, readid)

        if alignment is None or ((alignment and alignment['is_reverse'] == False and 
                                alignment['reference_end'] - alignment['reference_start'] > cw_len + 2 * primer_length) == False):
            log_failed_alignment(readid)
            cleanup_files([file_msg, file_seq])
            return f"Failed alignment check {i}"

        ref_start, ref_end = alignment['reference_start'], alignment['reference_end']
        
        if ref_start>= ref_end - cw_len - 2 * primer_length:
            cleanup_files([file_msg, file_seq])
            return f"Too short read"
        
        fwd_primer_start = rm.randint(ref_start, ref_end - cw_len - 2 * primer_length)

        # Load reference (each process needs its own handle)
        reference = pysam.Fastafile(REF_FILE)
        ref_sequence = reference.fetch(alignment["reference_name"],
                                     fwd_primer_start,
                                     fwd_primer_start + cw_len + 2 * primer_length)

        START_BARCODE, END_BARCODE = ref_sequence[:primer_length], ref_sequence[-primer_length:]
        payload = ref_sequence[primer_length : -primer_length]
        START_BARCODE_RC = helper.reverse_complement(END_BARCODE)
        END_BARCODE_RC = helper.reverse_complement(START_BARCODE)

        # STEP 4: decide random offset
        assert cw_len == len(payload)
        offset = [(base_dict[payload[j]] - base_dict[cw[j]]) % 4 for j in range(cw_len)]
        offset_str = ''.join([str(num) for num in offset])

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

        # truncate post according to barcode
        (start_pos, end_pos, dist_start, dist_end) = helper.find_barcode_pos_in_post(
            trans_filename, fastq_filename, START_BARCODE, END_BARCODE, 
            barcode_search_extend_len, args.barcode_extend_penalty)
        (start_pos_RC, end_pos_RC, dist_start_RC, dist_end_RC) = helper.find_barcode_pos_in_post(
            trans_filename, fastq_filename, START_BARCODE_RC, END_BARCODE_RC, 
            barcode_search_extend_len, args.barcode_extend_penalty)
        
        rc = False
        if dist_start + dist_end > dist_start_RC + dist_end_RC:
            rc = True
            start_pos = start_pos_RC
            end_pos = end_pos_RC

        if start_pos == -1 or end_pos - start_pos + 1 < MEM_CONV + MSG_LEN + 1:
            print(f"Failure in barcode removing for read {i}")
            log_failed_barcode_removal(readid)
            cleanup_files([fast5_filename, post_filename, trans_filename, fastq_filename, 
                          file_msg, file_seq])
            shutil.rmtree(fast5_dir, ignore_errors=True)
            return f"Failed barcode removal {i}"

        new_post_filename = f'tmp.{rnd}_{i}_{readid}.post.new'
        helper.truncate_post_file(post_filename, new_post_filename, start_pos, end_pos)

        rc_flag = ''
        if rc:
            rc_flag = '--rc'
            cleanup_files([fast5_filename, post_filename, new_post_filename, trans_filename, 
                          fastq_filename, file_msg, file_seq])
            shutil.rmtree(fast5_dir, ignore_errors=True)
            return f"RC flag set, skipped {i}"

        subprocess.run([PATH_TO_CPP_EXEC,'-m','decode','-i',new_post_filename,'-o',decoded_filename,
                       '--msg-len',str(MSG_LEN),'-l',str(LIST_SIZE),'-t',str(NUM_THREADS),
                       '--mem-conv',str(MEM_CONV),rc_flag,'--max-deviation','20','-r',str(RATE_CONV), 
                       '--offset', offset_str, '--complexity_log_file', cfile])

        with open(decoded_filename, 'a') as f:
            f.write(msg+'\n')  # write correct message on the last line

        # Cleanup temporary files
        cleanup_files([fast5_filename, post_filename, new_post_filename, trans_filename, 
                      fastq_filename, file_msg, file_seq])
        shutil.rmtree(fast5_dir, ignore_errors=True)
        
        reference.close()
        return f"Successfully processed {i}"

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

def main():
    # Load all reads first
    print("Loading all reads...")
    f_raw = h5py.File(FAST5_FILE, 'r')
    read_data_list = []

    with get_fast5_file(FAST5_FILE) as f5:
        for i, read in enumerate(f5.get_reads()):
            readid = read.read_id
            raw_data = read.get_raw_data()
            read_data_list.append((i, readid, raw_data))
    
    print(f"Found {len(read_data_list)} reads to process")

    num_reads = len(read_data_list)
    # Option 1: ProcessPoolExecutor (recommended)
    print(f"Starting parallel processing with {args.process_workers} workers...")
    with ProcessPoolExecutor(max_workers=args.process_workers) as executor:
        # Submit all jobs
        future_to_read = {executor.submit(process_single_read, read_data, num_reads, CFILE): read_data[0]
                         for read_data in read_data_list}
        
        # Process completed jobs
        for future in as_completed(future_to_read):
            read_idx = future_to_read[future]
            try:
                result = future.result()
                print(f"Read {read_idx}: {result}")
            except Exception as exc:
                print(f"Read {read_idx} generated an exception: {exc}")

if __name__ == "__main__":
    main()
