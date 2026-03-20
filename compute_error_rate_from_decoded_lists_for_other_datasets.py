import helper
import os
import argparse

def return_filtered_FERs(DECODED_LISTS_DIR, MAX_LIST_SIZE, pad):
    list_size_arr = [2**i for i in range(6)]
    list_size_arr = [el for el in list_size_arr if el <= MAX_LIST_SIZE]

    percent_correct_CRC = [0] * len(list_size_arr)
    percent_erased_CRC = [0] * len(list_size_arr)
    percent_incorrect_CRC = [0] * len(list_size_arr)
    FER_filtered = [0] * len(list_size_arr)

    for i, LIST_SIZE in enumerate(list_size_arr):
        num_reads, num_correct, num_erasure_CRC_index, num_error_CRC_index = process_directory(DECODED_LISTS_DIR, pad)

        percent_correct_CRC[i] = num_correct * 100 / num_reads
        percent_erased_CRC[i] = num_erasure_CRC_index * 100 / num_reads
        percent_incorrect_CRC[i] = num_error_CRC_index * 100 / num_reads
        FER_filtered[i] = percent_incorrect_CRC / (percent_correct_CRC + percent_incorrect_CRC)

    return percent_incorrect_CRC, FER_filtered

def process_directory(DECODED_LISTS_DIR, pad):
    num_reads = 0
    num_correct = 0
    num_erasure_CRC_index = 0
    num_error_CRC_index = 0
    decoded_index_dict = {}

    for filename in os.listdir(DECODED_LISTS_DIR):
            if not filename.startswith("list_"):
                continue
            num_reads += 1

            with open(os.path.join(DECODED_LISTS_DIR,filename)) as f:
                decoded_msg_list = [l.rstrip('\n') for l in f.readlines()]
                correct_msg = decoded_msg_list[-1]
                decoded_msg_list = decoded_msg_list[:-1]
                decoded_msg_list = decoded_msg_list[:LIST_SIZE]

            crc_check, decoded_msg = helper.validate_CRC(decoded_msg_list, pad)
            if crc_check == False:
                num_erasure_CRC_index += 1
            else:
                if decoded_msg == correct_msg:
                    num_correct += 1
                else:
                    num_error_CRC_index += 1
    return num_reads, num_correct, num_erasure_CRC_index, num_error_CRC_index

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='compite error rate from decoded list')
    parser.add_argument('--list_size',type=int,required=True)
    parser.add_argument('--dec_dir',type=str,required=True)
    args = parser.parse_args()

    # warnings.filterwarnings("ignore", category=DeprecationWarning)

    # PARAMETERS TO BE SET BEFORE RUNNING THE CODE
    MAX_LIST_SIZE = args.list_size # 8
    DECODED_LISTS_DIR = os.path.expanduser(args.dec_dir) #'../nanopore_dna_storage_data/decoded_lists/exp_7/'
    pad = False

    print('max list size:', MAX_LIST_SIZE)

    list_size_arr = [2**i for i in range(6)]
    list_size_arr = [el for el in list_size_arr if el <= MAX_LIST_SIZE]


    for i, LIST_SIZE in enumerate(list_size_arr):
        num_reads, num_correct, num_erasure_CRC_index, num_error_CRC_index = process_directory(DECODED_LISTS_DIR, pad)
        print('list size:',LIST_SIZE)
        print('num_reads:',num_reads)
        print('num_correct:',num_correct)
        print('num_erasure_CRC_index:',num_erasure_CRC_index)
        print('num_error_CRC_index:',num_error_CRC_index)
        print('overall FER: ', (num_erasure_CRC_index + num_error_CRC_index) * 100 / num_reads, '%')
        print('filtered FER: ', (num_error_CRC_index) * 100 / (num_error_CRC_index + num_correct), '% after discarding ', num_erasure_CRC_index * 100/ num_reads, '%')





