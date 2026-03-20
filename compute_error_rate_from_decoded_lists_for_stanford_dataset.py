import helper
import os
import argparse


def return_filtered_FERs_stanford_dataset(DECODED_LISTS_DIR, CONV_INPUT_FILE, code_params, MAX_LIST_SIZE,pad,bytes_per_oligo):
    with open(CONV_INPUT_FILE) as f:
        conv_input_list = [s.rstrip('\n') for s in f.readlines()]
    num_oligos = len(conv_input_list)

    list_size_arr = [2**i for i in range(6)]
    list_size_arr = [el for el in list_size_arr if el <= MAX_LIST_SIZE]

    num_list_sizes = len(list_size_arr)
    max_list_size = list_size_arr[-1]

    percent_correct_CRC = [0] * num_list_sizes
    percent_erased_CRC = [0] * num_list_sizes
    percent_incorrect_CRC = [0] * num_list_sizes
    FER_filtered = [0] * num_list_sizes

    num_reads = 0
    num_correct = [0] * num_list_sizes
    num_erasure_CRC_index = [0] * num_list_sizes
    num_error_CRC_index = [0] * num_list_sizes
    decoded_index_dict = [{}] * num_list_sizes

    for filename in os.listdir(DECODED_LISTS_DIR):
        if not filename.startswith("list_"):
            continue
        num_reads += 1
        with open(os.path.join(DECODED_LISTS_DIR,filename)) as f:
            decoded_msg_list = [l.rstrip('\n') for l in f.readlines()]
            decoded_msg_list = decoded_msg_list[:MAX_LIST_SIZE]
        for i, list_size in enumerate(list_size_arr):
            decoded_msg_list_cropped = decoded_msg_list[:list_size]
            (index, payload_bytes, decoded_msg) = helper.decode_list_CRC_index(decoded_msg_list_cropped,bytes_per_oligo,num_oligos,pad)
            if index == None:
                num_erasure_CRC_index[i] += 1
            else:
                if index in decoded_index_dict[i]:
                    found = False
                    for tup in decoded_index_dict[i][index]:
                        if tup[0] == decoded_msg:
                            tup[1] += 1
                            found = True
                            break
                    if not found:
                        decoded_index_dict[i][index].append([decoded_msg,1])
                    decoded_index_dict[i][index] = sorted(decoded_index_dict[i][index],key=lambda x: -x[1])
                else:
                    decoded_index_dict[i][index] = [[decoded_msg,1]]
                if decoded_msg == conv_input_list[index]:
                    num_correct[i] += 1
                else:
                    num_error_CRC_index[i] += 1

    for i in range(num_list_sizes):
        percent_correct_CRC[i] = num_correct[i] * 100 / num_reads
        percent_erased_CRC[i] = num_erasure_CRC_index[i] * 100 / num_reads
        percent_incorrect_CRC[i] = num_error_CRC_index[i] * 100 / num_reads
        FER_filtered[i] = percent_incorrect_CRC[i] / (percent_correct_CRC[i] + percent_incorrect_CRC[i])

    return percent_erased_CRC, FER_filtered


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='compute error rate from decoded list')
    parser.add_argument('--DECODED_LISTS_DIR',type=str,required=True)
    parser.add_argument('--list_size',type=int,required=True)
    parser.add_argument('--CONV_INPUT_FILE',type=str,required=False, default=2020)
    parser.add_argument('--pad',type=bool,required=False, default=False)
    args = parser.parse_args()
    pad=args.pad

    # PARAMETERS TO BE SET BEFORE RUNNING THE CODE
    ROOT_PATH="/dss/dsshome1/03/ge73wor2"
    # ----------- [4,3,8] ---------------------
    LIST_SIZE = args.list_size
    DECODED_LISTS_DIR = os.path.expanduser(args.DECODED_LISTS_DIR)
    CONV_INPUT_FILE = os.path.expanduser(args.CONV_INPUT_FILE)
    bytes_per_oligo = 18

    #print('list size:',LIST_SIZE)
    with open(CONV_INPUT_FILE) as f:
        conv_input_list = [s.rstrip('\n') for s in f.readlines()]
    num_oligos = len(conv_input_list)
    #print('num_oligos',num_oligos)

    num_reads = 0
    num_correct = 0
    num_erasure_CRC_index = 0
    num_error_CRC_index = 0
    decoded_index_dict = {}
    # map from index to [[decoded_msg,count]]

    for filename in os.listdir(DECODED_LISTS_DIR):
        if not filename.startswith("list_"):
            continue
        num_reads += 1
        with open(os.path.join(DECODED_LISTS_DIR,filename)) as f:
            decoded_msg_list = [l.rstrip('\n') for l in f.readlines()]
            decoded_msg_list = decoded_msg_list[:LIST_SIZE]
        (index, payload_bytes, decoded_msg) = helper.decode_list_CRC_index(decoded_msg_list,bytes_per_oligo,num_oligos,pad)
        if index == None:
            num_erasure_CRC_index += 1
        else:
            if index in decoded_index_dict:
                found = False
                for tup in decoded_index_dict[index]:
                    if tup[0] == decoded_msg:
                        tup[1] += 1
                        found = True
                        break
                if not found:
                    decoded_index_dict[index].append([decoded_msg,1])
                decoded_index_dict[index] = sorted(decoded_index_dict[index],key=lambda x: -x[1])
            else:
                decoded_index_dict[index] = [[decoded_msg,1]]
            if decoded_msg == conv_input_list[index]:
                num_correct += 1
            else:
                num_error_CRC_index += 1

    #print('num_reads:',num_reads)
    #print('num_correct:',num_correct)
    #print('num_erasure_CRC_index:',num_erasure_CRC_index)
    #print('num_error_CRC_index:',num_error_CRC_index)
    percent_discarded=num_erasure_CRC_index*100/(num_reads)
    FER=num_error_CRC_index/(num_error_CRC_index+num_correct)
    print(f'At list size: {LIST_SIZE}, FER={FER} with percent discarded={percent_discarded}% for {num_reads} reads')
