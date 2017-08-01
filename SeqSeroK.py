#!/usr/bin/env python3


from sys import argv
import os
import gzip
import io
import pickle
from argparse import ArgumentParser


def parse_args():
    "Parse the input arguments, use '-h' for help."
    parser = ArgumentParser(description='SalmID - rapid Kmer based Salmonella identifier from raw data')
    # inputs
    parser.add_argument(
        '-i','--input_file', type=str, required=True, metavar = 'your_fastqgz or fasta',
        help='Single fastq.gz or fasta file input, include path to file if file is not in same directory ')
    parser.add_argument(
        '-t', '--type', type=str, required=False, default='fastq.gz', metavar = 'file_extension',
        help='Type of data, "fastq.gz" (default) or "assembly"')
    parser.add_argument(
        '-m', '--mode', type=str, required=False, default='normal', metavar = 'normal or debug',
        help='When "debug" is chosen, the full list of matches for all alleles will be given')
    return parser.parse_args()


def reverse_complement(sequence):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'M': 'K', 'R': 'Y', 'W': 'W',
                            'S': 'S', 'Y': 'R', 'K': 'M', 'V': 'B', 'H': 'D', 'D': 'H', 'B': 'V'}
    return "".join(complement[base] for base in reversed(sequence))




def createKmerDict_reads(list_of_strings, kmer):
    kmer_table = {}
    for string in list_of_strings:
        sequence = string.strip('\n')
        for i in range(len(sequence)-kmer+1):
            new_mer =sequence[i:i+kmer].upper()
            new_mer_rc = reverse_complement(new_mer)
            if new_mer in kmer_table:
                kmer_table[new_mer.upper()] += 1
            else:
                kmer_table[new_mer.upper()] = 1
            if new_mer_rc in kmer_table:
                kmer_table[new_mer_rc.upper()] += 1
            else:
                kmer_table[new_mer_rc.upper()] = 1
    return kmer_table

def multifasta_dict(multifasta):
    multifasta_list = [line.strip() for line in open(multifasta, 'r') if len(line.strip()) > 0]
    headers = [i for i in multifasta_list if i[0] == '>']
    multifasta_dict = {}
    for h in headers:
        start = multifasta_list.index(h)
        for element in multifasta_list[start + 1:]:
            if element[0] == '>':
                break
            else:
                if h[1:] in multifasta_dict:
                    multifasta_dict[h[1:]] += element
                else:
                    multifasta_dict[h[1:]] = element
    return multifasta_dict

def multifasta_single_string(multifasta):
    multifasta_list = [line.strip() for line in open(multifasta, 'r') if (len(line.strip()) > 0) and (line.strip()[0] != '>')]
    return ''.join(multifasta_list)

def target_read_kmerizer(file, k, kmerDict):
    i = 1
    n_reads = 0
    total_coverage = 0
    target_mers = []
    for line in io.BufferedReader(gzip.open(file)):
        start = int((len(line) - k) // 2)
        if i % 4 == 2:
            s1 = line[start:k + start].decode()
            if s1 in kmerDict:
                n_reads += 1
                total_coverage += len(line)
                target_mers += [k for k in createKmerDict_reads([str(line.decode())], k)]
        i += 1
        if total_coverage >= 40000:
            break
    return set(target_mers)

def multifasta_single_string2(multifasta):
    single_string = ''
    with open(multifasta, 'r') as f:
        for line in f:
            if line.strip()[0] == '>':
                pass
            else:
                single_string += line.strip()
    return single_string

def multifasta_to_kmers_dict(multifasta):
    multi_seq_dict = multifasta_dict(multifasta)
    lib_dict = {}
    for h in multi_seq_dict:
        lib_dict[h] = set([k for k in createKmerDict_reads([multi_seq_dict[h]], 27)])
    return lib_dict

def main():
    #todo clean up def main, write some functions
    args = parse_args()
    input_file = args.input_file
    data_type = args.type
    output_mode = args.mode
    #print(len(input_fasta_ks))
    ex_dir = os.path.dirname(os.path.realpath(__file__))
    #print(ex_dir)
    try:
        f = open(ex_dir +'/antigens.pickle', 'rb')
        lib_dict = pickle.load(f)
        f.close
    except FileNotFoundError:
        lib_dict = multifasta_to_kmers_dict(ex_dir + '/H_and_O_and_specific_genes.fasta')
        f = open(ex_dir +'/antigens.pickle', "wb")
        pickle.dump(lib_dict, f)
        f.close()
        #print('Created lib')
    #kmers = [lib_dict[h] for h in lib_dict]
    kmers=[]
    for h in lib_dict:
        kmers += lib_dict[h]
    if data_type == 'assembly':
        input_Ks = set([k for k in createKmerDict_reads([multifasta_single_string(input_file)], 27)])
    if data_type == 'fastq.gz':
        input_Ks = target_read_kmerizer(input_file, 27, set(kmers))
    #print(len(input_Ks))
    #keep lists of O, H and special genes; both list of headers and
    O_dict = {}
    H_dict = {}
    Special_dict = {}
    for h in lib_dict:
        score = (len(lib_dict[h] & input_Ks)/ len(lib_dict[h])) * 100
        if output_mode == 'debug':
            print(h, score)
        if score > 15: # Arbitrary cut-off for similarity score very low but seems necessary to detect O-3,10 in some cases
            if h.startswith('O') and score > 15:
                O_dict[h] = score
            if h.startswith('fl') and score > 70:
                H_dict[h] = score
            if (h[:2] != 'fl') or (h[0] != 'O'):
                Special_dict[h] = score

    #call O:
    highest_O = '-'
    if len(O_dict) == 0:
        pass
    else:
        if 'O-9,46_wbaV__1002' in O_dict:  # not sure should use and float(O9_wbaV)/float(num_1) > 0.1
            if 'O-9,46_wzy__1191' in O_dict:  # and float(O946_wzy)/float(num_1) > 0.1
                highest_O = "O-9,46"
                # print "$$$Most possilble Otype:  O-9,46"
            elif "O-9,46,27_partial_wzy__1019" in O_dict:  # and float(O94627)/float(num_1) > 0.1
                highest_O = "O-9,46,27"
                # print "$$$Most possilble Otype:  O-9,46,27"
            else:
                highest_O = "O-9"  # next, detect O9 vs O2?
                O2 = 0
                O9 = 0
                for z in Special_dict:
                    if "tyr-O-9" in z:
                        O9 = Special_dict[z]
                    if "tyr-O-2" in z:
                        O2 = Special_dict[z]
                if O2 > O9:
                    highest_O = "O-2"
                elif O2 < O9:
                    pass
                else:
                    pass
        elif ("O-3,10_wzx__1539" in O_dict) and (
            "O-9,46_wzy__1191" in O_dict):  # and float(O310_wzx)/float(num_1) > 0.1 and float(O946_wzy)/float(num_1) > 0.1
            if "O-3,10_not_in_1,3,19__1519" in O_dict:  # and float(O310_no_1319)/float(num_1) > 0.1
                highest_O = "O-3,10"
                # print "$$$Most possilble Otype:  O-3,10 (contain O-3,10_not_in_1,3,19)"
            else:
                highest_O = "O-1,3,19"
                # print "$$$Most possilble Otype:  O-1,3,19 (not contain O-3,10_not_in_1,3,19)"
        ### end of special test for O9,46 and O3,10 family
        else:
            try:
                max_score = 0
                for x in O_dict:
                    if float(O_dict[x]) >= max_score:
                        max_score = float(O_dict[x])
                        highest_O = x.split("_")[0]
                if highest_O == "O-1,3,19":
                    #print(O_dict)
                    # always the second ?
                    O_list = [h for h in O_dict]
                    highest_O = O_list[1].split("_")[0]
                    # print "$$$Most possilble Otype: ",O_choice
            except:
                pass
    #call_fliC:
    highest_fliC = '-'
    highest_Score = 0
    for s in H_dict:
        if s.startswith('fliC'):
            if float(H_dict[s]) > highest_Score:
                highest_fliC = s.split('_')[1]
                highest_Score = float(H_dict[s])
    #call_fljB
    highest_fljB = '-'
    highest_Score = 0
    for s in H_dict:
        if s.startswith('fljB'):
            if float(H_dict[s]) > highest_Score:
                highest_fljB = s.split('_')[1]
                highest_Score = float(H_dict[s])
    print(input_file+'\t'+highest_O +':'+ highest_fliC + ':' + highest_fljB)

if __name__ == '__main__':
    main()