#!/usr/bin/env python3

from sys import argv
import os
import gzip
import io
import pickle
from argparse import ArgumentParser
import itertools


def parse_args():
    "Parse the input arguments, use '-h' for help."
    parser = ArgumentParser(description='SalmID - rapid Kmer based Salmonella identifier from raw data')
    # inputs
    parser.add_argument(
        '-i','--input_file', type=str, required=True, metavar = 'your_fastqgz or fasta',
        help='Single fastq.gz or fasta file input, include path to file if file is not in same directory ')
    parser.add_argument(
        '-t', '--type', type=str, required=False, default='fastq.gz', metavar = 'fastq.gz or assembly',
        help='Type of data; "fastq.gz" (default), "assembly", "minion_fasta", "minion_fastq"')
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

def minion_fasta_kmerizer(file, k, kmerDict):
    i = 1
    n_reads = 0
    total_coverage = 0
    target_mers = {}
    for line in open(file):
        if i % 2 == 0:
            for kmer, rc_kmer in kmers(line.strip().upper(), k):
                if (kmer in kmerDict) or (rc_kmer in kmerDict):
                    if kmer in target_mers:
                        target_mers[kmer] += 1
                    else:
                        target_mers[kmer] = 1
                    if rc_kmer in target_mers:
                        target_mers[rc_kmer] += 1
                    else:
                        target_mers[rc_kmer] = 1
        #if i == 20:
         #   break
        i += 1
    return set([h for h in target_mers])

def minion_fastq_kmerizer(file, k, kmerDict):
    i = 1
    n_reads = 0
    total_coverage = 0
    target_mers = {}
    for line in open(file):
        if i % 4 == 2:
            for kmer, rc_kmer in kmers(line.strip().upper(), k):
                if (kmer in kmerDict) or (rc_kmer in kmerDict):
                    if kmer in target_mers:
                        target_mers[kmer] += 1
                    else:
                        target_mers[kmer] = 1
                    if rc_kmer in target_mers:
                        target_mers[rc_kmer] += 1
                    else:
                        target_mers[rc_kmer] = 1
        #if i == 20:
         #   break
        i += 1
    return set([h for h in target_mers])

def multifasta_single_string2(multifasta):
    single_string = ''
    with open(multifasta, 'r') as f:
        for line in f:
            if line.strip()[0] == '>':
                pass
            else:
                single_string += line.strip()
    return single_string

def kmers(seq, k):
    rev_comp = reverse_complement(seq)
    for start in range(1,len(seq) - k + 1):
        yield seq[start:start + k], rev_comp[-(start + k):-start]

def multifasta_to_kmers_dict(multifasta):
    multi_seq_dict = multifasta_dict(multifasta)
    lib_dict = {}
    for h in multi_seq_dict:
        lib_dict[h] = set([k for k in createKmerDict_reads([multi_seq_dict[h]], 27)])
    return lib_dict

def Combine(b,c):
  fliC_combinations=[]
  fliC_combinations.append(",".join(c))
  temp_combinations=[]
  for i in range(len(b)):
    for x in itertools.combinations(b,i+1):
      temp_combinations.append(",".join(x))
  for x in temp_combinations:
    temp=[]
    for y in c:
      temp.append(y)
    temp.append(x)
    temp=",".join(temp)
    temp=temp.split(",")
    temp.sort()
    temp=",".join(temp)
    fliC_combinations.append(temp)
  return fliC_combinations


def seqsero_from_formula_to_serotypes(Otype,fliC,fljB,special_gene_list):
  #like test_output_06012017.txt
  #can add more varialbles like sdf-type, sub-species-type in future (we can conclude it into a special-gene-list)
  from Initial_Conditions import phase1
  from Initial_Conditions import phase2
  from Initial_Conditions import phaseO
  from Initial_Conditions import sero
  seronames=[]
  for i in range(len(phase1)):
    fliC_combine=[]
    fljB_combine=[]
    if phaseO[i]==Otype:
      ### for fliC, detect every possible combinations to avoid the effect of "["
      if phase1[i].count("[")==0:
        fliC_combine.append(phase1[i])
      elif phase1[i].count("[")>=1:
        c=[]
        b=[]
        if phase1[i][0]=="[" and phase1[i][-1]=="]" and phase1[i].count("[")==1:
          content=phase1[i].replace("[","").replace("]","")
          fliC_combine.append(content)
          fliC_combine.append("-")
        else:
          for x in phase1[i].split(","):
            if "[" in x:
              b.append(x.replace("[","").replace("]",""))
            else:
              c.append(x)
          fliC_combine=Combine(b,c) #Combine will offer every possible combinations of the formula, like f,[g],t: f,t  f,g,t
      ### end of fliC "[" detect
      ### for fljB, detect every possible combinations to avoid the effect of "["
      if phase2[i].count("[")==0:
        fljB_combine.append(phase2[i])
      elif phase2[i].count("[")>=1:
        d=[]
        e=[]
        if phase2[i][0]=="[" and phase2[i][-1]=="]" and phase2[i].count("[")==1:
          content=phase2[i].replace("[","").replace("]","")
          fljB_combine.append(content)
          fljB_combine.append("-")
        else:
          for x in phase2[i].split(","):
            if "[" in x:
              d.append(x.replace("[","").replace("]",""))
            else:
              e.append(x)
          fljB_combine=Combine(d,e)
      ### end of fljB "[" detect
      new_fliC=fliC.split(",") #because some antigen like r,[i] not follow alphabetical order, so use this one to judge and can avoid missings
      new_fliC.sort()
      new_fliC=",".join(new_fliC)
      new_fljB=fljB.split(",")
      new_fljB.sort()
      new_fljB=",".join(new_fljB)
      if (new_fliC in fliC_combine or fliC in fliC_combine) and (new_fljB in fljB_combine or fljB in fljB_combine):
        seronames.append(sero[i])
  #analyze seronames
  if len(seronames)==0:
    seronames=["N/A (The predicted antigenic profile does not exist in the White-Kauffmann-Le Minor scheme)"]
  star=""
  star_line=""
  if len(seronames)>1:#there are two possible predictions for serotypes
    star="*"
    star_line="The predicted serotypes share the same general formula:\t"+Otype+":"+fliC+":"+fljB+"\n"##
  #print ("\n")
  predict_form=Otype+":"+fliC+":"+fljB#
  predict_sero=(" or ").join(seronames)
  ###special test for Enteritidis
  if predict_form=="9:g,m:-":
    sdf="-"
    for x in special_gene_list:
      if x.startswith("sdf"):
        sdf="+"
    predict_form=predict_form+" Sdf prediction:"+sdf
    if sdf=="-":
      star="*"
      star_line="Additional characterization is necessary to assign a serotype to this strain.  Commonly circulating strains of serotype Enteritidis are sdf+, although sdf- strains of serotype Enteritidis are known to exist. Serotype Gallinarum is typically sdf- but should be quite rare. Sdf- strains of serotype Enteritidis and serotype Gallinarum can be differentiated by phenotypic profile or genetic criteria.\n"#+##
      predict_sero="Gallinarum/Enteritidis sdf -"
  ###end of special test for Enteritidis
  elif predict_form=="4:i:-":
    predict_sero="potential monophasic variant of Typhimurium"
  elif predict_form=="4:r:-":
    predict_sero="potential monophasic variant of Heidelberg"
  elif predict_form=="4:b:-":
    predict_sero="potential monophasic variant of Paratyphi B"
  elif predict_form=="8:e,h:1,2":
    predict_sero="Newport"
    star="*"
    star_line="Serotype Bardo shares the same antigenic profile with Newport, but Bardo is exceedingly rare."
  claim="The serotype(s) is/are the only serotype(s) with the indicated antigenic profile currently recognized in the Kauffmann White Scheme.  New serotypes can emerge and the possibility exists that this antigenic profile may emerge in a different subspecies.  Identification of strains to the subspecies level should accompany serotype determination; the same antigenic profile in different subspecies is considered different serotypes."##
  if "N/A" in predict_sero:
    claim=""
  if "Typhimurium" in predict_sero or predict_form=="4:i:-":
    normal=0
    mutation=0
    for x in special_gene_list:
      if "oafA-O-4_full" in x:
        normal = float(special_gene_list[x])
      elif "oafA-O-4_5-" in x:
        mutation = float(special_gene_list[x])
    if normal>mutation:
      #print "$$$Typhimurium"
      pass
    elif normal<mutation:
      predict_sero=predict_sero.strip()+"(O5-)"
      star="*"#
      star_line="Detected the deletion of O5-."
      #print "$$$Typhimurium_O5-"
    else:
      #print "$$$Typhimurium, even no 7 bases difference"
      pass
  return predict_form,predict_sero,star,star_line,claim

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
    if data_type == 'minion_2d_fasta':
        input_Ks = minion_fasta_kmerizer(input_file, 27, set(kmers))
    if data_type == 'minion_2d_fastq':
        input_Ks = minion_fastq_kmerizer(input_file, 27, set(kmers))
    #print(len(input_Ks))
    #keep lists of O, H and special genes; both list of headers and
    O_dict = {}
    H_dict = {}
    Special_dict = {}
    for h in lib_dict:
        score = (len(lib_dict[h] & input_Ks)/ len(lib_dict[h])) * 100
        if output_mode == 'debug':
            print(h, score)
        if score > 1: # Arbitrary cut-off for similarity score very low but seems necessary to detect O-3,10 in some cases
            if h.startswith('O') and score > 15:
                O_dict[h] = score
            if h.startswith('fl') and score > 40:
                H_dict[h] = score
            if (h[:2] != 'fl') and (h[0] != 'O'):
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
                        O9 = float(Special_dict[z])
                    if "tyr-O-2" in z:
                        O2 = float(Special_dict[z])
                if O2 > O9:
                    highest_O = "O-2"
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
                        #print('highest O:', highest_O)
                if highest_O == "O-1,3,19":
                    highest_O = '-'
                    max_score = 0
                    for x in O_dict:
                        if x == 'O-1,3,19_not_in_3,10__130':
                            pass
                        else:
                            if float(O_dict[x]) >= max_score:
                                max_score = float(O_dict[x])
                                highest_O = x.split("_")[0]
            except:
                pass
    #call_fliC:
    highest_fliC = '-'
    highest_fliC_raw = '-'
    highest_Score = 0
    for s in H_dict:
        if s.startswith('fliC'):
            if float(H_dict[s]) > highest_Score:
                highest_fliC = s.split('_')[1]
                highest_fliC_raw = s
                highest_Score = float(H_dict[s])
    #call_fljB
    highest_fljB = '-'
    highest_fljB_raw = '-'
    highest_Score = 0
    for s in H_dict:
        if s.startswith('fljB'):
            if float(H_dict[s]) > highest_Score:
                highest_fljB = s.split('_')[1]
                highest_fljB_raw = s
                highest_Score = float(H_dict[s])
    if output_mode == 'debug':
        print(input_file+'\t'+highest_O.split('-')[1] +':'+ highest_fliC + ':' + highest_fljB)
    result = seqsero_from_formula_to_serotypes(highest_O.split('-')[1],highest_fliC,highest_fljB,Special_dict)
    print(input_file+'\t' + result[0] + '\t' +result[1])


if __name__ == '__main__':
    main()
