#!/usr/bin/env python

# ./Mapping_miRNA_with_bowtie.py -i Input_folder/ -r miRNA_sequence.fa -o Output_folder/

# bowtie parameter: --norc -l 7 -n 2 -a -f
# All scripts must be present in same directory!
# Scripts used: FastQcollapse.py, parse_bowtie_result.py, count_miRNA.py, extract_bowtie_unmapped_reads.py

import subprocess
import os
import sys
import argparse
import pandas as pd
from multiprocessing import Pool

def get_input_file(file_list, dir, format=['.fq', '.fastq']):
    out_list = []
    if file_list:
        dir = ''
    else:
        if dir[-1] != '/': dir += '/'
        file_list = os.listdir(dir)
    for file in file_list:
        suffix = os.path.splitext(file)[-1]
        if suffix in format:
            fq_file = dir + file
            out_list.append(fq_file)
    return out_list

def run_bowtie_and_parse(args):
    fq_file, miR_file, out_name = args
    try:
        fq_col = out_name + '.col.fa'
        print(f"Collapsing FASTQ file: {fq_file}")
        run = subprocess.check_output(f"python {script_loc}/FastQcollapse.py {fq_file} > {fq_col}", shell=True)

        print(f"Building Bowtie index for: {fq_col}")
        discard = open(os.devnull, 'w')
        run = subprocess.check_output(f'bowtie-build -q {fq_col} {out_name}', stderr=discard, shell=True)

        fq_mapped = out_name + '.bowtie_mapped.txt'
        print(f"Running Bowtie for: {fq_col}")
        run = subprocess.check_output(f'bowtie {out_name} --norc -l 7 -n 2 -a -q -f {miR_file} > {fq_mapped}', stderr=discard, shell=True)

        if os.path.getsize(fq_mapped) != 0:
            out = fq_mapped.replace('.txt', '.parsed.txt')
            print(f"Parsing Bowtie result for: {fq_mapped}")
            run = subprocess.check_output(f'python {script_loc}/parse_bowtie_result.py -i {fq_mapped} -fq {fq_file} -o {out}', shell=True)

            if p['un']:
                fq_unmapped = out_name + '.bowtie_unmapped.fastq'
                print(f"Extracting unmapped reads for: {fq_mapped}")
                run = subprocess.check_output(f'python {script_loc}/extract_bowtie_unmapped_reads.py -mapped {fq_mapped} -fq {fq_file} -out {fq_unmapped}', shell=True)
        else:
            print(f"No mapped reads found in: {fq_mapped}")
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e.cmd}")
        print(f"Return code: {e.returncode}")
        print(f"Output: {e.output.decode()}")

def get_miRNA_count(mapped_list, Qscore, dir):
    if not mapped_list:
        print("No mapped files found.")
        return {}
    
    count_file = dir + 'miRNA_count.Q' + str(Qscore) + '.txt'
    print(f"Counting miRNAs in: {mapped_list}")
    run = subprocess.check_output(f'python {script_loc}/count_miRNA.py -i {" ".join(mapped_list)} -q {Qscore} -o {count_file}', shell=True)
    return count_mapped_read(count_file)

def count_total_read(fq_file):
    return int(subprocess.check_output(f'wc -l {fq_file} | awk \'{{print $1/4}}\'', shell=True))

def count_mapped_read(count_file):
    df = pd.read_table(count_file)
    df_col = [c for c in df.columns if c.split()[-1] != '(PM+1MM+2MM)']
    df = df[df_col]
    df['mm_numb'] = df['pos:mut'].str.count(':')
    df = df.groupby(['mm_numb']).sum()
    df = df.rename(index={0: 'PM', 1: '1MM', 2: '2MM'})
    return df.to_dict()

parser = argparse.ArgumentParser(usage='./Mapping_miRNA_with_bowtie.py [options] -r miRNA_sequence.fa', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-r', metavar='miR_seq.fa', required=True, help='\t\t\t(REQUIRED)')
parser.add_argument('-i', metavar='input.fq', nargs='*', help='\tdefault: None \t(optional)')
parser.add_argument('-d', metavar='input_dir/', default='./', help='\tdefault: ./ \t(optional)')
parser.add_argument('-o', metavar='output_dir/', default='./', help='\tdefault: ./ \t(optional)')
parser.add_argument('-p', metavar='int', choices=range(10), default=4, type=int, help='\tdefault: 4 \t(no. of worker threads)')
parser.add_argument('-q', metavar='int', choices=range(42), default=28, type=int, help='\tdefault: 38 \t(mismatch\'s Q-score must be higher or equal (>=) to this value)')
parser.add_argument('-un', action='store_true', help='\tdefault: False \t(save unmapped reads to fastq file)')
p = vars(parser.parse_args())

script_loc = '/'.join(os.path.realpath(__file__).split('/')[:-1])

fq_file_list = get_input_file(p['i'], p['d'])
out_dir = p['o']
miR_file = p['r']
q_threshold = p['q']
threads = p['p']

if out_dir[-1] != '/': out_dir += '/'
out_file_list = [out_dir + os.path.basename(f).split('.')[0] for f in fq_file_list]

pool = Pool(processes=threads)
if __name__ == '__main__':
    processes = pool.map(run_bowtie_and_parse, [(f, miR_file, out_file_list[i]) for i, f in enumerate(fq_file_list)])

bowtie_file_list = [f+'.bowtie_mapped.parsed.txt' for f in out_file_list if os.path.getsize(f+'.bowtie_mapped.txt') != 0]

if not bowtie_file_list:
    print("No bowtie mapped files found.")
else:
    mapped_count = get_miRNA_count(bowtie_file_list, q_threshold, out_dir)
    summary_file = out_dir + 'Mapping_summary.Q' + str(q_threshold) + '.txt'
    header = ['File name', 'Total read', 'PM read', '%', '1MM read', '%', '2MM read', '%']
    with open(summary_file, 'w') as f:
        f.writelines('\t'.join(header) + '\n')
        for file in fq_file_list:
            total = count_total_read(file)
            file_name = os.path.basename(file).split('.')[0]
            print_out = [file_name, total]
            if file_name in mapped_count:
                for mm in ['PM', '1MM', '2MM']:
                    mapped = mapped_count[file_name].get(mm, 0)  # Use .get to handle missing keys
                    print_out.append(str(mapped))
                    print_out.append(str(round(mapped / float(total) * 100, 2)))
            else:
                print_out += ['0', '0', '0', '0', '0', '0']
            print_out = list(map(str, print_out))  # Convert all items to strings
            f.writelines('\t'.join(print_out) + '\n')


