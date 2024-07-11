# python count_miRNA.py -i input1.txt input2.txt ... -o output.txt

import argparse
import os
import pandas as pd

def read_parsed_bowtie(file, qscore=0):
    bowtie = {}
    with open(file, 'r') as f:
        l = f.readline()
        if l.startswith('#'):
            l = f.readline()
        while l:
            data = l.rstrip('\n').split('\t')
            if MM_Qscore_pass(data[1], qscore):
                mir_name = data[0] + '#' + remove_Qscore(data[1])
                if mir_name in bowtie:
                    bowtie[mir_name][0] += 1
                else:
                    bowtie[mir_name] = [1]
            l = f.readline()
    return bowtie

def MM_Qscore_pass(mut, qscore=0):
    if mut == 'PM':
        return True
    else:
        for mm in mut.split(','):
            if not int(mm.split(':')[-1]) >= qscore:
                return False
        return True

def remove_Qscore(mut):
    if mut == 'PM':
        return '0'
    else:
        return ','.join([':'.join(x.split(':')[:-1]) for x in mut.split(',')])

parser = argparse.ArgumentParser(usage='python count_miRNA.py [ -i input1.txt input2.txt ... / -d input_dir ] -o output.txt')
parser.add_argument('-i', metavar='Input files', nargs='*', help='Parsed bowtie result files.')
parser.add_argument('-d', metavar='Input directory', help='Directory with parsed bowtie result files. (suffix= ".parsed.txt")')
parser.add_argument('-o', metavar='Output file', required=True, help='Name for output file.')
parser.add_argument('-q', metavar='Q-score', choices=range(42), type=int, default=38, help='[optional] Q-score filter threshold for mismatch base. (higher or equal, Q>=38)')
p = vars(parser.parse_args())

in_files = p['i']
if p['d']:
    in_files = []
    if p['d'][-1] != '/':
        p['d'] += '/'
    for file in os.listdir(p['d']):
        if '.parsed.txt' in file:
            in_files.append(p['d'] + file)

if len(in_files) != len(set(in_files)):
    print("Warning: duplicate files given. Please check the file names.\n")
    quit()

final_df = pd.DataFrame({})
for file in in_files:
    bowtie_dict = read_parsed_bowtie(file, p['q'])
    miR_list = sorted(bowtie_dict.keys())
    bowtie_df = pd.DataFrame({file: [bowtie_dict[x][0] for x in miR_list]}, index=miR_list)
    final_df = pd.concat([final_df, bowtie_df], axis=1).fillna(0)

final_df['miRNA name#pos:mut'] = final_df.index.astype(str)
final_df[['miRNA name', 'pos:mut']] = final_df.pop('miRNA name#pos:mut').str.split("#", expand=True)
tmp_df = final_df['pos:mut'].str.split(",", expand=True).fillna(0)

# Debugging print statements
print("tmp_df.head():", tmp_df.head())
print("tmp_df.columns:", tmp_df.columns)
print("tmp_df.shape:", tmp_df.shape)

# Ensure the correct number of columns
if tmp_df.shape[1] == 1:
    tmp_df.columns = ['1']
elif tmp_df.shape[1] == 2:
    tmp_df.columns = ['1', '2']
else:
    raise ValueError(f"Unexpected number of columns in 'pos:mut' split: {tmp_df.shape[1]}")

final_df = final_df.join(tmp_df)
for i in list(tmp_df.columns):
    tmp_df = final_df.pop(i).str.split(":", expand=True).fillna(0)
    # Debugging print statement
    print(f"tmp_df.columns before renaming ({i}):", tmp_df.columns)
    print(f"tmp_df.shape before renaming ({i}):", tmp_df.shape)
    if tmp_df.shape[1] == 1:
        tmp_df.columns = ['pos' + i]
        tmp_df['mut' + i] = 0
    elif tmp_df.shape[1] == 2:
        tmp_df.columns = ['pos' + i, 'mut' + i]
    else:
        raise ValueError(f"Unexpected number of columns in 'pos:mut' split for {i}: {tmp_df.shape[1]}")
    tmp_df['pos' + i] = tmp_df['pos' + i].astype(int)
    final_df = final_df.join(tmp_df)

# Create 'pos2' and 'mut2' columns if they do not exist
if 'pos2' not in final_df.columns:
    final_df['pos2'] = 0
if 'mut2' not in final_df.columns:
    final_df['mut2'] = 0

final_df['mm_numb'] = final_df['pos:mut'].str.count(':')
final_df = final_df.sort_values(['miRNA name', 'mm_numb', 'pos1', 'mut1', 'pos2', 'mut2'])
final_df = final_df[list(final_df.columns)[-7:-5] + list(final_df.columns)[:-7]]
final_df['pos:mut'] = final_df['pos:mut'].replace('0', 'PM')
final_df.columns = [os.path.basename(x).split('.')[0] if i > 1 else x for i, x in enumerate(final_df.columns)]

tmp_df = final_df.groupby(['miRNA name']).sum()
tmp_df.columns = [x + ' (PM+1MM+2MM)' for x in tmp_df.columns]
final_df = pd.merge(final_df, tmp_df, left_on='miRNA name', right_index=True)

final_df.to_csv(p['o'], sep='\t', index=False)


