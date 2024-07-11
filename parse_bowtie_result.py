# python parse_bowtie_result.py -i bowtie_result.txt -fq bowtie_input.fastq -o output.txt

""" * Alignment criteria
	1) PM > 1MM > 2MM
 	2) Counts "concatemer" as 1. (where, "concatemer" = >=2 identical miRNA sequence (including mismatch) present in a read
    - In case of concatemer, only the first-encountered miRNA is considered. (for Q-score)
"""
import argparse

def read_bowtie_result(file):
    mapped = {}
    with open(file, 'r') as f:
        for line in f:
            data = line.rstrip('\n').split('\t')
            if data[1] == '+':
                miR_name = data[0].split()[0]
                fq_seq = data[2]
                fq_count = 1
                if '#' in data[2]:
                    fq_seq = data[2].split('#')[0]
                    fq_count = int(data[2].split('#')[-1])
                mm_count = data[7].count('>')
                start = int(data[3])
                end = str(start + len(data[4]))
                data = [start, end, miR_name, data[7], fq_count]

                if fq_seq in mapped:
                    if mm_count in mapped[fq_seq]:
                        mapped[fq_seq][mm_count].append(data)
                    else:
                        mapped[fq_seq][mm_count] = [data]
                else:
                    mapped[fq_seq] = {mm_count: [data]}
    return mapped

def select_uniquely_mapped_read(bowtie_dict):
    uniq_mapped = {}
    multi_mapped = {}
    for fq_seq, score_dict in bowtie_dict.items():
        mm_count = min(score_dict.keys())

        if len(score_dict[mm_count]) == 1:
            uniq_mapped[fq_seq] = score_dict[mm_count][0]
        elif len(set([tuple(x[2:-1]) for x in score_dict[mm_count]])) == 1:
            start_index_list = [int(x[0]) for x in score_dict[mm_count]]
            start_index = start_index_list.index(min(start_index_list))
            uniq_mapped[fq_seq] = score_dict[mm_count][start_index]
        else:
            multi_mapped[fq_seq] = score_dict[mm_count]
    return uniq_mapped, multi_mapped

def get_qscore(file, seq_list=[]):
    res = {}
    with open(file, 'r') as r:
        l = r.readline()
        if l[0] == '@':
            while l:
                l = r.readline()  # Read the sequence line
                seq = l.rstrip('\n\r\a')
                r.readline()  # Skip the '+' line
                l = r.readline()  # Read the quality score line
                qscore = l.rstrip('\n\r\a')
                if seq in res:
                    res[seq].append(qscore)
                else:
                    res[seq] = [qscore]
                l = r.readline()  # Move to the next sequence identifier line
        else:
            print(f'Warning! "{file}" invalid fastq format. (does not start with @)\n')
            quit()
    if seq_list:
        return {seq: res[seq] for seq in seq_list}
    else:
        return res

q_ref = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"
q_dict = {x: i for i, x in enumerate(q_ref)}

def correct_bowtie_mut_Qscore(mut, start, qscore):
    if mut == '':
        return 'PM', qscore
    else:
        mismatch = []
        qscore = list(qscore)  # Convert to list for mutability

        for mm in mut.split(','):
            tmp = mm.split(':')
            mm_pos = int(tmp[0])
            qscore_pos = mm_pos + int(start)

            if qscore_pos < len(qscore):  # Ensure the position is within bounds
                mm_qscore = q_dict[qscore[qscore_pos]]
                qscore[qscore_pos] = ' '  # Replace the Q-score character with a space

                # Increment position by one
                new_pos = mm_pos + 1
                from_base = tmp[1][0]
                to_base = tmp[1][-1]

                mutation = f"{new_pos}:{from_base}>{to_base}:{mm_qscore}"
                mismatch.append(mutation)

        new_qscore = ''.join(qscore)
        return ','.join(mismatch), new_qscore

parser = argparse.ArgumentParser(usage='python parse_bowtie_result.py [options] -i bowtie_result.txt -fq bowtie_input.fastq -o output.txt')
parser.add_argument('-i', metavar='bowtie result', required=True, help='bowtie_result file in its default table format')
parser.add_argument('-o', metavar='file name', help='By default, ".parsed" will be added to the input file\'s name.')
parser.add_argument('-fq', metavar='fastq file', required=True, help='fastq file used for bowtie mapping')
p = vars(parser.parse_args())

if p['o'] is None:
    p['o'] = p['i'] + '.parsed'

mapped_reads, discard = select_uniquely_mapped_read(read_bowtie_result(p['i']))

q_score_dict = get_qscore(p['fq'], mapped_reads.keys())

heading = ['# miRNA_name', 'pos:mut:Q-score', 'Q-score', '5p_sequence', '5p_Q-score', '3p_sequence', '3p_Q-score']
with open(p['o'], 'w') as f:
    f.write('\t'.join(heading) + '\n')

    for fq_read, data in mapped_reads.items():
        start = int(data[0])
        end = int(data[1])
        for q_score in q_score_dict[fq_read]:
            mut, q_score = correct_bowtie_mut_Qscore(data[-2], start, q_score)
            f.write('\t'.join([data[2], mut, q_score[start:end], fq_read[:start], q_score[:start], fq_read[end:], q_score[end:]]) + '\n')

