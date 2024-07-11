# python extract_bowtie_unmapped_reads.py -fq bowtie_input.fastq -mapped bowtie_mapped.txt -out bowtie_unmapped.fastq


#### get mapped read ####
import argparse

def bowtie_mapped(file):
    mapped = set()
    with open(file, 'r') as f:
        l = f.readline()
        while l:
            mapped.add(l.split('\t')[2].split('#')[0])
            l = f.readline()
    return list(mapped)

def extract_unmapped_reads(file, mapped):
    unmapped = {}
    with open(file, 'r') as f:
        l = f.readline()
        while l:
            if l[0] == '@':
                data = l
                seq = f.readline().rstrip('\n\r\a')
                data += seq + '\n'
                data += f.readline()
                data += f.readline()

                if seq in unmapped:
                    unmapped[seq].append(data)
                else:
                    unmapped[seq] = [data]
            l = f.readline()

    for s in mapped:
        if s in unmapped:
            del unmapped[s]

    return unmapped

parser = argparse.ArgumentParser(usage='python extract_bowtie_unmapped_reads.py [options] -fq bowtie_input.fastq -mapped bowtie_mapped.txt -out bowtie_unmapped.fastq')
parser.add_argument('-fq', metavar='in.fq', required=True, help='fastq format file used for bowtie mapping.')
parser.add_argument('-mapped', metavar='mapped.txt', required=True, help='Default output file from bowtie.')
parser.add_argument('-out', metavar='unmapped.fq', required=True, help='Output file name.')
p = vars(parser.parse_args())

mapped_reads = bowtie_mapped(p['mapped'])
unmapped_reads = extract_unmapped_reads(p['fq'], mapped_reads)

with open(p['out'], 'w') as f:
    for data in unmapped_reads.values():
        f.write(''.join(data))
