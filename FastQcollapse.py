# python FastQcollapse.py "input.fastq" > "output.fasta"

import sys
import operator

def collapse_fastq(infile):
    with open(infile, 'r') as r:
        res = {}
        while True:
            l = r.readline()
            if not l:
                break
            if l[0] == '@':
                seq = r.readline().strip()
                if '@' not in seq:
                    r.readline()  # Skip the '+' line
                    r.readline()  # Skip the quality score line
                    if seq in res:
                        res[seq] += 1
                    else:
                        res[seq] = 1
                else:
                    continue

    for data in sorted(res.items(), key=operator.itemgetter(1), reverse=True):
        print(f'>{data[0]}#{data[1]}')
        print(data[0])

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python FastQcollapse.py <input.fastq>")
        sys.exit(1)
    input_file = sys.argv[1]
    collapse_fastq(input_file)

