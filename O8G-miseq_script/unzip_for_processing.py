import os
import gzip

def unzip_and_reformat_fastq(gzip_path, output_dir):
    """Unzips a .fastq.gz file, removes the UMI part from headers, and saves to the specified output directory."""
    filename = os.path.basename(gzip_path).replace('.gz', '')
    output_path = os.path.join(output_dir, filename)
    with gzip.open(gzip_path, 'rt') as f_in, open(output_path, 'wt') as f_out:
        while True:
            header = f_in.readline().strip()
            if not header:
                break
            sequence = f_in.readline().strip()
            plus = f_in.readline().strip()
            quality = f_in.readline().strip()
            
            # Remove the UMI part from the header
            header_parts = header.split()
            header_main = header_parts[0].split('_')[0]
            header = f"{header_main} {' '.join(header_parts[1:])}"
            
            f_out.write(f"{header}\n{sequence}\n{plus}\n{quality}\n")
    return output_path

def process_files(sample_file, output_dir):
    """Processes each file listed in the sample file."""
    with open(sample_file, 'r') as f:
        fastq_files = f.read().splitlines()
    
    for fastq_gz in fastq_files:
        unzipped_fastq = unzip_and_reformat_fastq(fastq_gz, output_dir)
        print(f"Unzipped and reformatted {fastq_gz} to {unzipped_fastq}")

if __name__ == "__main__":
    sample_file = "samples.txt"
    output_dir = "/path/to/output_dir"
    os.makedirs(output_dir, exist_ok=True)
    process_files(sample_file, output_dir)


