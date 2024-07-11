# Dependency
1) Anaconda3/2024.02-1
# create environment and download other packages:
1) python3/3.9.6 (with package: pandas)
2) bowtie version 1.3.1
4) samtools 1.20

# Note:
1) All 5 scripts must be in same directory!
   (scripts: Mapping_miRNA_with_bowtie.py, FastQcollapse.py, parse_bowtie_result.py, count_miRNA.py & extract_bowtie_unmapped_reads.py)
2) "mature.fa" must be converted from RNA to DNA sequence!
   (if it's downloaded from miRBase, it is in RNA sequence)
3) Adjust quality threshold (default = 33) in Mapping_miRNA_with_bowtie.py to represent data proerly

# Steps:
1. "bash pre-processing_QIAseq_mirna.sh"
2. "python unzip_to_process.py" (samples.txt file must be provided with paths to _deduped.fq.gz files)
3. "python Converting_RNAtoDNA.py" (change paths in script, filtering option for species (ex. mmu-) and name of mature sequencing file)
4. "python Mapping_miRNA_with_bowtie.py -i /path/to/files/*.fastq -r /path/to/reference/mature_miRNA_sequence_dna.fa -o path/to/output_dir"

# FINAL OUTPUT: 
1) Mapping_summary.txt	(format: "file_name" | "total_read" | "PM_read" | "PM%" | "1MM_read" | "1MM%" | "2MM_read" | "2MM%")
2) miRNA_count.txt	(format: "miR_name" | "pos:mut" | "read_count" | "total_miR_count (PM+1MM+2MM)")
+
3) bowtie_unmapped.fastq (optional, if option "-un" is used)
