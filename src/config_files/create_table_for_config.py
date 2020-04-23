# In bash do: (this creates a file with three columns: sample_specific_folder, read1fastq, read2fastq (NOTE: there is no header in this file)
# The following line prints in separated files the information about each sample. Taking advantage of the pattern in the folder names (starting as Sample_..), do:
# for i in Sample_*; do ls $i/*_R1_001.fastq.gz >> /home/jlanillos/Disco4tb/Projects/NGS_panel_custom_CNIC_181127/src/R1_fastq.txt; ls $i/*_R2_001.fastq.gz >> /home/jlanillos/Disco4tb/Projects/[ANALYSIS_DIR]/src/R2_fastq.txt; echo $i >> /home/jlanillos/Disco4tb/Projects/[ANALYSIS_DIR]/src/sample_names.txt; done # Go to directory where the final table is stored:
#cd /home/jlanillos/Disco4tb/Projects/[ANALYSIS_DIR]/src/
#paste sample_names.txt R1_fastq.txt R2_fastq.txt > samples_dir.csv


import pandas as pd
import sys
import re
import argparse

parser=argparse.ArgumentParser(description='create_table_for_config.py takes "samples_dir.csv" file and makes some modifications to be used as input for create_config_table.py script')
parser.add_argument('--file', required=True)
parser.add_argument('--datadir', required=True)


args=parser.parse_args()
file=args.file
data_dir = args.datadir

cols_names = ['sample_folder','read1','read2']
df = pd.read_csv(file, sep = '\t', header = None)
df.columns = cols_names

#This column will be used to fill the config["path"]["fastq"] field in the json file
df['path2fastq'] = data_dir + '/' + df.sample_folder + '/'
# Extract the samples-file-prefix name to be used in the config too. Example:
df['file_prefix'] = df['read1'].str.split('/').str[1].str.rstrip('1_001.fastq.gz')
# Add a column only with the sample name:
df['sample'] = df['file_prefix'].str.split('__').str[0]
# Write to output .csv (I overwrite the same file "samples_dir.csv")
df.to_csv(file, sep='\t',index=None)
