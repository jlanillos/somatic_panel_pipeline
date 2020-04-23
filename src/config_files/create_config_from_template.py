# This Python script reads a template json file and fills the necessary fields with the sample information from the file create by create_table_for_config.py
import json
import pandas as pd
import sys
import re
import argparse

parser=argparse.ArgumentParser(description='create_config_from_table.py takes "samples_dir.csv info and creates confi files per smple using template.json as template file')
parser.add_argument('--file', required=True)
parser.add_argument('--template', required=True)
parser.add_argument('--outdir', required=True)

args=parser.parse_args()
file=args.file
template_json = args.template
outdir = args.outdir

# Open the table with sample-specific info:
df = pd.read_csv(file, sep = '\t')
# Open the template.json and fill the fields needed using the info from the table created
with open(template_json) as f:
    data_aux = json.load(f)
    f.close()

for index, row in df.iterrows():
    data = data_aux
    # Fill some fields in the analysis section with sample-specific information and where the sample-specfic results will be stored
    data['analysis']['sample_id'] = row['sample']
    data['analysis']['result'] = row['sample'] + '_results/'
    # Path to sample-specific fastq files
    data['path']['fastq'] = row['path2fastq']
    # Path to save sample-specific results
    data['samples'] = {row['file_prefix']: {'file_prefix': row['file_prefix'], 'type': 'tumor', 'readpair_suffix': ['1', '2']}}
    config_file = outdir + 'config_' + row['sample'] + '.json'
    with open(config_file, 'w') as outfile:
        json.dump(data, outfile, indent = 4, sort_keys=True)
    outfile.close()

