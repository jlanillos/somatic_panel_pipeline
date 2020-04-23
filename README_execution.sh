#!/bin/bash
#########################################################################
# Create Config files (this approach might vary depending how your fastq files are organized)
#########################################################################
# Define some directories
data_dir=/mnt/64716603-5b56-4f9a-b195-c11560647a3a/Projects/MariaGarcia/data/MariaGarcia_DNA_Seq
analysis_dir=/home/jlanillos/Disco4tb/Projects/MariaGarcia
cd $data_dir
for i in Sample_*; do ls $i/*1_001.fastq.gz >> $analysis_dir/src/R1_fast.txt; ls $i/*2_001.fastq.gz >> $analysis_dir/src/R2_fast.txt; echo $i >> $analysis_dir/src/sample_names.txt; done
cd $analysis_dir
paste -d '\t' $analysis_dir/src/sample_names.txt $analysis_dir/src/R1_fast.txt $analysis_dir/src/R2_fast.txt > $analysis_dir/src/samples_dir.csv
python $analysis_dir/src/config_files/create_table_for_config.py --file $analysis_dir/src/samples_dir.csv --datadir $data_dir
# Before executing the next script, make sure you modify /src/config_files/template.json according to your analysis
python $analysis_dir/src/config_files/create_config_from_template.py --outdir $analysis_dir/src/config_files/ --file $analysis_dir/src/samples_dir.csv --template $analysis_dir/src/config_files/template.json


#########################################################################
# Run Snakemake over your config files (remove dryrun when ready for execution)
#########################################################################
# GATK: this pipeline makes use of Mutect2 and Haplotypecaller from GATK4. Download your preferred version make sure to modify the path in "template.json"
conda activate py36 # Activate conda environment (mine is py36) where all tools in the pipeline are installed (bwa, samtools, picard tools, bedtools, and SNAKEMAKE)
snakemake --snakefile $analysis_dir/src/workflows/allrules.py --configfile $analysis_dir/src/config_files/config_sampleID.json --dryrun



#########################################################################
# ANALYSIS POST-PROCESSING (generating matrix by merging all somatic VCF)
#########################################################################
# Merging all Mutect2 VCF files annotate them with VEP
analysis_dir=/home/jlanillos/Disco4tb/Projects/MariaGarcia/results
mkdir MERGED
for i in *results; do realpath $i/vcf/*_mutect.vcf.gz | grep -v 'filt' >> MERGED/mutect_vcf_list.txt; done
for i in *results; do bgzip $i/vcf/*haplo.tumor.vcf; tabix -p vcf $i/vcf/*haplo.tumor.vcf.gz realpath $i/vcf/*haplo.tumorq.vcf.gz>> MERGED/haplo_vcf_list.txt; done
grep ^'#' MERGED_haplotypecaller.vcf > NON_REF_MERGED_haplotypecaller.vcf; grep -v ^'#' MERGED_haplotypecaller.vcf | grep 'NON_REF' >> NON_REF_MERGED_haplotypecaller.vcf

cd MERGED/
# Copied matrix_scripts folders to use those script as ofr creting the matrix:
bcftools merge -o MERGED_all.vcf.gz -Oz -l mutect_vcf_list.txt
bcftools merge -o MERGED_haplotypecaller.vcf -Ov -l haplo_vcf_list.txt
# VEP annotation
conda deactivate
/home/jlanillos/Disco4tb/usr/ensembl-vep/vep -i MERGED_all.vcf.gz -o VEP_MERGED.vcf --fork 16 --vcf --poly p --sift p --variant_class -format vcf --offline -v --force_overwrite --assembly GRCh38 --everything
conda activate py36
bgzip VEP_MERGED.vcf; tabix -p vcf VEP_MERGED.vcf.gz
# Split multiallelic sites
bcftools norm -m -both -o VEP_MERGED_multiallelic.vcf.gz -Oz VEP_MERGED.vcf.gz
zcat VEP_MERGED_multiallelic.vcf.gz > VEP_MERGED_multiallelic.vcf # This output file will be used to create the matrix
# Parse VEP_MERGED VCF into a dataframe and perform some filtering
python MATRIXfromMERGEDVCF.py --file VEP_MERGED_multiallelic.vcf # It creates MUTECT_APPRIS_MATRIX.csv
python MATRIX_filtering.py --VCFfile VEP_MERGED_multiallelic.vcf --MATRIXfile MUTECT_APPRIS_MATRIX.csv # It creates stats_MUTECT_AF_IMPACT_APPRIS_MATRIX.csv
# Parse Haplotypecaller and create a matrix
python mergingVCFs_haplotypecaller.py --file NON_REF_MERGED_haplotypecaller.vcf # It creates HAPLOTYPECALLER_NON_REF_only.MATRIX.csv
# Adding Haplotypecaller information and getting the final matrix
python REANNOTATE_mutect_haplotypecaller.py --mutectfile stats_MUTECT_AF_IMPACT_APPRIS_MATRIX.csv --haplotypecallerfile HAPLOTYPECALLER_NON_REF_only.MATRIX.csv --VCFfile VEP_MERGED_multiallelic.vcf


# Move the Final matrix to a separate folder
mkdir MATRIX
mv MATRIX_MUTECT_HAPLO_AF_IMPACT_APPRIS.csv MATRIX
cd MATRIX
