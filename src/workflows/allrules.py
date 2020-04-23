#!python
# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#source activate py36
import json
# Directories:
bam_dir = config["analysis"]["analysis_dir"]  + config["analysis"]["result"]  + "bam/"
result_dir = config["analysis"]["analysis_dir"] + config["analysis"]["result"]
vcf_dir = config["analysis"]["analysis_dir"] + config["analysis"]["result"]  + "vcf/"
vep_dir =  config["analysis"]["analysis_dir"] + config["analysis"]["result"]  + "vep/"

######## RULE_ALL #######

rule all:
  input:
    expand(bam_dir + "{sample}.bam", sample=config["samples"]),
    expand(bam_dir + "{sample}.sorted.bam", sample=config["samples"]),
    expand(bam_dir + "{sample}.sorted.mrkdup.bam", sample = config["samples"]),
    expand(bam_dir + "{sample}.metrics.txt", sample = config["samples"]),
    expand(result_dir + "{sample}" + ".sorted.mrkdup.bam.picard.bedintervals", sample = config["samples"]),
    expand(result_dir + "{sample}.sorted.mrkdup.hsmetric", sample = config["samples"]),
    expand(result_dir + "{sample}.sorted.alignmetric", sample = config["samples"]),
    expand(result_dir + "{sample}.bedtools.coverage.bed", sample = config["samples"]),
    expand(vcf_dir + "{sample}" + config["vcf"]["mutect"]["default"], sample = config["samples"]),
    expand(vcf_dir + "filt_split_" + "{sample}" + config["vcf"]["mutect"]["default"], sample = config["samples"]),
    expand(vcf_dir + "filt_split_" + "{sample}" + config["vcf"]["mutect"]["default"], sample = config["samples"]),
    expand(vcf_dir + "{sample}" + '_haplo.tumor.vcf', sample = config["samples"]),
    expand(vep_dir + "vep_" + "{sample}" + config["vcf"]["mutect"]["default"], sample = config["samples"]),
    expand(bam_dir + "{sample}.insert_size_metrics.txt", sample = config["samples"]),
    expand(bam_dir + "{sample}.insert_size_histogram.pdf", sample = config["samples"]),

########## BWA MEM ############
rule bwa_mem:
  input:
    fa = config["path"]["genomefa"] + config["references"]["genomefa_GATK"],
    #fa = config["path"]["genomefa"] + config["references"]["genomefa"],
    read1 = config["path"]["fastq"] + "{sample}" + "1_001.fastq.gz",
    read2 = config["path"]["fastq"] + "{sample}" + "2_001.fastq.gz",
    #refidx = expand(config["path"]["genomefa"] + config["references"]["genomefa"] + ".{prefix}", prefix=["amb","ann","bwt","pac","sa"])
    refidx = expand(config["path"]["genomefa"] + config["references"]["genomefa_GATK"] + ".{prefix}", prefix=["amb","ann","bwt","pac","sa"]),
  output:
    bamout = bam_dir + "{sample}.bam",
  params:
    header_1 = "'@RG\\tID:" +  "{sample}" + "\\tSM:" + "{sample}" + "\\tPL:ILLUMINAi'",
  threads: 4
  shell:
    "bwa mem "
        "-t 4 "
        "-R  {params.header_1} "
        "-M "
        "-v 1 "
        "{input.fa} {input.read1} {input.read2} "
        "| samtools view -Sb - > {output.bamout}; "


######## SAMTOOLS_SORT #########
rule samtools_sort_index:
    input:
        bam_dir + "{sample}.bam",
    output:
        bam_dir + "{sample}.sorted.bam",
    threads: 4
    shell:
        "samtools sort -@ {threads} -o {output} {input}"


########## PICARD TOOLS MARK DUPLICATES ############

rule picard_mark_duplicate:
  input:
    bam_dir + "{sample}" + ".sorted.bam",
  output:
    mrkdup = bam_dir + "{sample}" + ".sorted.mrkdup.bam",
    metrics= bam_dir + "{sample}.metrics.txt",
  log:
    stats = bam_dir + "{sample}.sorted.mrkdup.txt",
  shell:
    "picard MarkDuplicates "
        "INPUT={input} "
        "OUTPUT={output.mrkdup} "
        "VALIDATION_STRINGENCY=SILENT "
        "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 "
        "REMOVE_DUPLICATES=FALSE "
        "METRICS_FILE='{output.metrics}'; "
    "samtools index {output.mrkdup}; "


########## PICARD COLLECTHSMETRICS ###########

rule picard_collectHSmetric:
  input:
    fadict = (config["path"]["genomefa"] + config["references"]["genomefa_GATK"]).replace(".fasta",".dict"),
    bed = config["path"]["panel"] + config["bed"]["variant_panel"],
    bam = bam_dir + "{sample}" + ".sorted.mrkdup.bam",
  output:
    hsmetric = result_dir + "{sample}.sorted.mrkdup.hsmetric",
    bed2interval = result_dir + "{sample}" + ".sorted.mrkdup.bam.picard.bedintervals",
  shell:
    "picard "
      "BedToIntervalList "
      "I={input.bed} "
      "O={output.bed2interval} "
      "SD={input.fadict}; "
    "picard "
      "CollectHsMetrics "
      "BI={output.bed2interval} "
      "TI={output.bed2interval} "
      "I={input.bam} "
      "O={output.hsmetric} "
      "COVERAGE_CAP=500 "
      "METRIC_ACCUMULATION_LEVEL=ALL_READS "
      "METRIC_ACCUMULATION_LEVEL=LIBRARY"

######## CALCULATE COVERAGE PER BED REGION ###################
rule bedtools_coverage:
  input:
    bed = config["path"]["panel"] + config["bed"]["variant_panel_nopad"],
    bam = expand(bam_dir + "{sample}" + ".sorted.mrkdup.bam", sample = config["samples"]),
  output:
    bedcov = result_dir + "{sample}.bedtools.coverage.bed",
  shell:
    "bedtools "
      "coverage "
      "-a {input.bed} "
      "-b {input.bam} "
      "-mean "
      "> {output.bedcov}"

######### PICARD COLLECTALIGNMENTSUMMARYMETRICS ##########

rule picard_CollectAlignmentSummaryMetrics:
  input:
    bam = bam_dir + "{sample}" + ".sorted.mrkdup.bam",
    fa = config["path"]["genomefa"] + config["references"]["genomefa_GATK"],
  output:
    metric = result_dir + "{sample}.sorted.alignmetric",
  params:
    adapter = config["QC"]["adapter"],
  shell:
    "picard "
      "CollectAlignmentSummaryMetrics "
      "R={input.fa} "
      "I={input.bam} "
      "O={output.metric} "
      "ADAPTER_SEQUENCE={params.adapter} "
      "METRIC_ACCUMULATION_LEVEL=ALL_READS "
      "METRIC_ACCUMULATION_LEVEL=LIBRARY;"

############ Picard CollectInsertSizeMetrics ################

rule picard_CollectInsertSizeMetrics:
  input:
    bam = bam_dir + "{sample}" + ".sorted.mrkdup.bam",
  output:
    metric = bam_dir + "{sample}.insert_size_metrics.txt",
    histogram = bam_dir + "{sample}.insert_size_histogram.pdf"
  params:
    gatk4 = config["analysis"]["gatk4"]
  shell:
    "conda activate py27; java -jar {params.gatk4} "
      "CollectInsertSizeMetrics "
      "--INPUT {input.bam} "
      "--OUTPUT {output.metric} "
      "--Histogram_FILE {output.histogram}; conda deactivate"

##############################################################
#####################  VARIANT CALLING #######################
##############################################################

#####################  Mutect 2  #############################

rule mutect2_somatic_tumor_only:
  input:
    fa = config["path"]["genomefa"] + config["references"]["genomefa_GATK"],
    dbsnpALL = config["path"]["genomefa"] + config["references"]["dbsnpALL"],
    COSMICcodMut = config["path"]["genomefa"] + config["references"]["cosmic"],
    bamTumor =  bam_dir + "{sample}" + ".sorted.mrkdup.bam",
    bed = config["path"]["panel"] + config["bed"]["variant_panel"],
    germline = config["path"]["genomefa"] + config["references"]["germline"]
  output:
    m2_vcf = vcf_dir + "{sample}" + config["vcf"]["mutect"]["default"],
  params:
    tumor_sample_name = list(config["samples"].keys())[0],
    gatk4 = config["analysis"]["gatk4"]
  threads: 6
  shell:
    "java -jar {params.gatk4} Mutect2 "
        "-R {input.fa} "
        "-I {input.bamTumor} "
        "-tumor {params.tumor_sample_name} "
        "-L {input.bed} "
        "-O {output.m2_vcf};"
#        "--germline-resource {input.germline} "  # Add this only if totally sure that some variants of interest will not be filtered out

rule splitmultiallelic:
  input:
    m2_vcf = vcf_dir + "{sample}" + config["vcf"]["mutect"]["default"]
  output:
    m2_vcf_gz = vcf_dir + "{sample}" + config["vcf"]["mutect"]["merged"],  # Do not take much into account the name "merged". Nothing was merged, but only a VCF compressed
    m2_vcf_splitalleles = vcf_dir + "split_" + "{sample}" + config["vcf"]["mutect"]["default"]
  shell:
    "bgzip -c {input.m2_vcf} > {output.m2_vcf_gz}; tabix -p vcf {output.m2_vcf_gz};" # Extra step in this rule: Be aware of the use of bcftools norm to split multiallelic locations
    "bcftools norm -m -both -o {output.m2_vcf_splitalleles} -Ov {output.m2_vcf_gz}" # Extra step in this rule: splitting multiallelic sites in the VCF


rule mutect2_filterMutectCalls:
  input:
   m2_vcf_splitalleles = vcf_dir + "split_" + "{sample}" + config["vcf"]["mutect"]["default"]
  output:
    filt_m2_vcf = vcf_dir + "filt_split_" + "{sample}" + config["vcf"]["mutect"]["default"]
  params:
    gatk4 = config["analysis"]["gatk4"]
  threads: 6
  shell:
    "java -jar {params.gatk4} FilterMutectCalls "
        "-V {input.m2_vcf_splitalleles} "
        "-O {output.filt_m2_vcf};"


##############################################################
################## HAPLOTYPECALLER ###########################
##############################################################

rule haplotypecaller:
  input:
    bamTumor = bam_dir + "{sample}" + ".sorted.mrkdup.bam",
    fa = config["path"]["genomefa"] + config["references"]["genomefa_GATK"],
    bed = config["path"]["panel"] + config["bed"]["variant_panel"],
  output:
    hp_vcf = vcf_dir + "{sample}" + '_haplo.tumor.vcf'
  params:
    gatk4 = config["analysis"]["gatk4"]
  threads: 6
  shell:
    "java -jar {params.gatk4} HaplotypeCaller "
        "-R {input.fa} "
        "-I {input.bamTumor} "
        "-O {output.hp_vcf} "
        "-L {input.bed} "
        "-ERC BP_RESOLUTION"


#############################################################
#####################  VCF ANNOTATION #######################
#############################################################

##############  VARIANT-EFFECT-PREDICTOR  ###################


rule vep_annotation_Mutect2:
  input:
    filt_m2_vcf = vcf_dir + "filt_split_" + "{sample}" + config["vcf"]["mutect"]["default"],
  output:
    vep_vcf = vep_dir + "vep_" + "{sample}" + config["vcf"]["mutect"]["default"],
  params:
    vep_exec = config["path"]["vep"],
    #vep_refseq = config["references"]["vep_cache_refseq"]
  threads: 4
  shell:
    "source deactivate; {params.vep_exec}vep "
        "-i {input.filt_m2_vcf} "
        "-o {output.vep_vcf} "
        "--fork 2 "
        "--vcf "
        "--poly p "
        "--sift p "
        "--variant_class "
        "-format vcf "
        "--offline "
        "-v "
        "--force_overwrite "
        "--assembly GRCh37 "
        "--af_gnomad; source activate py36;"
