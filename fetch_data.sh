#!/bin/bash

# Reference genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai

# ONT data
aws s3 sync --no-sign-request --exclude='*' --include='PAO89685.pass.cram' s3://ont-open-data/giab_2023.05/analysis/hg002/sup/ .
aws s3 sync --no-sign-request --exclude='*' --include='PAO89685.pass.cram.crai' s3://ont-open-data/giab_2023.05/analysis/hg002/sup/ .

# PacBio data
wget https://downloads.pacbcloud.com/public/revio/2022Q4/HG002-rep1/analysis/HG002.m84011_220902_175841_s1.GRCh38.bam
wget https://downloads.pacbcloud.com/public/revio/2022Q4/HG002-rep1/analysis/HG002.m84011_220902_175841_s1.GRCh38.bam.bai


# PacBio newer:
wget https://downloads.pacbcloud.com/public/2024Q4/Vega/HG002/Human-WGS-variant-pipeline/HiFi-human-WGS-WDL_v2.0.0-rc4/m21009_241011_231051.GRCh38.haplotagged.bam
wget https://downloads.pacbcloud.com/public/2024Q4/Vega/HG002/Human-WGS-variant-pipeline/HiFi-human-WGS-WDL_v2.0.0-rc4/m21009_241011_231051.GRCh38.haplotagged.bam.bai


# CMRG truth set
wget -r -np -nH --cut-dirs=7 ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_HG002_medical_genes_SV_benchmark_v0.01/HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.*

# GIAB truth set
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/GRCh38_HG2-T2TQ100-V1.1.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/GRCh38_HG2-T2TQ100-V1.1.vcf.gz.tbi
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed




# Make some test data
samtools view -bh --reference GCA_000001405.15_GRCh38_no_alt_analysis_set.fna m21009_241011_231051.GRCh38.haplotagged.bam chr1:1-100000 > chr1.pacbio_test.bam; samtools index chr1.pacbio_test.bam

samtools view -bh --reference GCA_000001405.15_GRCh38_no_alt_analysis_set.fna PAO89685.pass.cram  chr1:1-100000 > chr1.ont_test.bam; samtools index chr1.ont_test.bam
