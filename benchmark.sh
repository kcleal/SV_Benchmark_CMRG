#!/bin/bash

# Optionally run this in a Docker container. Set up first using:
# mkdir benchmark && cd benchmark
# docker run -it --memory="32g" --mount src="${PWD}",target=/results,type=bind condaforge/mambaforge
# mamba update conda -y && cd results

mamba create -c bioconda -c conda-forge -n bench python=3.9 awscli sniffles=2.0.7 cuteSV=2.0.3 truvari=4.0.0 delly=1.1.6 -y
conda activate bench
pip install dysgu==1.5.0


# Reference genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai
ref=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

# Get data, see ONT website https://labs.epi2me.io/giab-2023.05/
aws s3 sync --no-sign-request --include='PAO89685.pass.cram*' --exclude="*fail*" s3://ont-open-data/giab_2023.05/analysis/hg002/sup/ .

# PacBio data 
wget https://downloads.pacbcloud.com/public/revio/2022Q4/HG002-rep1/analysis/HG002.m84011_220902_175841_s1.GRCh38.bam
wget https://downloads.pacbcloud.com/public/revio/2022Q4/HG002-rep1/analysis/HG002.m84011_220902_175841_s1.GRCh38.bam.bai

# Truth set
wget -r -np -nH --cut-dirs=7 ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_HG002_medical_genes_SV_benchmark_v0.01/HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.*


# Run callers ONT
sniffles --input PAO89685.pass.cram --vcf HG002.PAO89685.sniffles.vcf

dysgu run --mode nanopore --procs 4 -x --clean $ref wd PAO89685.pass.cram > HG002.PAO89685.dysgu.vcf

mkdir wd_cuteSV
cuteSV -t 4 -s 3 --genotype PAO89685.pass.cram $ref HG002.PAO89685.cuteSV.vcf wd_cuteSV

delly lr -g $ref PAO89685.pass.cram > HG002.PAO89685.delly.vcf


# Run callers Revio
sniffles --input HG002.m84011_220902_175841_s1.GRCh38.bam --vcf HG002.pacbio.sniffles.vcf

dysgu call --mode pacbio --procs 4 -x --clean $ref wd HG002.m84011_220902_175841_s1.GRCh38.bam > HG002.pacbio.dysgu.vcf

mkdir wd_cuteSV
cuteSV -t 4 -s 3 --genotype HG002.m84011_220902_175841_s1.GRCh38.bam $ref HG002.pacbio.cuteSV.vcf wd_cuteSV

delly lr -g $ref  HG002.m84011_220902_175841_s1.GRCh38.bam > HG002.pacbio.delly.vcf


# Benchmark
callers=( "sniffles" "cuteSV" "delly" "dysgu" )

for name in "${callers[@]}"
do
  bgzip HG002.PAO89685.${name}.vcf
  tabix HG002.PAO89685.${name}.vcf.gz
  truvari bench -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
                -b HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.vcf.gz \
                --includebed HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.bed \
                -c HG002.PAO89685.${name}.vcf.gz \
                --passonly -r 1000 \
                -o truvari_${name}
done


for name in "${callers[@]}"
do
  bgzip HG002.pacbio.${name}.vcf
  tabix HG002.pacbio.${name}.vcf.gz
  truvari bench -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
                -b HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.vcf.gz \
                --includebed HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.bed \
                -c HG002.pacbio.${name}.vcf.gz \
                --passonly -r 1000 \
                -o truvari_${name}_pacbio
done

# Optionally run:
pip install matplotlib pandas seaborn
python3 plot_benchmark.py
