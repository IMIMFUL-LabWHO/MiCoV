#!/bin/sh

# micov

# This script is a workhorse of hCoV project born out of collaboration between IMI MF, University of Ljubljana
# and Hiscox Lab, University of Liverpool.

# This procedure performs read trimming,  mapping on the provided reference genome, primer trimming,
# calculation of the mapping statistics, consensus calling and variant analysis.

# Author: Alen Suljič (alen.suljic@mf.uni-lj.si)
# Inspired by: Martin Bosilj, Hannah Goldswain
# First wielding: 09.07.2024

# Prerequisites:
###################################################################################################################
# 1. micov.sh
# 2. Reference directory (reference_'amplicon-length'.fasta, reference_'amplicon-length'.gff3, reference_'amplicon-length'.primer.bed)
# 3. Singulaity container micov.sif
###################################################################################################################

# Usage:
# 1. Copy the script to the working directory (this will be the output directory). 
# 2. Make sure that "data" directory with renamed barcodes (use rename_barcodes.sh).
# 3. Reference directory must contain "reference_'amplicon-length'.fasta", "reference_'amplicon-length'.gff3" and "reference_'amplicon-length'.primer.bed" files!
# 4. Take care to include amplicon length when inputing reference file name (e.g. hcov_229e_1000)!
# 4. The script is run inside dedicated singularity container "micov.sif".
# 5. Before running the script, adjust the parameters according to ONT output summary!

# Command example:
# bash micov.sh reference_'amplicon-length' /grid/userdata/asuljic/Runs/Run_20240306/data/

# Parameters:
###################################################################################################################
## thread count
thr=24
## minimum phread quality score
qqp=12
## length required
lr=60
## cut mean quality score
cmq=12
## minimum base quality score
mbq=12
## minimum mapping quality score
mmq=10
###################################################################################################################

# initiate log files
exec > >(tee -a "experiment.log") 2>&1

# create sample list
ls "${2}" | cut -d "." -f 1 | sort -n | uniq > samples

## Fastp

echo "Trimming reads"

for i in $(cat samples); do

	echo $i;

	mkdir -p trimmed/ qc/

	fastp \
	--in1 "${2}"/${i}.fastq.gz \
	--out1 trimmed/${i}_trim.fastq.gz \
	--cut_front \
	--cut_tail \
	--cut_window_size 4 \
	--cut_mean_quality $cmq \
	--qualified_quality_phred $qqp \
	--unqualified_percent_limit 40 \
	--length_required $lr \
	--trim_poly_x --trim_poly_g \
	--html qc/${i}_fastp.html \
	--json qc/${i}_fastp.json \
	--thread $thr
	
done

# Mapping reads to reference with Minimap2
## Primer trimming with iVar
echo "Mapping reads"

for i in $(cat samples); do

	echo $i;

	mkdir -p mappings/${i}

	minimap2 \
	-a \
	-t $thr \
	-x map-ont \
	--secondary=no \
	reference/"${1}".fasta trimmed/${i}_trim.fastq.gz | \
		samtools view -bS -F 4 -@ $thr | \
		samtools sort -@ $thr | \
		ivar trim -b reference/"${1}".primer.bed -m 0 -q 0 -e | \
		samtools sort -@ $thr > mappings/${i}/${i}.bam

	samtools index mappings/${i}/${i}.bam > mappings/${i}/${i}.bam.bai

done

# Mapping statistics
echo "Calculating mapping statistics"

for i in $(cat samples); do

	echo $i;
	mkdir -p stats

	echo -e "${i}" > stats/${i}_mapstats_name.log
	samtools flagstat mappings/${i}/${i}.bam > stats/${i}_allstats.log
	samtools coverage -q $mbq -Q $mmq mappings/${i}/${i}.bam > stats/${i}_coverage.log
	samtools depth -aH -q $mbq -Q $mmq mappings/${i}/${i}.bam -o stats/${i}.covdepth
	
	cat stats/${i}_mapstats_name.log stats/${i}_allstats.log stats/${i}_coverage.log > stats/${i}_stats.log
	rm  stats/${i}_allstats.log stats/${i}_mapstats_name.log stats/${i}_coverage.log

done

## Generate mapping stats report

echo "Generating mapping statistics report"

mkdir -p results

echo -e "sample\trname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\tprimary_mapped\tnreads_raw" > results/mapstats.tsv

for i in $(cat samples); do

  sample="${i}"

  first_line=$(head -n 1 stats/${sample}_stats.log)
  
  nineteenth_line=$(sed -n '19p' stats/${sample}_stats.log)

  primary_mapped=$(grep 'primary mapped' stats/${sample}_stats.log | awk '{print $1}')

  nreads_raw=$(gunzip -c trimmed/${sample}_trim.fastq.gz | awk '{s++}END{print s/4}')
  
  echo -e "$first_line\t$nineteenth_line\t$primary_mapped\t$nreads_raw" >> results/mapstats.tsv
  
done


# Consensus generation with iVar
echo "Generating consensus sequences"

for i in $(cat samples); do

	echo $i;
	mkdir -p consensus/

	samtools mpileup -aa -A -d 0 -E -f reference/"${1}".fasta -q $mmq -Q $mbq mappings/${i}/${i}.bam | \
	ivar consensus -p consensus/${i} -t 0.5 -q $mmq -m 10

	seqtk seq consensus/${i}.fa | \
		tr "?RYSWKMBDHVN.ryswkmbdhvn" "N" | \
		tr "-" "N" | \
		tr [:lower:] [:upper:] | grep -iv ">" - > consensus/tmp.seq

	echo ">"$i > consensus/${i}.fa
	cat consensus/${i}.fa consensus/tmp.seq  > consensus/${i}.fasta
	find ./consensus -type f ! -name '*.fasta' -delete

done

# Variant calling with iVar
echo "Variant calling"

for i in $(cat samples); do

	echo $i;

	mkdir -p variants

	samtools mpileup -aa -A -d 0 -E -f reference/"${1}".fasta -q $mmq -Q $mbq mappings/${i}/${i}.bam | \
	ivar variants -p variants/${i} -t 0.01 -q $mmq -m 10 -r reference/"${1}".fasta -g reference/"${1}".gff3

done

# Results generation and beautification
## Create results directory and concatenate consensus sequences
echo "Results consolidation and beautification"

cat consensus/*.fasta > results/consensus_sequences.fasta

## Add sample names as a column to the tsv files 

for file in variants/*.tsv; do 

    filename_noext=$(basename -- "$file")
    filename_noext="${filename_noext%.*}"

    awk -v filename="$filename_noext" 'BEGIN{OFS="\t"} {$1=filename; print}' "$file" > "${file}.tmp"
    mv "${file}.tmp" "$file"

done

## Merge variants to a single tsv file

awk 'FNR==1 && NR!=1{next;}{print}' variants/*.tsv > results/sleek_variants.tsv


echo "Beautifying coverage data"

for file in stats/*.covdepth; do 

    filename_noext=$(basename -- "$file")
    filename_noext="${filename_noext%.*}"

    awk -v filename="$filename_noext" 'BEGIN{OFS="\t"} {$1=filename; print}' "$file" > "${file}.tmp"
    mv "${file}.tmp" "$file"

done

# Merge coverage files to single dataframe

awk 'FNR==1 && NR!=1{next;}{print}' stats/*.covdepth > results/coverage.tsv

tr '\t' ',' < results/coverage.tsv > results/coverage.csv

rm results/coverage.tsv

# Cleanup
mkdir -p logs
mv experiment.log logs/
mv samples logs/


echo "Analysis complete!"

