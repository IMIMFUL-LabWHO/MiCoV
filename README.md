# MiCov

MiCov is a Bash-based bioinformatics pipeline for processing **Oxford Nanopore Technologies (ONT)** sequencing data from human coronavirus (**hCoV**) samples.

The pipeline performs read trimming and quality control, reference mapping, primer trimming, mapping statistics, consensus sequence generation, variant calling, and result consolidation.

MiCov was developed through collaboration between **IMI MF, University of Ljubljana** and the **Hiscox Lab, University of Liverpool**.

**Author:** Alen Suljič ([alen.suljic@mf.uni-lj.si](mailto:alen.suljic@mf.uni-lj.si))
**Inspired by:** Martin Bosilj, Hannah Goldswain
**Initial version:** 9 July 2024

---

## Prerequisites

Before running MiCov, install:

* **Singularity**
  Installation instructions: https://singularity-tutorial.github.io/01-installation/

The workflow requires:

* `micov.sh` — main MiCov analysis pipeline
* MiCov Singularity definition file — used to build the analysis container
* `reference/` — reference genomes, genome annotations, and primer coordinates
* gzip-compressed ONT FASTQ files

The `reference/` directory must contain matching:

```text
<reference>_<amplicon-length>.fasta
<reference>_<amplicon-length>.gff3
<reference>_<amplicon-length>.primer.bed
```

files used for reference mapping, variant annotation, and primer trimming.

For example:

```text
hcov_229e_1000.fasta
hcov_229e_1000.gff3
hcov_229e_1000.primer.bed
```

The amplicon length must be included in the reference basename supplied to `micov.sh`.

The Singularity container must provide the bioinformatics software used by the pipeline, including:

* `fastp`
* `minimap2`
* `samtools`
* `ivar`
* `seqtk`

Standard Unix utilities such as `awk`, `grep`, `sed`, `tr`, `sort`, `cut`, `find`, and `gunzip` are also used.

---

## Repository overview

MiCov consists primarily of a Bash pipeline that converts ONT FASTQ reads into processed alignments, consensus sequences, coverage statistics, and variant tables.

```text
ONT FASTQ files
      │
      ▼
   micov.sh
      │
      ▼
Read trimming + QC
      │
      ▼
ONT reference mapping
      │
      ▼
Primer trimming
      │
      ├──────────────► Mapping statistics
      │
      ├──────────────► Coverage data
      │
      ├──────────────► Consensus sequences
      │
      └──────────────► Variant calls
                              │
                              ▼
                       Consolidated results
```

### Scripts

| Script               | Purpose                                                                                                                              |
| -------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| `micov.sh`           | Main ONT sequencing pipeline: read trimming, mapping, primer trimming, mapping statistics, consensus generation, and variant calling |
| `rename_barcodes.sh` | Barcode-renaming helper referenced by the MiCov workflow for preparing ONT input files                                               |

---

# Usage

## 1. Build the Singularity container

Build the `micov.sif` container from the supplied MiCov Singularity definition file.

Replace `<micov_definition_file>.def` below with the actual definition filename included in the repository:

```bash
sudo singularity build micov.sif <micov_definition_file>.def
```

If you do not have `sudo` privileges, build the container using `--fakeroot`:

```bash
singularity build --fakeroot micov.sif <micov_definition_file>.def
```

After a successful build, the working directory should contain approximately:

```text
MiCov/
├── micov.sh
├── micov.sif
├── <micov_definition_file>.def
├── reference/
│   ├── <reference>_<amplicon-length>.fasta
│   ├── <reference>_<amplicon-length>.gff3
│   └── <reference>_<amplicon-length>.primer.bed
└── ...
```

---

## 2. Prepare the input FASTQ files

MiCov processes gzip-compressed **single-end ONT FASTQ files**.

The script expects files following the naming convention:

```text
<sample>.fastq.gz
```

For example:

```text
1.fastq.gz
2.fastq.gz
3.fastq.gz
4.fastq.gz
```

or:

```text
sample01.fastq.gz
sample02.fastq.gz
sample03.fastq.gz
```

The workflow assumes that ONT barcode files have already been renamed appropriately before running `micov.sh`.

The pipeline source refers to:

```text
rename_barcodes.sh
```

for this preprocessing step.

Sample identifiers are generated from the portion of each filename before the first `.`.

For example:

```text
sample01.fastq.gz
```

is interpreted as:

```text
sample01
```

> **Important:** Avoid additional periods in sample filenames because the current sample-detection logic uses the portion of the filename preceding the first `.` as the sample identifier.

---

## 3. Select the appropriate reference

MiCov expects the reference basename to include the corresponding **amplicon length**.

For example:

```text
hcov_229e_1000
```

corresponds to:

```text
reference/hcov_229e_1000.fasta
reference/hcov_229e_1000.gff3
reference/hcov_229e_1000.primer.bed
```

All three files must use the same basename.

The files are used as follows:

* `.fasta` — reference genome used for mapping, consensus generation, and variant calling
* `.gff3` — genome annotation used for variant annotation
* `.primer.bed` — primer coordinates used by `ivar trim`

---

## 4. Check pipeline parameters

Before running MiCov, review the quality and length thresholds near the beginning of `micov.sh`.

These parameters may need to be adjusted according to the characteristics and quality summary of the ONT sequencing run.

Default values are:

```bash
# thread count
thr=24

# minimum Phred quality score
qqp=12

# minimum read length
lr=60

# cut mean quality score
cmq=12

# minimum base quality score
mbq=12

# minimum mapping quality score
mmq=10
```

| Parameter | Default | Description                                                |
| --------- | ------: | ---------------------------------------------------------- |
| `thr`     |    `24` | Number of processing threads                               |
| `qqp`     |    `12` | Minimum Phred quality threshold used during read filtering |
| `lr`      |    `60` | Minimum read length retained after trimming                |
| `cmq`     |    `12` | Mean quality threshold used during trimming                |
| `mbq`     |    `12` | Minimum base-quality parameter                             |
| `mmq`     |    `10` | Minimum mapping-quality parameter                          |

These values can be modified directly in `micov.sh` before execution.

---

## 5. Run MiCov

Run the pipeline inside the Singularity container using Bash:

```bash
singularity exec micov.sif bash micov.sh <reference_name> <fastq_directory>
```

For example:

```bash
singularity exec micov.sif bash micov.sh hcov_229e_1000 /path/to/run/data/
```

The first positional argument specifies the reference basename:

```text
hcov_229e_1000
```

MiCov therefore expects:

```text
reference/hcov_229e_1000.fasta
reference/hcov_229e_1000.gff3
reference/hcov_229e_1000.primer.bed
```

The second positional argument specifies the directory containing the ONT `.fastq.gz` files.

For example:

```text
/path/to/run/data/
```

---

# Main script overview

The main `micov.sh` script automates the processing of ONT sequencing reads from FASTQ files to consensus sequences, variant calls, sequencing statistics, and per-position coverage data.

For each sample, MiCov performs the following main steps:

1. **Sample detection**
   Identifies samples from gzip-compressed FASTQ filenames in the supplied input directory.

2. **Read trimming and QC**
   Uses `fastp` to perform quality trimming, remove poly-X/poly-G sequence, filter short reads, and generate HTML and JSON QC reports.

3. **ONT reference mapping**
   Maps trimmed reads to the selected reference genome using `minimap2` with the ONT-specific `map-ont` preset.

4. **Primer trimming**
   Removes primer-derived sequence from mapped reads using `ivar trim` and the corresponding primer BED file.

5. **BAM processing**
   Uses `samtools` to remove unmapped reads, sort alignments, and index the final BAM file.

6. **Mapping and coverage statistics**
   Calculates mapping statistics, reference coverage, and per-position sequencing depth using `samtools`.

7. **Consensus generation**
   Generates a consensus sequence using `samtools mpileup` and `ivar consensus`. Ambiguous nucleotide calls and gap characters are converted to `N`.

8. **Variant calling and annotation**
   Detects variants using `ivar variants` and annotates them using the supplied reference FASTA and GFF3 files.

9. **Result consolidation**
   Combines per-sample consensus sequences, variants, coverage data, and mapping statistics into consolidated result files.

The workflow can be summarized as:

```text
ONT FASTQ
    │
    ▼
fastp
Read trimming + QC
    │
    ▼
minimap2
ONT reference mapping
    │
    ▼
iVar primer trimming
    │
    ▼
Final BAM
    │
    ├──────────────► Mapping statistics
    │
    ├──────────────► Per-position coverage
    │
    ├──────────────► Consensus sequence
    │
    └──────────────► Variant calling
                            │
                            ▼
                     Consolidated results
```

---

# Pipeline details

## Read trimming and quality control

Raw ONT reads are processed using `fastp`.

The pipeline performs:

* front-end trimming
* tail-end trimming
* sliding-window quality filtering
* minimum quality filtering
* minimum read-length filtering
* poly-X trimming
* poly-G trimming

Trimmed reads are written to:

```text
trimmed/<sample>_trim.fastq.gz
```

For example:

```text
trimmed/sample01_trim.fastq.gz
```

QC reports are generated as:

```text
qc/<sample>_fastp.html
qc/<sample>_fastp.json
```

---

## Reference mapping

Trimmed reads are aligned to the reference genome using:

```bash
minimap2 -x map-ont
```

The `map-ont` preset is used for Oxford Nanopore reads.

Secondary alignments are disabled using:

```bash
--secondary=no
```

Unmapped reads are removed before downstream processing.

---

## Primer trimming

Mapped reads are processed with:

```text
ivar trim
```

using:

```text
reference/<reference>.primer.bed
```

Primer trimming is performed before the final BAM file is generated.

---

## Final BAM files

The final mapped and primer-trimmed alignment for each sample is written to:

```text
mappings/<sample>/<sample>.bam
```

with the corresponding index:

```text
mappings/<sample>/<sample>.bam.bai
```

For example:

```text
mappings/sample01/sample01.bam
mappings/sample01/sample01.bam.bai
```

---

# Mapping and coverage statistics

MiCov calculates several alignment and sequencing statistics using `samtools`.

For each sample, the pipeline runs:

* `samtools flagstat`
* `samtools coverage`
* `samtools depth`

Per-sample statistics are stored in:

```text
stats/
```

including:

```text
stats/<sample>_stats.log
stats/<sample>.covdepth
```

---

## Mapping statistics report

Per-sample mapping information is consolidated into:

```text
results/mapstats.tsv
```

The table contains:

```text
sample
rname
startpos
endpos
numreads
covbases
coverage
meandepth
meanbaseq
meanmapq
primary_mapped
nreads_raw
```

The fields include:

* `sample` — sample identifier
* `rname` — reference sequence
* `startpos` — first reference position
* `endpos` — final reference position
* `numreads` — number of reads included in the coverage calculation
* `covbases` — number of covered reference bases
* `coverage` — percentage of the reference sequence covered
* `meandepth` — mean sequencing depth
* `meanbaseq` — mean base quality
* `meanmapq` — mean mapping quality
* `primary_mapped` — number of primary mapped reads
* `nreads_raw` — number of reads counted from the trimmed FASTQ file

> **Note:** Despite the column name `nreads_raw`, the current implementation calculates this value from `trimmed/<sample>_trim.fastq.gz`, so it represents the number of reads remaining after `fastp` processing.

---

# Consensus sequence generation

Consensus sequences are generated using:

```text
samtools mpileup
```

followed by:

```text
ivar consensus
```

The current pipeline uses a consensus allele-frequency threshold of:

```text
0.5
```

and a minimum sequencing depth of:

```text
10
```

Ambiguous IUPAC nucleotide codes are converted to:

```text
N
```

Gap characters are also converted to `N`, and the resulting sequence is converted to uppercase.

Individual consensus sequences are written to:

```text
consensus/<sample>.fasta
```

All sample consensus sequences are subsequently concatenated into:

```text
results/consensus_sequences.fasta
```

---

# Variant calling

Variants are detected using:

```text
samtools mpileup
```

followed by:

```text
ivar variants
```

The current variant-frequency threshold is:

```text
0.01
```

corresponding to a minimum alternate allele frequency of 1%.

The minimum sequencing depth supplied to iVar is:

```text
10
```

The reference genome:

```text
reference/<reference>.fasta
```

and genome annotation:

```text
reference/<reference>.gff3
```

are used during variant calling and annotation.

Individual sample variant tables are written to:

```text
variants/<sample>.tsv
```

---

# Result consolidation

After processing all samples, MiCov consolidates the individual outputs into analysis-wide result files.

Sample names are added to variant and coverage tables before the per-sample files are merged.

The main consolidated outputs are written to:

```text
results/
```

---

# Output files

The main result files generated by `micov.sh` are:

```text
results/
├── mapstats.tsv
├── consensus_sequences.fasta
├── sleek_variants.tsv
└── coverage.csv
```

---

## `mapstats.tsv`

Contains consolidated mapping and sequencing statistics for all samples.

---

## `consensus_sequences.fasta`

Multi-FASTA file containing the final consensus sequence generated for each sample.

---

## `sleek_variants.tsv`

Combined iVar variant table containing variants detected across all samples.

Individual sample variant tables from:

```text
variants/
```

are consolidated into this file.

---

## `coverage.csv`

Combined per-position sequencing-depth data generated from:

```text
samtools depth
```

The first column contains the sample identifier followed by the corresponding reference position and sequencing depth.

---

# Complete output structure

After a successful MiCov analysis, the working directory will contain approximately:

```text
working_directory/
├── micov.sh
├── micov.sif
│
├── reference/
│   ├── <reference>_<amplicon-length>.fasta
│   ├── <reference>_<amplicon-length>.gff3
│   └── <reference>_<amplicon-length>.primer.bed
│
├── trimmed/
│   ├── <sample>_trim.fastq.gz
│   └── ...
│
├── qc/
│   ├── <sample>_fastp.html
│   ├── <sample>_fastp.json
│   └── ...
│
├── mappings/
│   ├── <sample>/
│   │   ├── <sample>.bam
│   │   └── <sample>.bam.bai
│   └── ...
│
├── stats/
│   ├── <sample>_stats.log
│   ├── <sample>.covdepth
│   └── ...
│
├── consensus/
│   ├── <sample>.fasta
│   └── ...
│
├── variants/
│   ├── <sample>.tsv
│   └── ...
│
├── results/
│   ├── mapstats.tsv
│   ├── consensus_sequences.fasta
│   ├── sleek_variants.tsv
│   └── coverage.csv
│
└── logs/
    ├── experiment.log
    └── samples
```

---

# Logging

The main pipeline writes standard output and standard error to:

```text
experiment.log
```

using:

```bash
exec > >(tee -a "experiment.log") 2>&1
```

At the end of the analysis, the log file and generated sample list are moved to:

```text
logs/
├── experiment.log
└── samples
```

The log file can be used to inspect individual processing steps and troubleshoot failed or incomplete analyses.

---

# Important assumptions

The current implementation assumes that:

1. sequencing data consist of Oxford Nanopore reads
2. input reads are gzip-compressed FASTQ files
3. FASTQ files follow the `<sample>.fastq.gz` naming convention
4. barcode files have been appropriately renamed before running MiCov
5. sample identifiers do not contain additional periods
6. the appropriate reference basename includes the amplicon length
7. matching FASTA, GFF3, and primer BED files are available for the selected reference
8. the three reference files use the same basename
9. the `reference/` directory is located directly inside the working directory
10. the script is executed using Bash
11. required bioinformatics software is available through the MiCov Singularity container
12. pipeline quality thresholds are reviewed and adjusted when appropriate for the ONT sequencing run
13. sufficient CPU, memory, and disk space are available for mapping and intermediate files

When repeating an analysis, using a new working directory is recommended to avoid mixing files generated by different sequencing runs or parameter configurations.

---

# Citation

If MiCov is used in published work, please cite this repository together with relevant publications describing the sequencing protocol, amplicon schemes, and analysis workflow.

For the seasonal human coronavirus amplicon sequencing workflow:

Ošep, A., Goldswain, H., Suljič, A. *et al.* Development of type-specific amplicon schemes for whole-genome sequencing of seasonal human coronaviruses from clinical samples. **Scientific Reports** (2026).
https://doi.org/10.1038/s41598-026-63549-1

---

