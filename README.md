# MiCov

MiCov is a Bash-based bioinformatics pipeline for processing **Oxford Nanopore Technologies (ONT)** sequencing data from human coronavirus (**hCoV**) samples.

The pipeline performs read trimming and quality control, reference mapping, primer trimming, mapping and coverage statistics, consensus sequence generation, variant calling, and result consolidation. Additional R scripts are included for downstream coverage analysis, visualization, and variant-data preparation.

MiCov was developed through collaboration between **IMI MF, University of Ljubljana** and the **Hiscox Lab, University of Liverpool**.

**Author:** Alen Suljič (alen.suljic@mf.uni-lj.si)

---

## Prerequisites

Before running MiCov, install:

* **Singularity**
  Installation instructions: https://singularity-tutorial.github.io/01-installation/

The workflow requires:

* `micov.sh` — main MiCov analysis pipeline
* Singularity definition file used to build the MiCov container
* `reference/` — reference genomes, genome annotations, and primer coordinates
* gzip-compressed ONT FASTQ files
* R scripts for downstream coverage and variant analysis

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

The Singularity container provides the bioinformatics software used by the main pipeline, including:

* `fastp`
* `minimap2`
* `samtools`
* `ivar`
* `seqtk`

The downstream R scripts require:

* R
* `tidyverse`

---

## Repository overview

MiCov consists of a primary Bash pipeline for ONT sequencing-data processing and R scripts for downstream coverage and variant analysis.

```text
ONT FASTQ files
      │
      ▼
   micov.sh
      │
      ├────────► QC reports
      │
      ├────────► Final BAM files
      │
      ├────────► Consensus sequences
      │
      ├────────► Mapping statistics
      │
      ├────────► coverage.csv
      │                │
      │                ▼
      │      coverage_visualisation.R
      │                │
      │                ▼
      │      Genome coverage figures
      │
      └────────► sleek_variants.tsv
                       │
                       ▼
             variants_data_prep.R
                       │
                       ▼
             variants_enhanced.csv
```

### Scripts

| Script                     | Purpose                                                                                                                              |
| -------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| `micov.sh`                 | Main ONT sequencing pipeline: read trimming, mapping, primer trimming, mapping statistics, consensus generation, and variant calling |
| `coverage_visualisation.R` | Calculates coverage metrics and generates genome-wide coverage visualizations for hCoV-229E, NL63, OC43, and HKU1                    |
| `variants_data_prep.R`     | Cleans, classifies, and enriches the consolidated iVar variant table for downstream analysis                                         |

---

## Usage

### 1. Build the Singularity container

Build the `micov.sif` container from the supplied Singularity definition file.

Replace `<definition_file>.def` with the name of the MiCov definition file provided in the repository:

```bash
sudo singularity build micov.sif <definition_file>.def
```

If you do not have `sudo` privileges, build the container using `--fakeroot`:

```bash
singularity build --fakeroot micov.sif <definition_file>.def
```

After a successful build, the repository should contain approximately:

```text
MiCov/
├── micov.sh
├── micov.sif
├── <definition_file>.def
├── rename_barcodes.sh
├── coverage_visualisation.R
├── variants_data_prep.R
└── reference/
    ├── <reference>_<amplicon-length>.fasta
    ├── <reference>_<amplicon-length>.gff3
    └── <reference>_<amplicon-length>.primer.bed
```

---

### 2. Prepare the input FASTQ files

MiCov processes gzip-compressed **single-end ONT FASTQ files**.

The main script expects input files following the naming convention:

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

ONT barcode files should be renamed appropriately before running MiCov. The workflow references:

Sample identifiers are inferred from the portion of each filename before the first `.`.

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

### 3. Select the appropriate reference

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
* `.primer.bed` — primer coordinates used for primer trimming with `ivar trim`

---

### 4. Review pipeline parameters

Before running MiCov, review the quality and length thresholds near the beginning of `micov.sh`.

These parameters can be adjusted according to the characteristics and quality of the ONT sequencing run.

The default values are:

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
| `cmq`     |    `12` | Mean quality threshold used during read trimming           |
| `mbq`     |    `12` | Minimum base-quality threshold                             |
| `mmq`     |    `10` | Minimum mapping-quality threshold                          |

These values can be modified directly in `micov.sh` before execution.

---

### 5. Run MiCov

Run the pipeline inside the Singularity container:

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

The second positional argument specifies the directory containing the ONT `.fastq.gz` files:

```text
/path/to/run/data/
```

The current working directory is used as the analysis output directory.

---

## Main script overview

The main `micov.sh` script automates processing of ONT sequencing reads from raw FASTQ files to consensus sequences, variant calls, sequencing statistics, and per-position coverage data.

For each sample, the pipeline performs the following steps:

1. **Sample detection**
   Identifies samples from gzip-compressed FASTQ filenames in the supplied input directory.

2. **Read trimming and QC**
   Uses `fastp` to perform quality trimming, minimum-length filtering, and poly-X/poly-G trimming. HTML and JSON QC reports are generated for each sample.

3. **ONT reference mapping**
   Maps trimmed reads to the selected reference genome using `minimap2` with the ONT-specific `map-ont` preset.

4. **Primer trimming**
   Removes primer-derived sequences using `ivar trim` and the corresponding primer BED file.

5. **BAM processing**
   Removes unmapped reads, sorts the alignments, and generates an indexed BAM file for each sample.

6. **Mapping and coverage statistics**
   Calculates alignment statistics, reference coverage, and per-position sequencing depth using `samtools`.

7. **Consensus generation**
   Generates a consensus sequence for each sample using `samtools mpileup` and `ivar consensus`. Ambiguous nucleotide calls and gaps are converted to `N`.

8. **Variant calling and annotation**
   Detects variants using `ivar variants` and annotates them using the selected reference FASTA and GFF3 files.

9. **Result consolidation**
   Combines individual sample outputs into analysis-wide mapping statistics, consensus sequences, variant tables, and coverage data.

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
    ├──────────────► Consensus sequences
    │
    └──────────────► Variant calling
                            │
                            ▼
                     Consolidated results
```

---

## Pipeline details

### Read trimming and quality control

Raw ONT reads are processed using `fastp`.

The current pipeline performs:

* front-end quality trimming
* tail-end quality trimming
* sliding-window quality filtering
* minimum Phred-quality filtering
* filtering based on the proportion of unqualified bases
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

Quality-control reports are generated as:

```text
qc/<sample>_fastp.html
qc/<sample>_fastp.json
```

---

### Reference mapping

Trimmed reads are aligned to the reference genome using `minimap2`.

The pipeline uses:

```bash
minimap2 -x map-ont --secondary=no
```

The `map-ont` preset is designed for Oxford Nanopore sequencing reads.

Secondary alignments are disabled, and unmapped reads are removed before downstream processing.

---

### Primer trimming

Primer-derived sequences are removed from mapped reads using:

```text
ivar trim
```

with the corresponding primer BED file:

```text
reference/<reference>.primer.bed
```

Primer-trimmed alignments are then sorted to produce the final BAM file.

---

### Final BAM files

The final mapped and primer-trimmed alignment for each sample is written to:

```text
mappings/<sample>/<sample>.bam
```

with the corresponding BAM index:

```text
mappings/<sample>/<sample>.bam.bai
```

For example:

```text
mappings/sample01/sample01.bam
mappings/sample01/sample01.bam.bai
```

---

## Mapping and coverage statistics

MiCov calculates alignment and sequencing statistics using:

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

### Mapping statistics report

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
* `startpos` — first reported reference position
* `endpos` — final reported reference position
* `numreads` — number of reads reported by `samtools coverage`
* `covbases` — number of covered reference bases
* `coverage` — percentage of reference bases covered
* `meandepth` — mean sequencing depth
* `meanbaseq` — mean base quality
* `meanmapq` — mean mapping quality
* `primary_mapped` — number of primary mapped reads reported by `samtools flagstat`
* `nreads_raw` — number of reads counted in the trimmed FASTQ file

> **Note:** Despite the current column name `nreads_raw`, the value is calculated from `trimmed/<sample>_trim.fastq.gz`. It therefore represents the number of reads remaining after `fastp` processing rather than the number of reads in the original untrimmed FASTQ file.

---

## Consensus sequence generation

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

Gap characters are also converted to `N`, and the resulting sequences are converted to uppercase.

Individual consensus sequences are written to:

```text
consensus/<sample>.fasta
```

All sample consensus sequences are subsequently concatenated into:

```text
results/consensus_sequences.fasta
```

---

## Variant calling

Variants are detected using:

```text
samtools mpileup
```

followed by:

```text
ivar variants
```

The current minimum alternate allele-frequency threshold is:

```text
0.01
```

corresponding to 1%.

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

are used for variant calling and annotation.

Individual sample variant tables are written to:

```text
variants/<sample>.tsv
```

---

## Result consolidation

After all samples have been processed, MiCov consolidates the individual outputs into analysis-wide result files.

Sample identifiers are added to the variant and coverage tables before the individual files are merged.

The main consolidated outputs are written to:

```text
results/
```

---

## Output files

The main result files generated by `micov.sh` are:

```text
results/
├── mapstats.tsv
├── consensus_sequences.fasta
├── sleek_variants.tsv
└── coverage.csv
```

### `mapstats.tsv`

Contains consolidated mapping and sequencing statistics for all samples.

### `consensus_sequences.fasta`

Multi-FASTA file containing the final consensus sequence generated for each sample.

### `sleek_variants.tsv`

Combined iVar variant table containing variants detected across all samples.

Individual sample variant tables generated by `ivar variants` are consolidated into this file for downstream analysis.

### `coverage.csv`

Combined per-position sequencing-depth data generated from `samtools depth`.

The sample identifier is added to the per-position coverage data before individual sample files are merged.

---

## Complete output structure

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

## Downstream analysis

MiCov produces two consolidated datasets used by the downstream R workflows:

```text
results/
├── coverage.csv
└── sleek_variants.tsv
```

The repository includes:

```text
coverage_visualisation.R
variants_data_prep.R
```

Both scripts require:

```r
library(tidyverse)
```

---

### Coverage analysis

The `coverage_visualisation.R` script processes `coverage.csv` files generated by MiCov and summarizes sequencing coverage across the viral reference genomes.

The script currently includes virus-specific analyses for:

* **hCoV-229E**
* **hCoV-NL63**
* **hCoV-OC43**
* **hCoV-HKU1**

For each coronavirus, the script:

1. loads the corresponding `coverage.csv` file
2. renames the input columns to sample, genomic position, and sequencing depth
3. calculates mean and median sequencing depth for each genomic position
4. calculates overall sequencing-depth statistics
5. determines whether individual positions have sufficient coverage for consensus generation
6. calculates consensus completeness for each sample
7. assigns genome positions to virus-specific annotated genes and genomic regions
8. calculates region-level coverage statistics
9. generates genome-wide sequencing coverage plots
10. exports publication-quality TIFF figures

The downstream workflow can be summarized as:

```text
results/coverage.csv
        │
        ▼
Load coverage data
        │
        ▼
Calculate per-position coverage
        │
        ├────► Mean / median sequencing depth
        │
        ├────► Consensus completeness
        │
        └────► Per-sample coverage statistics
        │
        ▼
Assign genomic regions
        │
        ▼
Calculate region-level coverage
        │
        ▼
Generate genome-wide coverage plot
        │
        ▼
coverage_<virus>.tiff
```

The script uses virus-specific reference genome coordinates to annotate regions including:

```text
5'UTR
ORF1ab
S
E
M
N
3'UTR
```

along with additional virus-specific genes and ORFs.

A genomic position is considered sufficiently covered for consensus-completeness calculations when:

```text
coverage > 9
```

corresponding to a minimum sequencing depth of 10 reads.

#### Input

The script uses the `coverage.csv` output generated by MiCov.

An example organization for an analysis containing all four seasonal human coronaviruses is:

```text
coverage_data/
├── e229/
│   └── coverage.csv
├── nl63/
│   └── coverage.csv
├── oc43/
│   └── coverage.csv
└── HKU1/
    └── coverage.csv
```

Before running the script, replace:

```r
setwd("/path/to/coverage/data")
```

with the directory containing the MiCov coverage results.

#### Output

The script generates:

```text
coverage_229e.tiff
coverage_nl63.tiff
coverage_oc43.tiff
coverage_hku1.tiff
```

The output location is controlled by the `path` argument supplied to `ggsave()` and should be adjusted for the local analysis environment.

---

### Variant data preparation

The `variants_data_prep.R` script performs downstream processing of the consolidated iVar variant table generated by MiCov:

```text
results/sleek_variants.tsv
```

The script requires:

```r
library(tidyverse)
```

#### 1. Load variants

Set the working directory to the directory containing `sleek_variants.tsv`:

```r
setwd("/path/to/results")
```

and load the consolidated table:

```r
d <- read_tsv("sleek_variants.tsv")
```

#### 2. Filter variants

By default, only variants satisfying:

```r
PASS == TRUE
```

are retained.

The filtering line can be removed if all detected variants should be analysed regardless of the iVar Fisher's exact test result.

#### 3. Classify variants by allele frequency

Variants are divided into two frequency classes:

```text
ALT_FREQ >= 0.5  → major
ALT_FREQ < 0.5   → minor
```

#### 4. Extract gene annotations

Gene information is extracted from the iVar:

```text
GFF_FEATURE
```

field.

#### 5. Generate nucleotide mutation identifiers

Reference nucleotide, genomic position, and alternate allele information are combined to generate a compact nucleotide mutation identifier.

#### 6. Classify mutation type

Variants are categorized as:

* `Deletion`
* `Insertion`
* `SNP`
* `Non-coding`

#### 7. Classify coding substitutions

Coding SNPs are further classified as:

* `Synonymous`
* `Non-synonymous`

#### 8. Remove duplicate variants

Duplicate sample/mutation combinations are removed from the dataset.

#### 9. Export the enhanced variant table

The processed dataset is exported as:

```text
variants_enhanced.csv
```

The output path in the script should be adjusted for the local analysis environment.

The workflow can be summarized as:

```text
results/sleek_variants.tsv
           │
           ▼
      Load iVar results
           │
           ▼
    Filter PASS variants
           │
           ▼
  Classify allele frequency
      │             │
    major          minor
           │
           ▼
    Extract gene annotation
           │
           ▼
    Identify mutation type
     │       │       │
    SNP   Insertion Deletion
     │
     ├── Synonymous
     └── Non-synonymous
           │
           ▼
     Remove duplicates
           │
           ▼
 variants_enhanced.csv
```

---

## Logging

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

## Important assumptions

The current implementation assumes that:

1. sequencing data consist of Oxford Nanopore reads
2. input reads are gzip-compressed FASTQ files
3. FASTQ files follow the `<sample>.fastq.gz` naming convention
4. ONT barcode FASTQ files have been appropriately renamed before running MiCov
5. sample identifiers do not contain additional periods
6. the appropriate reference basename includes the amplicon length
7. matching FASTA, GFF3, and primer BED files are available for the selected reference
8. all three reference files use the same basename
9. the `reference/` directory is located directly inside the working directory
10. `micov.sh` is executed using Bash
11. the required bioinformatics software is available through the MiCov Singularity container
12. pipeline quality thresholds are reviewed and adjusted when appropriate for the ONT sequencing run
13. sufficient CPU, memory, and disk space are available for mapping and intermediate files

When repeating an analysis, using a new working directory is recommended to avoid mixing files generated by different sequencing runs or parameter configurations.

---

## Citation

If MiCov is used in published work, please cite this repository together with the relevant publication describing the seasonal human coronavirus amplicon sequencing workflow:

Ošep, A., Goldswain, H., Suljič, A. *et al.* Development of type-specific amplicon schemes for whole-genome sequencing of seasonal human coronaviruses from clinical samples. **Scientific Reports** (2026).
https://doi.org/10.1038/s41598-026-63549-1

---
