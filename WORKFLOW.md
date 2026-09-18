# Illumina Influenza A Variant Calling and Annotation Workflow

## Overview

This workflow identifies single-nucleotide variants (SNVs) in Illumina paired-end sequencing data from Influenza A virus using Galaxy.

The workflow is designed for analysis of viral population variation relative to a sample-specific or laboratory Influenza A reference genome consisting of the eight viral genome segments:

- PB2
- PB1
- PA
- HA
- NP
- NA
- M
- NS

The pipeline performs:

1. FASTQ organization
2. Raw-read quality control
3. Read preprocessing with `fastp`
4. Post-trimming quality control
5. Alignment with `BWA-MEM`
6. PCR/optical duplicate removal
7. Alignment quality assessment
8. Per-segment coverage assessment
9. Optional consensus construction
10. Local realignment using LoFreq Viterbi
11. BAQ/IDAQ calculation
12. SNV calling using LoFreq
13. Variant filtering
14. Custom SnpEff database construction
15. Functional variant annotation

The analysis is intentionally restricted to **SNVs** for the variant-analysis component.

---

# Workflow Summary

```text
Paired-end FASTQ
       │
       ▼
   FastQC #1
       │
       ▼
     fastp
       │
       ▼
   FastQC #2
       │
       ▼
    BWA-MEM
       │
       ▼
 MarkDuplicates
       │
       ├──────────────► Samtools flagstat
       │
       ├──────────────► Samtools idxstats
       │
       ├──────────────► Samtools depth
       │
       │
       ├──────────────► iVar consensus [optional]
       │
       ▼
 LoFreq Viterbi
       │
       ▼
LoFreq BAQ + IDAQ
       │
       ▼
 LoFreq Call Variants
       │
       ▼
  LoFreq Filter
       │
       ▼
   Filtered VCF
       │
       ▼
  SnpEff Annotation
       │
       ▼
Annotated VCF + CSV/report outputs
````

---

# 1. Input Data Organization

## Galaxy collection type

Import sequencing reads as:

**Collection: `List of Pairs`**

Each biological sample must therefore contain:

```text
Sample_01
├── Sample_01_R1.fastq.gz
└── Sample_01_R2.fastq.gz

Sample_02
├── Sample_02_R1.fastq.gz
└── Sample_02_R2.fastq.gz
```

Using a paired collection allows the same Galaxy workflow to be executed independently across multiple samples while preserving the R1/R2 pairing.

---

# 2. Raw Read Quality Control — FastQC

**Tool:** FastQC
**Galaxy version:** `0.74+galaxy1`

Run FastQC on the raw paired-end FASTQ datasets before trimming.

## Parameters

| Parameter       | Setting |
| --------------- | ------- |
| Tool parameters | Default |

## Purpose

Inspect:

* per-base sequence quality
* per-sequence quality
* per-base sequence content
* GC distribution
* sequence length distribution
* duplication levels
* adapter content
* overrepresented sequences

Do not remove reads simply because a FastQC module reports `WARN` or `FAIL`; use the reports to determine whether preprocessing is necessary.

---

# 3. Read Preprocessing — fastp

**Tool:** fastp
**Galaxy version:** `1.3.6+galaxy0`

## Parameters

| Parameter                                  | Setting               |
| ------------------------------------------ | --------------------- |
| Single-end or paired reads                 | **Paired Collection** |
| Merge forward and reverse reads            | **No**                |
| Disable adapter trimming                   | **No**                |
| Adapter sequence R1                        | Blank                 |
| Adapter sequence R2                        | Blank                 |
| Adapter auto-detection for paired-end data | **Yes**               |
| Enable overrepresented sequence analysis   | **Yes**               |
| Overrepresentation sampling                | `20` / default        |
| Disable quality filtering                  | **No**                |
| Qualified quality Phred                    | `20`                  |
| Unqualified percent limit                  | `20`                  |
| N-base limit                               | `5`                   |
| Disable length filtering                   | **No**                |
| Minimum required length                    | `50 bp`               |
| Maximum length                             | `0` / no limit        |
| Enable low-complexity filter               | **No**                |
| Enable duplicate-read analysis             | **Yes**               |
| Drop duplicate reads/pairs                 | **No**                |
| PolyG trimming                             | Automatic             |
| PolyG minimum length                       | `10` / default        |
| PolyX trimming                             | **No**                |
| Enable UMI                                 | **No**                |
| Cut by quality at 5′ end                   | **No**                |
| Cut by quality at 3′ end                   | **No**                |
| Cut right                                  | **No**                |
| Enable base correction                     | **Yes**               |
| Output HTML report                         | **Yes**               |
| Output JSON report                         | **Yes**               |

## Important distinction

Duplicate reads are **analyzed but not removed by fastp**.

```text
Enable duplicated-read analysis = YES
Drop duplicate reads/pairs       = NO
```

Duplicate handling is performed later after alignment using Picard `MarkDuplicates`.

## Output

For every sample:

```text
filtered_R1.fastq.gz
filtered_R2.fastq.gz
fastp.html
fastp.json
```

---

# 4. Post-Preprocessing QC — FastQC

**Tool:** FastQC
**Galaxy version:** `0.74+galaxy1`

Run FastQC again on the fastp-filtered reads.

## Parameters

| Parameter       | Setting |
| --------------- | ------- |
| Tool parameters | Default |

Compare the post-fastp results with the original FastQC reports.

Particular attention should be paid to:

* base quality
* adapter contamination
* sequence length
* overrepresented sequences
* sequence-content bias

---

# 5. Reference Genome

Use one consistent Influenza A reference FASTA throughout the complete workflow.

The reference must contain the eight Influenza A genome segments as separate FASTA records.

Example:

```fasta
>PB2
ATG...

>PB1
ATG...

>PA
ATG...

>HA
ATG...

>NP
ATG...

>NA
ATG...

>M
ATG...

>NS
ATG...
```

The same reference sequence must be used for:

* BWA-MEM
* LoFreq Viterbi
* LoFreq alignment-quality calculation
* LoFreq variant calling
* SnpEff database construction

Consistency of chromosome/segment names is essential.

For example:

```text
PB2
PB1
PA
HA
NP
NA
M
NS
```

must be used consistently between the FASTA, GFF3, BAM and VCF files.

---

# 6. Read Mapping — BWA-MEM

**Tool:** Map with BWA-MEM
**Galaxy version:** `0.7.19+galaxy1`

## Parameters

| Parameter          | Setting                                         |
| ------------------ | ----------------------------------------------- |
| Reference source   | Use a genome from history and build index       |
| Reference sequence | Exact eight-segment Influenza A reference FASTA |
| BWT algorithm      | Auto                                            |
| Reads              | **Paired**                                      |
| First reads        | fastp-filtered R1                               |
| Second reads       | fastp-filtered R2                               |
| Insert length      | Blank                                           |
| Analysis mode      | **Simple Illumina mode**                        |
| BAM sorting        | **Sort by chromosomal coordinates**             |

## Output

Coordinate-sorted BAM containing reads mapped against all eight influenza segments.

---

# 7. Duplicate Removal — MarkDuplicates

**Tool:** MarkDuplicates
**Galaxy version:** `3.1.1.0`

Input:

```text
BWA-MEM coordinate-sorted BAM
```

## Parameters

| Parameter                          | Setting                 |
| ---------------------------------- | ----------------------- |
| Select SAM/BAM dataset             | BWA-MEM BAM             |
| Comment                            | Blank                   |
| `REMOVE_DUPLICATES`                | **True**                |
| `ASSUME_SORTED`                    | **True**                |
| `DUPLICATE_SCORING_STRATEGY`       | `SUM_OF_BASE_QUALITIES` |
| `READ_NAME_REGEX`                  | Blank / default         |
| `OPTICAL_DUPLICATE_PIXEL_DISTANCE` | `100`                   |
| Barcode tag                        | Blank                   |
| Validation stringency              | **Lenient**             |

## Output

```text
deduplicated.bam
```

This BAM becomes the primary alignment dataset for downstream analysis.

---

# 8. Mapping Statistics — Samtools flagstat

**Tool:** Samtools flagstat
**Galaxy version:** `2.0.8`

## Parameters

| Parameter     | Setting          |
| ------------- | ---------------- |
| Input         | Deduplicated BAM |
| Output format | `txt`            |

## Evaluate

Important statistics include:

```text
total reads
mapped reads
mapping percentage
properly paired reads
singletons
secondary alignments
supplementary alignments
duplicates
```

The mapping statistics should be evaluated for every sample before variant calling.

---

# 9. Segment Representation — Samtools idxstats

**Tool:** Samtools idxstats
**Galaxy version:** `2.0.8`

## Parameters

| Parameter | Setting          |
| --------- | ---------------- |
| Input     | Deduplicated BAM |

## Purpose

Use `idxstats` to determine whether reads mapped to each Influenza A segment.

Expected segments:

```text
PB2
PB1
PA
HA
NP
NA
M
NS
```

The output is particularly useful for detecting:

* poorly represented segments
* unexpectedly high/low segment abundance
* completely missing segments
* unmapped reads

---

# 10. Position-Level Coverage — Samtools depth

**Tool:** Samtools depth
**Galaxy version:** `1.22+galaxy1`

## Parameters

| Parameter                                     | Setting                                  |
| --------------------------------------------- | ---------------------------------------- |
| Filter by regions                             | **No**                                   |
| Output all positions                          | **Yes — including zero-depth positions** |
| Ignore reads shorter than                     | Blank                                    |
| Maximum alignments beginning at each position | `0`                                      |
| Minimum base quality                          | `20`                                     |
| Minimum mapping quality                       | `20`                                     |
| Exclude unmapped reads                        | **Yes**                                  |
| Exclude non-primary alignments                | **Yes**                                  |
| Exclude reads failing platform/vendor QC      | **Yes**                                  |
| Exclude PCR/optical duplicates                | **Yes**                                  |
| Require flags (`-g`)                          | None / blank                             |
| Include deletions (`-J`)                      | **No**                                   |
| Count overlapping paired reads once (`-s`)    | **Yes**                                  |
| Print header (`-H`)                           | **Yes**                                  |

## Purpose

This step produces per-nucleotide coverage and allows examination of coverage continuity across all eight genome segments.

Because zero-depth positions are included, this output can identify:

```text
coverage gaps
segment ends with poor coverage
completely missing positions
regions below the variant-calling threshold
```

---

# 11. Optional Consensus Genome — iVar

**Tool:** iVar consensus
**Galaxy version:** `1.4.4+galaxy0`

This step is optional for variant calling but can be used to reconstruct a majority-rule consensus genome.

## Parameters

| Parameter                      | Setting                   |
| ------------------------------ | ------------------------- |
| BAM file                       | MarkDuplicates output BAM |
| Minimum quality score (`-q`)   | `20`                      |
| Minimum frequency (`-t`)       | `0.50`                    |
| Minimum indel frequency (`-c`) | `0.50`                    |
| Minimum depth (`-m`)           | `10`                      |
| Low-coverage representation    | `N`                       |

Thus, a nucleotide becomes part of the consensus when:

```text
frequency ≥ 50%
depth ≥ 10
base quality ≥ Q20
```

Positions failing the depth requirement are represented as:

```text
N
```

---

# 12. LoFreq Viterbi Realignment

**Tool:** Realign reads with LoFreq Viterbi
**Galaxy version:** `2.1.5+galaxy0`

## Parameters

| Parameter                        | Setting                                    |
| -------------------------------- | ------------------------------------------ |
| Reads to realign                 | MarkDuplicates output BAM                  |
| Reference source                 | History                                    |
| Reference                        | Exact same eight-segment Influenza A FASTA |
| Keep MC, MD, NM and A flags      | **No / unchecked**                         |
| Handle base qualities equal to 2 | **Keep unchanged**                         |

## Purpose

LoFreq Viterbi performs local realignment to improve the placement of reads around regions where alignment ambiguity may occur.

Output:

```text
LoFreq-Viterbi-realigned BAM
```

---

# 13. Add LoFreq Alignment Quality Scores

**Tool:** Add LoFreq alignment quality scores
**Galaxy version:** `2.1.5+galaxy1`

Input:

```text
LoFreq Viterbi BAM
```

## Parameters

| Parameter                        | Setting                                             |
| -------------------------------- | --------------------------------------------------- |
| Reads                            | LoFreq Viterbi output BAM                           |
| Reference source                 | History                                             |
| Reference                        | Exact same Influenza A reference                    |
| Alignment quality scores         | **Base and indel alignment qualities (BAQ + IDAQ)** |
| Extended BAQ (`-e`)              | **Yes**                                             |
| Overwrite existing values (`-r`) | **Yes**                                             |

## Output

A BAM containing alignment-aware base and indel quality information for LoFreq variant calling.

---

# 14. SNV Calling — LoFreq

**Tool:** Call variants with LoFreq
**Galaxy version:** `2.1.5+galaxy3`

Input:

```text
BAM from Add LoFreq alignment quality scores
```

## General Parameters

| Parameter                  | Setting                                |
| -------------------------- | -------------------------------------- |
| Reference source           | History                                |
| Reference                  | Exact same 8-segment Influenza A FASTA |
| Call variants across       | **Whole reference**                    |
| Variant types              | **Only SNVs**                          |
| Variant calling parameters | Configure settings                     |

## Coverage

| Parameter                      |     Setting |
| ------------------------------ | ----------: |
| Minimal coverage (`--min-cov`) |        `10` |
| Coverage cap (`--max-depth`)   | `1,000,000` |

## Read Pair Handling

| Parameter                                            | Setting      |
| ---------------------------------------------------- | ------------ |
| Use anomalously mapped/orphan pairs (`--use-orphan`) | **No / OFF** |

## Base Quality

| Parameter                                       |                 Setting |
| ----------------------------------------------- | ----------------------: |
| Minimum base quality (`--min-bq`)               |                    `20` |
| Minimum alternate-base quality (`--min-alt-bq`) |                    `20` |
| Base quality used for alternate bases           | Original base qualities |

## Alignment Quality

| Parameter                               | Setting             |
| --------------------------------------- | ------------------- |
| Consider base/indel alignment qualities | **Yes**             |
| Alignment qualities                     | Existing BAQ + IDAQ |
| Extended BAQ (`-e`)                     | **Yes / ON**        |

## Mapping Quality

| Parameter                               | Setting |
| --------------------------------------- | ------: |
| Minimum mapping quality (`--min-mq`)    |    `20` |
| Consider mapping quality during calling | **Yes** |
| Maximum mapping quality (`--max-mq`)    |   `255` |

Mapping quality is incorporated into LoFreq's joint quality score.

## Source / Joint Quality

| Parameter                                                    | Setting |
| ------------------------------------------------------------ | ------: |
| Compute source quality                                       |  **No** |
| Minimum joined quality (`--min-jq`)                          |     `0` |
| Minimum joined quality for alternate bases (`--min-alt-jq`)  |     `0` |
| Overwrite joined quality of alternate bases (`--def-alt-jq`) |     `0` |

## Internal Variant Filtering

Use:

```text
Preset filtering on QUAL + coverage + strand bias
(LoFreq default)
```

## Output

```text
raw_lofreq_snvs.vcf
```

---

# 15. Post-Calling Variant Filtering — LoFreq Filter

**Tool:** LoFreq filter
**Galaxy version:** `2.1.5+galaxy0`

Input:

```text
LoFreq Call VCF
```

## Variant Type

| Parameter                 | Setting       |
| ------------------------- | ------------- |
| Types of variants to keep | **SNVs only** |

## Call Quality

| Parameter                         | Setting |
| --------------------------------- | ------- |
| Filter SNVs based on call quality | **Yes** |
| Minimum QUAL                      | `0`     |

## Coverage

| Parameter        |              Setting |
| ---------------- | -------------------: |
| Minimum coverage |              **500** |
| Maximum coverage | `0` / no upper limit |

The final analysis therefore retains variants supported by:

```text
DP ≥ 500
```

## Variant Allele Frequency

| Parameter                |              Setting |
| ------------------------ | -------------------: |
| Minimum allele frequency |             **0.05** |
| Maximum allele frequency | `0` / no upper limit |

The final analysis therefore detects variants at:

```text
AF ≥ 5%
```

subject to the other quality criteria.

## Strand Bias

| Parameter                        | Setting                                        |
| -------------------------------- | ---------------------------------------------- |
| Strand-bias filtering            | **Yes**                                        |
| Test                             | Multiple-testing corrected strand-bias p-value |
| Corrected p-value threshold      | `0.001`                                        |
| Multiple-testing correction      | **False Discovery Rate (FDR)**                 |
| Compound strand-bias filter      | **No / OFF**                                   |
| Strand-bias filtering for indels | **No / OFF**                                   |

## Failed Variants

| Parameter | Setting                                       |
| --------- | --------------------------------------------- |
| Action    | **Drop variants failing one or more filters** |

## Final SNV Criteria

The principal explicit downstream thresholds are:

```text
Variant type = SNV
DP ≥ 500
AF ≥ 0.05
strand-bias adjusted P ≥ required LoFreq filter criterion
```

with upstream read/base filtering requiring:

```text
Base quality ≥ Q20
Mapping quality ≥ 20
```

Output:

```text
filtered_lofreq_snvs.vcf
```

---

# 16. Build Custom Influenza SnpEff Database

Because the analysis uses a custom Influenza A reference, build a SnpEff database corresponding exactly to the reference genome.

**Tool:** SnpEff build
**Galaxy version:** `5.4+galaxy0`

## Required Inputs

```text
Reference FASTA
+
matching GFF3 annotation
```

## Parameters

| Parameter                | Setting      |
| ------------------------ | ------------ |
| Input annotations are in | **GFF**      |
| Reference genome source  | **History**  |
| Genetic code             | **Standard** |

## Critical Requirement

The sequence identifiers in the GFF3 must exactly match the FASTA headers.

Example:

```text
FASTA:
>PB2

GFF3:
PB2    RefSeq    gene ...
```

Do not mix identifiers such as:

```text
PB2
segment_1
NC_026438
```

unless the same naming system is used everywhere.

---

# 17. Functional Annotation — SnpEff eff

**Tool:** SnpEff eff
**Galaxy version:** `5.2+galaxy1`

Input:

```text
filtered_lofreq_snvs.vcf
```

## Parameters

| Parameter                       | Setting                                        |
| ------------------------------- | ---------------------------------------------- |
| Sequence changes                | Filtered sample VCF                            |
| Input format                    | **VCF**                                        |
| Output format                   | **VCF**                                        |
| Create CSV report               | **Yes**                                        |
| Produce Summary Statistics      | **Yes**                                        |
| Produce Gene Statistics         | **Yes**                                        |
| Genome source                   | **Custom SnpEff database from Galaxy history** |
| Genetic code                    | **Standard**                                   |
| Upstream/downstream length      | `0 bp`                                         |
| Splice-site donor/acceptor size | `2 bp`                                         |
| Splice-region settings          | Default                                        |
| Effect terminology              | **Sequence Ontology**                          |
| Amino-acid annotation           | **HGVS**                                       |
| Add LOF/NMD annotations         | **No**                                         |
| Filter specific effects         | **No**                                         |
| Suppress usage statistics       | **Yes**                                        |

## Example Annotation Categories

SnpEff may classify SNVs as:

```text
synonymous_variant
missense_variant
stop_gained
stop_lost
start_lost
intergenic_region
upstream_gene_variant
downstream_gene_variant
```

Because upstream/downstream length is set to zero in this workflow, annotation is focused primarily on sequence features defined explicitly in the supplied GFF3.

---

# 18. Recommended Final Variant Table

For downstream analyses, convert the annotated VCF into a tabular dataset containing at least:

| Field   | Description                     |
| ------- | ------------------------------- |
| Sample  | Sample ID                       |
| Segment | PB2/PB1/PA/HA/NP/NA/M/NS        |
| POS     | Reference nucleotide position   |
| REF     | Reference allele                |
| ALT     | Alternate allele                |
| DP      | Read depth                      |
| AF      | Variant allele frequency        |
| QUAL    | Variant call quality            |
| Gene    | Influenza gene                  |
| Effect  | SnpEff Sequence Ontology effect |
| HGVS.c  | Coding nucleotide change        |
| HGVS.p  | Amino-acid change               |
| Impact  | LOW/MODERATE/HIGH/MODIFIER      |

Example:

```text
Sample    Segment  POS   REF ALT   DP    AF     Effect             HGVS.p
Pig01     HA       548   G   A     1821  0.23   missense_variant   p.Asp173Asn
Pig01     PB2      2026  C   A     2150  0.67   missense_variant   p.Thr676Asn
Pig02     NP       771   T   C     1940  0.08   synonymous_variant p.Leu257=
```

---

# 19. QC Checkpoints

Do not interpret variants until the following QC stages have been reviewed.

## FASTQ QC

```text
FastQC raw
   ↓
fastp
   ↓
FastQC filtered
```

Confirm that trimming did not cause excessive read loss.

## Alignment QC

Review:

```text
Samtools flagstat
Samtools idxstats
Samtools depth
```

Look for:

* high overall mapping rate
* proper pairing
* representation of all expected segments
* coverage gaps
* abnormal segment-specific coverage
* unexpectedly high unmapped-read counts

## Variant QC

For each retained variant inspect:

```text
DP
AF
QUAL
strand-bias status
functional annotation
```

For biologically important variants, manual inspection of the BAM in IGV/JBrowse is recommended.

---

# 20. Workflow Parameter Summary

| Step | Tool           | Key parameters                                                            |
| ---- | -------------- | ------------------------------------------------------------------------- |
| 1    | Import         | List of paired datasets                                                   |
| 2    | FastQC         | Default                                                                   |
| 3    | fastp          | Q20; ≤20% low-quality bases; ≤5 Ns; min length 50; adapters auto-detected |
| 4    | FastQC         | Default                                                                   |
| 5    | BWA-MEM        | Paired; Simple Illumina; coordinate sorted                                |
| 6    | MarkDuplicates | Remove duplicates = TRUE                                                  |
| 7    | flagstat       | Deduplicated BAM                                                          |
| 8    | idxstats       | Deduplicated BAM                                                          |
| 9    | depth          | BQ≥20; MQ≥20; overlapping mates counted once                              |
| 10   | iVar           | Q20; frequency 0.5; depth 10; low coverage=N                              |
| 11   | LoFreq Viterbi | Same influenza reference                                                  |
| 12   | LoFreq alnqual | BAQ + IDAQ; extended BAQ                                                  |
| 13   | LoFreq Call    | SNVs; minCov=10; BQ≥20; MQ≥20                                             |
| 14   | LoFreq Filter  | **DP≥500; AF≥0.05; FDR strand-bias filter**                               |
| 15   | SnpEff build   | Custom GFF + reference                                                    |
| 16   | SnpEff eff     | SO effects + HGVS; standard genetic code                                  |

---

# 21. Recommended Repository Structure

```text
influenza-illumina-variant-analysis/
│
├── README.md
├── WORKFLOW.md
│
├── workflow/
│   ├── influenza_variant_calling.ga
│   └── tool_parameters.xlsx
│
├── reference/
│   ├── README.md
│   ├── influenza_reference.fasta
│   └── influenza_reference.gff3
│
├── docs/
│   ├── pipeline_overview.md
│   ├── QC_guidelines.md
│   └── variant_interpretation.md
│
├── examples/
│   ├── example_variants.tsv
│   └── example_annotation.tsv
│
└── results/
    └── README.md
```

Large FASTQ, BAM and intermediate sequencing datasets should generally not be committed directly to GitHub.

---

# 22. Minimal Methods Description

Illumina paired-end reads were initially assessed using FastQC and processed with fastp. Adapter trimming and paired-end adapter autodetection were enabled, and reads were filtered using a minimum Phred quality score of 20, a maximum of 20% unqualified bases, a maximum of five ambiguous bases, and a minimum read length of 50 bp. Filtered paired-end reads were mapped against an eight-segment Influenza A reference genome using BWA-MEM in simple Illumina mode, and alignments were coordinate sorted. PCR and optical duplicates were subsequently removed using Picard MarkDuplicates.

Mapping quality and genome coverage were assessed using Samtools flagstat, idxstats, and depth. For depth calculations, bases and reads with quality values below Q20 and MAPQ 20, respectively, were excluded and overlapping paired-end reads were counted once.

Deduplicated reads were realigned using LoFreq Viterbi, followed by calculation of base and indel alignment qualities using extended BAQ. LoFreq was used to call SNVs across the complete eight-segment reference genome with minimum base and mapping quality thresholds of 20. Variants were subsequently filtered to retain SNVs with a minimum depth of 500 reads and a variant allele frequency of at least 5%. Multiple-testing-corrected strand-bias filtering was performed using a false-discovery-rate correction with a threshold of 0.001.

A custom SnpEff database was generated from the reference genome FASTA and corresponding GFF3 annotation. Filtered variants were annotated using Sequence Ontology effect terms and HGVS nomenclature.

---

# References / Workflow Basis

This workflow was developed using the accompanying Galaxy parameter record and Galaxy Training Network documentation for viral and non-diploid variant analysis.

Relevant concepts include:

* paired-end FASTQ processing and Galaxy dataset collections
* BWA-MEM reference mapping
* duplicate handling
* LoFreq Viterbi realignment
* BAQ/IDAQ calculation
* LoFreq variant discovery
* variant filtering
* SnpEff functional annotation
* viral/non-diploid allele-frequency analysis
* segmented Influenza A genome analysis
