# Influenza Variant Analysis

## A Practical Guide from Raw Illumina Reads to Interpreted Viral Variants

This guide describes a reference-based workflow for **influenza A variant analysis**, with emphasis on paired-end Illumina sequencing, viral/non-diploid variant calling, and implementation in **Galaxy**.

It covers:

- Core variant-analysis concepts
- Influenza-specific considerations
- FASTQ quality control
- Adapter and quality trimming
- Reference preparation
- Read mapping
- BAM processing
- Coverage and mapping QC
- Duplicate handling
- LoFreq variant calling
- Variant filtering
- VCF interpretation
- Functional annotation
- Manual validation
- Consensus construction
- Recommended reporting standards

---

# Table of Contents

1. [What Variant Analysis Measures](#1-what-variant-analysis-measures)
2. [Core Terminology](#2-core-terminology)
3. [Why Influenza Is a Non-Diploid Variant-Calling Problem](#3-why-influenza-is-a-non-diploid-variant-calling-problem)
4. [Influenza Genome Architecture](#4-influenza-genome-architecture)
5. [Reference Selection](#5-reference-selection)
6. [FASTQ and Phred Quality](#6-fastq-and-phred-quality)
7. [Complete Workflow](#7-complete-workflow)
8. [Step 1: Organize Paired-End Reads](#8-step-1-organize-paired-end-reads)
9. [Step 2: Raw Read Quality Control](#9-step-2-raw-read-quality-control)
10. [Step 3: Adapter and Quality Processing](#10-step-3-adapter-and-quality-processing)
11. [Step 4: Post-Trim QC](#11-step-4-post-trim-qc)
12. [Step 5: Prepare the Influenza Reference](#12-step-5-prepare-the-influenza-reference)
13. [Step 6: Read Mapping](#13-step-6-read-mapping)
14. [SAM/BAM and Mapping Quality](#14-sambam-and-mapping-quality)
15. [Step 7: Sort the BAM](#15-step-7-sort-the-bam)
16. [Step 8: Alignment QC](#16-step-8-alignment-qc)
17. [Duplicate Handling](#17-duplicate-handling)
18. [LoFreq Variant Calling](#18-lofreq-variant-calling)
19. [Variant Filtering](#19-variant-filtering)
20. [VCF Interpretation](#20-vcf-interpretation)
21. [Functional Annotation](#21-functional-annotation)
22. [Visual Validation](#22-visual-validation)
23. [Consensus Sequence Construction](#23-consensus-sequence-construction)
24. [PCR and Amplicon Bias](#24-pcr-and-amplicon-bias)
25. [Tool Summary](#25-tool-summary)
26. [QC Checkpoints](#26-qc-checkpoints)
27. [Recommended Final Variant Table](#27-recommended-final-variant-table)
28. [Worked Interpretation Example](#28-worked-interpretation-example)
29. [Practical Galaxy Workflow](#29-practical-galaxy-workflow)
30. [References](#30-references)

---

# 1. What Variant Analysis Measures

Variant analysis asks:

> At each nucleotide position in the sequenced viral population, is there reliable evidence that one or more nucleotides differ from the selected reference sequence?

Example:

```text
Reference:  ... A C T G A A C ...
Read 1:     ... A C T G A A C ...
Read 2:     ... A C T A A A C ...
Read 3:     ... A C T G A A C ...
Read 4:     ... A C T A A A C ...
````

A possible variant call could be:

```text
REF = G
ALT = A
AF ≈ 0.50
```

A variant is therefore always defined **relative to a reference sequence**.

The Galaxy microbial variant-calling tutorial defines variant calling as the process of identifying differences between sequencing reads and a known reference genome, commonly focusing on SNPs/SNVs and small insertions or deletions.

---

# 2. Core Terminology

| Term                   | Definition                                                                                                 |
| ---------------------- | ---------------------------------------------------------------------------------------------------------- |
| **Variant**            | A sequence state that differs from the selected reference.                                                 |
| **SNV**                | Single-nucleotide variant; one nucleotide differs from the reference.                                      |
| **SNP**                | Single-nucleotide polymorphism. Frequently used interchangeably with SNV in software and older literature. |
| **Indel**              | An insertion or deletion relative to the reference.                                                        |
| **REF**                | Reference allele at a genomic position.                                                                    |
| **ALT**                | Alternate allele supported by sequencing reads.                                                            |
| **AF / VAF**           | Alternate or variant allele frequency.                                                                     |
| **DP**                 | Read depth or coverage at a site.                                                                          |
| **QUAL**               | Variant-caller confidence score.                                                                           |
| **MAPQ**               | Mapping quality; confidence that a read is aligned to the correct location.                                |
| **Consensus sequence** | One representative nucleotide selected for every genomic position.                                         |

---

## Variant Allele Frequency

Variant allele frequency is approximately:

```text
AF = reads supporting ALT / total usable reads covering the site
```

Example:

```text
Total usable reads = 1,000
Reads supporting G = 240

AF = 240 / 1000
AF = 0.24
AF = 24%
```

AF should be interpreted as an estimate based on sequencing observations.

It is **not** a direct count of intact virions.

---

## Minority and Majority Variants

A practical descriptive distinction is:

```text
AF < 0.50  → minority variant
AF > 0.50  → majority variant
AF ≈ 1.00  → near-fixed relative to the reference
```

Example:

```text
A = 55%
G = 45%
```

The consensus sequence may contain:

```text
A
```

while the underlying viral population still contains a substantial:

```text
G = 45%
```

Therefore, consensus genomes can hide within-host viral diversity.

---

# 3. Why Influenza Is a Non-Diploid Variant-Calling Problem

Human germline variant calling is often based on diploid genotypes such as:

```text
0/0
0/1
1/1
```

which approximately correspond to:

```text
0%
50%
100%
```

Influenza samples contain populations of viral genomes.

Therefore, an alternate allele can theoretically occur at any frequency:

```text
1%
7%
23%
52%
89%
100%
```

For viral analysis, the important question is not simply:

```text
What genotype does this sample have?
```

Instead:

```text
At what frequency is the alternate allele observed,
and is that frequency distinguishable from sequencing error?
```

The Galaxy tutorial on non-diploid variant calling identifies the distinction between real low-frequency variation and sequencing noise as a major analytical challenge.

---

# 4. Influenza Genome Architecture

Influenza A has eight negative-sense RNA genome segments:

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

A reference FASTA should therefore normally contain eight separate sequence records:

```fasta
>PB2
SEQUENCE...

>PB1
SEQUENCE...

>PA
SEQUENCE...

>HA
SEQUENCE...

>NP
SEQUENCE...

>NA
SEQUENCE...

>M
SEQUENCE...

>NS
SEQUENCE...
```

Influenza A evolves through both:

* Mutation
* Reassortment

Reassortment allows complete genome segments to be exchanged between viruses.

Example:

```text
PB2 → strain A
PB1 → strain A
PA  → strain B
HA  → strain C
NP  → strain A
NA  → strain D
M   → strain A
NS  → strain B
```

This makes influenza reference selection more complicated than for many non-segmented viruses.

---

# 5. Reference Selection

Reference choice directly determines which differences will be reported as variants.

Example:

```text
Sample nucleotide = G
Reference 1 = A
```

Result:

```text
A → G variant
```

But if:

```text
Sample nucleotide = G
Reference 2 = G
```

Result:

```text
No variant
```

The sample did not change.

The reference changed.

Therefore, variant calls should always be interpreted as:

> Differences relative to reference X.

---

## Influenza-Specific Reference Strategies

### Strategy A: Known Challenge/Inoculum Reference

Useful for controlled experimental infections.

```text
Challenge virus
     ↓
reference consensus
     ↓
experimental samples
```

This answers:

> What changed relative to the virus introduced into the experiment?

---

### Strategy B: Public Reference Strain

Useful when comparing many samples to a standardized external baseline.

This answers:

> How does each sample differ from the selected public reference strain?

---

### Strategy C: Segment-Specific Hybrid Reference

The Galaxy avian influenza workflow describes selecting an appropriate reference separately for each influenza segment.

Conceptually:

```text
PB2 best reference
PB1 best reference
PA best reference
HA best reference
NP best reference
NA best reference
M best reference
NS best reference
        ↓
combined hybrid reference
```

This approach can be useful when reassortment or substantial HA/NA divergence makes one fixed strain a poor mapping reference.

---

# 6. FASTQ and Phred Quality

Illumina sequencing generally produces paired FASTQ files:

```text
Sample01_R1.fastq.gz
Sample01_R2.fastq.gz
```

FASTQ stores both:

1. The nucleotide sequence
2. A quality score for every base

Example:

```text
@READ_ID
ACGTACGTACGT
+
IIIIIIIIIIII
```

A FASTQ record contains four lines:

```text
1. Read identifier
2. Sequence
3. "+" separator
4. Quality string
```

The sequence and quality strings must have the same length.

---

## Phred Quality Scores

Phred quality is defined as:

```text
Q = -10 × log10(Perror)
```

Typical values:

| Phred score | Approximate error probability |
| ----------: | ----------------------------: |
|         Q10 |                           10% |
|         Q20 |                            1% |
|         Q30 |                          0.1% |
|         Q40 |                         0.01% |

This is especially important for low-frequency variant analysis.

For example, if a variant is observed at:

```text
AF = 0.3%
```

but the sequencing error rate is also close to that range, distinguishing a true biological variant from sequencing noise becomes difficult.

---

# 7. Complete Workflow

```text
Paired Illumina FASTQ
        │
        ▼
Raw read QC
        │
        ▼
Adapter / quality trimming
        │
        ▼
Post-trim QC
        │
        ▼
Prepare influenza reference
        │
        ▼
Read mapping
        │
        ▼
BAM sorting
        │
        ▼
Alignment / coverage QC
        │
        ▼
Appropriate duplicate treatment
        │
        ▼
LoFreq preprocessing
        │
        ▼
Variant calling
        │
        ▼
Variant filtering
        │
        ▼
Functional annotation
        │
        ▼
Manual validation
        │
        ├──────────────► Within-host variant table
        │
        ▼
Consensus construction
        │
        ▼
Comparative / evolutionary analyses
```

A closely related Galaxy viral workflow uses:

```text
fastp
  ↓
BWA-MEM
  ↓
MarkDuplicates
  ↓
Samtools statistics
  ↓
LoFreq Viterbi
  ↓
LoFreq quality processing
  ↓
LoFreq variant calling
  ↓
SnpEff
  ↓
SnpSift
  ↓
MultiQC
```

---

# 8. Step 1: Organize Paired-End Reads

Each sample contains a forward and reverse FASTQ:

```text
Sample01_R1.fastq.gz
Sample01_R2.fastq.gz
```

For multiple samples in Galaxy, use a paired collection:

```text
Sample01
├── forward
└── reverse

Sample02
├── forward
└── reverse

Sample03
├── forward
└── reverse
```

This allows Galaxy tools to process all samples in parallel while preserving the relationship between R1 and R2.

---

# 9. Step 2: Raw Read Quality Control

Use:

```text
FastQC
```

Inspect:

* Per-base quality
* Per-sequence quality
* Sequence length
* Adapter content
* Overrepresented sequences
* GC distribution
* N content
* Per-base sequence composition
* Unexpected sequence at read starts or ends

The goal is **not** simply to make every FastQC category green.

Instead ask:

* Are adapters present?
* Does read quality decline strongly toward the ends?
* Is there unexpected technical sequence?
* Is the observed sequence composition expected for the library protocol?
* Are there signs of a failed sequencing run?

---

# 10. Step 3: Adapter and Quality Processing

A commonly used tool is:

```text
fastp
```

fastp can perform:

* Adapter removal
* Quality trimming
* Read filtering
* Paired-end processing
* Basic QC reporting

Workflow:

```text
RAW FASTQ
    ↓
  fastp
    ↓
CLEAN FASTQ
```

The goal is to remove technical sequence and clearly unreliable bases while retaining as much valid biological information as possible.

---

> **Do not over-trim**
>
> Do not remove sequence solely to make FastQC plots look better. Trimming should have a technical or quality-based justification.

---

# 11. Step 4: Post-Trim QC

Run FastQC again after fastp.

Compare:

```text
RAW READS
    vs.
TRIMMED READS
```

Confirm:

* Adapter contamination decreased
* Poor-quality tails improved
* Read lengths remain reasonable
* Excessive numbers of reads were not removed
* No unexpected new problems appeared

---

# 12. Step 5: Prepare the Influenza Reference

The reference should contain separate sequences for:

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

Sequence identifiers must match across:

```text
FASTA
GFF/GTF
BAM
VCF
```

Example of a correct relationship:

```text
FASTA:
>HA

GFF:
HA    source    gene    ...
```

A mismatch such as:

```text
FASTA:
>segment_4_HA

GFF:
HA
```

can prevent annotation software from recognizing that the variants belong to the HA segment.

---

# 13. Step 6: Read Mapping

For paired Illumina influenza sequencing, use:

```text
BWA-MEM
```

Input:

```text
R1 FASTQ
R2 FASTQ
Reference FASTA
```

Output:

```text
SAM/BAM
```

Conceptually:

```text
FASTQ reads
      +
Reference FASTA
      ↓
   BWA-MEM
      ↓
     BAM
```

BWA-MEM determines where each sequencing read most plausibly aligns to the reference genome.

---

# 14. SAM/BAM and Mapping Quality

## SAM

Sequence Alignment/Map format.

Text-based alignment representation.

---

## BAM

Binary compressed form of SAM.

A BAM file contains information such as:

```text
Read ID
Reference segment
Mapping coordinate
Mapping quality
CIGAR string
Strand
Mate information
Read sequence
Base quality
```

A BAM therefore describes:

> Where each read aligns and how it aligns.

---

## Mapping Quality

MAPQ describes confidence in read placement.

Low MAPQ:

```text
This read may align equally well elsewhere.
```

High MAPQ:

```text
The aligner has much stronger confidence in this placement.
```

Low-quality mappings can create false variant evidence.

---

# 15. Step 7: Sort the BAM

Use:

```text
Samtools sort
```

Workflow:

```text
Unsorted BAM
     ↓
Samtools sort
     ↓
Coordinate-sorted BAM
```

Many downstream tools require or strongly prefer coordinate-sorted BAM files.

---

# 16. Step 8: Alignment QC

Do not proceed directly from mapping to variant calling.

Inspect the alignment first.

Important metrics include:

```text
Total reads
Mapped reads
Mapping percentage
Properly paired reads
Unmapped reads
Coverage depth
Coverage uniformity
Zero-coverage positions
Per-segment coverage
```

---

## Influenza-Specific Requirement

Inspect all eight segments independently.

A sample might show:

```text
PB2    excellent coverage
PB1    excellent coverage
PA     good coverage
HA     excellent coverage
NP     good coverage
NA     very low coverage
M      moderate coverage
NS     poor coverage
```

A good whole-genome mapping percentage does not mean every influenza segment is sufficiently covered.

---

> **Important distinction**
>
> "No variant detected" is not equivalent to "strong evidence for the reference allele" if the site has little or no coverage.

---

# 17. Duplicate Handling

Tools such as:

```text
MarkDuplicates
```

identify reads that appear to represent the same original DNA fragment.

The rationale is:

```text
One original molecule
        ↓
       PCR
        ↓
Many copies
        ↓
Sequencing
        ↓
Artificially inflated read support
```

However, duplicate removal is **not automatically appropriate for every dataset**.

---

## Random-Fragment Libraries

Duplicate marking/removal may be reasonable when the library was generated by random fragmentation.

---

## Very Deep Sequencing

At extremely high depth, independent molecules may naturally share the same start and end positions.

A study by Zhou et al. showed that aggressive duplicate removal in ultra-deep sequencing can overcorrect read counts and bias allele-frequency estimates.

---

## Amplicon Sequencing

This is particularly important for influenza amplicon sequencing.

Amplicons have predefined boundaries.

Therefore many legitimate reads can naturally begin and end at the same coordinates.

Conceptually:

```text
Amplicon 1
|-----------------------|

Read A
|-----------------------|

Read B
|-----------------------|

Read C
|-----------------------|
```

These reads may not represent unwanted PCR duplication.

They may simply reflect the assay design.

Therefore:

```text
Random-fragment WGS
→ duplicate removal may be appropriate

Amplicon sequencing
→ do not blindly remove duplicates
```

Duplicate handling must be determined by the library protocol.

---

# 18. LoFreq Variant Calling

LoFreq is commonly used for sensitive viral variant calling.

The key analytical question is:

```text
Is the observed ALT allele statistically credible
given sequencing and alignment uncertainty?
```

rather than:

```text
Is the sample genotype 0/0, 0/1, or 1/1?
```

LoFreq is suitable for detecting viral variants that occur below majority frequency.

---

## LoFreq Viterbi Realignment

Local alignments around indels can be ambiguous.

Example:

```text
Reference: AAAACCCCCGGGG
Read:      AAAA-CCCCGGGG
```

A gap may be represented in multiple nearby positions.

Poor alignment around an indel can create artificial SNP-like mismatches.

LoFreq Viterbi attempts to improve local read alignment before variant calling.

Workflow:

```text
Mapped BAM
   ↓
LoFreq Viterbi
   ↓
Realigned BAM
```

---

## Indel Quality Processing

A typical LoFreq workflow includes indel-quality processing before variant calling.

Conceptually:

```text
BAM
 ↓
Viterbi realignment
 ↓
Indel-quality processing
 ↓
Variant calling
```

Even when the final analysis is restricted to SNVs, preprocessing should normally follow the caller's intended workflow.

---

## Variant Calling

At each genomic site, the caller evaluates evidence such as:

```text
Reference allele
ALT-supporting reads
Base quality
Mapping quality
Coverage
Sequencing error probability
Strand support
```

Example:

```text
Reference nucleotide: A

Observed:
A = 7,895
G = 1,932
T = 11
C = 4
```

The caller determines whether:

```text
A → G
```

is sufficiently supported to be emitted as a variant.

---

# 19. Variant Filtering

Raw variant calls should generally not be interpreted directly.

Important filtering dimensions include:

```text
QUAL
DP
AF
Base quality
Mapping quality
Strand evidence
Variant type
Local alignment context
```

---

## Allele Frequency Threshold

Example:

```text
AF ≥ 0.20
```

would retain:

```text
20%
47%
92%
```

and remove:

```text
4%
12%
19%
```

The AF threshold determines which viral subpopulations your study is capable of analyzing.

---

### Higher AF Thresholds

Example:

```text
AF ≥ 0.50
```

Primarily captures majority variants.

---

### Intermediate AF Thresholds

Example:

```text
AF ≥ 0.20
```

Includes substantial minority viral populations.

---

### Very Low AF Thresholds

Example:

```text
AF ≥ 0.01
```

moves into a range where:

```text
true biological variation
vs.
sequencing error
```

becomes increasingly difficult to resolve without specialized validation.

---

## Depth Threshold

Consider two variants:

```text
Variant A:
AF = 0.20
DP = 10
ALT ≈ 2 reads
```

versus:

```text
Variant B:
AF = 0.20
DP = 1000
ALT ≈ 200 reads
```

The AF is identical.

The evidential support is very different.

Therefore, DP and AF should be evaluated together.

---

## Strand Bias

A suspicious site might look like:

```text
Forward ALT reads = 57
Reverse ALT reads = 0
```

A more convincing site might show:

```text
Forward ALT reads = 31
Reverse ALT reads = 29
```

Strong directional imbalance can indicate a technical artifact.

However, strand bias must be interpreted in the context of the library design and local alignment structure.

---

# 20. VCF Interpretation

Variant callers commonly produce:

```text
VCF
```

Example:

```text
#CHROM   POS   REF   ALT   QUAL   FILTER   INFO
HA       521   A     G     126    PASS     DP=932;AF=0.28
```

Interpretation:

```text
Segment = HA
Position = 521
Reference = A
Alternate = G
Depth = 932
AF = 28%
```

---

## Important VCF Fields

### CHROM

For influenza:

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

---

### POS

Position within the segment.

---

### REF

Reference nucleotide.

---

### ALT

Alternate nucleotide.

---

### DP

Depth at the site.

---

### AF

Alternate allele frequency.

---

### QUAL

Variant-caller confidence score.

Do not confuse:

```text
Base quality
Mapping quality
Variant QUAL
```

These represent different uncertainty levels.

---

# 21. Functional Annotation

A raw VCF may tell you:

```text
HA
521
A → G
```

But biological interpretation requires:

```text
Which gene?
Which codon?
Which amino acid?
Synonymous?
Missense?
Stop gained?
Non-coding?
```

Tools commonly used:

```text
SnpEff
SnpSift
```

---

## SnpEff

SnpEff predicts the biological consequence of variants using genome annotation.

Example output:

```text
missense_variant
p.Glu172Ala
```

---

## SnpSift

SnpSift can extract selected fields from the annotated VCF into a tabular format.

Useful fields include:

```text
CHROM
POS
REF
ALT
QUAL
DP
AF
ANN.EFFECT
ANN.IMPACT
ANN.GENE
ANN.AA_POS
ANN.HGVS_C
ANN.HGVS_P
```

---

## Synonymous Variant

A nucleotide changes without altering the encoded amino acid.

Example:

```text
GAA → GAG
```

Both encode:

```text
Glu
```

Therefore:

```text
Nucleotide changed
Protein sequence unchanged
```

---

## Missense Variant

A nucleotide substitution alters the amino acid.

Example:

```text
GAA → GCA
```

which changes:

```text
Glu → Ala
```

A missense mutation does **not** automatically mean that the mutation is functionally important.

---

## Nonsense / Stop-Gained Variant

A nucleotide substitution generates a premature stop codon.

Example:

```text
TGG → TAG
```

Result:

```text
Trp → STOP
```

This may have a substantial predicted effect on the protein, but functional significance still requires biological interpretation.

---

# 22. Visual Validation

Important variants should be inspected directly against mapped reads using a genome browser such as:

```text
IGV
JBrowse
```

Inspect:

* Read depth
* REF/ALT balance
* Strand support
* Read-end clustering
* Soft clipping
* Nearby indels
* Local alignment quality
* Mapping ambiguity

---

## More Convincing Pattern

ALT evidence is:

```text
present in many reads
distributed across both orientations
not restricted to read ends
supported by good-quality alignments
```

---

## Suspicious Pattern

ALT evidence occurs:

```text
only at read ends
only on one strand
only in a few poor-quality reads
next to a large indel
inside a badly aligned region
```

---

# 23. Consensus Sequence Construction

Variant calling and consensus construction answer different questions.

Variant analysis asks:

> What nucleotide diversity exists in the sample?

Consensus analysis asks:

> Which single nucleotide should represent the sample at each position?

Example:

```text
Reference = A

A = 62%
G = 38%
```

Variant output:

```text
A → G
AF = 0.38
```

Consensus:

```text
A
```

Therefore:

```text
Variant table
```

and:

```text
Consensus genome
```

must both be retained for within-host viral evolution studies.

---

## Uncertain Positions

Consensus pipelines may represent positions with insufficient support as:

```text
N
```

Example:

```text
ACGTACNNNNGTAC
```

This is preferable to assigning a nucleotide where the data do not provide adequate evidence.

---

# 24. PCR and Amplicon Bias

Amplicon-based sequencing introduces additional considerations.

PCR can produce:

* Unequal template amplification
* Primer-dependent bias
* Primer mismatch effects
* Chimeric products
* Late-cycle artifacts
* Distorted template abundance

Kanagawa's review of multitemplate PCR describes how amplification bias and artifact formation can distort apparent sequence abundance.

Therefore:

```text
Observed VAF
```

should be interpreted as an estimate produced after:

```text
Sample
 ↓
RNA extraction
 ↓
cDNA synthesis
 ↓
PCR amplification
 ↓
Library preparation
 ↓
Sequencing
 ↓
Bioinformatic processing
```

It is not a perfectly unbiased measurement of the original virion population.

---

# 25. Tool Summary

| Tool                    | Input                             | Output              | Purpose                                      |
| ----------------------- | --------------------------------- | ------------------- | -------------------------------------------- |
| **FastQC**              | FASTQ                             | QC report           | Diagnose read quality                        |
| **fastp**               | FASTQ                             | Clean FASTQ         | Adapter and quality processing               |
| **BWA-MEM**             | FASTQ + FASTA                     | SAM/BAM             | Align reads to the influenza reference       |
| **Samtools sort**       | BAM                               | Sorted BAM          | Coordinate-sort alignments                   |
| **Samtools / Qualimap** | BAM                               | QC metrics          | Assess mapping and coverage                  |
| **MarkDuplicates**      | BAM                               | Marked/filtered BAM | Handle duplicate-like reads when appropriate |
| **LoFreq Viterbi**      | BAM + FASTA                       | Realigned BAM       | Improve local alignment                      |
| **LoFreq**              | BAM + FASTA                       | VCF                 | Sensitive SNV/indel calling                  |
| **SnpEff**              | VCF + annotation                  | Annotated VCF       | Predict gene/protein consequences            |
| **SnpSift**             | Annotated VCF                     | Table               | Extract useful annotation fields             |
| **IGV/JBrowse**         | BAM + VCF + FASTA                 | Visualization       | Inspect candidate variants                   |
| **Consensus tool**      | Reference + variant/read evidence | FASTA               | Build representative viral genome            |

---

# 26. QC Checkpoints

A pipeline should not be treated as:

```text
Input
 ↓
Run every tool
 ↓
Accept VCF
```

There should be decision points throughout the workflow.

---

## FASTQ Checkpoint

Ask:

```text
Are the reads usable?
Are adapters present?
Are quality scores acceptable?
Is sequence composition compatible with the protocol?
```

---

## Mapping Checkpoint

Ask:

```text
What percentage mapped?
Are reads properly paired?
Did every influenza segment map?
```

---

## Coverage Checkpoint

Ask:

```text
Is coverage sufficient across every segment?
Are there zero-coverage regions?
Are there extreme coverage spikes?
```

---

## Variant Checkpoint

Ask:

```text
Is DP sufficient?
Is AF credible?
Is QUAL acceptable?
Is support present on both orientations?
Is the site next to an alignment artifact?
```

---

## Annotation Checkpoint

Ask:

```text
Do FASTA and GFF segment names match?
Does the variant fall in the expected gene?
Does the amino-acid position make biological sense?
```

---

## Biological Checkpoint

Ask:

```text
Minority or majority?
Synonymous or nonsynonymous?
Persistent or transient?
Shared or unique?
Present at multiple time points?
Reference-dependent?
```

---

# 27. Recommended Final Variant Table

A useful final table should contain at least:

| Field      | Description                                        |
| ---------- | -------------------------------------------------- |
| Sample     | Sample identifier                                  |
| Segment    | PB2, PB1, PA, HA, NP, NA, M, or NS                 |
| Position   | Genomic coordinate                                 |
| REF        | Reference nucleotide                               |
| ALT        | Alternate nucleotide                               |
| DP         | Read depth                                         |
| ALT count  | Number of alternate-supporting reads, if available |
| AF         | Alternate allele frequency                         |
| QUAL       | Variant confidence                                 |
| Filter     | PASS/failure information                           |
| Gene       | Gene name                                          |
| Effect     | Synonymous, missense, etc.                         |
| AA change  | Predicted protein-level change                     |
| Time point | Experimental sampling point                        |
| Group      | Experimental treatment group                       |
| Validation | Manual QC note if needed                           |

Example:

```text
Sample  Segment  Position  REF  ALT  DP    AF    Effect       AA_change
B01     HA       540       G    A    1250  0.37  missense     p.Asp173Asn
B02     PB2      1909      A    G    2130  0.74  missense     p.Xxx###Yyy
B03     NP       812       C    T    684   0.33  synonymous   -
```

---

# 28. Worked Interpretation Example

Suppose the final variant record is:

```text
Sample: B01
Segment: HA
Position: 540
REF: G
ALT: A
DP: 1250
AF: 0.37
Effect: missense
AA change: p.Asp173Asn
```

A defensible interpretation is:

> Sequencing reads from B01 contain an HA SNV at position 540 relative to the selected reference. The alternate A allele represents approximately 37% of qualifying sequencing evidence at a read depth of 1,250. Functional annotation predicts an Asp-to-Asn amino-acid substitution.

This result alone does **not** establish:

```text
Adaptive evolution
Immune escape
Increased fitness
Increased transmission
Causality
```

Those require additional statistical, experimental, structural, or literature evidence.

---

# 29. Practical Galaxy Workflow

For paired-end Illumina influenza SNV analysis:

```text
1. Upload FASTQ R1/R2
      ↓
2. Create paired collection
      ↓
3. FastQC — raw reads
      ↓
4. fastp
      ↓
5. FastQC — processed reads
      ↓
6. Prepare 8-segment influenza FASTA
      ↓
7. BWA-MEM
      ↓
8. Samtools sort
      ↓
9. Mapping and coverage QC
      ↓
10. Appropriate duplicate handling
      ↓
11. LoFreq Viterbi
      ↓
12. LoFreq quality preprocessing
      ↓
13. LoFreq variant calling
      ↓
14. Restrict to SNVs if required
      ↓
15. Filter by:
      - QUAL
      - DP
      - AF
      - artifact criteria
      ↓
16. SnpEff annotation
      ↓
17. SnpSift field extraction
      ↓
18. IGV/JBrowse inspection
      ↓
19. Final SNV table
      ↓
20. Consensus genome generation
```

---

# Practical Checklist

* [ ] Create a dedicated Galaxy history.
* [ ] Upload paired R1/R2 FASTQ files.
* [ ] Build a paired collection.
* [ ] Run raw FastQC.
* [ ] Process reads with fastp.
* [ ] Run post-trim FastQC.
* [ ] Prepare the eight-segment reference FASTA.
* [ ] Verify FASTA/GFF segment IDs.
* [ ] Map reads with BWA-MEM.
* [ ] Sort BAM files.
* [ ] Inspect mapping rate.
* [ ] Inspect properly paired reads.
* [ ] Calculate coverage per segment.
* [ ] Identify zero-coverage positions.
* [ ] Decide duplicate strategy based on library design.
* [ ] Run LoFreq Viterbi.
* [ ] Perform required LoFreq quality processing.
* [ ] Call variants with LoFreq.
* [ ] Restrict to SNVs if required.
* [ ] Apply justified AF, DP, and QUAL filters.
* [ ] Annotate variants with SnpEff.
* [ ] Extract useful fields with SnpSift.
* [ ] Inspect important calls in IGV/JBrowse.
* [ ] Export a final variant table.
* [ ] Construct consensus sequences separately.
* [ ] Retain both consensus and within-host AF information.

---

# Key Conceptual Summary

Three distinct biological layers should always be kept separate:

## 1. Reads

```text
Direct sequencing observations
```

## 2. Variants

```text
Estimated nucleotide differences and frequencies relative to a reference
```

## 3. Consensus Genome

```text
One representative nucleotide per genomic coordinate
```

Example:

```text
HA position 800

A reads = 620
G reads = 380
```

Variant analysis:

```text
A → G
AF = 0.38
```

Consensus analysis:

```text
A
```

If only the consensus is retained, the 38% minority population disappears.

For within-host influenza evolution studies, **variant-frequency information should therefore be retained alongside consensus genomes**.

---

# References

1. Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM.
   **The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants.**
   *Nucleic Acids Research.* 2010;38(6):1767-1771.
   DOI: `10.1093/nar/gkp1137`

2. Blankenberg D, Gordon A, Von Kuster G, et al.
   **Manipulation of FASTQ data with Galaxy.**
   *Bioinformatics.* 2010;26(14):1783-1785.
   DOI: `10.1093/bioinformatics/btq281`

3. Galaxy Training Network.
   **Microbial Variant Calling.**
   Updated January 23, 2026.
   [https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/microbial-variants/tutorial.html](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/microbial-variants/tutorial.html)

4. Galaxy Training Network.
   **Calling variants in non-diploid systems.**
   [https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/non-dip/tutorial.html](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/non-dip/tutorial.html)

5. Galaxy Training Network.
   **Avian influenza viral strain analysis from gene segment sequencing data.**
   [https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/aiv-analysis/tutorial.html](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/aiv-analysis/tutorial.html)

6. Galaxy Training Network.
   **Calling very rare variants.**
   [https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/dunovo/tutorial.html](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/dunovo/tutorial.html)

7. Galaxy Training Network.
   **NGS data logistics.**
   [https://training.galaxyproject.org/training-material/topics/introduction/tutorials/galaxy-intro-ngs-data-managment/tutorial.html](https://training.galaxyproject.org/training-material/topics/introduction/tutorials/galaxy-intro-ngs-data-managment/tutorial.html)

8. Galaxy Training Network.
   **From NCBI's Sequence Read Archive to Galaxy: SARS-CoV-2 variant analysis.**
   [https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/sars-cov-2/tutorial.html](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/sars-cov-2/tutorial.html)

9. Zhou W, Chen T, Zhao H, et al.
   **Bias from removing read duplication in ultra-deep sequencing experiments.**
   *Bioinformatics.* 2014;30(8):1073-1080.
   DOI: `10.1093/bioinformatics/btt771`

10. Galaxy Training Network.
    **Mutation calling, viral genome reconstruction and lineage/clade assignment from SARS-CoV-2 sequencing data.**
    [https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/sars-cov-2-variant-discovery/tutorial.html](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/sars-cov-2-variant-discovery/tutorial.html)

11. Kanagawa T.
    **Bias and Artifacts in Multitemplate Polymerase Chain Reactions (PCR).**
    *Journal of Bioscience and Bioengineering.* 2003;96(4):317-323.

12. Nielsen R, Paul JS, Albrechtsen A, Song YS.
    **Genotype and SNP calling from next-generation sequencing data.**
    *Nature Reviews Genetics.* 2011;12(6):443-451.
    DOI: `10.1038/nrg2986`

13. Seemann T.
    **Snippy: Rapid haploid variant calling and core genome alignment.**
    [https://github.com/tseemann/snippy](https://github.com/tseemann/snippy)

---

# Scope

This guide is intended as a conceptual and workflow reference.

Parameters such as:

```text
Minimum AF
Minimum DP
Minimum QUAL
Duplicate handling
Consensus threshold
Mapping-quality cutoff
```

should not be treated as universal constants.

They should be chosen based on:

* Sequencing platform
* Library preparation method
* Amplicon vs random-fragment sequencing
* Expected sequencing error
* Read depth
* Reference choice
* Experimental design
* Biological question
* Desired sensitivity for minority variants
