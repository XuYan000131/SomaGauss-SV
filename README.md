# SomaGauss-SV: Detection of Somatic Structural Variations

## Overview

SomaGauss-SV is a pipeline designed to identify somatic structural variations (SVs) with high confidence, particularly focusing on large-scale insertions (INSs) and deletions (DELs). It integrates a series of tools and custom scripts to detect, refine, and validate somatic SVs, improving the accuracy and reliability of somatic SV identification.

## Pipeline Workflow

The pipeline consists of the following key steps:

### Step 1: Merge BAM Files and Detect Structural Variations
Merge sorted tag BAM files of paired samples using Samtools.
Use Sniffles2 to detect structural variations from the merged BAM file.
Filter SVs based on read alignment quality, supporting read count, SV length, and minimum alignment length.

### Step 2: Screening for Candidate Somatic Structural Variations
Preliminary filter SVs reported by Sniffles2, retaining only those supported by reads with the "tumor" tag.
For DEL types, further refine the screening process:
- Amplify regions around candidate DELs and generate small BAM files.
- Detect DELs in these regions using pysam and apply a mixed Gaussian distribution model.
- Resolve overfitting issues using a student.test method to compare DEL lengths between blood and tumor samples.

For INS types, adopt a screening strategy using Straglr for re-detection and preliminary screening.

### Step 3: Normal Population Structural Variation Filtering
Integrate VCF files using Jasmine software, comparing somatic SVs with those in the Chinese normal population.
Remove structural variations identical to those in the normal population, retaining the rest as somatic structural variations.

## Installation

### Prerequisites
Ensure the following tools are installed:

- Samtools (version 1.9)
- Sniffles2 (version 2.0.6)
- Pysam
- Straglr (version 1.4.1)
- Jasmine (version 1.1.5)

## Usage

### Input Data
- Tumor sample BAM files.
- Paired blood sample BAM files.

### Output
- `{Sample}_GaussianMixture_Judgement.bed`: This file contains the filtered structural variations detected in the pipeline.
- `somatic.{Sample}.reslut.AF.PON.vcf`: This file contains the final list of somatic structural variations after PON steps.

## Directory Structure

```
SomaGauss-SV/
│
├── 01.Base.Merge.Seek.Candidate.Somatic.SV/   # Merging BAM files and seeking candidate somatic SVs
│   ├── ONT_clean_fq_minmap2.2.26_blood.sh
│   ├── ONT_clean_fq_minmap2.2.26_tumor.sh
│   ├── config.py
│   ├── function.py
│   └── pipeline.py
│
├── 02.Straglr.Seek.Somatic.INS/             # Straglr for seeking somatic INS
│   ├── config.py
│   ├── function.py
│   └── pipeline.py
│
├── 03.Gaussian.Mixture.Seek.Somatic.DEL/     # Gaussian mixture model for seeking somatic DEL
│   ├── config.py
│   ├── function.py
│   └── pipeline.py
│
└── 04.PON.and.Bulit.New.VCF/                # PON and building new VCF
│   ├── config.py
│   ├── function.py
│   └── pipeline.py
```
