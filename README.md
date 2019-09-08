# 🦉 ddOWL: Detection of DNM Origin With Long Reads
This is a collection of scripts used for phasing putative pathogenic de novo mutations (DNMs) to the parental allele of origin.

## General idea

- Sequence a region of a few kb around the mutation with long read sequencing technologies (eg. Oxford Nanopore)
- Identify common variants in the vicinity that happen to differ between father and mother
- Find reads that overlap both the mutation and the common variant
- Count how often which alleles are seen together on the same read to phase the mutation

The script requires a list of known variants to test. This can be obtained from a database such as dbSNP or by directly calling variants in the sequencing data.

ddOWL is usually run with a full trio, in which case uninformative common variants can be filtered out automatically and alleles will be labelled as paternal/maternal.
However, it can also be run with the proband sample only, in which case the phase of the mutation with all listed variants will be output.

# phaser.py
Python script to find informative variants and calculate counts. Takes a sample table, variant tables and BAM files and generates allele counts in CSV format. Depending on the amount of samples and data this may take a few hours to run.

Requires Python 3 with pandas + pysam.

## Input
To run the phaser script you need to provide three parameters:

### Path to samples table CSV
This comma-separated file specifies the families that should be analysed. Each family should consist of mother, father and proband. If only proband data is available the mother and father can be left out, but then the origin of SNP alleles cannot be determined automatically.

Example: `examples/samples.csv`

#### Columns:
- BC (sample barcode, must match barcode in BAM filename)
- FamilyID (must be the same for all members of a family)
- relationship (mother/father/proband)

### Filename mask for variant TSV files
For each family there must be a tab-separated TSV file that provides the location and alleles of the variants/SNPs in the sequenced region. File names need to contain the family ID but must otherwise be exactly the same. The filename mask has to contain the characters `%s`, which will be replaced by the family ID (column FamilyID in sample table).

Example: `examples/snps_%s.tsv`

In each of these files there should also be exactly one variant that represents the mutation in that family.
For this entry, the `ID` column has to say "mutation", the `REF` allele is the wild-type allele and the `ALT` allele is the mutation allele.

(untested:) Deletions can be specified like in VCF format - as a variant at the preceding position where REF is the preceding base plus the deleted sequence and ALT is just the preceding base.
Eg. a "CT" deletion at position 10000029-10000030 preceded by an A at 10000028 would be specified as a variant at position 10000028 with alleles ref="ACT" and alt="A".

#### Columns:
- `#CHROM` (chromosome, must match naming scheme in BAM)
- `POS` (one-based coordinates)
- `ID` (any ID)
- `REF` (ref allele)
- `ALT` (alt allele)

### Filename mask for BAM files
For each individual there must be a sorted BAM file with their aligned reads (for example generated with minimap2).
The filename mask has to contain the characters `%s`, which will be replaced by the sample barcode (column BC in sample table).

Example: `/insert/path/to/bam/files/here/%s.sorted.bam`

### Additional parameters
Run `python phaser.py --help` to see additional parameters.

## Output
Several CSV files starting with the sample table filename (and in same directory):
Informative SNPs, allele counts, read stats, phase evidence, etc

## Example command

    module purge && module load python3-cbrg/201703
    python phaser.py \
        'examples/samples.csv' \
        'examples/snps_%s.tsv' \
        '/insert/path/to/bam/files/here/%s.sorted.bam' \
    ;

# phaser.r
R script to load CSV files from phaser.py and make heatmap plots for each SNP etc.

## Input
To run the phaser script you need to provide two parameters - the sample CSV and the variant TSV filename mask, as described above.
Run `Rscript phaser.r --help` to see additional parameters.

## Output
PDF files with various plots (such as the heatmaps) and CSV files with summary data. All output file names will start with sample table filename (and in same directory).

## Example command

    module purge && module load R/3.3.1-newgcc
    Rscript phaser.r \
        'examples/samples.csv' \
        'examples/snps_%s.tsv' \
    ;

# phaser_call.r
R script to load CSV files from phaser.py and do advanced statistical analysis (combining data from multiple SNPs, p-values, etc)

## Input
To run the phaser script you need to provide two parameters - the sample CSV and the variant TSV filename mask, as described above.
Run `Rscript phaser_call.r --help` to see additional parameters.

## Output
PDF file with phasing plots, named `SAMPLES_PATH.phasing.pdf`

## Example command

    module purge && module load R/3.3.1-newgcc
    Rscript phaser_call.r \
        'examples/samples.csv' \
        'examples/snps_%s.tsv' \
    ;

# Citation and License
Licensed under the Apache License, version 2.0.
Copyright 2019 Nils Koelling.
When you use ddOWL, please cite it in your work. For example:

> Koelling, Nils. "ddOWL: Detection of DNM Origin With Long Reads". https://github.com/koelling/ddowl/
