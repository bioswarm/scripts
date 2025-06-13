# Genome Data Processing and Protein Sequence Analysis Tools

This repository contains a collection of tools for genome data processing and protein sequence analysis. These tools are primarily used for processing GFF files, protein sequences, CDS sequences, and performing sequence alignment analysis.

## Tool List

### 1. GFF File Processing Tools

- **miniprot_GFF_2_EVM_GFF3.py**: Convert miniprot GFF output to EVM specification format
- **gff_unique_for_flybase.pl**: Process GFF files from FlyBase database to generate unique protein and CDS sequences

### 2. Protein Sequence Generation Tools

- **generate_protein_from_refseq1.pl**: Generate protein sequences from RefSeq data (Version 1)
- **generate_protein_from_refseq2.pl**: Generate protein sequences from RefSeq data (Version 2)
- **generate_protein_from_genbank.pl**: Generate protein sequences from GenBank data
- **generate_protein_from_egapx.pl**: Generate protein sequences from EGAPX data

### 3. Sequence Position Conversion Tool

- **convert_position.py**: Convert positions between aligned sequences and original sequences
  - Convert alignment positions to original sequence positions
  - Convert original sequence positions to alignment positions

### 4. Protein Position Importance Analysis Tool

- **calculate_position_importance.py**: Analyze the importance of positions in protein sequences
  - Calculate information gain
  - Calculate chi-square test scores
  - Calculate perfect discrimination scores
  - Calculate biochemical property difference scores
  - Generate visualization results

## Usage Instructions

### 1. GFF File Processing

```bash
# Convert miniprot GFF to EVM format
python miniprot_GFF_2_EVM_GFF3.py <miniprot_GFF_file>

# Process FlyBase GFF files
perl gff_unique_for_flybase.pl
```

### 2. Protein Sequence Generation

```bash
# Generate protein sequences from RefSeq data
perl generate_protein_from_refseq1.pl <output_prefix>
perl generate_protein_from_refseq2.pl <output_prefix>

# Generate protein sequences from GenBank data
perl generate_protein_from_genbank.pl <output_prefix>

# Generate protein sequences from EGAPX data
perl generate_protein_from_egapx.pl <output_prefix>
```

### 3. Sequence Position Conversion

```bash
# Convert alignment positions to original sequence positions
python convert_position.py <fasta_file> <gene_name>

# Convert original sequence positions to alignment positions
python convert_position.py <fasta_file> <gene_name> -v
```

### 4. Protein Position Importance Analysis

```bash
python calculate_position_importance.py
```

## Dependencies

- Python 3.x
- Perl 5.x
- Required Python packages:
  - numpy
  - pandas
  - scikit-learn
  - scipy
  - biopython
  - matplotlib
  - seaborn

## Notes

1. Ensure correct input file formats when using protein sequence generation tools
2. Position conversion tool requires properly formatted FASTA input files
3. Protein position importance analysis tool requires label files and alignment files
