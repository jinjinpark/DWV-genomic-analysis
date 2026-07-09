# DWV Genomic Analysis Scripts

Python scripts used in:

**"Whole-genome based Phylogenetic Analysis of Deformed Wing Virus from Korea and Vietnam"**  
Jinyoung Park, Bui Thi Thuy Duong, Pham Thi Lanh, Seil Kim, Dong Van Quyen

---

## Scripts

### Genome Assembly

| Script | Description |
|--------|-------------|
| `fill_gap_sequences.py` | Fills gap regions in ONT consensus sequences using Sanger sequencing amplicons aligned to the consensus |

### Ancestral Sequence Reconstruction

| Script | Description |
|--------|-------------|
| `ancestral_sequence_extractor.py` | Extracts ancestral sequences from RAxML-NG `.ancestralStates` output for specified nodes |
| `resolve_ambiguous_nodes.py` | Resolves IUPAC ambiguous bases by selecting the highest posterior probability base from `.ancestralProbs` |

### Variant Analysis

| Script | Description |
|--------|-------------|
| `fasta_to_vcf.py` | Converts multi-FASTA alignment to VCF format using a specified reference sequence |
| `convert_protein_to_dna_gff.py` | Converts protein-coordinate GFF3 annotations to DNA coordinates for SnpEff database construction |
| `analyze_variants_by_sample.py` | Counts synonymous and missense variants per sample from annotated VCF |
| `analyze_all_variants_by_gene.py` | Summarizes all variant types (synonymous, missense) by gene and sample |
| `analyze_signature_variants.py` | Identifies lineage-unique and rare variants; calculates unique variant proportions per sample |
| `analyze_variant_groups.py` | Analyzes variant-sharing patterns between samples; computes pairwise Jaccard similarity |

### Selection Analysis

| Script | Description |
|--------|-------------|
| `dnds_analysis.py` | Computes pairwise dN/dS (Nei-Gojobori, Jukes-Cantor corrected) for each RAS lineage versus the MRCA reference, per protein-coding region |

### Primer Analysis

| Script | Description |
|--------|-------------|
| `primer_coverage_analyzer.py` | Evaluates in silico primer coverage against target genome sequences, accounting for degenerate bases and mismatches |

---

## Pipeline Overview

```
ONT consensus sequences
        ↓
fill_gap_sequences.py  (Sanger gap filling)
        ↓
Complete genome sequences
        ↓
RAxML-NG ancestralStates/ancestralProbs
        ↓
ancestral_sequence_extractor.py + resolve_ambiguous_nodes.py
        ↓
RAS_nucleotide_sequence.fas  +  convert_protein_to_dna_gff.py (SnpEff DB)
        ↓
fasta_to_vcf.py → variants.vcf → SnpEff annotation → variants_annotated.vcf
        ↓
analyze_variants_by_sample.py
analyze_all_variants_by_gene.py
analyze_signature_variants.py
analyze_variant_groups.py

RAS_nucleotide_sequence.fas
        ↓
dnds_analysis.py  (pairwise dN/dS per protein, lineage vs MRCA)
```

---

## Requirements

```
biopython
pandas
```

Install with:
```bash
pip install biopython pandas
```

---

## Usage

See individual script docstrings and `--help` flags for usage details.

```bash
# Gap filling
python fill_gap_sequences.py

# Ancestral sequence extraction
python ancestral_sequence_extractor.py --all
python resolve_ambiguous_nodes.py -s *.ancestralStates -p *.ancestralProbs --all

# Variant analysis
python fasta_to_vcf.py
python analyze_variants_by_sample.py
python analyze_all_variants_by_gene.py
python analyze_signature_variants.py
python analyze_variant_groups.py

# Selection analysis (dN/dS)
python dnds_analysis.py --fasta RAS_nucleotide_sequence.fas --reference Node167 --outdir results

# Primer coverage
python primer_coverage_analyzer.py --help
```
