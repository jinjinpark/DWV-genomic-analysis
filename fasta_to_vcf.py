#!/usr/bin/env python3
"""
Generate VCF file from multi-FASTA alignment
Compares all sequences to a reference sequence
"""

from Bio import SeqIO
from datetime import datetime

def fasta_to_vcf(fasta_file, reference_name, output_vcf):
    """Convert aligned FASTA to VCF format"""

    # Read all sequences
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_name = record.id
        sequences[seq_name] = str(record.seq)

    if reference_name not in sequences:
        raise ValueError(f"Reference '{reference_name}' not found in FASTA file")

    ref_seq = sequences[reference_name]
    ref_len = len(ref_seq)

    # Get all sample names (excluding reference)
    samples = [name for name in sequences.keys() if name != reference_name]

    # Open VCF file for writing
    with open(output_vcf, 'w') as vcf:
        # Write VCF header
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
        vcf.write(f"##source=fasta_to_vcf.py\n")
        vcf.write(f"##reference={reference_name}\n")
        vcf.write(f"##contig=<ID={reference_name},length={ref_len}>\n")
        vcf.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n')
        vcf.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
        vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')

        # Write column headers
        vcf.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{chr(9).join(samples)}\n")

        # Find variants
        variants_found = 0
        for pos in range(ref_len):
            ref_base = ref_seq[pos]

            # Skip if reference is gap
            if ref_base in ['-', 'N', 'n']:
                continue

            # Collect alleles at this position
            alleles = {}
            for sample in samples:
                if pos < len(sequences[sample]):
                    alt_base = sequences[sample][pos]
                    if alt_base != ref_base and alt_base not in ['-', 'N', 'n']:
                        alleles[alt_base] = alleles.get(alt_base, 0) + 1

            # If there are variants, write VCF line
            if alleles:
                alt_bases = sorted(alleles.keys())
                alt_str = ','.join(alt_bases)

                # Calculate allele frequency
                total_alleles = len(samples)
                alt_count = sum(alleles.values())
                af = alt_count / total_alleles

                # Write VCF record
                chrom = reference_name
                pos_1based = pos + 1
                vcf_id = '.'
                ref = ref_base
                alt = alt_str
                qual = '.'
                filter_val = 'PASS'
                info = f"DP={total_alleles};AF={af:.4f}"
                format_str = 'GT'

                # Genotypes
                genotypes = []
                for sample in samples:
                    if pos < len(sequences[sample]):
                        sample_base = sequences[sample][pos]
                        if sample_base == ref_base:
                            gt = '0'  # Reference
                        elif sample_base in alt_bases:
                            gt = str(alt_bases.index(sample_base) + 1)  # Alt allele
                        else:
                            gt = '.'  # Missing
                    else:
                        gt = '.'
                    genotypes.append(gt)

                genotype_str = '\t'.join(genotypes)

                vcf.write(f"{chrom}\t{pos_1based}\t{vcf_id}\t{ref}\t{alt}\t{qual}\t{filter_val}\t{info}\t{format_str}\t{genotype_str}\n")
                variants_found += 1

    print(f"VCF file created: {output_vcf}")
    print(f"Reference: {reference_name}")
    print(f"Samples: {len(samples)}")
    print(f"Variants found: {variants_found}")

# Run the conversion
fasta_file = 'RAS_nucleotide_sequence.fas'
reference_name = 'Node167'
output_vcf = 'variants.vcf'

fasta_to_vcf(fasta_file, reference_name, output_vcf)
