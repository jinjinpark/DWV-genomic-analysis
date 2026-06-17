#!/usr/bin/env python3
"""
Analyze variants by sample from annotated VCF
"""

import re
from collections import defaultdict

def parse_vcf_by_sample(vcf_file):
    """Parse VCF and extract variants per sample"""

    samples = []
    sample_variants = defaultdict(lambda: {
        'total': 0,
        'synonymous': 0,
        'missense': 0,
        'positions': [],
        'missense_details': []
    })

    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                continue

            if line.startswith('#CHROM'):
                # Extract sample names
                parts = line.strip().split('\t')
                samples = parts[9:]  # Sample names start from column 10
                continue

            # Parse variant line
            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue

            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            info = parts[7]

            # Extract genotypes (starting from column 9)
            genotypes = parts[9:]

            # Parse ANN field for variant effects
            ann_match = re.search(r'ANN=([^;]+)', info)
            if not ann_match:
                continue

            ann_field = ann_match.group(1)

            # Get first annotation (most severe)
            first_ann = ann_field.split(',')[0]
            ann_parts = first_ann.split('|')

            if len(ann_parts) < 4:
                continue

            alt_allele = ann_parts[0]
            effect = ann_parts[1]
            gene = ann_parts[3]

            # DNA and protein changes
            dna_change = ann_parts[9] if len(ann_parts) > 9 else ''
            protein_change = ann_parts[10] if len(ann_parts) > 10 else ''

            # Check each sample
            for i, (sample, gt) in enumerate(zip(samples, genotypes)):
                # Genotype is first field (GT)
                gt_value = gt.split(':')[0] if ':' in gt else gt

                # If genotype is not reference (0) and not missing (.)
                if gt_value not in ['0', '.', './.']:
                    sample_variants[sample]['total'] += 1
                    sample_variants[sample]['positions'].append(pos)

                    # Count variant types
                    if 'synonymous' in effect:
                        sample_variants[sample]['synonymous'] += 1
                    elif 'missense' in effect:
                        sample_variants[sample]['missense'] += 1
                        sample_variants[sample]['missense_details'].append({
                            'pos': pos,
                            'ref': ref,
                            'alt': alt_allele,
                            'dna': dna_change,
                            'protein': protein_change,
                            'gene': gene
                        })

    return samples, sample_variants

# Parse VCF
vcf_file = 'variants_annotated.vcf'
samples, sample_variants = parse_vcf_by_sample(vcf_file)

# Write summary report
output_file = 'variants_by_sample.txt'
with open(output_file, 'w') as out:
    out.write("=" * 80 + "\n")
    out.write("VARIANT ANALYSIS BY SAMPLE (Node167 as Reference)\n")
    out.write("=" * 80 + "\n\n")

    # Summary table
    out.write("SUMMARY TABLE:\n")
    out.write("-" * 80 + "\n")
    out.write(f"{'Sample':<20} {'Total Variants':<15} {'Synonymous':<15} {'Missense':<15}\n")
    out.write("-" * 80 + "\n")

    for sample in samples:
        stats = sample_variants[sample]
        out.write(f"{sample:<20} {stats['total']:<15} {stats['synonymous']:<15} {stats['missense']:<15}\n")

    out.write("\n" + "=" * 80 + "\n\n")

    # Detailed missense variants for each sample
    for sample in samples:
        stats = sample_variants[sample]

        out.write(f"\n{'=' * 80}\n")
        out.write(f"SAMPLE: {sample}\n")
        out.write(f"{'=' * 80}\n")
        out.write(f"Total variants: {stats['total']}\n")
        out.write(f"Synonymous: {stats['synonymous']}\n")
        out.write(f"Missense: {stats['missense']}\n")

        if stats['missense'] > 0:
            out.write(f"\nMISSENSE VARIANTS (n={stats['missense']}):\n")
            out.write("-" * 80 + "\n")
            out.write(f"{'Pos':<8} {'Ref>Alt':<10} {'DNA Change':<15} {'Protein Change':<20} {'Gene':<20}\n")
            out.write("-" * 80 + "\n")

            for variant in stats['missense_details']:
                ref_alt = f"{variant['ref']}>{variant['alt']}"
                out.write(f"{variant['pos']:<8} {ref_alt:<10} {variant['dna']:<15} {variant['protein']:<20} {variant['gene']:<20}\n")

        out.write("\n")

print(f"Analysis complete!")
print(f"Report saved to: {output_file}")
print(f"\nSummary:")
for sample in samples:
    stats = sample_variants[sample]
    print(f"{sample}: {stats['total']} variants ({stats['missense']} missense)")
