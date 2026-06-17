#!/usr/bin/env python3
"""
Analyze ALL variants (including synonymous) by gene and sample
"""

import re
from collections import defaultdict

def parse_vcf_all_variants(vcf_file):
    """Parse VCF and extract all variants per gene"""

    samples = []
    gene_variants = defaultdict(lambda: defaultdict(list))

    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                continue

            if line.startswith('#CHROM'):
                parts = line.strip().split('\t')
                samples = parts[9:]
                continue

            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue

            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            info = parts[7]
            genotypes = parts[9:]

            # Parse ANN field
            ann_match = re.search(r'ANN=([^;]+)', info)
            if not ann_match:
                continue

            ann_field = ann_match.group(1)
            first_ann = ann_field.split(',')[0]
            ann_parts = first_ann.split('|')

            if len(ann_parts) < 4:
                continue

            alt_allele = ann_parts[0]
            effect = ann_parts[1]
            gene = ann_parts[3]
            dna_change = ann_parts[9] if len(ann_parts) > 9 else ''
            protein_change = ann_parts[10] if len(ann_parts) > 10 else ''

            # Only interested in CDS variants
            if gene.startswith('PRO_'):
                # Check each sample
                for i, (sample, gt) in enumerate(zip(samples, genotypes)):
                    gt_value = gt.split(':')[0] if ':' in gt else gt

                    if gt_value not in ['0', '.', './.']:
                        variant_type = 'synonymous' if 'synonymous' in effect else \
                                      'missense' if 'missense' in effect else 'other'

                        gene_variants[gene][sample].append({
                            'pos': pos,
                            'ref': ref,
                            'alt': alt_allele,
                            'type': variant_type,
                            'dna': dna_change,
                            'protein': protein_change
                        })

    return samples, gene_variants

# Parse VCF
vcf_file = 'variants_annotated.vcf'
samples, gene_variants = parse_vcf_all_variants(vcf_file)

# Gene names mapping
gene_names = {
    'PRO_0000460700': 'Genome polyprotein',
    'PRO_0000460701': 'Leader protein',
    'PRO_0000460702': 'Capsid VP2',
    'PRO_0000460703': 'Capsid VP4',
    'PRO_0000460704': 'Capsid VP3',
    'PRO_0000460705': 'Capsid VP1',
    'PRO_0000460706': 'Protein 2A-like',
    'PRO_0000460707': 'Helicase',
    'PRO_0000460708': '3C-like protease',
    'PRO_0000460709': 'RNA polymerase'
}

# Write report
output_file = 'all_variants_by_gene.txt'
with open(output_file, 'w') as out:
    out.write("=" * 100 + "\n")
    out.write("ALL VARIANTS BY GENE AND SAMPLE (INCLUDING SYNONYMOUS)\n")
    out.write("=" * 100 + "\n\n")

    for gene in sorted(gene_variants.keys()):
        gene_name = gene_names.get(gene, gene)

        out.write(f"\n{'=' * 100}\n")
        out.write(f"GENE: {gene} - {gene_name}\n")
        out.write(f"{'=' * 100}\n\n")

        # Count variants per sample
        out.write(f"{'Sample':<20} {'Total':<10} {'Synonymous':<15} {'Missense':<15}\n")
        out.write("-" * 100 + "\n")

        for sample in samples:
            variants = gene_variants[gene][sample]
            syn_count = sum(1 for v in variants if v['type'] == 'synonymous')
            mis_count = sum(1 for v in variants if v['type'] == 'missense')

            if len(variants) > 0:
                out.write(f"{sample:<20} {len(variants):<10} {syn_count:<15} {mis_count:<15}\n")

        out.write("\n")

print(f"Report saved to: {output_file}")
