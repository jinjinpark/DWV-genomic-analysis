#!/usr/bin/env python3
"""
Identify signature variants for each sample
특히 참조서열과 거리가 먼 샘플의 특징적 변이 분석
"""

import re
from collections import defaultdict

def parse_vcf_for_signatures(vcf_file):
    """Parse VCF and identify sample-specific variants"""

    samples = []
    variant_data = []

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

            # Determine variant type
            if 'synonymous' in effect:
                variant_type = 'synonymous'
            elif 'missense' in effect:
                variant_type = 'missense'
            else:
                variant_type = 'other'

            # Check which samples have this variant
            samples_with_variant = []
            for i, (sample, gt) in enumerate(zip(samples, genotypes)):
                gt_value = gt.split(':')[0] if ':' in gt else gt
                if gt_value not in ['0', '.', './.']:
                    samples_with_variant.append(sample)

            variant_data.append({
                'pos': pos,
                'ref': ref,
                'alt': alt_allele,
                'type': variant_type,
                'gene': gene,
                'dna': dna_change,
                'protein': protein_change,
                'samples': samples_with_variant,
                'count': len(samples_with_variant)
            })

    return samples, variant_data

def analyze_signatures(samples, variant_data):
    """Identify unique and group-specific variants"""

    # 1. Unique variants (found in only one sample)
    unique_variants = defaultdict(list)

    # 2. Rare variants (found in 1-2 samples)
    rare_variants = defaultdict(list)

    # 3. Group-specific variants
    group_variants = defaultdict(lambda: defaultdict(list))

    for var in variant_data:
        sample_count = var['count']

        if sample_count == 1:
            sample = var['samples'][0]
            unique_variants[sample].append(var)

        if sample_count <= 2:
            for sample in var['samples']:
                rare_variants[sample].append(var)

        # Group by sample combination
        sample_key = tuple(sorted(var['samples']))
        for sample in var['samples']:
            group_variants[sample][sample_key].append(var)

    return unique_variants, rare_variants, group_variants

# Parse VCF
print("Parsing VCF file...")
vcf_file = 'variants_annotated.vcf'
samples, variant_data = parse_vcf_for_signatures(vcf_file)

print("Analyzing signature variants...")
unique_variants, rare_variants, group_variants = analyze_signatures(samples, variant_data)

# Gene names
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

# Calculate total variants per sample
sample_total = defaultdict(int)
sample_missense = defaultdict(int)
for var in variant_data:
    for sample in var['samples']:
        sample_total[sample] += 1
        if var['type'] == 'missense':
            sample_missense[sample] += 1

# Write report
output_file = 'signature_variants.txt'
with open(output_file, 'w') as out:
    out.write("=" * 100 + "\n")
    out.write("SIGNATURE VARIANT ANALYSIS\n")
    out.write("특정 샘플의 특징적 변이 분석\n")
    out.write("=" * 100 + "\n\n")

    # Summary statistics
    out.write("SUMMARY: Unique Variants per Sample\n")
    out.write("-" * 100 + "\n")
    out.write(f"{'Sample':<20} {'Total Variants':<15} {'Unique Variants':<18} {'Unique %':<12}\n")
    out.write("-" * 100 + "\n")

    # Sort by genetic distance (total variants)
    sorted_samples = sorted(samples, key=lambda s: sample_total[s], reverse=True)

    for sample in sorted_samples:
        total = sample_total[sample]
        unique = len(unique_variants.get(sample, []))
        unique_pct = (unique / total * 100) if total > 0 else 0
        out.write(f"{sample:<20} {total:<15} {unique:<18} {unique_pct:<12.1f}\n")

    out.write("\n" + "=" * 100 + "\n\n")

    # Detailed analysis for each sample
    for sample in sorted_samples:
        out.write("\n" + "=" * 100 + "\n")
        out.write(f"SAMPLE: {sample}\n")
        out.write(f"Total variants: {sample_total[sample]}, Missense: {sample_missense[sample]}\n")
        out.write("=" * 100 + "\n\n")

        # Unique variants
        uniq_vars = unique_variants.get(sample, [])
        if uniq_vars:
            uniq_missense = [v for v in uniq_vars if v['type'] == 'missense']

            out.write(f"UNIQUE VARIANTS: {len(uniq_vars)} total ({len(uniq_missense)} missense)\n")
            out.write("샘플 고유 변이 (다른 샘플에는 없음)\n")
            out.write("-" * 100 + "\n")

            if uniq_missense:
                out.write("Missense variants:\n")
                out.write(f"{'Pos':<8} {'Ref>Alt':<10} {'Protein Change':<20} {'Gene':<30}\n")
                out.write("-" * 100 + "\n")

                for var in sorted(uniq_missense, key=lambda x: x['pos']):
                    ref_alt = f"{var['ref']}>{var['alt']}"
                    gene_name = gene_names.get(var['gene'], var['gene'])
                    out.write(f"{var['pos']:<8} {ref_alt:<10} {var['protein']:<20} {gene_name:<30}\n")
                out.write("\n")

            # Unique synonymous
            uniq_syn = [v for v in uniq_vars if v['type'] == 'synonymous']
            if uniq_syn:
                out.write(f"Synonymous variants: {len(uniq_syn)}\n")
                out.write(f"Positions: {', '.join(str(v['pos']) for v in sorted(uniq_syn, key=lambda x: x['pos'])[:20])}")
                if len(uniq_syn) > 20:
                    out.write(f" ... and {len(uniq_syn)-20} more")
                out.write("\n\n")
        else:
            out.write("UNIQUE VARIANTS: None\n")
            out.write("이 샘플은 고유 변이가 없음 (모든 변이가 다른 샘플과 공유됨)\n\n")

        # Rare variants (1-2 samples)
        rare_vars = [v for v in rare_variants.get(sample, []) if v['type'] == 'missense' and v['count'] == 2]
        if rare_vars:
            out.write(f"RARE MISSENSE VARIANTS (shared with only 1 other sample): {len(rare_vars)}\n")
            out.write("-" * 100 + "\n")
            out.write(f"{'Pos':<8} {'Ref>Alt':<10} {'Protein Change':<20} {'Shared with':<20} {'Gene':<30}\n")
            out.write("-" * 100 + "\n")

            for var in sorted(rare_vars, key=lambda x: x['pos']):
                ref_alt = f"{var['ref']}>{var['alt']}"
                other_samples = [s for s in var['samples'] if s != sample]
                shared_with = ', '.join(other_samples)
                gene_name = gene_names.get(var['gene'], var['gene'])
                out.write(f"{var['pos']:<8} {ref_alt:<10} {var['protein']:<20} {shared_with:<20} {gene_name:<30}\n")
            out.write("\n")

        out.write("\n")

print(f"\nAnalysis complete!")
print(f"Report saved to: {output_file}")

# Print summary to console
print("\n" + "=" * 80)
print("UNIQUE VARIANTS SUMMARY (샘플별 고유 변이)")
print("=" * 80)
print(f"{'Sample':<20} {'Distance':<12} {'Unique Vars':<15} {'% Unique'}")
print("-" * 80)

for sample in sorted_samples:
    total = sample_total[sample]
    unique = len(unique_variants.get(sample, []))
    unique_pct = (unique / total * 100) if total > 0 else 0
    print(f"{sample:<20} {total:<12} {unique:<15} {unique_pct:.1f}%")
