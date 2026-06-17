#!/usr/bin/env python3
"""
Analyze variant sharing patterns to identify phylogenetic groups
계통 그룹 식별을 위한 변이 공유 패턴 분석
"""

import re
from collections import defaultdict

def parse_vcf_for_groups(vcf_file):
    """Parse VCF and collect variant sharing patterns"""

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

            effect = ann_parts[1]
            protein_change = ann_parts[10] if len(ann_parts) > 10 else ''

            # Determine variant type
            if 'synonymous' in effect:
                variant_type = 'synonymous'
            elif 'missense' in effect:
                variant_type = 'missense'
            else:
                variant_type = 'other'

            # Check which samples have this variant
            samples_with_variant = set()
            for i, (sample, gt) in enumerate(zip(samples, genotypes)):
                gt_value = gt.split(':')[0] if ':' in gt else gt
                if gt_value not in ['0', '.', './.']:
                    samples_with_variant.add(sample)

            if samples_with_variant:
                variant_data.append({
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'type': variant_type,
                    'protein': protein_change,
                    'samples': frozenset(samples_with_variant)
                })

    return samples, variant_data

def analyze_groups(samples, variant_data):
    """Identify groups based on shared variants"""

    # Count variants by sample group
    group_counts = defaultdict(lambda: {'total': 0, 'missense': 0, 'variants': []})

    for var in variant_data:
        sample_set = var['samples']
        group_counts[sample_set]['total'] += 1
        if var['type'] == 'missense':
            group_counts[sample_set]['missense'] += 1
            group_counts[sample_set]['variants'].append(var)

    # Sort groups by size (number of samples) and variant count
    sorted_groups = sorted(group_counts.items(),
                          key=lambda x: (len(x[0]), x[1]['total']),
                          reverse=True)

    return sorted_groups

# Parse VCF
print("Parsing VCF file...")
vcf_file = 'variants_annotated.vcf'
samples, variant_data = parse_vcf_for_groups(vcf_file)

print("Analyzing variant groups...")
sorted_groups = analyze_groups(samples, variant_data)

# Calculate total variants per sample
sample_totals = defaultdict(int)
for var in variant_data:
    for sample in var['samples']:
        sample_totals[sample] += 1

# Write report
output_file = 'variant_groups.txt'
with open(output_file, 'w') as out:
    out.write("=" * 120 + "\n")
    out.write("VARIANT SHARING PATTERN ANALYSIS\n")
    out.write("변이 공유 패턴 분석 - 계통 그룹 식별\n")
    out.write("=" * 120 + "\n\n")

    out.write("MAJOR GROUPS (공유 변이 패턴으로 식별된 주요 그룹)\n")
    out.write("=" * 120 + "\n\n")

    # Focus on groups with multiple samples and significant variant counts
    major_groups = [(g, c) for g, c in sorted_groups if len(g) >= 2 and c['total'] >= 5]

    for i, (sample_group, counts) in enumerate(major_groups[:30], 1):
        sample_list = sorted(list(sample_group))
        out.write(f"\nGROUP {i}: {len(sample_list)} samples\n")
        out.write("-" * 120 + "\n")
        out.write(f"Samples: {', '.join(sample_list)}\n")
        out.write(f"Shared variants: {counts['total']} total ({counts['missense']} missense)\n")

        if counts['missense'] > 0:
            out.write(f"\nShared missense variants:\n")
            out.write(f"{'Pos':<10} {'Ref>Alt':<10} {'Protein Change':<30}\n")
            out.write("-" * 120 + "\n")

            for var in sorted(counts['variants'][:20], key=lambda x: x['pos']):
                ref_alt = f"{var['ref']}>{var['alt']}"
                out.write(f"{var['pos']:<10} {ref_alt:<10} {var['protein']:<30}\n")

            if len(counts['variants']) > 20:
                out.write(f"... and {len(counts['variants']) - 20} more missense variants\n")

        out.write("\n")

    # Pairwise similarity analysis
    out.write("\n" + "=" * 120 + "\n")
    out.write("PAIRWISE SIMILARITY MATRIX\n")
    out.write("샘플 간 유사도 매트릭스 (공유 변이 개수)\n")
    out.write("=" * 120 + "\n\n")

    # Build pairwise matrix
    sample_variants = defaultdict(set)
    for var in variant_data:
        var_id = (var['pos'], var['ref'], var['alt'])
        for sample in var['samples']:
            sample_variants[sample].add(var_id)

    sorted_samples = sorted(samples, key=lambda s: sample_totals[s], reverse=True)

    # Write matrix
    out.write(f"{'Sample':<15}")
    for s in sorted_samples:
        out.write(f"{s:<10}")
    out.write("\n" + "-" * 120 + "\n")

    for s1 in sorted_samples:
        out.write(f"{s1:<15}")
        for s2 in sorted_samples:
            if s1 == s2:
                out.write(f"{'-':<10}")
            else:
                shared = len(sample_variants[s1] & sample_variants[s2])
                out.write(f"{shared:<10}")
        out.write("\n")

    # Jaccard similarity
    out.write("\n" + "=" * 120 + "\n")
    out.write("JACCARD SIMILARITY INDEX\n")
    out.write("자카드 유사도 지수 (0-1, 높을수록 유사)\n")
    out.write("=" * 120 + "\n\n")

    out.write(f"{'Sample':<15}")
    for s in sorted_samples:
        out.write(f"{s:<10}")
    out.write("\n" + "-" * 120 + "\n")

    for s1 in sorted_samples:
        out.write(f"{s1:<15}")
        for s2 in sorted_samples:
            if s1 == s2:
                out.write(f"{'-':<10}")
            else:
                intersection = len(sample_variants[s1] & sample_variants[s2])
                union = len(sample_variants[s1] | sample_variants[s2])
                jaccard = intersection / union if union > 0 else 0
                out.write(f"{jaccard:<10.3f}")
        out.write("\n")

print(f"\nAnalysis complete!")
print(f"Report saved to: {output_file}")

# Print key findings
print("\n" + "=" * 80)
print("KEY FINDINGS:")
print("=" * 80)

# Find closest pairs
print("\nMost similar pairs (highest Jaccard similarity):")
similarities = []
for s1 in samples:
    for s2 in samples:
        if s1 < s2:  # Avoid duplicates
            intersection = len(sample_variants[s1] & sample_variants[s2])
            union = len(sample_variants[s1] | sample_variants[s2])
            jaccard = intersection / union if union > 0 else 0
            similarities.append((s1, s2, jaccard, intersection))

similarities.sort(key=lambda x: x[2], reverse=True)
for s1, s2, jaccard, shared in similarities[:5]:
    print(f"{s1} <-> {s2}: Jaccard={jaccard:.3f}, Shared={shared} variants")

# Find large groups
print("\nLargest variant-sharing groups:")
large_groups = [(g, c) for g, c in sorted_groups if len(g) >= 3 and c['total'] >= 10][:5]
for i, (group, counts) in enumerate(large_groups, 1):
    print(f"Group {i}: {len(group)} samples, {counts['total']} shared variants ({counts['missense']} missense)")
    print(f"  Samples: {', '.join(sorted(list(group)))}")
