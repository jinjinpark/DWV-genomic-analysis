#!/usr/bin/env python3
"""
Convert protein coordinate GFF3 to DNA coordinate GFF3
Protein position n -> DNA position (n-1)*3+1 to n*3
"""

def protein_to_dna(protein_start, protein_end):
    """Convert protein coordinates to DNA coordinates"""
    dna_start = (protein_start - 1) * 3 + 1
    dna_end = protein_end * 3
    return dna_start, dna_end

def convert_gff_line(line):
    """Convert a single GFF3 line from protein to DNA coordinates"""
    if line.startswith('#'):
        return line

    parts = line.strip().split('\t')
    if len(parts) < 9:
        return line

    # GFF3 format: seqid, source, type, start, end, score, strand, phase, attributes
    try:
        protein_start = int(parts[3])
        protein_end = int(parts[4])

        dna_start, dna_end = protein_to_dna(protein_start, protein_end)

        # Update coordinates
        parts[3] = str(dna_start)
        parts[4] = str(dna_end)

        # For CDS features, we need to add phase
        if parts[2] == 'CDS' or parts[2] == 'Chain':
            parts[2] = 'CDS'  # Convert Chain to CDS
            parts[6] = '+'    # strand
            parts[7] = '0'    # phase

        return '\t'.join(parts)
    except (ValueError, IndexError):
        return line

# Read and convert
input_file = 'reference.tsv'
output_file = 'reference_dna.gff3'

with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
    for line in f_in:
        converted_line = convert_gff_line(line)
        f_out.write(converted_line + '\n')

print(f"Converted {input_file} to {output_file}")
print("Protein coordinates -> DNA coordinates (x3)")
