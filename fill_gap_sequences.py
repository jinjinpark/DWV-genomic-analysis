import os
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Create 'result' subfolder
result_dir = "result"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
    
# Find all files with .fas extension in the current directory
input_files = [f for f in os.listdir('.') if f.endswith('.fas')]

def find_gap_positions(sequence):
    """Find positions of '-' or 'N' in the sequence"""
    return [i for i, base in enumerate(sequence) if base in ['-', 'N', 'n']]
    
def fill_gaps(consensus_seq, amplicon_seqs):
    """Fill gaps in consensus sequence using amplicon sequences"""
    filled_seq = list(consensus_seq)
    gap_positions = find_gap_positions(consensus_seq)
    
    for pos in gap_positions:
        for amplicon in amplicon_seqs:
            if pos < len(amplicon) and amplicon[pos] not in ['-', 'N', 'n']:
                filled_seq[pos] = amplicon[pos]
                break
        else:
            print(f"Could not fill gap at position {pos}")
    
    return ''.join(filled_seq)
    
# Iterate over each file
for input_file in input_files:
    # Set file path
    file_path = os.path.join('.', input_file)
    
    # Read aligned sequences
    alignment = AlignIO.read(file_path, "fasta")
    aligned_consensus_seq = str(alignment[0].seq)  # Extract consensus sequence
    aligned_amplicon_seqs = [str(record.seq) for record in alignment[1:]]  # Extract amplicon sequences
    
    # Fill gaps
    filled_sequence = fill_gaps(aligned_consensus_seq, aligned_amplicon_seqs)
    
    # Set output file path
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_filename = f"{base_name}_gapfilled.fasta"
    output_path = os.path.join(result_dir, output_filename)
    
    # Save results
    filled_record = SeqRecord(Seq(filled_sequence), id=f"filled_consensus_S{base_name}")
    SeqIO.write(filled_record, output_path, "fasta")
    
    # Output message
    print(f"Gap-filled sequence for '{input_file}' has been saved to '{output_path}'")
    
    # First calculate statistics
    original_gaps = aligned_consensus_seq.count('-') + aligned_consensus_seq.count('N') + aligned_consensus_seq.count('n')
    filled_gaps = filled_sequence.count('-') + filled_sequence.count('N') + filled_sequence.count('n')
    
    # Then write log file
    log_filename = f"{base_name}_gap_statistics.log"
    log_path = os.path.join(result_dir, log_filename)
    with open(log_path, 'w') as log_file:
        log_file.write(f"Original number of gaps in '{input_file}': {original_gaps}\n")
        log_file.write(f"Number of gaps after filling in '{input_file}': {filled_gaps}\n")
        log_file.write(f"Number of gaps filled in '{input_file}': {original_gaps - filled_gaps}\n")