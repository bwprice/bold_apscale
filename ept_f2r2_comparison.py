import time
from Bio import SeqIO
import csv
from collections import defaultdict
import os

def main():
    # Start timing
    start_time = time.time()
    
    # File paths (using absolute paths)
    ept_file = "D:\\_claude_files\\projects\\geneflow\\geneflow_EPT_ESVs.fasta"
    f2r2_file = "D:\\_claude_files\\projects\\geneflow\\geneflow_F2R2_ESVs.fasta"
    output_csv = "D:\\_claude_files\\projects\\geneflow\\sequence_matches.csv"
    
    print(f"Reading EPT sequences from {ept_file}...")
    # Create a dictionary to store EPT sequences and their IDs
    ept_sequences = {}
    for record in SeqIO.parse(ept_file, "fasta"):
        # Convert sequence to string and store in dictionary with ID as key
        ept_sequences[str(record.seq).upper()] = record.id
    
    print(f"Loaded {len(ept_sequences)} EPT sequences")
    
    print(f"Reading F2R2 sequences from {f2r2_file}...")
    # Dictionary to store F2R2 sequences' first 142bp and their IDs
    f2r2_sequences = defaultdict(list)
    for record in SeqIO.parse(f2r2_file, "fasta"):
        # Take only the first 142bp and convert to uppercase
        truncated_seq = str(record.seq[:142]).upper()
        # Store ID in a list associated with this sequence
        # (multiple F2R2 sequences might have identical first 142bp)
        f2r2_sequences[truncated_seq].append(record.id)
    
    print(f"Loaded {sum(len(ids) for ids in f2r2_sequences.values())} F2R2 sequences")
    
    # Find matches
    print("Finding matches...")
    matches = []
    
    # For each EPT sequence, check if it exists in F2R2 sequences
    for ept_seq, ept_id in ept_sequences.items():
        if ept_seq in f2r2_sequences:
            # If match found, add all matching F2R2 IDs to results
            for f2r2_id in f2r2_sequences[ept_seq]:
                matches.append((ept_id, f2r2_id))
    
    print(f"Found {len(matches)} matches")
    
    # Write matches to CSV
    print(f"Writing results to {output_csv}...")
    with open(output_csv, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        # Write header
        csv_writer.writerow(['EPT_ID', 'F2R2_ID'])
        # Write matches
        csv_writer.writerows(matches)
    
    # Calculate and print elapsed time
    elapsed_time = time.time() - start_time
    print(f"Completed in {elapsed_time:.2f} seconds")
    print(f"Results saved to {output_csv}")

if __name__ == "__main__":
    main()