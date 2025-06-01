import sys
import argparse

def read_fasta(file_path):
    """Read FASTA file and return a dictionary of gene names and sequences."""
    sequences = {}
    current_gene = None
    current_sequence = []
    
    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if current_gene is not None:
                        sequences[current_gene] = ''.join(current_sequence)
                    current_gene = line[1:].split()[0]  # Extract gene name (ignore description)
                    current_sequence = []
                else:
                    current_sequence.append(line)
            if current_gene is not None:
                sequences[current_gene] = ''.join(current_sequence)
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        sys.exit(1)
    
    return sequences

def get_original_positions(gene_name, positions, sequences):
    """Convert multiple alignment positions to original sequence positions."""
    if gene_name not in sequences:
        return [f"Error: Gene {gene_name} not found in the FASTA file."] * len(positions)
    
    sequence = sequences[gene_name]
    results = []
    
    for position in positions:
        if position < 1 or position > len(sequence):
            results.append(f"Error: Position {position} is out of range for sequence {gene_name} (length: {len(sequence)}).")
        else:
            # Extract subsequence from start to the given position (1-based indexing)
            sub_sequence = sequence[:position]
            # Remove gaps ('-') and count non-gap characters
            original_position = len(sub_sequence.replace('-', ''))
            results.append(original_position)
    
    return results

def get_alignment_positions(gene_name, original_positions, sequences):
    """Convert multiple original sequence positions to alignment positions."""
    if gene_name not in sequences:
        return [f"Error: Gene {gene_name} not found in the FASTA file."] * len(original_positions)
    
    sequence = sequences[gene_name]
    max_non_gap = len(sequence.replace('-', ''))  # Total non-gap characters
    
    results = []
    for orig_pos in original_positions:
        if orig_pos < 1 or orig_pos > max_non_gap:
            results.append(f"Error: Original position {orig_pos} is out of range (max non-gap length: {max_non_gap}).")
            continue
        
        non_gap_count = 0
        for align_pos, char in enumerate(sequence, 1):  # 1-based indexing
            if char != '-':
                non_gap_count += 1
            if non_gap_count == orig_pos:
                results.append(align_pos)
                break
        else:
            results.append(f"Error: Could not find alignment position for original position {orig_pos}.")
    
    return results

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Convert positions between alignment and original sequence.")
    parser.add_argument("file_path", help="Path to the FASTA file")
    parser.add_argument("gene_name", help="Name of the gene (e.g., gene-AT1G12940)")
    parser.add_argument("-v", "--invert", action="store_true", help="Convert original sequence positions to alignment positions")
    
    try:
        args = parser.parse_args()
    except SystemExit:
        print("Error: Invalid command-line arguments.")
        print("Usage: python convert_position.py <file_path> <gene_name> [-v]")
        sys.exit(1)
    
    # Debug: Print parsed arguments
    print(f"Debug: Parsed arguments - file_path: {args.file_path}, gene_name: {args.gene_name}, invert: {args.invert}")
    
    # Read FASTA file
    sequences = read_fasta(args.file_path)
    
    if args.invert:
        # Reverse conversion: original positions to alignment positions
        positions_input = input("Enter original sequence positions (comma-separated, e.g., 1,2,3): ")
        try:
            original_positions = [int(p.strip()) for p in positions_input.split(',')]
        except ValueError:
            print("Error: Invalid input. Please enter comma-separated integers (e.g., 1,2,3).")
            sys.exit(1)
        
        results = get_alignment_positions(args.gene_name, original_positions, sequences)
        
        for orig_pos, result in zip(original_positions, results):
            if isinstance(result, int):
                print(f"Original position {orig_pos} in gene {args.gene_name} corresponds to alignment position: {result}")
            else:
                print(result)
    else:
        # Forward conversion: alignment positions to original positions
        positions_input = input("Enter alignment positions (comma-separated, e.g., 1,2,3): ")
        try:
            positions = [int(p.strip()) for p in positions_input.split(',')]
        except ValueError:
            print("Error: Invalid input. Please enter comma-separated integers (e.g., 1,2,3).")
            sys.exit(1)
        
        results = get_original_positions(args.gene_name, positions, sequences)
        
        for pos, result in zip(positions, results):
            if isinstance(result, int):
                print(f"Gene {args.gene_name} at alignment position {pos} corresponds to original sequence position: {result}")
            else:
                print(result)

if __name__ == '__main__':
    main()
