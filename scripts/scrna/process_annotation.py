#!/usr/bin/env python3
import re
import sys

def process_file(input_file, output_file):
    """
    Process an annotation file and write a five-column output file.
    Columns: Symbol, Mutation ID, Gene ID, Genome region, Function prediction
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Write header
        outfile.write("Symbol\tMutation ID\tGene ID\tGenome region\tFunction prediction\n")
        
        for line in infile:
            line = line.strip()
            if not line:  # skip empty lines
                continue
            
            # Split each line; three columns are expected
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            
            # Parse the first column: chr2_3576510_G_A
            mutation_parts = parts[0].split('_')
            if len(mutation_parts) == 4:
                chr_part, pos_part, ref_part, alt_part = mutation_parts
                # Convert to chr1:39034563,T>A format
                mutation_id = f"{chr_part}:{pos_part},{ref_part}>{alt_part}"
            else:
                # If it is not the expected 4-part format, keep the original value
                mutation_id = parts[0]
            
            # Parse the second column: exonic(synonymous_SNV) or ncRNA_exonic
            second_col = parts[1]
            genome_region = ""
            function_pred = "unknown"
            
            # Check whether parentheses are present
            match = re.search(r'^(.*?)\((.*?)\)$', second_col)
            if match:
                # Parentheses present
                genome_region = match.group(1)  # text before parentheses
                function_pred = match.group(2)  # text inside parentheses
            else:
                # No parentheses
                genome_region = second_col
                # function_pred remains "unknown"
            
            # Third column: Symbol
            symbol = parts[0]
            
            # Write the output line
            outfile.write(f"{symbol}\t{mutation_id}\t{symbol}\t{genome_region}\t{function_pred}\n")

def main():
    # Check command-line arguments
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
        print("Example: python script.py Org4S15D63.benchmark_identifier.anno_info.txt output.txt")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    try:
        process_file(input_file, output_file)
        print(f"Processing complete! Output written to {output_file}")
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
