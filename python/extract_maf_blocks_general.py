# General-purpose utility for extracting MAF blocks containing a specified
# number of species or sequences.

import sys


def write_block_with_required_species(block, output_file, required_species_count):
    species_count = sum(1 for line in block if line.startswith('s'))

    if species_count == required_species_count:
        output_file.writelines(block)
        output_file.write('\n')  # Separate MAF blocks with a blank line


def extract_maf_blocks(
    input_filename,
    output_filename,
    required_species_count
):
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        current_block = None

        for line in infile:
            # Detect the start of a new MAF block
            if line.startswith('a'):
                if current_block is not None:
                    write_block_with_required_species(
                        current_block,
                        outfile,
                        required_species_count
                    )
                current_block = [line]
            elif current_block is not None:
                current_block.append(line)

        # Process the final MAF block
        if current_block is not None:
            write_block_with_required_species(
                current_block,
                outfile,
                required_species_count
            )


if __name__ == "__main__":
    input_file = sys.argv[1]  # Input MAF file path
    output_file = sys.argv[2]  # Output MAF file path
    required_species_count = int(sys.argv[3])  # Required number of species

    extract_maf_blocks(
        input_file,
        output_file,
        required_species_count
    )

    print(f"Processing complete. Results saved to {output_file}")