# General-purpose utility for identifying fixed differences between
# an ingroup and the remaining sequences in MAF alignment blocks.

import sys


def process_block(
    block,
    output_file,
    ingroup_indices,
    outgroup_indices,
    required_sequence_count
):
    sequence_count = sum(1 for line in block if line.startswith('s'))

    # Analyze blocks containing the required number of sequences
    if sequence_count != required_sequence_count:
        return

    sequence_list = []

    # Extract chromosome information and alignment start position
    reference_info = block[1].split()
    chromosome_parts = reference_info[1].split('.')
    chromosome_id = '.'.join(chromosome_parts[1:])
    start_position = int(reference_info[2])

    # Extract aligned sequences from MAF sequence lines
    for line in block[1:]:
        if line.startswith('s'):
            sequence_info = line.split()
            sequence_list.append(sequence_info[6])

    # Identify fixed differences between the ingroup and outgroup
    for site_index in range(len(sequence_list[0])):
        ingroup_sites = ''.join(
            sequence_list[index][site_index].upper()
            for index in ingroup_indices
        )
        outgroup_sites = ''.join(
            sequence_list[index][site_index].upper()
            for index in outgroup_indices
        )

        # Require identical bases within each group
        if len(set(ingroup_sites)) == 1 and len(set(outgroup_sites)) == 1:
            ingroup_base = ingroup_sites[0]
            outgroup_base = outgroup_sites[0]

            # Require a fixed difference and exclude outgroup gaps
            if ingroup_base != outgroup_base and outgroup_base != '-':
                gap_count = sequence_list[0][:site_index].count('-')
                genomic_position = start_position + site_index - gap_count

                # Write the variant as a BED record
                output_file.write(
                    '\t'.join([
                        chromosome_id,
                        str(genomic_position),
                        str(genomic_position + 1),
                        f'{ingroup_base},{outgroup_base}'
                    ]) + '\n'
                )


def extract_ingroup_specific_sites_from_maf(
    input_filename,
    output_filename,
    ingroup_indices,
    outgroup_indices,
    required_sequence_count
):
    with open(input_filename, 'r') as input_file, \
            open(output_filename, 'w') as output_file:

        current_block = None

        for line in input_file:
            # Detect the beginning of a new MAF block
            if line.startswith('a'):
                if current_block is not None:
                    process_block(
                        current_block,
                        output_file,
                        ingroup_indices,
                        outgroup_indices,
                        required_sequence_count
                    )
                current_block = [line]
            elif current_block is not None:
                current_block.append(line)

        # Process the final MAF block
        if current_block is not None:
            process_block(
                current_block,
                output_file,
                ingroup_indices,
                outgroup_indices,
                required_sequence_count
            )


if __name__ == "__main__":
    input_filename = sys.argv[1]  # Input MAF file path
    output_filename = sys.argv[2]  # Output BED file path
    required_sequence_count = int(sys.argv[3])  # Total number of sequences

    ingroup_indices = [0]
    outgroup_indices = list(range(1, required_sequence_count))

    extract_ingroup_specific_sites_from_maf(
        input_filename,
        output_filename,
        ingroup_indices,
        outgroup_indices,
        required_sequence_count
    )

    print(f'Processing complete. Results saved to {output_filename}')