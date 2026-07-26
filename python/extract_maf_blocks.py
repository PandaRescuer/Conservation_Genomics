import sys

def process_block(block, output_file):
    s_count = 0
    for line in block:
        if line.startswith('s'):
            s_count += 1
    if s_count == 6:
        output_file.writelines(block)
        output_file.write('\n')  

def filter_maf(input_filename, output_filename):
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        current_block = None
        for line in infile:
            if line.startswith('a'):
                if current_block is not None:
                    process_block(current_block, outfile)
                current_block = [line]
            else:
                if current_block is not None:
                    current_block.append(line)
        if current_block is not None:
            process_block(current_block, outfile)

if __name__ == "__main__":
    input_file = sys.argv[1]   
    output_file = sys.argv[2]  
    filter_maf(input_file, output_file)

