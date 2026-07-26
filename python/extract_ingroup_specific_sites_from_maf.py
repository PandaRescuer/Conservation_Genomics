import sys

def process_block(block, output_file):
    s_count = 0
    output=False
    for line in block:
        if line.startswith('s'):
            s_count += 1
    if s_count == 6:
        seq_list=[]
        ref_info=block[1].split()
        Chr_id_info=ref_info[1].split('.')
        Chr_id='.'.join(Chr_id_info[1:])
        start_index = int(ref_info[2])
        for align_item in block[1:]:
            if len(align_item)>6:
                align_item_info=align_item.split()
                seq_list.append(align_item_info[6])
                
        for site_index in range(len(seq_list[0])):
            
            ingroup_sites=''
            for ingroup_index_item in ingroup_index:
                ingroup_site=seq_list[ingroup_index_item][site_index].upper()
                ingroup_sites+=ingroup_site

            outgroup_sites=''
            for outgroup_index_item in outgroup_index:
                outgroup_site=seq_list[outgroup_index_item][site_index].upper()
                outgroup_sites+=outgroup_site
                
            if len(set(outgroup_sites))==1 and len(set(ingroup_sites))==1:
                if ingroup_site!=outgroup_site and outgroup_site!='-':
                    output=True
                    site_index_with_line=seq_list[0][:site_index+1].count('-')
                    output_file.write('\t'.join([Chr_id,str(start_index+site_index-site_index_with_line),str(start_index+site_index-site_index_with_line+1),ingroup_site+','+outgroup_site])+'\n')
                   
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
    ingroup_index=[0]
    outgroup_index=[1,2,3,4,5]
    filter_maf(input_file, output_file)

