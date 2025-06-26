#!/bin/bash

# Extract the longest protein
python ${python_script} ${gff_file} ${output_file}
gffread ${output_file} -g ${genome_folder}/${species_name}.fa -y ${pep_folder}/${species_name}_longest_pep.fa

# Orthofinder
python ${orthofinder.py} -f ${longest_pep} -t 8 -a 1 

# Run cactus using docker
# docker run cactus
# docker run cactus-hal2maf

# Extract MAF
python ${python_script} cactus_bear_no_anc.maf 6_cactus_bear_no_anc.maf

# Identify variant sites between giant panda and other species
python ${python_script} 6_cactus_bear_no_anc.maf PSV_add_site.bed

# Replace RMA and perform deleterious mutation prediction
python ${python_script} ${bed} ${vcf} ${ref_fa} ${RMA_fa} ${RMA_vcf}

# SIFT
# Build database
perl make-SIFT-db-all.pl -config ../sift_out/bears_V3_RMA.txt
#docker run sift4g

# SIFT4G annotation
java -jar SIFT4G_Annotator.jar -c -i ${missense_snp_RMA.vcf} -d ${bears_V3_RMA} -r ${sift_result} -t

# find deleterious and tolerated mutations
python ${python_script} 


