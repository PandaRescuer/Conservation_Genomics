#!/bin/bash

# Extract the longest protein
python ${python_script} ${gff_file} ${output_file}
gffread ${output_file} -g ${genome_folder}/${species_name}.fa -y ${pep_folder}/${species_name}_longest_pep.fa

# Orthofinder to get guide tree
python ${orthofinder.py} -f ${longest_pep} -t 8 -a 1 

# Run cactus using docker
# docker run cactus
# docker run cactus-hal2maf

# Extract MAF
python extract_maf_blocks.py ${cactus.no_anc.maf} ${cactus.all_species.maf}

# Identify variant sites between giant panda and other species
python extract_ingroup_specific_sites_from_maf.py ${cactus.all_species.maf} ${site.bed}

# Replace RMA
python correct_RMA_sites_in_vcf_and_reference.py ${site.bed} ${vcf} ${ref_fa} ${RMA_fa} ${RMA_vcf}

# SIFT
# Build database
perl make-SIFT-db-all.pl -config ${config_file}
# docker run sift4g

# SIFT4G annotation
java -jar SIFT4G_Annotator.jar -c -i ${missense_snp_RMA.vcf} -d ${sift_database} -r ${sift_result} -t

# find deleterious and tolerated mutations
python classify_vcf_by_sift_prediction.py


