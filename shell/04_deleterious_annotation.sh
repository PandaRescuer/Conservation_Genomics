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
python extract_sites_from_maf.py ${cactus.all_species.maf} ${site.bed}

# Replace RMA
python correct_RMA_sites_in_vcf.py ${site.diff.bed} ${site.same.bed}  ${mis.vcf} ${mis.RMA.vcf} ${syn.vcf} ${syn.RMA.vcf} ${CNE.vcf} ${CNE.RMA.vcf} ${LOF.vcf} ${LOF.RMA.vcf}

# SIFT
# Build database
perl make-SIFT-db-all.pl -config ${config_file}
# docker run sift4g

# SIFT4G annotation
java -jar SIFT4G_Annotator.jar -c -i ${mis.RMA.vcf} -d ${sift_database} -r ${sift_result} -t

# find deleterious and tolerated mutations
python classify_vcf_by_sift_prediction.py

# Extract the strict subset from the SIFT prediction set
python filter_vcf_by_RMA.py ${mis.RMA.vcf} ${DEL.vcf} ${TOL.vcf} ${DEL.filtered.vcf} ${TOL.filtered.vcf} 2> filter_RMA.stats.log

#Correct the reference genome using the RMA VCF.
correct_RMA_sites_in_ref_genome.py ${site.diff.bed} ${site.same.bed} Ailuropoda_melanoleuca.fasta ASM200744v3_AncestralNonLOF_V1.fa ${mis.vcf} ${mis.RMA.vcf} ${syn.vcf} ${syn.RMA.vcf} ${CNE.vcf} ${CNE.RMA.vcf}

#Correct the reference genome using the RMA VCF, including LOF.
correct_RMA_sites_in_ref_genome.py ${site.diff.bed} ${site.same.bed} Ailuropoda_melanoleuca.fasta ASM200744v3_AncestralAllSNP_V1.fa ${mis.vcf} ${mis.RMA.vcf} ${syn.vcf} ${syn.RMA.vcf} ${CNE.vcf} ${CNE.RMA.vcf} ${LOF.vcf} ${LOF.RMA.vcf}
