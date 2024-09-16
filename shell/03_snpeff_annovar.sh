#!/bin/bash


#snpeff
#snpEff databases > snpEff.databases.list.txt
#1. make database
snpEff -Xmx40g build -c snpEff.config  -gtf22 -v pandaV3
#2. annotation snp and indel
snpEff -Xmx40g -c snpEff.config -v pandaV3 ${vcf_file} > ${snpeff_out}/anno.vcf -s ${snpeff_out}/summary.html  -csvStats ${snpeff_out}/summary.csv
cat ${snpeff_out}/anno.vcf | SnpSift filter "(exists LOF[*].PERC)  & (LOF[*].PERC > 0)" > ${snpeff_out}/LOF.anno.vcf
cat ${snpeff_out}/anno.vcf | SnpSift filter "ANN[*].EFFECT has 'missense_variant'" > ${snpeff_out}/missense.anno.vcf
cat ${snpeff_out}/anno.vcf | SnpSift filter "ANN[*].EFFECT has 'synonymous_variant'" > ${snpeff_out}/synonymous.anno.vcf

#annovar
#1.
gtfToGenePred -genePredExt  ${gtf} ${annovar_dir}/reference_refGene.txt -ignoreGroupsWithoutExons
#2.
refasta_pl='retrieve_seq_from_fasta.pl'
perl ${refasta_pl} --format refGene --seqfile ${refgenomefasta} ${annovar_dir}/reference_refGene.txt --out ${annovar_dir}/reference_refGeneMrna.fa
#3.
convert2annovar_pl='convert2annovar.pl'
perl ${convert2annovar_pl} --allsample -withfreq -format vcf4 ${vcf} > ${annovar_dir}/annovar.avinput 
#4.
annotate_variation_pl='annotate_variation.pl'
perl ${annotate_variation_pl} --dbtype refGene -build reference -geneanno -outfile ${annovar_dir}/anno ${annovar_dir}/annovar.avinput reference --separate
