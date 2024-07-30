#!/bin/bash

samtools faidx ${ref_fasta}

gatk CreateSequenceDictionary -R ${ref_fasta} -O ${ref_dict}

bwa index -p ${refgenome} -a bwtsw  ${ref_fasta}

# fastp to filter low quality raw reads
fastp -i ${fastq_1} -I ${fastq_2} -o ${fsg_fastq_dir}/${sample_array[$i+j]}/${sample_array[$i+j]}_1.clean.fq.gz -O ${fsg_fastq_dir}/${sample_array[$i+j]}/${sample_array[$i+j]}_2.clean.fq.gz -z 6 -W 5 -r -M 20 -g -q 5 -u 50 -n 15 -l 50 --overlap_diff_limit 1 --overlap_diff_percent_limit 10 -h ${fastqc_dir}/${sample_array[$i+j]}.fastp.html

#bwa
bwa mem -M -t 10 -R "@RG\tID:${sample_array[$i+j]}\tPL:UNKNOWN\tLB:${sample_array[$i+j]}\tSM:${sample_array[$i+j]}" ${refgenome} ${fastq_1} ${fastq_2} | samtools view -S -b - > ${bwa_out_dir}/${sample_array[$i+j]}.bam

#sort
samtools sort -@ 10 -m 40G -O bam -o ${sort_out_dir}/${sample_name}.sorted.bam ${sort_out_dir}/${sample_array_filtered[$i+j]}

#markdup
picard MarkDuplicates -Xmx20g REMOVE_DUPLICATES=true TMP_DIR=${markdup_out_dir}/tmp I=${bam_dir}/${sample_array_filtered[$i+j]} O=${markdup_out_dir}/${sample_name}.bam M=${markdup_out_dir}/${sample_name}.markdup_metrics.txt

#samtools remove bad reads
samtools view -hb -f 2 -F 256 -q 30 ${markdup_out_dir}/${sample_array_filtered[$i+j]} | samtools view -hb -F 1024 > ${samtools_dir}/${sample_name}.removebadreads.bam

#index
samtools index -@ 4 ${samtools_dir}/${sample_array_filtered[$i+j]}

#mosdepth
mosdepth -t 4 -b ${bed} -T 1,2,3,4,5,10,15,20,30,40 ${coverage_out_dir}/${sample_array_filtered[$i+j]} ${samtools_dir}/${sample_array_filtered[$i+j]}

# calculate 99% max depth of individual
python cal_99depth.py

# merge 99% depth file
python cal_and_merge_depth.py

#gatk
gatk HaplotypeCaller --java-options "-Xmx40g" -mbq 20  --emit-ref-confidence GVCF -R ${ref_fasta} -I ${samtools_dir}/${sample_array_filtered[$i+j]} -O ${gatk_out_dir}/${sample_name}.HC.g.vcf.gz

gatk CombineGVCFs --java-options "-Xmx40g" -R ${ref_fasta} ${gvcfs} -L NC_0482${chr_num}.1 -O ${gatk_combine_out_dir}/group1.NC_0482${chr_num}.1.HC.g.vcf.gz
gatk CombineGVCFs --java-options "-Xmx40g" -R ${ref_fasta} ${gvcfs} -L ${scaffold_bed} -O ${gatk_combine_out_dir}/group1.scaffold.HC.g.vcf.gz
gatk CombineGVCFs --java-options "-Xmx20g" -R ${ref_fasta} ${gvcfs} -O ${gatk_combine_out_dir}/all.NC_0482${chr_num}.1.HC.g.vcf.gz
gatk CombineGVCFs --java-options "-Xmx20g" -R ${ref_fasta} ${gvcfs} -O ${gatk_combine_out_dir}/all.scaffold.HC.g.vcf.gz

gatk GenotypeGVCFs --java-options "-Xmx100g" -stand-call-conf 0 -R ${ref_fasta} -V ${gatk_combine_out_dir}/${combine_name}.NC_0482${i}.1.HC.g.vcf.gz -O ${gatk_combine_out_dir}/${combine_name}.NC_0482${i}.1.HC.vcf.gz
gatk GenotypeGVCFs --java-options "-Xmx100g" -stand-call-conf 0 -R ${ref_fasta} -V ${gatk_combine_out_dir}/${combine_name}.scaffold.HC.g.vcf.gz -O ${gatk_combine_out_dir}/${combine_name}.scaffold.HC.vcf.gz

picard MergeVcfs ${vcfs} -O ${gatk_combine_out_dir}/merge.vcf.gz

gatk SelectVariants -select-type SNP -V ${gatk_combine_out_dir}/test10_merge.vcf.gz -O ${gatk_combine_out_dir}/${combine_name}.HC.snp.vcf.gz
gatk SelectVariants -select-type INDEL -V ${gatk_final_vcf_out_dir}/merge.vcf.gz -O ${gatk_final_vcf_out_dir}/merge.indel.vcf.gz

#hard filter snp
gatk VariantFiltration -V ${gatk_combine_out_dir}/${combine_name}.HC.snp.vcf.gz --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "snp_filter" -O ${gatk_combine_out_dir}/${combine_name}.HC.snp.filter.vcf.gz

#simple filter
vcftools --gzvcf ${raw_vcf_file} --minDP 3 --minQ 30 --minGQ 20 --recode-INFO-all --recode --out ${QC_vcf_file}

#keep autosome
python filter_chr.py ${QC_vcf_file}.recode.vcf ${chr_vcf_file}

#final filter
python filter_snp.py ${chr_vcf_file} ${final_vcf_file} ${filtered_vcf_file} ${depth99_file}

#hard filter indel
gatk VariantFiltration -V ${gatk_final_vcf_out_dir}/merge.indel.vcf.gz --filter-expression "QD < 2.0 || FS > 200.0 || MQ < 40.0 || SOR > 3.0|| ReadPosRankSum < -20.0  " --filter-name "indel_filter" -O ${gatk_final_vcf_out_dir}/merge.indel.filter.vcf.gz
vcftools --gzvcf ${gatk_final_vcf_out_dir}/merge.indel.vcf.gz --minDP 4 --minGQ 20 --recode-INFO-all --recode --out ${gatk_final_vcf_out_dir}/merge.indel.QC
python filter_chr.py ${gatk_final_vcf_out_dir}/merge.indel.QC.recode.vcf ${gatk_final_vcf_out_dir}/merge.chr.indel.QC.recode.vcf
python filter_indel.py ${gatk_final_vcf_out_dir}/merge.chr.indel.QC.recode.vcf ${final_vcf_file} ${filtered_vcf_file} ${depth99_file}