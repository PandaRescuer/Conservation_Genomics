#!/bin/bash


# -all-sites to call more than monomorphic sites
gatk GenotypeGVCFs --java-options "-Xmx40g" -all-sites -stand-call-conf 0 -R ${ref_fasta} -V ${gatk_combine_out_dir}/all.NC_0482${i}.1.HC.g.vcf.gz -O ${gatk_vcf_out_dir}/all.NC_0482${i}.1.HC.vcf.gz

# NO_VARIATION
gatk SelectVariants -select-type NO_VARIATION -V ${gatk_final_vcf_out_dir}/merge.vcf.gz -O ${gatk_final_vcf_out_dir}/merge.novariant.vcf.gz

# only use autosome 
python filter_chr.py ${gatk_final_vcf_out_dir}/merge.novariant.QC.recode.vcf ${novariant_vcf}

# merge snp and novariant
picard MergeVcfs -I ${novariant_vcf} -I ${snp_vcf} -O ${het_vcf}

# calculate genome-wide heterozygosity
plink -bfile ${bfile} --out ${out_dir}/all_heterozygosity --allow-extra-chr --noweb --het

# calculate per-site heterozygosity in non-overlapping 1-Mb windows across called autosome
python cal_het_sliding_window_single_sample.py ${het_vcf} ${bed_file} ${var} ${window_dir}

# vcf 2 bfile
plink -vcf ${vcf_file} --recode --out ${bfile} --allow-extra-chr --noweb --make-bed

# missing
plink -bfile ${bfile} --missing --out ${out_dir}/miss  --allow-extra-chr --noweb

# ROH
plink -bfile ${bfile} --homozyg --homozyg-density 50 --homozyg-gap 1000 --homozyg-kb 100 --homozyg-snp 25 --homozyg-window-het 1 --homozyg-window-snp 100 --homozyg-window-threshold 0.05 --out ${ROH_dir}/ROH

# PCA
gcta --bfile ${bfile} --make-grm --make-grm-alg 1 --out  ${pca_out}/all_grm
gcta --grm ${pca_out}/all_grm --pca 5 --out ${pca_out}/all_pca

# change chr id to number
python rename_chr.py ${final_vcf_file} ${rename_vcf}

# admixture
plink --vcf ${rename_snp_vcf} --allow-extra-chr --recode --out ${admixture_file}
plink -file ${admixture_file} --out ${admixture_file} --allow-extra-chr --noweb --make-bed
plink -bfile ${admixture_file} --missing --out ${admixture_dir}/admixture  --allow-extra-chr --noweb
plink --noweb --file ${admixture_file} --geno 0.05 --maf 0.05 --hwe 0.0001 --indep-pairwise 50 10 0.1 --make-bed --out ${QC_admixture} --allow-extra-chr
plink -bfile ${QC_admixture} --missing --out ${admixture_dir}/QC.admixture  --allow-extra-chr --noweb
plink  --bfile ${QC_admixture} --extract admixture.QC.prune.in --make-bed --out ${pruned_admixture}

for i in $(seq 2 5);do admixture --cv ${pruned_admixture}.bed ${i}|tee ${admixture_file}.log${i}.out;done
grep CV ${admixture_file}.log*.out|awk '{print NR","${NF}' >${admixture_dir}/admixture.ce.csv

# tree
python vcf2phylip.py -i ${rename_vcf}

phylip seqboot < seqboot.par
mv outfile seqboot.out
phylip dnadist < dnadist.par
mv outfile dnadist.out
phylip neighbor < neighbor.par
mv outfile neighbor.out
mv outtree neighbor.tree
phylip consense < consense.par
mv outfile consense.out
mv outtree consense.tree

