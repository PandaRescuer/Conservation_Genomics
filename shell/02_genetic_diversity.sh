#!/bin/bash

# Heterozygosity

# -all-sites to call more than monomorphic sites

gatk GenotypeGVCFs --java-options "-Xmx40g" -all-sites -stand-call-conf 0 -R ${ref_fasta} -V ${gatk_combine_out_dir}/all.NC_0482${i}.1.HC.g.vcf.gz -O ${gatk_vcf_out_dir}/all.NC_0482${i}.1.HC.vcf.gz
gatk SelectVariants -select-type NO_VARIATION -V ${gatk_final_vcf_out_dir}/merge.vcf.gz -O ${gatk_final_vcf_out_dir}/merge.novariant.vcf.gz
vcftools --gzvcf ${gatk_final_vcf_out_dir}/merge.novariant.vcf.gz --minDP 3 --minGQ 20 --minQ 30 --recode-INFO-all --recode --out ${gatk_final_vcf_out_dir}/merge.novariant.QC
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
bcftools +fill-tags ${vcf_file} -Oz -o ${vcf_file.filltags} -- -t AF,MAF,F_MISSING
bcftools roh -G30 --AF-tag AF --rec-rate 1e-8 -o ${bcftools_roh} ${vcf_file.filltags}

# PCA
gcta --bfile ${bfile} --make-grm --make-grm-alg 1 --out  ${pca_out}/all_grm
gcta --grm ${pca_out}/all_grm --pca 5 --out ${pca_out}/all_pca

# change chr id to number
python rename_chr.py ${final_vcf_file} ${rename_vcf}

# Admixture
plink --vcf "$rename_vcf" --allow-extra-chr --double-id --make-bed --out "$raw_prefix" --noweb
plink --bfile "$raw_prefix" --geno 0.05 --maf 0.05 --allow-extra-chr --make-bed --out "$qc_prefix" --noweb
plink --bfile "$qc_prefix" --indep-pairwise 50 5 0.2 --allow-extra-chr --out "$prune_prefix" --noweb
plink --bfile "$qc_prefix" --extract "${prune_prefix}.prune.in" --allow-extra-chr --make-bed --out "$pruned_prefix" --noweb

local_prefix=$(basename "$pruned_prefix")

for k in $(seq 2 5); do
  for seed in 1 2 3; do
    admixture --cv=10 -s "$seed" "${pruned_prefix}.bed" "$k" -j8 \
      | tee "${out_dir}/admixture.K${k}.seed${seed}.log"

    mv "${out_dir}/${local_prefix}.${k}.Q" "${out_dir}/admixture.K${k}.seed${seed}.Q"
    mv "${out_dir}/${local_prefix}.${k}.P" "${out_dir}/admixture.K${k}.seed${seed}.P"
  done
done

# Extract CV errors
printf "K,seed,CV_error\n" > "${out_dir}/admixture.cv_error.csv"

for log in "${out_dir}"/admixture.K*.seed*.log; do
  k=$(basename "$log" | sed -E 's/admixture\.K([0-9]+)\.seed([0-9]+)\.log/\1/')
  seed=$(basename "$log" | sed -E 's/admixture\.K([0-9]+)\.seed([0-9]+)\.log/\2/')
  cv=$(grep "CV error" "$log" | sed -E 's/.*CV error \(K=[0-9]+\): ([0-9.eE+-]+)/\1/')
  printf "%s,%s,%s\n" "$k" "$seed" "$cv" >> "${out_dir}/admixture.cv_error.csv"
done

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

