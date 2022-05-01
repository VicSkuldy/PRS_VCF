#! /bin/bash
gunzip -c Height.gwas.txt.gz | head
md5sum Height.gwas.txt.gz

# Standard GWAS QC
gunzip -c Height.gwas.txt.gz |\
awk 'NR==1 || ($11 > 0.01) && ($10 > 0.8) {print}' |\
gzip  > Height.gz

# Duplicate SNPs
gunzip -c Height.gz |\
awk '{seen[$3]++; if(seen[$3]==1){ print}}' |\
gzip - > Height.nodup.gz

# Ambiguous SNPs
gunzip -c Height.nodup.gz |\
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' |\
    gzip > Height.QC.gz

unzip EUR.zip

# Standard GWAS QC
plink \
    --bfile EUR \
    --maf 0.01 \
    --hwe 1e-6 \
    --geno 0.01 \
    --mind 0.01 \
    --write-snplist \
    --make-just-fam \
    --out EUR.QC

plink \
    --bfile EUR \
    --keep EUR.QC.fam \
    --extract EUR.QC.snplist \
    --indep-pairwise 200 50 0.25 \
    --out EUR.QC

plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.fam \
    --het \
    --out EUR.QC

Rscript 1.r

Rscript 2.r

plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.valid.sample \
    --check-sex \
    --out EUR.QC

Rscript 3.r

plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.valid \
    --rel-cutoff 0.125 \
    --out EUR.QC

plink \
    --bfile EUR \
    --make-bed \
    --keep EUR.QC.rel.id \
    --out EUR.QC \
    --extract EUR.QC.snplist \
    --exclude EUR.mismatch \
    --a1-allele EUR.a1

Rscript 4.r

plink \
    --bfile EUR.QC \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump Height.QC.Transformed \
    --clump-snp-field SNP \
    --clump-field P \
    --out EUR

awk 'NR!=1{print $3}' EUR.clumped >  EUR.valid.snp

awk '{print $3,$8}' Height.QC.Transformed > SNP.pvalue

echo "0.001 0 0.001" > range_list
echo "0.05 0 0.05" >> range_list
echo "0.1 0 0.1" >> range_list
echo "0.2 0 0.2" >> range_list
echo "0.3 0 0.3" >> range_list
echo "0.4 0 0.4" >> range_list
echo "0.5 0 0.5" >> range_list

plink \
    --bfile EUR.QC \
    --score Height.QC.Transformed 3 4 12 header \
    --q-score-range range_list SNP.pvalue \
    --extract EUR.valid.snp \
    --out EUR

# First, we need to perform prunning
plink \
    --bfile EUR.QC \
    --indep-pairwise 200 50 0.25 \
    --out EUR
# Then we calculate the first 6 PCs
plink \
    --bfile EUR.QC \
    --extract EUR.prune.in \
    --pca 6 \
    --out EUR

Rscript 5.r

Rscript 6.r

Rscript PRSice.R \
    --prsice PRSice_linux \
    --base Height.QC.gz \
    --target EUR.QC \
    --binary-target F \
    --pheno EUR.height \
    --cov EUR.covariate \
    --base-maf MAF:0.01 \
    --base-info INFO:0.8 \
    --stat OR \
    --or \
    --out EUR
