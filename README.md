# Carrier statistic
Carrier statistic is a statistical framework to prioritize disease-related rare variants by integrating gene expression data.

### Step1. Computing carrier statistic

```
Rscript /step1_carrier_stat.R \
--genotype=GENOTYPE_PREFIX \
--variants=VARIANTS_PREFIX \
--rna=RNA_PREFIX \
--gene=GENE_FILE \
--variants_gene_pair=VARIANTS_GENE_PAIR_FILE \
--outfile=OUTFILE_PREFIX
```
where the inputs are

* `GENOTYPE_PREFIX` (required): The prefix for genotype files. This prefix should correspond to `GENOTYPE_PREFIX_case.txt` for case group and `GENOTYPE_PREFIX_ctrl.txt` for control group.
* `VARIANTS_PREFIX` (required): The prefix for variant information files accompanying the genotype files. This prefix should correspond to `VARIANTS_PREFIX_case.txt` for case group and `VARIANTS_PREFIX_ctrl.txt` for control group.
* `RNA_PREFIX` (required): The prefix for gene expression data files. This prefix should correspond to `RNA_PREFIX_case.txt` for case group and `RNA_PREFIX_ctrl.txt` for control group.
* `GENE_FILE` (required): The full path to the gene information file accompanying the gene expression data files.
* `VARIANTS_GENE_PAIR_FILE` (required): The full path to the variant-gene-pair information file.
* `OUTFILE_PREFIX` (required): The prefix for output carrier statistic files. Two files will be generated, `OUTFILE_PREFIX_case.txt` for case group and `OUTFILE_PREFIX_ctrl.txt` for control group.

### A concrete example
```bash
#!/bin/bash
## ---- Download the required reference panel and example data for LOGODetect ---- ##
wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/LOGODetect/LOGODetect_data.tar.gz
tar -zxvf LOGODetect_data.tar.gz
rm -rf LOGODetect_data.tar.gz


## ---- Applying LOGODetect ---- ##
cd LOGODetect_data

mkdir ./results

conda activate ldsc

Rscript /LOGODetect/LOGODetect.R \
--sumstats ./sumstats/BIP.txt,./sumstats/SCZ.txt \
--n_gwas 51710,105318 \
--ref_dir ./LOGODetect_1kg_ref \
--pop EUR \
--ldsc_dir /LOGODetect/ldsc \
--block_partition /LOGODetect/block_partition.txt \
--out_dir ./results \
--n_cores 25

conda deactivate

```

### Output
LOGODetect outputs a whitespace-delimited text file `LOGODetect_regions.txt` in `PATH_TO_OUTFILE` specified by the user, with each row representing one small segment and the columns as such:
* `chr`: The chromosome. 
* `begin_pos`: The starting position of this detected small region.
* `stop_pos`: The stopping position of this detected small region.
* `stat`: The scan statistic value. Positive value means positive local genetic correlation between two traits. 
* `pval`: The p-value of this detected small region.
* `pval_adj`: The adjusted p-value of this detected small region.

We have prepared the example output file for you in `/LOGODetect_data/results/LOGODetect_regions.txt`. 

