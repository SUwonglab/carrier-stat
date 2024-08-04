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
```
cd carrier-stat
Rscript ./step1_carrier_stat.R \
--genotype=./example/genotype \
--variants=./example/variants \
--rna=./example/rna \
--gene=./example/gene.txt \
--variants_gene_pair=./example/variants_gene_pair.txt \
--outfile=./example/carrier_stat
```

### Input format




### Output format
LOGODetect outputs a whitespace-delimited text file `LOGODetect_regions.txt` in `PATH_TO_OUTFILE` specified by the user, with each row representing one small segment and the columns as such:
* `chr`: The chromosome. 
* `begin_pos`: The starting position of this detected small region.
* `stop_pos`: The stopping position of this detected small region.
* `stat`: The scan statistic value. Positive value means positive local genetic correlation between two traits. 
* `pval`: The p-value of this detected small region.
* `pval_adj`: The adjusted p-value of this detected small region.

We have prepared the example output file for you in `/LOGODetect_data/results/LOGODetect_regions.txt`. 

