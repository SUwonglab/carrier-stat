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
Genotype file (`GENOTYPE_PREFIX_case.txt` and `GENOTYPE_PREFIX_ctrl.txt`): Allelic dosage file (number of ALT alleles, only 0/1/2 are supported) without a header line, one row per sample and one column per variant. The number of columns (i.e., the number of variants) must be equal to the number of rows in the variant information file. The number of rows (i.e., the number of samples) must be equal to the number of columns in the gene expression data file. 

Variant information file (`VARIANTS_PREFIX_case.txt` and `VARIANTS_PREFIX_case.txt`): A text file with a header line (CHROM: chromosome; POS: position; ID: variant name; REF: reference allele; ALT: alternative allele). The number of rows (i.e., the number of variants) must be equal to the number of columns in the genotype file. 

Gene expression data file (`RNA_PREFIX_case.txt` and `RNA_PREFIX_ctrl.txt`): RNA reads count file without a header line, one row per gene and one column per sample. The number of columns (i.e., the number of samples) must be equal to the number of rows in the genotype file. The number of rows (i.e., the number of genes) must be equal to the number of rows in the gene information file. 

Gene information file (`GENE_FILE`): A text file with a header line (CHROM: chromosome; MINBP: start position of the gene; MAXBP: end position of the gene; GENE: gene name). The number of rows (i.e., the number of genes) must be equal to the number of rows in the gene expression data file.

Variant-gene-pair information file (`VARIANTS_GENE_PAIR_FILE`): A text file with a header line (CHROM: chromosome; POS: position; ID: variant name; REF: reference allele; ALT: alternative allele; GENE: gene name). 

### Output format
Output carrier statistic file (`OUTFILE_PREFIX_case.txt` and `OUTFILE_PREFIX_ctrl.txt`): A text file with a header line (CHROM: chromosome; POS: position; ID: variant name; REF: reference allele; ALT: alternative allele; GENE: gene name; n_carrier: number of samples carrying the variant; carrier_stat: carrier statistic value). 


