# Simulate data

```
Rscript simulate_data.R \
--gene=GENE_FILE \
--variants=VARIANTS_FILE \
--variants_gene_pair=VARIANTS_GENE_PAIR_FILE \
--n_gene=N_GENE \
--n_variant=N_VARIANT \
--p_carrier=P_CARRIER \
--p_noncarrier=P_NONCARRIER \
--z=Z \
--n_case=N_CASE \
--n_ctrl=N_CTRL \
--outfile=OUTFILE_PREFIX

```
where the inputs are

* `GENOTYPE_PREFIX` (required): The prefix for genotype files. This prefix should correspond to `GENOTYPE_PREFIX_case.txt` for case group and `GENOTYPE_PREFIX_ctrl.txt` for control group.
* `VARIANTS_PREFIX` (required): The prefix for variant information files accompanying the genotype files. This prefix should correspond to `VARIANTS_PREFIX_case.txt` for case group and `VARIANTS_PREFIX_ctrl.txt` for control group.
* `RNA_PREFIX` (required): The prefix for gene expression data files. This prefix should correspond to `RNA_PREFIX_case.txt` for case group and `RNA_PREFIX_ctrl.txt` for control group.

* `GENE_FILE` (required): The full path to the gene information file.
* `VARIANTS_FILE` (required): The full path to the variants information file.
* `VARIANTS_GENE_PAIR_FILE` (required): The full path to the variants-gene pair information file.
* `N_GENE` (required): The number of causal genes.
* `N_VARIANT` (required): The number of causal variants per causal gene.
* `P_CARRIER` (required): The penetrance of causal variants.
* `P_NONCARRIER` (required): The prevalence in causal variants noncarriers.
* `Z` (required): Gene expression perturbation level for causal variants.
* `N_CASE` (required): The number of simulated cases.
* `N_CTRL` (required): The number of simulated controls.
* `OUTFILE_PREFIX` (required): The folder name for simulated data. 

### A concrete example
```
cd carrier-stat/simulation

Rscript ./simulate_data.R \
--gene=./gene_all.txt \
--variants=./variants_all.txt \
--variants_gene_pair=./variants_gene_pair.txt \
--n_gene=10 \
--n_variant=5 \
--p_carrier=0.7 \
--p_noncarrier=0.01 \
--z=5 \
--n_case=500 \
--n_ctrl=500 \
--outfile=./sim
```

### Input format
Genotype file (`GENOTYPE_PREFIX_case.txt` and `GENOTYPE_PREFIX_ctrl.txt`): Allelic dosage file (number of ALT alleles, only 0/1/2 are supported) without a header line, one row per sample and one column per variant. The number of columns (i.e., the number of variants) must be equal to the number of rows in the variant information file. The number of rows (i.e., the number of samples) must be equal to the number of columns in the gene expression data file. 

Variant information file (`VARIANTS_PREFIX_case.txt` and `VARIANTS_PREFIX_case.txt`): A text file with a header line (CHROM: chromosome; POS: position; ID: variant name; REF: reference allele; ALT: alternative allele). The number of rows (i.e., the number of variants) must be equal to the number of columns in the genotype file. 

Gene expression data file (`RNA_PREFIX_case.txt` and `RNA_PREFIX_ctrl.txt`): RNA reads count file without a header line, one row per gene and one column per sample. The number of columns (i.e., the number of samples) must be equal to the number of rows in the genotype file. The number of rows (i.e., the number of genes) must be equal to the number of rows in the gene information file. 

Gene information file (`GENE_FILE`): A text file with a header line (CHROM: chromosome; MINBP: start position of the gene; MAXBP: end position of the gene; GENE: gene name). The number of rows (i.e., the number of genes) must be equal to the number of rows in the gene expression data file.

Variant-gene-pair information file (`VARIANTS_GENE_PAIR_FILE`): A text file with a header line (CHROM: chromosome; POS: position; ID: variant name; REF: reference allele; ALT: alternative allele; GENE: gene name). 

### Output format
Output carrier statistic file (`OUTFILE_PREFIX_case.txt` and `OUTFILE_PREFIX_ctrl.txt`): A text file with a header line (CHROM: chromosome; POS: position; ID: variant name; REF: reference allele; ALT: alternative allele; GENE: gene name; n_carrier: number of samples carrying the variant; carrier_stat: carrier statistic value). 



## Step2. Prioritize rare variant-gene pairs with extreme carrier statistic
```
Rscript step2_analysis.R \
--carrier_stat=CARRIER_STAT_PREFIX \
--outfile=OUTFILE_PREFIX \
--fdr_thre=FDR_THRESHOLD
```
where the inputs are

* `CARRIER_STAT_PREFIX` (required): The prefix for carrier statistic files output from Step 1. This prefix should correspond to `CARRIER_STAT_PREFIX_case.txt` for case group and `CARRIER_STAT_PREFIX_ctrl.txt` for control group.
* `OUTFILE_PREFIX` (required): The prefix for files containing significant variant-gene pairs. Two files will be generated, `OUTFILE_PREFIX_downregulated_fdr_FDR_THRESHOLD.txt` for significant variant-gene pairs with negative carrier statistics and `OUTFILE_PREFIX_upregulated_fdr_FDR_THRESHOLD.txt` for significant variant-gene pairs with positive carrier statistics.
* `FDR_THRESHOLD` (optional): FDR cutoff. Default is 0.05.

### A concrete example
```
cd carrier-stat

Rscript ./step2_analysis.R \
--carrier_stat=./example/carrier_stat \
--outfile=./example/carrier_stat \
--fdr_thre=0.2
```

### Output format
Output files containing significant variant-gene pairs (`OUTFILE_PREFIX_downregulated_fdr_FDR_THRESHOLD.txt` and `OUTFILE_PREFIX_upregulated_fdr_FDR_THRESHOLD.txt`): A text file with a header line (CHROM: chromosome; POS: position; ID: variant name; REF: reference allele; ALT: alternative allele; GENE: gene name; n_carrier: number of samples carrying the variant; carrier_stat: carrier statistic value; fdr: FDR). 




