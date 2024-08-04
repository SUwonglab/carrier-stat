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

### Output
Gene information file: `OUTFILE_PREFIX/gene.txt`.
Gene information file for causal genes: `OUTFILE_PREFIX/gene_causal.txt`.
Genotype files: `OUTFILE_PREFIX/genotype_case.txt` and `OUTFILE_PREFIX/genotype_ctrl.txt`.
Variant information file: `OUTFILE_PREFIX/variants_case.txt` and `OUTFILE_PREFIX/variants_case.txt`.
Genotype files restricted to causal variants: `OUTFILE_PREFIX/genotype_case_causal.txt` and `OUTFILE_PREFIX/genotype_ctrl_causal.txt`.
Variant information file restricted to causal variants: `OUTFILE_PREFIX/variants_causal.txt`.
Gene expression data file: `OUTFILE_PREFIX/rna_case.txt` and `OUTFILE_PREFIX/rna_ctrl.txt`.
Variant-gene pair information file: `OUTFILE_PREFIX/variants_gene_pair.txt`.
