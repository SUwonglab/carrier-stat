library(data.table)
library(optparse)
options(stringsAsFactors = F)
option_list = list(
  make_option('--gene', action = 'store', default = NA, type = 'character'), 
  make_option('--variants', action = 'store', default = NA, type = 'character'), 
  make_option('--variants_gene_pair', action = 'store', default = NA, type = 'character'), 
  make_option('--n_gene', action = 'store', default = NA, type = 'integer'),
  make_option('--n_variant', action = 'store', default = NA, type = 'integer'),
  make_option('--p_carrier', action = 'store', default = NA, type = 'numeric'),
  make_option('--p_noncarrier', action = 'store', default = NA, type = 'numeric'),
  make_option('--z', action = 'store', default = NA, type = 'numeric'),
  make_option('--prop_pos', action = 'store', default = 1, type = 'numeric'),
  make_option('--n', action = 'store', default = 125748, type = 'integer'),
  make_option('--n_case', action = 'store', default = NA, type = 'integer'),
  make_option('--n_ctrl', action = 'store', default = NA, type = 'integer'),
  make_option('--outfile', action = 'store', default = NA, type = 'character')
)
opt = parse_args(OptionParser(option_list=option_list))
gene_info = data.frame(fread(opt$gene))
variants = data.frame(fread(opt$variants))
variants_gene_pair = data.frame(fread(opt$variants_gene_pair))
n_gene = opt$n_gene
n_variant = opt$n_variant
p_carrier = opt$p_carrier
p_noncarrier = opt$p_noncarrier
z = opt$z
prop_pos = opt$prop_pos
n = opt$n
n_case = opt$n_case
n_ctrl = opt$n_ctrl
outfile_prefix = opt$outfile
IID = 1:n
if(!dir.exists(outfile_prefix)){
  dir.create(outfile_prefix)
}
set.seed(123)

write.table(gene_info[, 1:4], paste0(outfile_prefix, '/gene.txt'), row.names = F, col.names = T, quote = F, sep = '\t')
gene_ind = sort(sample(1:nrow(gene_info), n_gene, replace = F))
gene_info_sim = gene_info[gene_ind, ]
gene_info_sim$effect = -1
pos_ind = sort(sample(1:n_gene, ceiling(n_gene*prop_pos), replace = F))
gene_info_sim$effect[pos_ind] = 1
write.table(gene_info_sim, paste0(outfile_prefix, '/gene_causal.txt'), row.names = F, col.names = T, quote = F, sep = '\t')

variant_info_sim = c()
for(j in 1:n_gene){
  gene = gene_info_sim$GENE[j]
  variants_gene_pair_tmp = variants_gene_pair[variants_gene_pair$GENE == gene, ]
  if(nrow(variants_gene_pair_tmp) > 0){
    ind = sort(sample(1:nrow(variants_gene_pair_tmp), min(nrow(variants_gene_pair_tmp), n_variant), replace = F))
    variant_info_sim_new = variants_gene_pair_tmp[ind, ]
    variant_info_sim = rbind(variant_info_sim, variant_info_sim_new)
  }
}

carrier_ind = list()
carrier_ind_all = c()
affected_ind_carrier = c()
for(k in 1:nrow(variant_info_sim)){
  id = variant_info_sim$ID[k]
  ac = variants$AC[variants$ID == id]
  carrier_ind[[k]] = sort(sample(1:n, ac, replace = F))
  carrier_ind_all = c(carrier_ind_all, carrier_ind[[k]])
  is.affected = rbinom(ac, 1, prob = p_carrier)
  affected_ind_carrier = c(affected_ind_carrier, carrier_ind[[k]][is.affected == 1])
}
carrier_ind_all = sort(unique(carrier_ind_all))
noncarrier_ind_all = setdiff(1:n, carrier_ind_all)
is.affected = rbinom(length(noncarrier_ind_all), 1, prob = p_noncarrier)
affected_ind_noncarrier = noncarrier_ind_all[is.affected == 1]
affected_ind = sort(c(affected_ind_carrier, affected_ind_noncarrier))
nonaffected_ind = setdiff(1:n, affected_ind)
case_ind = sort(sample(affected_ind, n_case, replace = F))
ctrl_ind = sort(sample(nonaffected_ind, n_ctrl, replace = F))

genotype_case_causal = matrix(0, nrow = n_case, ncol = nrow(variant_info_sim))
genotype_ctrl_causal = matrix(0, nrow = n_ctrl, ncol = nrow(variant_info_sim))
for(k in 1:nrow(variant_info_sim)){
  ind1 = case_ind %in% carrier_ind[[k]]
  ind2 = ctrl_ind %in% carrier_ind[[k]]
  if(sum(ind1) > 0){
    genotype_case_causal[ind1, k] = 1
  }
  if(sum(ind2) > 0){
    genotype_ctrl_causal[ind2, k] = 1
  }
}
ac = apply(genotype_case_causal, 2, sum) + apply(genotype_ctrl_causal, 2, sum)
variant_info_sim = variant_info_sim[ac > 0, c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'GENE')]
genotype_case_causal = genotype_case_causal[, ac > 0]
genotype_ctrl_causal = genotype_ctrl_causal[, ac > 0]
variant_info_sim$effect = 1
for(j in 1:n_gene){
  variant_info_sim$effect[variant_info_sim$GENE == gene_info_sim$GENE[j]] = gene_info_sim$effect[j]
}
write.table(variant_info_sim, paste0(outfile_prefix, '/variants_causal.txt'), row.names = F, col.names = T, quote = F, sep = '\t')
write.table(genotype_case_causal, paste0(outfile_prefix, '/genotype_case_causal.txt'), row.names = F, col.names = F, quote = F, sep = '\t')
write.table(genotype_ctrl_causal, paste0(outfile_prefix, '/genotype_ctrl_causal.txt'), row.names = F, col.names = F, quote = F, sep = '\t')



rna_case = matrix(0, nrow = nrow(gene_info), ncol = n_case)
rna_ctrl = matrix(0, nrow = nrow(gene_info), ncol = n_ctrl)
colnames(rna_case) = paste0('case', 1:n_case)
colnames(rna_ctrl) = paste0('ctrl', 1:n_ctrl)
for(j in 1:nrow(gene_info)){
  rna_case[j, ] = rnorm(n_case) * gene_info$sd_log2_expr[j] + gene_info$mean_log2_expr[j]
  rna_ctrl[j, ] = rnorm(n_ctrl) * gene_info$sd_log2_expr[j] + gene_info$mean_log2_expr[j]
}
for(j in 1:nrow(gene_info_sim)){
  gene_ind = which(gene_info$GENE == gene_info_sim$GENE[j])
  variant_ind = which(variant_info_sim$GENE == gene_info_sim$GENE[j])
  if(length(variant_ind) > 0){
    tmp_case = as.matrix(genotype_case_causal[, variant_ind])
    tmp_ctrl = as.matrix(genotype_ctrl_causal[, variant_ind])
    z_case = apply(tmp_case, 1, sum) * z * gene_info_sim$effect[j]
    z_ctrl = apply(tmp_ctrl, 1, sum) * z * gene_info_sim$effect[j]
    rna_case[gene_ind, ] = rna_case[gene_ind, ] + z_case * gene_info_sim$sd_log2_expr[j]
    rna_ctrl[gene_ind, ] = rna_ctrl[gene_ind, ] + z_ctrl * gene_info_sim$sd_log2_expr[j]
  }
}
rna_case = floor(2^rna_case)
rna_ctrl = floor(2^rna_ctrl)
write.table(rna_case, paste0(outfile_prefix, '/rna_case.txt'), row.names = F, col.names = F, quote = F, sep = '\t')
write.table(rna_ctrl, paste0(outfile_prefix, '/rna_ctrl.txt'), row.names = F, col.names = F, quote = F, sep = '\t')



genotype_case = matrix(0, nrow = n_case, ncol = nrow(variants))
genotype_ctrl = matrix(0, nrow = n_ctrl, ncol = nrow(variants))
for(k in 1:nrow(variants)){
  ind = sample(1:n, size = variants$AC[k], replace = F)
  genotype_case[case_ind %in% ind, k] = 1
  genotype_ctrl[ctrl_ind %in% ind, k] = 1
}

for(k in 1:nrow(variant_info_sim)){
  ind = which(variants$ID == variant_info_sim$ID[k])
  genotype_case[, ind] = genotype_case_causal[, k]
  genotype_ctrl[, ind] = genotype_ctrl_causal[, k]
}

ac = apply(genotype_case, 2, sum) + apply(genotype_ctrl, 2, sum)
variants = variants[ac > 0, c('CHROM', 'POS', 'ID', 'REF', 'ALT')]
variants_gene_pair = variants_gene_pair[variants_gene_pair$ID %in% variants$ID, c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'GENE')]
genotype_case = genotype_case[, ac > 0]
genotype_ctrl = genotype_ctrl[, ac > 0]
write.table(variants_gene_pair, paste0(outfile_prefix, '/variants_gene_pair.txt'), row.names = F, col.names = T, quote = F, sep = '\t')

ac_case = apply(genotype_case, 2, sum)
genotype_case = genotype_case[, ac_case > 0]
variants_case = variants[ac_case > 0, ]
write.table(genotype_case, paste0(outfile_prefix, '/genotype_case.txt'), row.names = F, col.names = F, quote = F, sep = '\t')
write.table(variants_case, paste0(outfile_prefix, '/variants_case.txt'), row.names = F, col.names = T, quote = F, sep = '\t')

ac_ctrl = apply(genotype_ctrl, 2, sum)
genotype_ctrl = genotype_ctrl[, ac_ctrl > 0]
variants_ctrl = variants[ac_ctrl > 0, ]
write.table(genotype_ctrl, paste0(outfile_prefix, '/genotype_ctrl.txt'), row.names = F, col.names = F, quote = F, sep = '\t')
write.table(variants_ctrl, paste0(outfile_prefix, '/variants_ctrl.txt'), row.names = F, col.names = T, quote = F, sep = '\t')






