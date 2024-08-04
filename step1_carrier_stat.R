library(data.table)
library(optparse)
options(stringsAsFactors = F)
option_list = list(
  make_option('--genotype', action = 'store', default = NA, type = 'character'), 
  make_option('--variants', action = 'store', default = NA, type = 'character'), 
  make_option('--rna', action = 'store', default = NA, type = 'character'), 
  make_option('--gene', action = 'store', default = NA, type = 'character'), 
  make_option('--variants_gene_pair', action = 'store', default = NA, type = 'character'), 
  make_option('--outfile', action = 'store', default = NA, type = 'character')
)
opt = parse_args(OptionParser(option_list=option_list))
genotype_prefix = opt$genotype
variants_prefix = opt$variants
rna_prefix = opt$rna
gene_file = opt$gene
variants_gene_pair_file = opt$variants_gene_pair
outfile_prefix = opt$outfile



pseudo_count = 1
AC_thre = 5
cal_z = function(noncarrier_vec, target_vec){
  mean0 = mean(noncarrier_vec)
  sd0 = sd(noncarrier_vec)
  if(sd0 > 0){
    target_vec = (target_vec - mean0)/sd0
  }else{
    target_vec = 0
  }
  return(mean(target_vec))
}

GROUP = c('case', 'ctrl')
for(group in GROUP){
  gene_info = data.frame(fread(gene_file))
  variants_gene_pair = data.frame(fread(variants_gene_pair_file))
  rna = as.matrix(data.frame(fread(paste0(rna_prefix, '_', group, '.txt'))))
  rna_log = log10(rna + pseudo_count)
  
  variants = data.frame(fread(paste0(variants_prefix, '_', group, '.txt')))
  genotype = as.matrix(data.frame(fread(paste0(genotype_prefix, '_', group, '.txt'))))
  colnames(genotype) = variants$ID
  variants$AC = apply(genotype, 2, function(x){return(sum(x > 0))})
  ind = (variants$AC <= AC_thre) & (variants$AC > 0)
  variants = variants[ind, ]
  genotype = genotype[, ind]
  variants_gene_pair = variants_gene_pair[variants_gene_pair$ID %in% variants$ID, ]
  
  result = variants_gene_pair
  result$n_carrier = 0
  result$carrier_stat = 0
  for(k in 1:nrow(result)){
    variant_id = result$ID[k]
    ind1 = genotype[, variant_id] >= 1
    ind2 = genotype[, variant_id] == 0
    result$n_carrier[k] = sum(ind1)
    
    rna_vec = rna_log[gene_info$GENE == result$GENE[k], ]
    carrier_vec = rna_vec[ind1]
    noncarrier_vec = rna_vec[ind2]
    result$carrier_stat[k] = cal_z(noncarrier_vec, carrier_vec)
  }
  write.table(result, paste0(outfile_prefix, '_', group, '.txt'), row.names = F, col.names = T, quote = F, sep = '\t')
}

