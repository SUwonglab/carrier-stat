library(data.table)
library(optparse)
options(stringsAsFactors = F)
option_list = list(
  make_option('--carrier_stat', action = 'store', default = NA, type = 'character'), 
  make_option('--outfile', action = 'store', default = NA, type = 'character')
)
opt = parse_args(OptionParser(option_list=option_list))
carrier_stat_prefix = opt$carrier_stat
outfile_prefix = opt$outfile
fdr_thre = 0.2

stat_case = data.frame(fread(paste0(carrier_stat_prefix, '_case.txt')))
stat_case = stat_case[order(stat_case$carrier_stat), ]
stat_ctrl = data.frame(fread(paste0(carrier_stat_prefix, '_ctrl.txt')))
stat_ctrl = stat_ctrl[order(stat_ctrl$carrier_stat), ]
dat1 = data.frame(unique(stat_case$carrier_stat))
dat2 = data.frame(unique(stat_ctrl$carrier_stat))
colnames(dat1) = colnames(dat2) = 'z'
dat1$type = 'case'
dat2$type = 'ctrl'
dat_pl = rbind(dat1, dat2)
z1 = sort(dat1$z)
z2 = sort(dat2$z)
z1_d = rev(z1)
z2_d = rev(z2)
n1 = nrow(dat1)
n2 = nrow(dat2)

z_thre_left = max(z1[1:min(length(z1), 1000)], z2[1:min(length(z2), 1000)])
z_thre_right = min(z1_d[1:min(length(z1), 1000)], z2_d[1:min(length(z2), 1000)])
z1_left = z1[z1 <= z_thre_left]
z2_left = z2[z2 <= z_thre_left]
z1_right = z1_d[z1_d >= z_thre_right]
z2_right = z2_d[z2_d >= z_thre_right]
re_left = data.frame(sort(unique(c(z1_left, z2_left))))
re_right = data.frame(sort(unique(c(z1_right, z2_right)), decreasing = T))
colnames(re_left) = colnames(re_right) = 'z0'
re_left$fdr = Inf
re_right$fdr = Inf
for(i in 1:nrow(re_left)){
  z0 = re_left$z0[i]
  re_left$fdr[i] = sum(z2_left <= z0)/sum(z1_left <= z0)/n2*n1
}
for(i in 1:nrow(re_right)){
  z0 = re_right$z0[i]
  re_right$fdr[i] = sum(z2_right >= z0)/sum(z1_right >= z0)/n2*n1
}

if(sum(re_left$fdr <= fdr_thre) > 0){
  z0_thre = re_left$z0[max(which(re_left$fdr <= fdr_thre))]
  tmp = stat_case[stat_case$carrier_stat <= z0_thre, ]
  tmp$fdr = 1
  for(i in 1:nrow(tmp)){
    tmp$fdr[i] = re_left$fdr[min(which(re_left$z0 >= tmp$carrier_stat[i]))]
  }
  for(i in nrow(tmp):1){
    tmp$fdr[i] = min(tmp$fdr[i:nrow(tmp)])
  }
  write.table(tmp, paste0(outfile_prefix, '_downregulated_fdr_', fdr_thre, '.txt'), row.names = F, col.names = T, quote = F, sep = '\t')
}
if(sum(re_right$fdr <= fdr_thre) > 0){
  z0_thre = re_right$z0[max(which(re_right$fdr <= fdr_thre))]
  tmp = stat_case[stat_case$carrier_stat >= z0_thre, ]
  tmp = tmp[order(tmp$carrier_stat, decreasing = T), ]
  tmp$fdr = 1
  for(i in 1:nrow(tmp)){
    tmp$fdr[i] = re_right$fdr[min(which(re_right$z0 <= tmp$carrier_stat[i]))]
  }
  for(i in nrow(tmp):1){
    tmp$fdr[i] = min(tmp$fdr[i:nrow(tmp)])
  }
  write.table(tmp, paste0(outfile_prefix, '_upregulated_fdr_', fdr_thre, '.txt'), row.names = F, col.names = T, quote = F, sep = '\t')
}





