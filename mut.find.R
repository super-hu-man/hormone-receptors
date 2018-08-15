mut = read.csv('mut.tag',header = F,stringsAsFactors = F)[,1]

files = paste('./cli',list.files('./cli','na'),sep ='/')
put = read.csv(files[1],header = T,sep='\t',stringsAsFactors = F)
share = intersect(colnames(put),mut)
dat = put[,share]
for (f in files[-1]){
  put = read.csv(f,header = T,sep='\t',stringsAsFactors = F)
  share = intersect(colnames(put),mut)
  da = put[,share]
  dat = merge(dat, da ,by= intersect(colnames(dat),colnames(da)),all = T)
}

write.table(dat, 'mut.find.txt',col.names = T,row.names = F,sep='\t',quote = F)
