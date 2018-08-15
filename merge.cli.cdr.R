files = paste('./cli_filted/',list.files('./cli_filted/','txt'),sep='')
cdr = read.csv('TCGA-CDR.txt',header = T,sep='\t',stringsAsFactors = F)
mut = read.csv('mut.find.adjust.txt',header = T,row.names=1,sep='\t',stringsAsFactors = F)
rownames(cdr) = cdr[,1]
projs = sort(unique(cdr$type))

for(i in projs){
  p.cdr = cdr[cdr$type==i,]
  p.mut = mut[rownames(p.cdr),]
  p.tag = read.csv(paste('./cli_filted/',tolower(i),'.txt',sep=''),header = T,sep='\t',stringsAsFactors = F)$Tag
  p.cli = read.csv(paste('./cli/nationwidechildrens.org_clinical_patient_',tolower(i),'.txt',sep=''),header=T,sep='\t',stringsAsFactors=F)
  p.cli = p.cli[-c(1:2),-1] ; rownames(p.cli) = p.cli$bcr_patient_barcode
  p.cli = p.cli[,c('bcr_patient_barcode',p.tag)]
  
  for (c in colnames(p.cdr)){
    if (length(table(p.cdr[,c]))<2){p.cdr = p.cdr[,-which(c==colnames(p.cdr))]}
  }
  for (t in colnames(p.mut)){
    if (length(table(p.mut[,t]))<2){p.mut = p.mut[,-which(t==colnames(p.mut))]}
    if(length(colnames(p.mut))>1){p.mut = p.mut[p.cdr$bcr_patient_barcode,]}
  }
  for (s in colnames(p.cli)){
    #if (length(table(p.cli[,s]))<2){print(table(p.cli[,s]));p.cli = p.cli[,-which(s==colnames(p.cli))]}
    if (length(colnames(p.cli))>1){p.cli = p.cli[p.cli$bcr_patient_barcode,]}
  }
  #one = merge(p.cdr, p.mut[p.cdr$bcr_patient_barcode,],p.cli[p.cdr$bcr_patient_barcode,], by='bcr_patient_barcode',all=T)
  one = cbind(p.cdr, p.mut, p.cli)
  write.table(one,paste('./cli_filted/',i,'.merge.txt',sep=''),col.names=T,row.names=F,sep='\t',quote=F)
}