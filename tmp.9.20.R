cli = read.csv('./files/clinical.literature.final_use.txt',header=T,sep='\t',stringsAsFactors=F)
clis= cli[,1:32]
cli = cli[,-c(10,12,13,15,16,17,19,20,21,33)]

proj = sort(unique(cli$type))

for (i in proj){
  print(i)
  rna = read.csv(paste('../download.cbioprotal/',tolower(i),'_tcga_pan_can_atlas_2018/data_expression.txt',sep=''),header=T,sep='\t',stringsAsFactors=F)
  samples = substr(gsub('[.]','-',colnames(rna)[3:ncol(rna)]),1,12)
  
  cancer = cli[cli$type==i & (cli$bcr_patient_barcode %in% samples),]
  sums = apply(cancer,2,function(x){sum(is.na(x))})
  cancer = cancer[,sums!=nrow(cancer)]
  out = data.frame()
  write.table(cancer,paste('./cli.from.TCGA-CDR/',i,'.filter-NA.txt',sep=''),row.names=F,col.names=T,sep='\t',quote=F)
}

x=read.csv('./result/two.group.survival.xls',header=T,row.names=1,sep='\t',stringsAsFactors=F)
