files = list.files('./result/Meth.correlation','cor.mRNA_Meth')

for ( f in files) {
  dat = read.csv(f,header=T,sep='\t',stringsAsFactors=F)
  
}
