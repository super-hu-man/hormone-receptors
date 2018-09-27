
files = list.files('./result/Meth.correlation', 'cor.mRNA_Meth.txt',full.names = T)
GPL13534 = read.csv('C:/Users/root/Desktop/Meth_ID_match/GPL13534-11288.txt',header = T,row.names=1,stringsAsFactors = F,sep='\t')
colnames(GPL13534) = GPL13534[37,] ; GPL13534 = GPL13534[-c(1:37), ]

for (f in files ){
  name = strsplit(strsplit(f,'/',fixed=T)[[1]][4],'.',fixed=T)[[1]][1]
  dat0 = read.csv(f,sep='\t',header = T,row.names = 1,stringsAsFactors = F)
  
  dat1 = data.frame(Gene_name=GPL13534[rownames(dat0),'UCSC_RefGene_Name'] , Group=GPL13534[rownames(dat0),'UCSC_RefGene_Group'] ,
                    dat0 , stringsAsFactors = F)
  dat = data.frame(Meth_IDs=rownames(dat1),dat1)
  write.table(dat,f,col.names = T,row.names = F,sep='\t',quote = F)
  
  dat2 = dat1[(dat1$logE_M_pearson_q<=0.05) & (abs(dat1$logE_M_pearson_cor)>=0.3),]
  genes = unique(unlist(strsplit(dat2$Gene_name,';',fixed=T))) ; genes = genes[genes!='']
  write.table(genes,paste('./result/Meth.correlation/',name,'.filter.genes.txt',sep=''),col.names = F,row.names = F,sep='\t',quote = F)
  
  #out = data.frame(Meth_IDs=rownames(dat3),dat3)
  ##out$Group = unlist(lapply(out$Group,function(x) {y=unlist(strsplit(x,';',fixed=T));z=unique(y);x=paste(z,collapse=';');return(x)}))
  ##out$Gene_name = unlist(lapply(out$Gene_name,function(x) {y=unlist(strsplit(x,';',fixed=T));z=unique(y);x=paste(z,collapse=';');return(x)}))
  #write.table(out,paste('./result/Meth.correlation/',name,'.filter.rna_meth.cor.txt',sep=''),col.names = T,row.names = F,sep='\t',quote = F)
}

##
"
quit()
a = list.files( 'C:/Users/root/Desktop/Meth_ID_match' , full.names=T )
GPL13534 = read.csv(a[1],header = T,stringsAsFactors = F,sep='\t')
GPL16304 = read.csv(a[2],header = T,stringsAsFactors = F,sep='\t')
GPL18809 = read.csv(a[3],header = T,stringsAsFactors = F,sep='\t')
IDmeth = read.csv(a[4],header = T,stringsAsFactors = F,sep='\t')
GPLTCGA = read.csv(a[5],header = T,stringsAsFactors = F,sep='\t')
ACC = read.csv(a[6],header = T,stringsAsFactors = F,sep='\t')
"