library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
argvs = commandArgs(trailingOnly = T)
print(argvs)
argv=argvs[1]
n=as.numeric(argvs[2])
#library(RDAVIDWebService)

files = list.files('./result/diff',paste('deseq',argv,sep='.'),full.names = T)
out.path = paste('./result/clusterprofiler',argv,sep='/')
if(! dir.exists(out.path)){dir.create(out.path,recursive=T)}
print(out.path)
analysis = function(filted,dat,name,out.path){
  print('BP start')
  bp = enrichGO(gene = filted$ID , OrgDb = org.Hs.eg.db , universe = dat$ID , keyType = 'ENTREZID' , readable = TRUE ,
                ont = 'BP' , pAdjustMethod = 'BH' , pvalueCutoff = 0.05 , qvalueCutoff  = 0.05 )
  bps=simplify(bp)
  b = sum(duplicated(unlist(strsplit(as.data.frame(bps)[,8],'/',fixed=T))))
  if(nrow(as.data.frame(bps))>2){
  pdf(paste(out.path,'/',name,'.BP.dot.pdf',sep=''),width=12,height=8) ; dot=dotplot(bps,showCategory=50) ;plot(dot); dev.off()
  if (b >= 1){ pdf(paste(out.path,'/',name,'.BP.emap.pdf',sep=''),width=12,height=8) ; emap=emapplot(bps,showCategory=50) ;plot(emap); dev.off() }
  pdf(paste(out.path,'/',name,'.BP.cnet.pdf',sep=''),width=12,height=8) ; cnet=cnetplot(bps,showCategory=10) ;plot(cnet); dev.off()
  }
  print('CC start')
  cc = enrichGO(gene = filted$ID , OrgDb = org.Hs.eg.db , universe = dat$ID , keyType = 'ENTREZID' , readable = TRUE ,
                ont = 'CC' , pAdjustMethod = 'BH' , pvalueCutoff = 0.05 , qvalueCutoff  = 0.05 )
  ccs=simplify(cc) ; c = sum(duplicated(unlist(strsplit(as.data.frame(ccs)[,8],'/',fixed=T))))
  if(nrow(as.data.frame(ccs))>2){
  pdf(paste(out.path,'/',name,'.CC.dot.pdf',sep=''),width=12,height=8) ; dot=dotplot(ccs,showCategory=50) ; plot(dot);dev.off()
  if (c >= 1){pdf(paste(out.path,'/',name,'.CC.emap.pdf',sep=''),width=12,height=8) ; emap=emapplot(ccs,showCategory=50) ;plot(emap); dev.off()}
  pdf(paste(out.path,'/',name,'.CC.cnet.pdf',sep=''),width=12,height=8) ; cnet=cnetplot(ccs,showCategory=10) ;plot(cnet); dev.off()
  }
  print('MF start')
  mf = enrichGO(gene = filted$ID , OrgDb = org.Hs.eg.db , universe = dat$ID , keyType = 'ENTREZID' , readable = TRUE ,
                ont = 'MF' , pAdjustMethod = 'BH' , pvalueCutoff = 0.05 , qvalueCutoff  = 0.05 )
  mfs=simplify(mf) ; m = sum(duplicated(unlist(strsplit(as.data.frame(mfs)[,8],'/',fixed=T))))
  if(nrow(as.data.frame(mfs))>2){
  pdf(paste(out.path,'/',name,'.MF.dot.pdf',sep=''),width=12,height=8) ; dot=dotplot(mfs,showCategory=50) ;plot(dot); dev.off()
  if (m >= 1){pdf(paste(out.path,'/',name,'.MF.emap.pdf',sep=''),width=12,height=8) ; emap=emapplot(mfs,showCategory=50) ;plot(emap); dev.off()}
  pdf(paste(out.path,'/',name,'.MF.cnet.pdf',sep=''),width=12,height=8) ; cnet=cnetplot(mfs,showCategory=10) ;plot(cnet); dev.off()
  }
  print('KEGG start')
  kegg = enrichKEGG(gene = filted$ID , organism = 'hsa' , universe = dat$ID ,
                    pvalueCutoff = 0.05 , qvalueCutoff = 0.05 , use_internal_data=F )
  k = sum(duplicated(unlist(strsplit(as.data.frame(kegg)[,8],'/',fixed=T))))
  if(nrow(as.data.frame(kegg))>2){
  pdf(paste(out.path,'/',name,'.KEGG.dot.pdf',sep=''),width=12,height=8) ; dot=dotplot(kegg,showCategory=50) ;plot(dot); dev.off()
  if (k >= 1){ pdf(paste(out.path,'/',name,'.KEGG.emap.pdf',sep=''),width=12,height=8) ; emap=emapplot(kegg,showCategory=50) ;plot(emap); dev.off()}
  pdf(paste(out.path,'/',name,'.KEGG.cnet.pdf',sep=''),width=12,height=8) ; cnet=cnetplot(kegg,showCategory=10) ;plot(cnet); dev.off()
  }
  write.table(as.data.frame(bp),paste(out.path,'/',paste(name,argv,sep='_'),'.GO.BP.txt',sep=''),sep='\t',quote=F,col.names=T,row.names=F)
  write.table(as.data.frame(cc),paste(out.path,'/',paste(name,argv,sep='_'),'.GO.CC.txt',sep=''),sep='\t',quote=F,col.names=T,row.names=F)
  write.table(as.data.frame(mf),paste(out.path,'/',paste(name,argv,sep='_'),'.GO.MF.txt',sep=''),sep='\t',quote=F,col.names=T,row.names=F)
  write.table(as.data.frame(kegg),paste(out.path,'/',paste(name,argv,sep='_'),'.KEGG.txt',sep=''),sep='\t',quote=F,col.names=T,row.names=F)
  write.table(as.data.frame(bps),paste(out.path,'/',paste(name,argv,sep='_'),'.GO.BP.simplify.txt',sep=''),sep='\t',quote=F,col.names=T,row.names=F)
  write.table(as.data.frame(ccs),paste(out.path,'/',paste(name,argv,sep='_'),'.GO.CC.simplify.txt',sep=''),sep='\t',quote=F,col.names=T,row.names=F)
  write.table(as.data.frame(mfs),paste(out.path,'/',paste(name,argv,sep='_'),'.GO.MF.simplify.txt',sep=''),sep='\t',quote=F,col.names=T,row.names=F)
  #write.table(as.data.frame(keggs),paste(out.path,'/',paste(name,argv,sep='_'),'.KEGG.simplify.txt',sep=''),sep='\t',quote=F,col.names=T,row.names=F)
}

for(f in files[n]){
  name = strsplit(strsplit(f,'/',fixed=T)[[1]][4],'.',fixed=T)[[1]][1]
  print(name)
  dat = read.csv(f,header = F,row.names = 1, sep='\t',stringsAsFactors = F) 
  dat = dat[! rownames(dat)%in%c('Lgroup','Hgroup'),] ; colnames(dat) = dat[1,] ; dat = dat[-1,]
  x = matrix(unlist(strsplit(rownames(dat),'|',fixed=T)),ncol=2,byrow = T,dimnames = list(rownames(dat),list('Gene','ID')))
  dat = data.frame(x,dat,stringsAsFactors = F) ; dat = dat[order(dat$padj),]
  filted = dat[intersect(which(as.numeric(dat$padj)<=0.05),which(abs(as.numeric(dat$log2FoldChange))>1)),]
  write.table(filted$Gene,paste('./result/diff/',paste(name,argv,sep='_'),'.genesymbol.txt',sep=''),sep='\t',col.names=F,row.names=F,quote = F)
  analysis(filted, dat ,name,out.path)
}
print ('work done')
