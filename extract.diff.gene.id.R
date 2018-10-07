library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
argvs = commandArgs(trailingOnly = T)

files = list.files('./result/diff',paste('deseq',argv,sep='.'),full.names = T)
out.path = paste('./result/clusterprofiler',argv,sep='/')
if(! dir.exists(out.path)){dir.create(out.path,recursive=T)}
print(out.path)
# defined enrich function
analysis = function(genes,type,name,out.path){
  print(c(name,'BP start'))
  bp = enrichGO(gene = genes , OrgDb = org.Hs.eg.db , keyType = 'ENTREZID' , readable = TRUE ,
                ont = 'BP' , pAdjustMethod = 'BH' , pvalueCutoff = 0.05 , qvalueCutoff  = 0.05 )
  bps=simplify(bp)
  print(c(name,'CC start'))
  cc = enrichGO(gene = genes , OrgDb = org.Hs.eg.db , keyType = 'ENTREZID' , readable = TRUE ,
                ont = 'CC' , pAdjustMethod = 'BH' , pvalueCutoff = 0.05 , qvalueCutoff  = 0.05 )
  ccs=simplify(cc) 
  print(c(name,'MF start'))
  mf = enrichGO(gene = genes , OrgDb = org.Hs.eg.db , keyType = 'ENTREZID' , readable = TRUE ,
                ont = 'MF' , pAdjustMethod = 'BH' , pvalueCutoff = 0.05 , qvalueCutoff  = 0.05 )
  mfs=simplify(mf) ; 
  print(c(name,'KEGG start'))
  kegg = enrichKEGG(gene = genes , organism = 'hsa' ,
                    pvalueCutoff = 0.05 , qvalueCutoff = 0.05 , use_internal_data=F )
  #
  write.table(as.data.frame(kegg),paste(out.path,'/',paste(name,type,sep='_'),'.KEGG.txt',sep=''),sep='\t',quote=F,col.names=T,row.names=F)
  write.table(as.data.frame(bps),paste(out.path,'/',paste(name,type,sep='_'),'.GO.BP.simplify.txt',sep=''),sep='\t',quote=F,col.names=T,row.names=F)
  write.table(as.data.frame(ccs),paste(out.path,'/',paste(name,type,sep='_'),'.GO.CC.simplify.txt',sep=''),sep='\t',quote=F,col.names=T,row.names=F)
  write.table(as.data.frame(mfs),paste(out.path,'/',paste(name,type,sep='_'),'.GO.MF.simplify.txt',sep=''),sep='\t',quote=F,col.names=T,row.names=F)
}

# genes summary
plot(c(-1,1)~c(0,0),type='l',xlim=c(-10,10),ylim=c(0,1), lty=2)
tmp = data.frame(setNames(replicate(1,numeric(0), simplify = F), c("Gene")))
for(f in files){
  name = strsplit(strsplit(f,'/',fixed=T)[[1]][4],'.',fixed=T)[[1]][1]
  print(name)
  dat = read.csv(f,header = T, sep='\t',stringsAsFactors = F) 
  rownames(dat) = dat$Gene 
  x = matrix(unlist(strsplit(rownames(dat),'|',fixed=T)),ncol=2,byrow = T,dimnames = list(rownames(dat),list('Gene_name','ID')))
  dat = data.frame(x,dat,stringsAsFactors = F) ; dat = dat[order(dat$padj),]
  filted = dat[intersect(which(as.numeric(dat$padj)<=0.05),which(abs(as.numeric(dat$log2FoldChange))>1)),]
  x = filted[c(3,5)] ; colnames(x)[2] = name
  tmp = merge(tmp, x , all=T) 
}
rownames(tmp) = tmp$Gene ; tmp=tmp[,-1]
share = apply(tmp,1,function(x){s=sum(!is.na(x));return(s)})
num = apply(tmp,2,function(x){s=sum(!is.na(x));return(s)})
out= cbind(Number=share,tmp) ; out = rbind(Numbers=c(NA,as.character(num)),out)
out = data.frame(Names=row.names(out),out)
write.table(out,'./HvsL.DEGenes/cancers.DEgenes.summary.txt',col.names=T,row.names=F,sep='\t',quote=F)

#gene enrich analysis
dat = read.csv('./HvsL.DEGenes/cancers.DEgenes.summary.txt',header=T,row.names=1,sep='\t',stringsAsFactors=F)
dat = dat[-1,-1] ; dat[is.na(dat)]=0
x = matrix(unlist(strsplit(rownames(dat),'|',fixed=T)),ncol=2,byrow = T,dimnames = list(rownames(dat),list('Gene','ID')))
for(n in 1:ncol(dat)){
  name = colnames(dat)[n]
  print(name)
  up = x[dat[,n]>0,2] ; 
  down = x[dat[,n]<0,2]
  analysis( up ,'up',name,out.path)
  analysis(down,'down',name,out.path)
}
print ('work done')
