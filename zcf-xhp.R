library(limma)
library(edgeR)
library(DESeq2)
library(parallel)
library(survival)
library(gridExtra)
library(survminer)
#argv = commandArgs(trailingOnly = T)
fiho = sort(list.files('../download.firehose','uncv2.mRNAseq_raw_counts.txt',full.names=TRUE))
cli = read.csv('./ucsc/clinical.literature.final_use.txt',header = T,sep='\t',stringsAsFactors = F)
rownames(cli) = cli[,1]
n=1
for(f in fiho){
  name = strsplit(strsplit(f,'/',fixed=T)[[1]][3],'.',fixed=T)[[1]][1]
  print(name)
  dat = read.csv(f,header=T,sep='\t',stringsAsFactors=F,row.names=1)
  colnames(dat) = gsub('[.]','-',colnames(dat))
  Tnames = colnames(dat)[substr(colnames(dat),14,14)==0]
  Nnames = colnames(dat)[substr(colnames(dat),14,14)==1]
  
  datuniq = dat[,Tnames[duplicated(substr(Tnames,1,12))==FALSE]]
  clidat = datuniq[,which(substr(colnames(datuniq),1,12) %in% cli$bcr_patient_barcode)]
  share = intersect(substr(Tnames,1,12),substr(Nnames,1,12))
  
  if (length(share)<=10) {print('sample is not enough') ; next}
  if (! dir.exists(paste('../tmp/',name,sep=''))){dir.create(paste('../tmp/',name,sep=''))}
  print(length(share))
  Td = dat[,Tnames[substr(Tnames,1,12)%in%share]] ; Tg = data.frame(sample=colnames(Td),Group='g2',type='Tumor')
  Nd = dat[,Nnames[substr(Nnames,1,12)%in%share]] ; Ng = data.frame(sample=colnames(Nd),Group='g1',type='Normal')
  TvsN = cbind(Nd,Td) ; group = rbind(Ng,Tg)
  TvsN = TvsN[rowSums(TvsN)>ncol(TvsN),]
  TvsN = TvsN[apply(TvsN,1,function(x) sum(x=='0')<length(x)*0.3),]
  
  des = DESeqDataSetFromMatrix(round(TvsN), group, design=~Group)
  des2= DESeq(des)
  res = results(des2)
  #exp =assay(vst(des,blind=FALSE))
  resOrdered = res[order(res$padj),]
  resOrdered=as.data.frame(resOrdered)
  out = data.frame(Gene=rownames(resOrdered),resOrdered,stringsAsFactors=F)
  out$FoldChange = 2^out$log2FoldChange
  sig = out[((out$FoldChange>5) | (out$FoldChange<0.2))&(out$padj<0.05),]
  
  path=paste('../download.cbioprotal/',tolower(name),'_tcga_pan_can_atlas_2018/data_expression_merged.txt',sep='')
  nexp = read.csv(path,header=T,sep='\t',stringsAsFactors=F)
  rownames(nexp) = paste(nexp$Hugo_Symbol,nexp$Entrez_Gene_Id,sep='|') ; nexp = nexp[,3:ncol(nexp)]
  colnames(nexp) = gsub('[.]','-',colnames(nexp))
  datuniq = nexp[,colnames(nexp) %in% Tnames[duplicated(substr(Tnames,1,12))==FALSE]]
  clidats = datuniq[,which(substr(colnames(datuniq),1,12) %in% cli$bcr_patient_barcode)]
  clidat = clidats[apply(clidats,1,function(x) sum(duplicated(x)))<ncol(clidats)*0.25,]
  
  
  x = list()
  for (gene in sig$Gene ){
    if (! gene %in% rownames(clidat)){next}
    clis = cli[substr(colnames(clidat),1,12),]
    clis$Exp = log(as.numeric(clidat[gene,])+1,2)
    clis$means = 'Low' ; clis[clis$Exp>mean(clis$Exp),'means']='High' 
    clis$median = 'Low' ; clis[clis$Exp>median(clis$Exp),'median']='High'
    #
    os = survfit(formula= Surv(clis$OS.time, clis$OS=='1')~means ,data=clis)
    pv1 = survdiff(formula= Surv(clis$OS.time, clis$OS=='1')~means ,data=clis)
    do = ggsurvplot(os,pval=TRUE,legend.title="Group: ",legend.labs=c('High','Low'),pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days')
    d1 = do$plot+labs(title=paste(name,' OS - mean',sep=''))+theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1)
    #
    os = survfit(formula= Surv(clis$OS.time, clis$OS=='1')~median ,data=clis)
    pv2 = survdiff(formula= Surv(clis$OS.time, clis$OS=='1')~median ,data=clis)
    do = ggsurvplot(os,pval=TRUE,legend.title="Group: ",legend.labs=c('High','Low'),pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days')
    d2 = do$plot+labs(title=paste(name,' OS - median',sep=''))+theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1)
    #
    os = survfit(formula= Surv(clis$PFI.time, clis$PFI=='1')~means ,data=clis)
    pv3 = survdiff(formula= Surv(clis$PFI.time, clis$PFI=='1')~means ,data=clis)
    do = ggsurvplot(os,pval=TRUE,legend.title="Group: ",legend.labs=c('High','Low'),pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days')
    d3 = do$plot+labs(title=paste(name,' PFI - mean',sep=''))+theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1)
    #
    os = survfit(formula= Surv(clis$PFI.time, clis$PFI=='1')~median ,data=clis)
    pv4 = survdiff(formula= Surv(clis$PFI.time, clis$PFI=='1')~median ,data=clis)
    do = ggsurvplot(os,pval=TRUE,legend.title="Group: ",legend.labs=c('High','Low'),pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days')
    d4 = do$plot+labs(title=paste(name,' PFI - median',sep=''))+theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1)
    #
    sig[gene,'OS.mean'] = 1-pchisq(pv1$chisq,1); 
    sig[gene,'OS.dedian'] = 1-pchisq(pv2$chisq, 1) ; 
    sig[gene,'PFI.mean'] = 1-pchisq(pv3$chisq, 1); 
    sig[gene,'PFI.median'] = 1-pchisq(pv4$chisq, 1)
    ge = strsplit(gene,'|',fixed=T)[[1]][1] ; path2 = paste('../tmp/',name,'/',name,'.survival.',sep='')
    tiff(paste(path2,ge,'.tiff',sep=''),width =16,height = 14,bg="white",units='cm',compression='lzw' ,res=500)
    grid.arrange(d1,d2,d3,d4,nrow=2,ncol=2)
    dev.off()
    }
  write.table(sig,paste('../tmp/',name,'.deg.txt',sep=''),col.names = TRUE,row.names=F,sep='\t',quote=F)
}
