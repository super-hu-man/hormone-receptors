##mRNA.
args='mRNA.RSEM' # commandArgs( trailingOnly = T)  ## 'RSEM' or 'Protein.RPPA' 
library('ggplot2')
library(ggbeeswarm)

pro.list = c('BLCA','BRCA','COAD','HNSC','KICH','KIRC','KIRP','LIHC','LUAD','LUSC','PRAD','STAD','THCA','UCEC')
cli = read.table('./ucsc/clinical.literature.final_use.txt',header = T,sep='\t',row.names = 1,stringsAsFactors = F)
pro.list = unique(cli$type)
table = read.table(paste('./ucsc/',args,'.combined.txt',sep=''),header = T,row.names=1,sep='\t',stringsAsFactors = F)
table$Cancer = cli[table$bcr_patient_barcode,1]
table = table[is.na(table$Cancer)=='FALSE',]
t.e = table[substr(rownames(table),14,14)==0,]
n.e = table[substr(rownames(table),14,14)!=0,]

colour = c("darkred","#3300FFFF")  ## '#f8766d'和'#00b0f6'
for (gop in colnames(table)[2:5]){
  df = NULL ; p.values = NULL
  df.all = NULL ;p.vall = NULL
  for(pro in pro.list){
    t.p = t.e[t.e[,'Cancer']==pro,]
    n.p = n.e[n.e[,'Cancer']==pro,]
    index = which(pro.list==pro)
    ## tumor vs normal
    if(nrow(n.p)<= 10){
      t.p$type = 'Tumor' ; t.p$point = index+0.1
      p.vall = rbind(p.vall,c(pro,1,1))
      data = rbind(t.p[,c(gop,'Cancer','type','point')],c(0,pro,'Normal',index-0.1))
      df.all = rbind(df.all , data)
    }else{
      t.p$type = 'Tumor' ; t.p$point = index+0.1
      n.p$type = 'Normal'; n.p$point = index-0.1
      wil = wilcox.test(as.numeric(t.p[,gop]) , as.numeric(n.p[,gop]))$p.value
      tup = t.test(as.numeric(t.p[,gop]) , as.numeric(n.p[,gop]))$p.value
      p.vall = rbind(p.vall,c(pro,tup,wil))
      data = rbind(t.p[,c(gop,'Cancer','type','point')],n.p[,c(gop,'Cancer','type','point')])
      df.all = rbind(df.all , data)
    }}
    colnames(df.all)[1] = 'exp' ; colnames(df.all)[3] = 'Type'
    df.all[,4]=as.numeric(df.all[,4]) ; df.all[,1]=as.numeric(df.all[,1])
    pva = as.data.frame(p.vall,stringsAsFactors=F)
    pva[,2] = as.numeric(pva[,2]) ; pva[,3] = as.numeric(pva[,3])
    p = ggplot(data=df.all ,aes(x=Cancer , y=exp,fill = Type))
    p = p + geom_boxplot(alpha = 0,outlier.alpha = F,width=0.42,notch=F,notchwidth=0.1,size=0.3)
    p = p + geom_quasirandom(aes(x=point),colour='grey90',size=1.1,shape=21,width=0.1)
    p = p + geom_boxplot(colour='grey33',alpha = 0.3,outlier.alpha = F,width=0.42,notch=F,notchwidth=0.1,size=0.3)
    p = p + annotate('text',x=pro.list[which(pva[,2]<0.05)], y=rep(max(as.numeric(df.all[,1]))+0.6,length(which(pva[,2]<0.05))), label=paste('t=',format(round(pva[,2],4))[which(pva[,2]<0.05)],sep=''),size=2,colour='red')
    p = p + annotate('text',x=pro.list[which(pva[,3]<0.05)], y=rep(max(as.numeric(df.all[,1]))+1,length(which(pva[,3]<0.05))), label=paste('w=',format(round(pva[,2],4))[which(pva[,3]<0.05)],sep=''),size=2,colour='red')
    p = p + labs(title=paste(args[1],"differential expression of tumor and normal",sep=' '),x="Cancer Type",y ="Expression of mRNA (log2 RSEM+1 )")
    p = p + theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(size=8))
    #tiff(filename = paste(gop,"mRNA.matched.q.tiff",sep='.'),width =30,height = 12,units ="cm",compression="lzw",bg="white",res=1000)
    ggsave(filename = paste(gop,args,"TvsN.p.pdf",sep='.'),p,width =36,height = 12,units ="cm",dpi=1000)
    
    #normal vs matched tumor
  pro.match = pro.list
  for(pro in pro.match){
    t.p = t.e[t.e[,'Cancer']==pro,]
    n.p = n.e[n.e[,'Cancer']==pro,]
    index = which(pro.match==pro)
    if (nrow(n.p)<=10){print(paste(pro,' sample is not enough !',sep='')) ; pro.match = pro.match[-index]}else{
      t.p=t.p[sort(rownames(t.p)),] ; t.p = t.p[duplicated(t.p[,1])=='FALSE',]
      n.p=n.p[sort(rownames(n.p)),] ; n.p = n.p[duplicated(n.p[,1])=='FALSE',]
      share = intersect(t.p[,1],n.p[,1])
      t.u = t.p[t.p[,1]%in%share,] ; t.u$type = 'Tumor' ; t.u$point = index+0.1
      n.u = n.p[n.p[,1]%in%share,] ; n.u$type = 'Normal'; n.u$point = index-0.1
      wpa = wilcox.test(as.numeric(t.u[,3]) , as.numeric(n.u[,3]),paired = T)$p.value
      tpa = t.test(as.numeric(t.u[,3]) , as.numeric(n.u[,3]), paired=T)$p.value
      p.values = rbind(p.values , c(pro, tpa, wpa))
      data = rbind(t.u[,c(gop,'Cancer','type','point')],n.u[,c(gop,'Cancer','type','point')])
      df = rbind(df , data)
  } }
    colnames(df)[1] = 'exp' ; colnames(df)[3] = 'Type'
    df[,4]=as.numeric(df[,4]) ; df[,1]=as.numeric(df[,1])
    pva = as.data.frame(p.values,stringsAsFactors=F)
    pva[,2] = as.numeric(pva[,2]) ; pva[,3] = as.numeric(pva[,3])
    p = ggplot(data=df ,aes(x=Cancer , y=exp,fill = Type))
    p = p + geom_boxplot(alpha = 0,outlier.alpha = F,width=0.38,notch=F,notchwidth=0.1,size=0.3)
    p = p + geom_quasirandom(aes(x=point),colour='grey90',size=1.1,shape=21,width=0.1)
    p = p + geom_boxplot(colour='grey33',alpha = 0.3,outlier.alpha = F,width=0.38,notch=F,notchwidth=0.1,size=0.3)
    p = p + annotate('text',x=pro.match[which(pva[,2]<0.05)], y=rep(max(as.numeric(df[,1]))+0.6,length(which(pva[,2]<0.05))), label=paste('t=',format(round(pva[,2],4))[which(pva[,2]<0.05)],sep=''),size=2,colour='red')
    p = p + annotate('text',x=pro.match[which(pva[,3]<0.05)], y=rep(max(as.numeric(df[,1]))+1,length(which(pva[,3]<0.05))), label=paste('w=',format(round(pva[,2],4))[which(pva[,3]<0.05)],sep=''),size=2,colour='red')
    p = p + labs(title=paste(args[1],"differential expression of matched samples",sep=' '),x="Cancer Type",y ="Expression of mRNA (log2 RSEM+1 )")
    p = p + theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(size=8))
    #tiff(filename = paste(gop,"mRNA.matched.q.tiff",sep='.'),width =30,height = 12,units ="cm",compression="lzw",bg="white",res=1000)
    ggsave(filename = paste(gop,args,"matched.p.pdf",sep='.'),p,width =36,height = 12,units ="cm",dpi=1000)
}