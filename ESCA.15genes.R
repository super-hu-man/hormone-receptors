cli=read.csv('TCGA-CDR.txt',header = T,sep='\t',stringsAsFactors = F)
cli.esca=cli[cli$type=='ESCA',]
rna=read.csv('gene.other.txt',header = T,sep='\t',stringsAsFactors = F,row.names = 1)
rna=as.data.frame(t(rna))
rownames(rna)=gsub('[.]','-',rownames(rna))
rna.t=rna[substr(rownames(rna),14,14)==0,]
rna.esca=rna.t[substr(rownames(rna.t),1,12) %in% cli.esca$bcr_patient_barcode,]
rna.x=rna.esca[substr(rownames(rna.esca),14,15)=='01',]
rownames(rna.x)=substr(rownames(rna.x),1,12)
rna.x$bcr_patient_barcode=rownames(rna.x)
comb=merge(cli.esca,rna.x,by='bcr_patient_barcode')
comb$histological_type = gsub('Esophagus Adenocarcinoma, NOS','NOS',comb$histological_type)
comb$histological_type = gsub('Esophagus Squamous Cell Carcinoma','ESCA',comb$histological_type)

library(survminer)
library(survival)
library(gridExtra)

x = comb
for(g in colnames(rna.esca)){m=median(x[,g]); x[x[,g]<=m,g]='L'; x[x[,g]!='L',g]='H'}
x1 = comb[y$histological_type!='NOS',]
for(g in colnames(rna.esca)){m=median(x1[,g]); x1[x1[,g]<=m,g]='L'; x1[x1[,g]!='L',g]='H'}
x2 = comb[y$histological_type=='NOS',]
for(g in colnames(rna.esca)){m=median(x2[,g]); x2[x2[,g]<=m,g]='L'; x2[x2[,g]!='L',g]='H'}

col = c("#2E9FDF","#E7B800","#E272AC","#82B541","#5B7E2D","#E34625")
for (g in colnames(rna.esca)){
  x$Group = x[,g]
  x1$Group = paste(x1$histological_type,x1[,g],sep=' ')
  x2$Group = paste(x2$histological_type,x2[,g],sep=' ')
  
  os = survfit(Surv(x$OS.time,x$OS==1)~x$Group,data=x)
  x.os = ggsurvplot(os,pval=T,conf.int=F,legend.title="",risk.table = TRUE,palette=col,
                    pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days',
                    ggtheme=theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1))
  x.os = x.os+labs(title='Histological Type OS')
  pf = survfit(Surv(x$PFI.time,x$PFI==1)~x$Group,data=x)
  x.pf = ggsurvplot(pf,pval=T,conf.int=F,legend.title="",risk.table = TRUE,palette=col,
                    pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days',
                    ggtheme=theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1))
  x.pf = x.pf+labs(title='Histological Type PFS')
  
  os = survfit(Surv(x1$OS.time,x1$OS==1)~x1$Group,data=x1)
  x1.os = ggsurvplot(os,pval=T,conf.int=T,legend.title="",risk.table = TRUE,palette=col[1:2],
                     pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days',
                     ggtheme=theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1))
  x1.os = x1.os+labs(title='ESCA OS')
  pf = survfit(Surv(x1$PFI.time,x1$PFI==1)~x1$Group,data=x1)
  x1.pf = ggsurvplot(pf,pval=T,conf.int=T,legend.title="",risk.table = TRUE,palette=col[1:2],
                     pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days',
                     ggtheme=theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1))
  x1.pf = x1.pf+labs(title='ESCA PFS')
  
  os = survfit(Surv(x2$OS.time,x2$OS==1)~x2$Group,data=x2)
  x2.os = ggsurvplot(os,pval=T,conf.int=T,legend.title="",risk.table = TRUE,palette=col[1:2],
                     pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days',
                     ggtheme=theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1))
  x2.os = x2.os+labs(title='NOS OS')
  pf = survfit(Surv(x2$PFI.time,x2$PFI==1)~x2$Group,data=x2)
  x2.pf = ggsurvplot(pf,pval=T,conf.int=T,legend.title="",risk.table = TRUE,palette=col[1:2],
                     pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days',
                     ggtheme=theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1))
  x2.pf = x2.pf+labs(title='NOS PFS')
  
  draw = arrange_ggsurvplots(list(x.os,x1.os,x2.os,x.pf,x1.pf,x2.pf),nrow=3,ncol=2,print=F)
  ggsave(paste(g,'Histological.pdf',sep='.'),draw,width =36,height = 54,units ='cm',dpi=1200)
}

