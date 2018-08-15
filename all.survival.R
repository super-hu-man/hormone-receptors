library(survminer)
library(survival)
mrna = read.csv('./ucsc/mRNA.RSEM.Tumors.matchinfor.txt',header = T,sep='\t',row.names = 1,stringsAsFactors = F)

#lis = list()
col = c("#2E9FDF","#E7B800","#E272AC","#82B541","#5B7E2D","#E34625")
for(g in c('AR','ESR1','ESR2','PGR')){
  print(g)
  x=mrna[,c(2:6,29:36)]
  proj = sort(unique(x$type))
  for(p in proj){
    m = median(as.numeric(x[x$type==p,g]))
    x[x$type==p,][as.numeric(x[x$type==p,g])<=m,g] = 'L'
    x[x$type==p,][x[x$type==p,g]!='L',g]='H'
  }
  xf=mrna[mrna$gender=='FEMALE',c(2:6,29:36)]
  proj = sort(unique(xf$type))
  for(p in proj){
    m = median(as.numeric(xf[xf$type==p,g]))
    xf[xf$type==p,][as.numeric(xf[xf$type==p,g])<=m,g] = 'L'
    xf[xf$type==p,][xf[xf$type==p,g]!='L',g]='H'
  }
  xm=mrna[mrna$gender=='MALE',c(2:6,29:36)]
  proj = sort(unique(xm$type))
  for(p in proj){
    m = median(as.numeric(xm[xm$type==p,g]))
    xm[xm$type==p,][as.numeric(xm[xm$type==p,g])<=m,g] = 'L'
    xm[xm$type==p,][xm[xm$type==p,g]!='L',g]='H'
  }
  os = survfit(Surv(x$OS.time,x$OS==1)~x[,g],data=x)
  x.os = ggsurvplot(os,pval=T,conf.int=F,legend.title="",risk.table = TRUE,palette=col,
                    pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days',
                    ggtheme=theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1))
  x.os = x.os+labs(title='Gene All Sample OS')
  pf = survfit(Surv(x$PFI.time,x$PFI==1)~x[,g],data=x)
  x.pf = ggsurvplot(pf,pval=T,conf.int=F,legend.title="",risk.table = TRUE,palette=col,
                    pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days',
                    ggtheme=theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1))
  x.pf = x.pf+labs(title='Gene All Sample OS')
  
  os = survfit(Surv(xf$OS.time,xf$OS==1)~xf[,g],data=xf)
  xf.os = ggsurvplot(os,pval=T,conf.int=F,legend.title="",risk.table = TRUE,palette=col,
                    pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days',
                    ggtheme=theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1))
  xf.os = xf.os+labs(title='Gene All Sample OS')
  pf = survfit(Surv(xf$PFI.time,xf$PFI==1)~xf[,g],data=xf)
  xf.pf = ggsurvplot(pf,pval=T,conf.int=F,legend.title="",risk.table = TRUE,palette=col,
                    pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days',
                    ggtheme=theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1))
  xf.pf = xf.pf+labs(title='Gene All Sample OS')
  
  os = survfit(Surv(xm$OS.time,xm$OS==1)~xm[,g],data=xm)
  xm.os = ggsurvplot(os,pval=T,conf.int=F,legend.title="",risk.table = TRUE,palette=col,
                    pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days',
                    ggtheme=theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1))
  xm.os = xm.os+labs(title='Gene All Sample OS')
  pf = survfit(Surv(xm$PFI.time,xm$PFI==1)~xm[,g],data=xm)
  xm.pf = ggsurvplot(pf,pval=T,conf.int=F,legend.title="",risk.table = TRUE,palette=col,
                    pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days',
                    ggtheme=theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1))
  xm.pf = xm.pf+labs(title='Gene All Sample OS')
  draw = arrange_ggsurvplots(list(x.os,xf.os,xm.os,x.pf,xf.pf,xm.pf),nrow=3,ncol=2,print=F)
  ggsave(paste(g,'Histological.pdf',sep='.'),draw,width =36,height = 54,units ='cm',dpi=1200)
  #lis[which(g==c('AR','ESR1','ESR2','PGR'))]=x.os
}