library(survival)
library(survminer)

## 整理临床数据信息 ##
df = c('Cancer','Characteristic','Type','Total','Percent','Events')
log = NULL
files = paste('./cli_filted/',list.files('./cli_filted/','merge'),sep='')
for (f in files){
  cli = read.csv(f,header = T,sep='\t',stringsAsFactors = F,row.names = 1)
  p = strsplit(strsplit(f,'/',fixed=T)[[1]][3],'.',fixed=T)[[1]][1]
  df = rbind(df, c(p,'','',nrow(cli),'',sum(cli$OS)))
  x = which(colnames(cli)=='age_at_initial_pathologic_diagnosis')
  y = which(colnames(cli)=='OS')-1
  con.tag = colnames(cli)[1:x]
  cat.tag = colnames(cli)[(x+1):y]
  for (i in con.tag){
    n.na = length(which(is.na(cli[,i])=='TRUE'))
    rang = range(cli[,i][which(is.na(cli[,i])!='TRUE')])
    med = median(cli[,i][which(is.na(cli[,i])!='TRUE')])
    mea = mean(cli[,i][which(is.na(cli[,i])!='TRUE')])
    sd = sd(cli[,i][which(is.na(cli[,i])!='TRUE')])
    line1 = c('',i,'Median(Range)',paste(med,'(',paste(rang,collapse='-'),')',sep=''),'',sum(cli[is.na(cli[,i])!='TRUE','OS']))
    line2 = c('','','Mean(SD)',paste(mea,'(±',sd,')',sep=''),'','')
    line3 = c('','','N/A',n.na,paste(round(100*n.na/nrow(cli),2),'%',sep=''),'')
    df = rbind(df,line1,line2,line3)
  }
  for (j in cat.tag){
    df = rbind(df, c('',j,'','','',''))
    t = as.data.frame(table(cli[,j]),stringsAsFactors=F)
    t = t[order(t[,2],decreasing=T),]
    n.na = length(which(is.na(cli[,j])==TRUE))
    if (n.na > 0 ){ t = rbind(t,c('N/A', n.na )) }
    t$perc = paste(round(100*as.numeric(as.character(t[,2]))/nrow(cli),2),'%',sep='')
    for (tag in t[,1]){t[t[,1]==tag,'Events'] = sum(cli[cli[,j]==tag,'OS'][is.na(cli[cli[,j]==tag,'OS'])=='FALSE'])}
    t = data.frame('','',t,stringsAsFactors = F)
    colnames(t) = colnames(df)
    df = rbind(df, as.matrix(t))
  }
}
write.table(df,'x.xls',sep='\t',col.names = F,row.names = F,quote = F)
######################################################

tags=read.csv('no.na.txt',header = T,sep='\t',stringsAsFactors = F)
mrna=read.csv('./ucsc/Protein.RPPA.Tumors.matchinfor.txt',header = T,row.names = 1,sep='\t',stringsAsFactors = F)
mrna=mrna[order(rownames(mrna)),1:6]
mrna = mrna[duplicated(mrna$bcr_patient_barcode)=='FALSE',]
gene.list = colnames(mrna)[2:5]
library(survival)
library(survminer)
rownames(mrna)=mrna[,1]
proj = c(which(tags$Cancer!=''),2274)
sigs = which(tags$Characteristic!='')
types = which(tags$Type!='')
write.table('Cancer\tSignature\tType\tHR\tLow CI95\tUp CI95\tp-value\tNumber_Event\tN/A_number\tLikelihood_ratio_test\tWald_test\tScore_test\tzph\tconcordance','y.txt',sep='\t',row.names=F,col.names=F,quote=F)
for (p in proj[1:33]){
  pname = tags[p,'Cancer']
  use = read.csv(paste('./cli_filted/',pname,'.merge.txt',sep=''),header=T,row.names=1,sep='\t',stringsAsFactors=F)
  psigs = c(sigs[(sigs>p) == (sigs<proj[which(proj%in%p)+1])],proj[which(proj%in%p)+1])
  if (pname != 'LAML'){
  for(g in gene.list ){
    exp = mrna[mrna$type==pname,]
    n2 = paste(g,'v2',sep='_')
    exp[as.numeric(exp[,g])<=median(exp[,g]),n2]='L' ; exp[as.numeric(exp[,g])>median(exp[,g]),n2]='H'
    dat = data.frame(use[,c('OS','OS.time')],exp[rownames(use),c(g,n2)])
    dat[,n2] = factor(dat[,n2],levels = c('L','H'))
    cox = coxph(formula = Surv(dat$OS.time,dat$OS=='1')~dat[,g],data=dat)
    t1 = as.data.frame(summary(cox)$conf.int, stringsAsFactors=F)[,-2]
    t1$p_value = summary(cox)$coefficients[,5]
    t1$Num = paste(cox$n,cox$nevent,sep='/')
    n.na = length(cox$na.action) ; t1$'N/A_num' = paste(n.na,'(',round(100*n.na/nrow(use),2),'%)',sep='')
    Likelihood_ratio_test = summary(cox)$logtest[3] ; t1$Likelihood_ratio_test=Likelihood_ratio_test
    Wald_test = summary(cox)$waldtest[3] ; t1$Wald_test = Wald_test
    Score_test = summary(cox)$sctest[3] ; t1$Score_test = Score_test
    zph = cox.zph(cox) ; t1$zph.test = zph$table[,3]
    concor = summary(cox)$concordance[1] ; t1$concordance = concor
    t1=data.frame(Cancer=tags$Cancer[p],Type=g,Signature='(continuous)',t1)
    write.table(t1,'y.txt',sep='\t',row.names=F,col.names=F,quote=F,append = T)
    #
    cox = coxph(formula = Surv(dat$OS.time,dat$OS=='1')~dat[,n2],data=dat)
    t2 = as.data.frame(summary(cox)$conf.int, stringsAsFactors=F)[,-2]
    t2$p_value = summary(cox)$coefficients[,5]
    t2$Num = paste(paste(nrow(na.omit(dat)[na.omit(dat)[,n2]=='L',]),sum(na.omit(dat)[na.omit(dat)[,n2]=='L',][,'OS']),sep='/'),
                   paste(nrow(na.omit(dat)[na.omit(dat)[,n2]=='H',]),sum(na.omit(dat)[na.omit(dat)[,n2]=='H',][,'OS']),sep='/'),sep='--')
    rownames(t2) = paste(gsub('dat.*]','',rownames(t2)),paste(' (vs ',cox$xlevels[[1]][1],')',sep=''),sep='')
    n.na = length(cox$na.action) ; t2$'N/A_num' = paste(n.na,'(',round(100*n.na/nrow(use),2),'%)',sep='')
    Likelihood_ratio_test = summary(cox)$logtest[3] ; t2$Likelihood_ratio_test=Likelihood_ratio_test
    Wald_test = summary(cox)$waldtest[3] ; t2$Wald_test = Wald_test
    Score_test = summary(cox)$sctest[3] ; t2$Score_test = Score_test
    zph = cox.zph(cox) ; t2$zph.test = zph$table[,3]
    concor = summary(cox)$concordance[1] ; t2$concordance = concor
    t2=data.frame(Cancer=tags$Cancer[p],Type=g,Signature=rownames(t2),t2)
    write.table(t2,'y.txt',sep='\t',row.names=F,col.names=F,quote=F,append = T)
  }}
  for(sig in psigs[-length(psigs)]){
    if((sig+2)%in%psigs){
      sname = tags[sig,'Characteristic']
      cox = coxph(formula = Surv(use$OS.time,use$OS=='1')~use[,sname],data=use)
      t1 = as.data.frame(summary(cox)$conf.int, stringsAsFactors=F)[,-2]
      t1$p_value = summary(cox)$coefficients[,5]
      t1$Num = paste(cox$n,cox$nevent,sep='/')
      rownames(t1) = tags$Characteristic[sig]
      n.na = length(cox$na.action) ; t1$'N/A_num' = paste(n.na,'(',round(100*n.na/nrow(use),2),'%)',sep='')
      Likelihood_ratio_test = summary(cox)$logtest[3] ; t1$Likelihood_ratio_test=Likelihood_ratio_test
      Wald_test = summary(cox)$waldtest[3] ; t1$Wald_test = Wald_test
      Score_test = summary(cox)$sctest[3] ; t1$Score_test = Score_test
      zph = cox.zph(cox) ; t1$zph.test = zph$table[,3]
      concor = summary(cox)$concordance[1] ; t1$concordance = concor
      t1=data.frame(Cancer=tags$Cancer[p],Type=tags$Characteristic[sig],Signature='(continuous)',t1)
      write.table(t1,'y.txt',sep='\t',row.names=F,col.names=F,quote=F,append = T)
      #sur = survdiff(formula = Surv(use$OS.time,use$OS=='1')~use[,1],data=use)
    }else{
      type = types[(types>psigs[which(psigs%in%sig)])==(types<psigs[which(psigs%in%sig)+1])]
      if('TRUE' %in% is.na(tags$X[type])){order = tags$Type[type]}else{order = tags$Type[type][order(tags$X[type])]}
      use[,tags$Characteristic[sig]] = factor(use[,tags$Characteristic[sig]],levels=order)
      cox = coxph(formula = Surv(use$OS.time,use$OS=='1')~use[,tags$Characteristic[sig]],data=use)
      t2 = as.data.frame(summary(cox)$conf.int, stringsAsFactors=F)[,-2]
      t2$p_value = summary(cox)$coefficients[,5]
      nnn = paste(tags$Total[type[which(order%in%tags$Type)]][order(tags$X[type])],tags$Events[type[which(order%in%tags$Type)]][order(tags$X[type])],sep='/')
      t2$Num = paste(nnn[2:length(nnn)],nnn[1],sep='--')
      rownames(t2) = paste(gsub('use.*]','',rownames(t2)),paste(' (vs ',cox$xlevels[[1]][1],')',sep=''),sep='')
      if(length(cox$xlevels[[1]])==2){
        n.na = length(cox$na.action) ; t2$'N/A_num' = paste(n.na,'(',round(100*n.na/nrow(use),2),'%)',sep='')
        Likelihood_ratio_test = summary(cox)$logtest[3] ; t2$Likelihood_ratio_test=Likelihood_ratio_test
        Wald_test = summary(cox)$waldtest[3] ; t2$Wald_test = Wald_test
        Score_test = summary(cox)$sctest[3] ; t2$Score_test = Score_test
        zph = cox.zph(cox) ; t2$zph.test = zph$table[,3]
        concor = summary(cox)$concordance[1] ; t2$concordance = concor
        t2=data.frame(Cancer=tags$Cancer[p],Type=tags$Characteristic[sig],Signature=rownames(t2),t2)
      }else{
        n.na = length(cox$na.action) ; t2['Global','N/A_num'] = paste(n.na,'(',round(100*n.na/nrow(use),2),'%)',sep='')
        Likelihood_ratio_test = summary(cox)$logtest[3] ; t2['Global','Likelihood_ratio_test'] = Likelihood_ratio_test
        Wald_test = summary(cox)$waldtest[3] ; t2['Global','Wald_test'] = Wald_test
        Score_test = summary(cox)$sctest[3] ; t2['Global','Score_test'] = Score_test
        zph = cox.zph(cox) ; t2$zph.test = zph$table[,3]
        concor = summary(cox)$concordance[1] ; t2['Global','concordance'] = concor
        t2=data.frame(Cancer=tags$Cancer[p],Type=tags$Characteristic[sig],Signature=rownames(t2),t2)
      }
      write.table(t2,'y.txt',sep='\t',row.names=F,col.names=F,quote=F,append = T)
      #sur = survdiff(formula = Surv(use$OS.time,use$OS=='1')~use[,1],data=use)
    }
  }
}


#########################################################################


## 分析 ##
projs= sort(unique(cli$type))
df = read.csv('statistical_description.txt',header=T,sep='\t',stringsAsFactors=F)
res=NULL
for (one in projs){
    n.tag = colnames(dat)[2:11]
    s.tag = df[df[,1]==one,2]
    s.tag = s.tag[is.na(s.tag)=='FALSE']
    sigs = c(n.tag, s.tag)
    #res =c('Cancer','Signature','HR','CI95.LOW','CI95.UP','P-value','Likelihood_ratio_test','Wald_test','Score_test','Number','Event','Log-Rank')
    #res=NULL
    for (sig in sigs) {
        use = dat[dat$type==one ,c(sig,'OS','OS.time')]
        if(sig %in% s.tag){use[,1] = as.factor(use[,1])}else{use[,1] = as.numeric(as.character(use[,1]))}
        if (length(table(use[,sig]))<=1){next}
        if (length(which(is.na(use[,1])))>0.2*nrow(use)){next}
        cox = coxph(formula = Surv(use$OS.time,use$OS=='1')~use[,1],data=use)
        t1 = as.data.frame(summary(cox)$conf.int, stringsAsFactors=F)[,-2]
        t1$p_value = summary(cox)$coefficients[,5]
        rownames(t1) = gsub('use[[], 1[]]',paste(cox$xlevels[[1]][1],'-',sep=''),rownames(t1))
        #t1=as.data.frame(t1,stringsAsFactors = F)
        Likelihood_ratio_test = summary(cox)$logtest[3] ; t1['Global','Likelihood_ratio_test'] = Likelihood_ratio_test
        Wald_test = summary(cox)$waldtest[3] ; t1['Global','Wald_test'] = Wald_test
        Score_test = summary(cox)$sctest[3] ; t1['Global','Score_test'] = Score_test
        zph = cox.zph(cox) ; t1$zph.test = zph$table[,3]
        sur = survdiff(formula = Surv(use$OS.time,use$OS=='1')~use[,1],data=use)
        if (sig %in% n.tag){
            num = data.frame(Number=c(summary(cox)$n,''),Event=c(summary(cox)$nevent,''),stringsAsFactors=F)
            num[,'log_rank.test']=NA;
        }else{
            num = data.frame(sur$n,sur$obs) ; colnames(num) = c('','Number','Event')
            num = rbind(num[-1,],num[1,])[,-1]
            num[nrow(num),'log_rank.test'] = 1-pchisq(sur$chisq, length(table(use[,1]))-1)}
        t1 = cbind(t1,num)
        t1 = cbind(Cancer=one,Signature=sig,'type'=c(rownames(t1)),t1) #; rownames(t1)=NULL
        res = rbind(res, t1)
    }
}
write.table(res,'cox.result.txt',col.names=T,row.names=F,quote=F,sep='\t')

b=read.csv('BRCA.merged_only_clinical_clin_format.txt',header = T,row.names = 1,sep='\t',stringsAsFactors = F,na.strings='N/A')
b=t(b)
for (i in colnames(b)){if(length(table(b[,i]))<=2 || length(which(b[,5]=='NA')) > 0.5*nrow(b)){print(i);b=b[,-which(i==colnames(b))]}}
