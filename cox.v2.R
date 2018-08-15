library(survival)
library(survminer)

tags=read.csv('no.na.txt',header = T,sep='\t',stringsAsFactors = F)
mrna=read.csv('mRNA.RSEM.Tumors.txt',header = T,row.names = 1,sep='\t',stringsAsFactors = F)
mrna=mrna[order(rownames(mrna)),]
mrna = mrna[duplicated(mrna$bcr_patient_barcode)=='FALSE',]
library(survival)
library(survminer)

numbers = function(n) {
  x = length(which(use[,tags$Characteristic[sig]]==n))
  y = sum(use[use[,tags$Characteristic[sig]]==n,'OS'])
  out =  paste(x,y,sep='/')
  return(out)
}

rownames(mrna)=mrna[,1]
proj = c(which(tags$Cancer!=''),2274)
sigs = which(tags$Characteristic!='')
types = which(tags$Type!='')
write.table('Cancer\tSignature\tType\tHR\tLow CI95\tUp CI95\tp-value\tNumber_Event\tN/A_number\tLikelihood_ratio_test\tWald_test\tScore_test\tzph\tconcordance','y.txt',sep='\t',row.names=F,col.names=F,quote=F)
for (p in proj[1:33]){
  pname = tags[p,'Cancer']
  use = read.csv(paste('./cli_filted/',pname,'.merge.txt',sep=''),header=T,row.names=1,sep='\t',stringsAsFactors=F)
  psigs = c(sigs[(sigs>p) == (sigs<proj[which(proj%in%p)+1])],proj[which(proj%in%p)+1])
  for(g in c('AR','ESR1','ESR2','PGR')){
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
  }
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
      if(('TRUE' %in% is.na(use[,tags$Characteristic[sig]]))&&
        ((length(which(is.na(use[,tags$Characteristic[sig]])))>5)||
        (length(which(is.na(use[,tags$Characteristic[sig]])))/nrow(use)>0.05))){
        use[,tags$Characteristic[sig]][is.na(use[,tags$Characteristic[sig]])]='Unknow'
        order = c(order,'Unknow')
      }
      use[,tags$Characteristic[sig]] = factor(use[,tags$Characteristic[sig]],levels=order)
      cox = coxph(formula = Surv(use$OS.time,use$OS=='1')~use[,tags$Characteristic[sig]],data=use)
      t2 = as.data.frame(summary(cox)$conf.int, stringsAsFactors=F)[,-2]
      t2$p_value = summary(cox)$coefficients[,5]
      nnn = paste(tags$Total[type[which(order%in%tags$Type)]][order(tags$X[type])],tags$Events[type[which(order%in%tags$Type)]][order(tags$X[type])],sep='/')
      nnn = sapply(order,numbers)
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

quit()

#########################################################################

write.table(df,'statistical_description.txt',col.names=F,row.names=F,quote=F,sep='\t')

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
