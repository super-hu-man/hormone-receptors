##
args=commandArgs( trailingOnly = T)
dat = read.table(args[1],header=T,row.names=1,sep='\t',stringsAsFactors=F)
list = c('AR','ESR1','ESR2','PGR')

test.res = NULL
for (one in rownames(dat) ){
    tmp = as.numeric(dat[one,])
    line = NULL
    for (g in list){
        gene = as.numeric(dat[g,])
        p = cor.test(tmp, gene, method = 'pearson')
        s = cor.test(tmp, gene, method = 'spearman')
        k = cor.test(tmp, gene, method = 'kendall')
        line = c(line , p$estimate, p$p.value ,'' , s$estimate, s$p.value ,'')
    }
    line = c(one, line)
  print(line)
    test.res = rbind(test.res, line)
}
tags = NULL
for(g in list){for(m in c('pearson','spearman')){for(x in c('cor','p-v','q-v')){tags=c(tags,paste(g,m,x,sep='.'))}}}
colnames(test.res) = c('Gene',tags)
for(i in seq(3,25,3)){test.res[,i+1]=p.adjust(test.res[,i])}

write.table('x',colnames=T,rownames=F,sep='\t',quote=F)
