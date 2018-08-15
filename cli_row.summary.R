if(file.exists('./cli_summary')=='FALSE'){dir.create('./cli_summary')}
files = paste('cli_raw/',list.files('cli_raw/','na'),sep='')
for (f in files){
  pro = toupper(stringr::str_extract_all(f,'[[:alpha:]]+')[[1]][7])
  write.table('Factor1\tFactor2\tTags...',paste('./cli_summary/',pro,'.tag.summary.txt',sep=''),sep='\t',quote = F,row.names = F,col.names = F)
  cli = read.csv(f,header = T,sep='\t',stringsAsFactors = F,na.strings = T)
  for(i in colnames(cli)){
    x = as.data.frame(table(cli[-c(1,2),][,i]))
    y = t(x[order(x[,2],decreasing=T),])
    y = data.frame(Factor1=i,Factor2=cli[1,which(colnames(cli)%in%i)],y)
    write.table(y,paste('./cli_summary/',pro,'.tag.summary.txt',sep=''),sep='\t',quote = F,row.names = F,col.names = F,append = T)
  }
}
