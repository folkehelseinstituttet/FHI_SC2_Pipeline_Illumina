
nc.path<-list.files("/home/docker/Fastq/", pattern = "_Nextclade.results.csv", full.names = TRUE, recursive = TRUE)
na.path<-list.files("/home/docker/Fastq/", pattern = "insertions.csv", full.names = TRUE, recursive = TRUE)

nc<-read.csv(nc.path,sep = ";")
na<-read.csv(na.path)
na$aaInsertions<-gsub("(.*):(.*):(.*)","\\1:ins\\2\\3", na$aaInsertions)
nc<-merge(nc, na[,c(1,3)], by="seqName", all.x=TRUE)
nc$aaSubstitutions<-paste(nc$aaSubstitutions, nc$aaInsertions,sep = ",")
nc$aaSubstitutions<-gsub(",$","",nc$aaSubstitutions)

write.table(nc, nc.path, sep = ";", row.names = FALSE, quote = TRUE, na = "")