#LW preparation Illumina
args=commandArgs(TRUE)
input<-list.files("/home/docker/Fastq/", pattern ="_summaries_and_Pangolin.*", full.names = TRUE, recursive=TRUE)


results<-read.csv(input, sep = "\t")

results$Pangolin_full<- gsub(",.*","",gsub("Alias of ","",results$description))
results$Pangolin_full[grep(" ",results$Pangolin_full)]<-NA

colnames(results)[which(colnames(results)=="pangolin_version")]<-"pangolin_SW_version"
colnames(results)[which(colnames(results)=="version")]<-"pangolin_version"
results$pangolin_nextclade<-gsub(".*/","",results$clade)
results$clade_nextclade<-gsub("/.*","",results$clade)
results$clade<-NULL
results$who.name <- gsub("\\).*","", gsub(".*\\(", "", results$clade_nextclade))
results$who.name[-grep("\\(",results$clade_nextclade)]  <- NA
results$nextclade_version<-args[1]
results$nextclade_SW_version<-args[2]
results$nextclade_version<-gsub("\"", "", results$nextclade_version )

write.table(results, gsub("\\.csv","_LW.csv",input ) , sep = "\t", quote = FALSE, row.names = FALSE)
