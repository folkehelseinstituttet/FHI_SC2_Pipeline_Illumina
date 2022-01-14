library(readxl)
starting.folder<-"/home/docker/Fastq/"
df<-read.table(paste(starting.folder,"results/nextclade_for_FHI.tsv",sep = ""),sep = "\t")

if(paste(rownames(df),collapse = "")!=paste(c(1:nrow(df)), collapse="")){
  df<-cbind(rownames(df), df)
  rownames(df)<-NULL
  colnames(df)<-colnames(df)[-1]
  df<-df[,c(1:15)]
}

noise.file<-list.files(paste(starting.folder,"results/2_bam/noiseextractor/",sep = ""), pattern = "ResultsNoisExtractor.*\\.xlsx", full.names = TRUE)

df2<-read_xlsx(noise.file)
df2<-df2[,c("Sample", "Class","Status")]
colnames(df2)<-c("name", "Class","Status")
df2$name<-paste(df2$name,"_ivar_masked",sep = "")

df<-merge(df, df2, by="name")
write.table(df, paste(starting.folder, "results/nextclade_and_noise_for_FHI.tsv",sep = ""), row.names = FALSE, quote = FALSE,sep = "\t")
  
