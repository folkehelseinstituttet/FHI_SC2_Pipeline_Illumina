library(readxl)
starting.folder<-"/home/docker/Fastq/"
df<-read.table(paste(starting.folder,"results/nextclade_for_FHI.tsv",sep = ""),sep = "\t")
summaries.file<-list.files(starting.folder, full.names = TRUE, recursive = TRUE, pattern = "report_v"  )

summaries<-read.csv(summaries.file,sep = "\t") 

if(paste(rownames(df),collapse = "")!=paste(c(1:nrow(df)), collapse="")){
  df<-cbind(rownames(df), df)
  rownames(df)<-NULL
  colnames(df)<-colnames(df)[-1]
  df<-df[,-ncol(df)]
}

noise.file<-list.files(paste(starting.folder,"results/2_bam/noiseextractor/",sep = ""), pattern = "ResultsNoisExtractor.*\\.xlsx", full.names = TRUE)

df2<-read_xlsx(noise.file)
df2<-df2[,c("Sample", "Class","Status")]
colnames(df2)<-c("name", "Class","Status")
df2$name<-paste(df2$name,"_ivar_masked",sep = "")

df<-merge(df, df2, by="name")

summaries$name<-paste(summaries$Name,"_ivar_masked",sep = "")
if(length(which(df$Class=="C" & df$name %in% summaries$name[which(summaries$WGS_median<2500)]))>0){
df$Class[which(df$Class=="C" & df$name %in% summaries$name[which(summaries$WGS_median<2500)])]<-"LQ"
}
if(length(which(df$Class=="C"))>0) df$Class[which(df$Class=="C")]<-"MIX"
if(length(which(df$Class=="MIX"))>0) df$Status[which(df$Class=="MIX")]<-"PASS"

write.table(df, paste(starting.folder, "results/nextclade_and_noise_for_FHI.tsv",sep = ""), row.names = FALSE, quote = FALSE,sep = "\t")
  
