input<-list.files("/home/docker/Fastq", pattern ="*_NextcladeAndPangolin.csv", full.names = TRUE, recursive=TRUE)

#input<-list.files("/media/nacho/Data/temp/pangoparser/", pattern ="*_NextcladeAndPangolin.csv", full.names = TRUE, recursive=TRUE)

input2<-list.files("/home/docker/Fastq", pattern ="*_Nextclade.new.results.csv", full.names = TRUE, recursive=TRUE)
#input2<-list.files("/media/nacho/Data/temp/pangoparser/", pattern ="*test_Nextclade.new.results.csv", full.names = TRUE, recursive=TRUE)

  
results<-read.csv(input, sep = "\t")

pango.vars<-read.csv("https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json", sep = ":", header = FALSE)

if(nrow(pango.vars)>0){
  
  pango.vars<-pango.vars[-c(1,nrow(pango.vars)),]
  colnames(pango.vars)<-c("Short","Long")
  pango.vars$Long<-gsub(",","",pango.vars$Long)
  pango.vars$Long<-gsub(" ","",pango.vars$Long)
  pango.vars$Short<-gsub(" ","",pango.vars$Short)
  
  if(length(which(pango.vars$Long==""))>0){
    pango.vars$Long[which(pango.vars$Long=="")]<-pango.vars$Short[which(pango.vars$Long=="")]
  }
  
  short.lineage<-results$lineage
  short.lineage<-gsub("\\..*","",short.lineage)
  coil<-gsub("[[:upper:]]","",results$lineage)
  
  pango.replacement<-pango.vars$Short
  names(pango.replacement)<-pango.vars$Long
  start<-names(pango.replacement)[match(short.lineage, pango.replacement)]
  results$Pangolin_full<-paste(start, coil, sep = "")
  if(length(grep("NA", results$Pangolin_full))>0) results$Pangolin_full[grep("NA", results$Pangolin_full)]<-NA

  }else{
  results$Pangolin_full<- gsub(",.*","",gsub("Alias of ","",results$description))
  results$Pangolin_full[grep(" ",results$Pangolin_full)]<-NA  
}


if(length(input2)>0){
new.nc<-read.csv(input2, sep = ";")
results$clade<-NULL
results<-merge(results, new.nc[,c("seqName","clade")], by="seqName")
}


write.table(results, input , sep = "\t", quote = FALSE, row.names = FALSE)


