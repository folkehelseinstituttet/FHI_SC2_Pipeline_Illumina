#!/usr/bin/env Rscript

file.vec<-list.files(pattern = "report.*\\.tsv")

nc<-list.files(pattern = "Nextclade.new.results.csv")
if(length(nc)>0){ 
  nc.table<-read.csv(nc, sep = ";",stringsAsFactors = FALSE)
  nc.table<-nc.table[,c("seqName","clade_who")]
  colnames(nc.table)<-c("Name", "who.name")
  nc.table$Name<-gsub("_.*","", nc.table$Name)
   
}
#pangolin_ivar_lineage

for (i in 1:length(file.vec) ) {
  input<-file.vec[i]
  results<-read.csv(input, sep = "\t")
  
  if(length(which(colnames(results)=="Pangolin_full"))==0){
  
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
    
    short.lineage<-results$pangolin_ivar_lineage
    short.lineage<-gsub("\\..*","",short.lineage)
    coil<-gsub("[[:upper:]]","",results$pangolin_ivar_lineage)
    
    pango.replacement<-pango.vars$Short
    names(pango.replacement)<-pango.vars$Long
    start<-names(pango.replacement)[match(short.lineage, pango.replacement)]
    results$Pangolin_full<-paste(start, coil, sep = "")
    if(length(grep("NA", results$Pangolin_full))>0) results$Pangolin_full[grep("NA", results$Pangolin_full)]<-NA
    
  }else{
    results$Pangolin_full<- gsub(",.*","",gsub("Alias of ","",results$description))
    results$Pangolin_full[grep(" ",results$Pangolin_full)]<-NA  
  }
  }
  if(length(nc)>0){
    results<-merge(results,nc.table, by="Name", all.x = TRUE)
  }
  
  
  write.table(results, input , sep = "\t", quote = FALSE, row.names = FALSE)
  print(paste("Adding Full Pangolin to",input))
}
