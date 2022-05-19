library(seqinr)
library(Biostrings)

input.nc<-"/home/nacho/Downloads/test/ivar_masked_consensus_Nremoved_v11.csv"
genemap<-"/media/nacho/Data/DockerImages/FHISC2/FHI_SC2_Pipeline_Nanopore/CommonFiles/corona genemap.csv"
reference.path<-"/media/nacho/Data/DockerImages/FHISC2/FHI_SC2_Pipeline_Nanopore/CommonFiles/reference_nc.fasta"

df<-read.csv(input.nc,sep = ";")

genes<-read.csv(genemap)
reference<-read.fasta(reference.path)
reference<-as.character(reference$MN908947)

df$insertionsClean<-NA

for (i in 1:nrow(df)) {
  ins.dumm<-df$insertions[i]
  try(rm(total.neo.tag))
  ins.total <- unlist(base::strsplit(df$insertions[i], ","))
  if(length(ins.total)>0){
  for (k in 1:length(ins.total)) {
    ins.dumm<-ins.total[k]

  if(ins.dumm!=""){
   position<-as.numeric(gsub(":.*", "", ins.dumm))+1 
   gene<-genes[which(position > genes$start & position < genes$end ),]
   if(nrow(gene==1)){
    
    seq<-gsub(".*:", "", ins.dumm) 
    len<-nchar(seq)
    
    aapos<- round((position-gene$start) / 3)
    rest<-(position-gene$start) %% 3
    if(rest == 1) aapos<-aapos 
    gene.oi<-reference[gene$start:gene$end] 
    neo.gene<- c(gene.oi[1:(position-gene$start-1)], unlist(base::strsplit(tolower(seq),"")), gene.oi[(position-gene$start):length(gene.oi)])
    neo.protein <- unlist(base::strsplit(as.character(translate(DNAString( base::paste(toupper(neo.spike), collapse = "")))), "")  )
    #old.protein<-unlist(base::strsplit(as.character(translate(DNAString( base::paste(toupper(gene.oi), collapse = "")))), "")  )
    neo.seq<-neo.protein[aapos: (aapos+round(len/3)-1) ]
    neo.tag<- paste(gene$gene,":ins",aapos,paste(neo.seq,collapse = ""), sep = "")
   }
  }
    if(!exists("total.neo.tag")){
      total.neo.tag<-neo.tag
    }else{
      total.neo.tag<-paste(total.neo.tag, neo.tag, sep = ",")
    }
    
  }
  df$insertionsClean[i]<-total.neo.tag
  df$insertionsClean[which(is.na(df$insertionsClean))]<-""
  }
} 
df$aaSubstitutions<-paste(df$aaSubstitutions, df$insertionsClean,sep = ",")
df$aaSubstitutions<-gsub(",$","",df$aaSubstitutions)

write.table(df, input.nc, sep = ";", row.names = FALSE, quote = TRUE, na = "")