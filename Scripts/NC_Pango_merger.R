nc.file<-list.files("/home/docker/Fastq/", pattern ="Nextclade.results2.csv", full.names = TRUE, recursive=TRUE)
nc.new.file<-list.files("/home/docker/Fastq/", pattern ="Nextclade.new.results.csv", full.names = TRUE, recursive=TRUE)
pango.file<-list.files("/home/docker/Fastq/", pattern ="_pangolin_out.csv", full.names = TRUE, recursive=TRUE)

nc<-read.csv(nc.file[1], sep = "\t")
nc.new<-read.csv(nc.new.file[1], sep = ";")
pango<-read.csv(pango.file[1])

nc$dummy<-rownames(nc)
rownames(nc)<-NULL
nc<-nc[,c(length(nc), 1:(ncol(nc)-1))]
names<-colnames(nc)
nc<-nc[,-ncol(nc)]
colnames(nc)<-names[-1]
colnames(pango)[1]<-"name"
pango$taxon<-pango$name
ncpango<-merge(nc, pango[,c("name","lineage","version","pangolin_version","taxon")], by="name", all=TRUE)

colnames(nc.new)[1]<-"name"
try(nc.new$clade<-paste(nc.new$clade, nc.new$Nextclade_pango,sep = "/"))
ncpango<-merge(ncpango, nc.new[,c("name", "clade")], by="name", all=TRUE)
#Including technology
ncpango$Technology<-"Illumina"
write.table(ncpango, "/home/docker/Fastq/NextcladeAndPangolin.out2.csv", sep = "\t", quote = FALSE, row.names = FALSE)