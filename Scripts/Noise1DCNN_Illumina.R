#Classification
library(keras)

tsvs<-list.files("/home/docker/Fastq/temp/", full.names = TRUE, pattern = "\\.tsv")

dum.arr<-array(data = NA, dim = c(nrow(tsvs),29903,2))
pb<-txtProgressBar(min = 1, max = nrow(tsvs),initial = 1)
for (i in 1:nrow(tsvs)) {
  setTxtProgressBar(pb,i)

    dummy<-read.csv(tsvs[i],sep = "\t", header = FALSE)
    dummy<-dummy[,c(1,2,3)]
    colnames(dummy)<-c("Base","Noise","Reads")
    
    genome.position<-as.data.frame(c(1:29903))
    colnames(genome.position)<-"Base"
    dummy<-merge(genome.position, dummy, by="Base", all=TRUE)
    dummy$Reads[which(is.na(dummy$Reads))]<-0
    
    dummy$Reads<-dummy$Reads/500
    dummy$Reads[which(dummy$Reads>1)]<-1
    dummy$Noise[which(is.na(dummy$Noise))]<-0
    
    dum.arr[i,,1]<-dummy$Noise
    dum.arr[i,,2]<-dummy$Reads
    
}


model<-load_model_hdf5("/home/docker/CommonFiles/NoiseModelIllumina.hdf5")
test.pred<-as.data.frame(predict(model, val.array))