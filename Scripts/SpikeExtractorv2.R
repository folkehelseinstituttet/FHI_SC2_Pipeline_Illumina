library(seqinr)

multifasta<-read.fasta("/home/docker/Fastq/MultifastaForSpike.fasta")
spike<-lapply(multifasta, function(x) x[21563:25384])
write.fasta(spike, "/home/docker/Fastq/Spike_Aligned.fa", names = paste(names(spike), "_Spike", sep=""))