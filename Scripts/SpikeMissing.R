#MIssing aminoacids which are missing in Spike
library(writexl)

genelocation<-read.csv("/home/docker/CommonFiles/genemap.gff",sep = "\t", skip = 4, header = FALSE)
genelocation$V9<-gsub(".*=","",genelocation$V9)

#Spike masking
Spike.wu1<-"MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFR
                     SSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIR
                     GWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVY
                     SSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQ
                     GFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFL
                     LKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITN
                     LCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCF
                     TNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYN
                     YLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPY
                     RVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFG
                     RDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAI
                     HADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPR
                     RARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTM
                     YICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFG
                     GFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFN
                     GLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQN
                     VLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGA
                     ISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMS
                     ECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAH
                     FPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELD
                     SFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELG
                     KYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSE
                     PVLKGVKLHYT"
Spike.wu1<-gsub("\n","",Spike.wu1)
Spike.wu1<-gsub(" ","",Spike.wu1)
Spike.wu1<-unlist(base::strsplit(Spike.wu1,""))
Spike.wu1<-paste("S:",Spike.wu1, c(1:length(Spike.wu1)), "X",sep = "" )
Spike.wu1<-as.data.frame(Spike.wu1)
Spike.wu1$nt1<-NA
Spike.wu1$nt2<-NA
Spike.wu1$nt3<-NA
gene<-"S"
for (i in 1:nrow(Spike.wu1)) {
  
  position<-as.numeric(gsub("\\*","",gsub("[A-Z]","",gsub(".*:","", Spike.wu1$Spike.wu1[i]))))
  
  mut.nt.st<- genelocation$V4[which(genelocation$V9==gene)] +(position*3) -3
  mut.nt.end<- genelocation$V4[which(genelocation$V9==gene)] +(position*3) -1
  
  Spike.wu1$nt1[i]<-(mut.nt.st:mut.nt.end)[1]
  Spike.wu1$nt2[i]<-(mut.nt.st:mut.nt.end)[2]
  Spike.wu1$nt3[i]<-(mut.nt.st:mut.nt.end)[3]
}

nc<-list.files("/home/docker/Fastq", pattern = "Nextclade.results.csv", recursive = TRUE, full.names = TRUE)
df<-read.csv(nc,sep = ";")

missing.S<-df[,c("seqName","missing")]
missing.S$Missing.S<-NA

pb<-txtProgressBar(min = 1, max = nrow(missing.S), initial = 1)
for (kk in 1:nrow(missing.S)) {
  setTxtProgressBar(pb,kk)
  missing.dummy<-unlist(strsplit(missing.S$missing[kk],",") )
  
  present<-c(1:29903)
  if(length(which(present %in% as.numeric(missing.dummy[-grep("-",missing.dummy)])) > 0)){
    present<-present[-which(present %in% as.numeric(missing.dummy[-grep("-",missing.dummy)]))]  
  }
  missing.dummy<-missing.dummy[grep("-", missing.dummy)]
  for (ii in 1:length(missing.dummy)) {
    gap.dummy<-as.numeric(unlist(strsplit(missing.dummy[ii],"-") ))
    present<-present[-which(present>=gap.dummy[1] & present<=gap.dummy[2])]
  }
   missing.dummy<-c(1:29903)[-which(c(1:29903) %in% present )]
   missing.aa<-c(Spike.wu1$Spike.wu1[which(Spike.wu1$nt1  %in% missing.dummy)],
                 Spike.wu1$Spike.wu1[which(Spike.wu1$nt2  %in% missing.dummy)],
                 Spike.wu1$Spike.wu1[which(Spike.wu1$nt3  %in% missing.dummy)])
   
   missing.aa<-unique(missing.aa)
   if(length(missing.aa)>0) missing.S$Missing.S[kk]<- paste(gsub("X","U",missing.aa), collapse = ";")
   
}

write_xlsx(missing.S[,c(1,3)], "/home/docker/Fastq/MissingAA.Spike.xlsx")