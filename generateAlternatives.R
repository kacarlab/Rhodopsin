#generate alternative sequence
probabilities=read.csv("CSRho27_2sequences_aligned.asr.probdists.csv")

generate_between_ML_and_Worst=function(probabilities,intNo,seedNo){
  probs=probabilities
  cutoff=0.7
  intNo=intNo
  seedNo=seedNo
  residues=c()
  nodes=c(570,571,572,811,895)
  for (i in nodes){
    prnode=probs[probs[,1]==i,]
    count=0
    lowProbSites=c()
    for(x in 1:nrow(prnode)){
      if(sort(prnode[x,3:23])[21]>=cutoff){
        count=count
      }else{
        lowProbSites=c(lowProbSites,x)
        count=count+1
      }
    }
    intervals=sapply(seq(from=5,to=count,length.out=50),floor)
    set.seed(seedNo)
    tobechanged=sort(sample(lowProbSites,intervals[intNo]))
    for (a in 1:nrow(prnode)){
      if(a %in% tobechanged){
        residues=c(residues, names(sort(prnode[a,3:23])[20]))
      }else{
        residues=c(residues, names(sort(prnode[a,3:23])[21]))
      }
    }
  }
  nodeNames=c(rep(570,nrow(prnode)),rep(571,nrow(prnode)),rep(572,nrow(prnode)),rep(811,nrow(prnode)),rep(895,nrow(prnode)))
  residues_mat=cbind(nodeNames,residues)

  altSeq=c()
  for (b in nodes){
    altSeq=c(altSeq,paste(">",b,sep=""),paste0(residues_mat[which(residues_mat[,1]==b),2],sep="",collapse=""))
  }
  altSeq=as.matrix(altSeq)

  #write new sequences into a file
  outname=paste("alternative_ancestral_sequences",intNo,".fasta", sep="")
  write.table(altSeq, outname, quote=F, row.names=F, col.names = F)
}

for(x in 1:50){
  generate_between_ML_and_Worst(probabilities,x,x)
}
