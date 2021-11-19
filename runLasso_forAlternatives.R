#RUN LASSO (GrplassoSeq)
library(GrplassoSeq)

resultsV3<- train_glasso("Alignment_alternative_ancestral_modern_train.fasta","WAVELENGTHS.txt")
waves3=resultsV3$predicted_wavelen
waves3=waves3[569:874]
waves3=matrix(waves3, nrow=6, ncol=51)
colnames(waves3)=c(1:50,"ML")
rownames(waves3)=c("Anc5-570","Anc4-571","Anc3-572",Anc2-811","Anc1-895")
waves3=t(waves3)
waves3=as.data.frame(waves3)
boxplot(waves3$Anc5,waves3$Anc4,waves3$Anc3,waves3$Anc2,waves3$Anc1, names=c("Anc5","Anc4","Anc3","Anc2","Anc1"), ylab="Peak Absorbtion (nm)")


anc570=cbind(rep("570",51),rep(1,51),waves3[,1])
anc571=cbind(rep("571",51),rep(2,51),waves3[,2])
anc572=cbind(rep("572",51),rep(3,51),waves3[,3])
anc811=cbind(rep("811",51),rep(5,1),waves3[,5])
anc895=cbind(rep("895",51),rep(6,1),waves3[,6])
wt=cbind(rep("H. salinarum NRC-1",1),rep(7,1),568)
bestWaves=as.data.frame(rbind(anc570,anc571,anc572,anc811,anc895))
colnames(bestWaves)=c("Ancestors","Node","Wavelength")
bestWaves$Node=as.factor(bestWaves$Node)
bestWaves

seqType=rep(c(rep("Alt",49),"Worst","ML"),5)
bestWaves=cbind(bestWaves,seqType)

#PLOT LASSO RESULTS
library(ggplot2)
ggplot(bestWaves,aes(x=Node, y=as.numeric(Wavelength),fill=seqType))+
  geom_boxplot(fill="gray")+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  scale_fill_manual(values=c("black","#00ba38","#619cff","black"))+
  labs(x=' ', y='Peak Absorption (nm)')+
  scale_x_discrete(breaks=bestWaves$Node, labels=bestWaves$Ancestors, guide = guide_axis(angle = 45))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x= element_text(size=11),axis.text.y= element_text(size=11),legend.position = "none")
