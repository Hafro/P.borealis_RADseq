setwd("C:/Users/aki/Desktop/P.borealis/P.borealis/")
require(ggplot2)
require(vcfR)
require(cowplot)
require(ggrepel)
require(seqinr)

#require(BiocManager)
#BiocManager::install("qvalue")
#require(devtools)
#devtools::install_github("whitlock/OutFLANK")
require(OutFLANK)

#BiocManager::install("SNPRelate")
require(SNPRelate)

#remotes::install_github("privefl/bigsnpr")
require(bigsnpr)

#install.packages("robust")
require(robust)

#install.packages("bigstatsr")
require(bigstatsr)

###############################################################
# Pre-filtration (population and SNP level) data visualization 
###############################################################
# x<-read.vcfR("P.borealis_stacks.vcf")
# 
# queryMETA(x)
# queryMETA(x, element = 'DP')
# 
# ad <- extract.gt(x, element = "AD", as.numeric=TRUE)
# 
# dp <- extract.gt(x, element = "DP", as.numeric=TRUE)
# dp2<-data.frame(colnames(dp),colMeans(dp,na.rm = T))
# colnames(dp2)<-c("Sample","MeanDP")
# dp2$Max <- apply(dp, 2, function(x) max(x, na.rm = TRUE))
# dp2$Min <- apply(dp, 2, function(x) min(x, na.rm = TRUE))
# dp2$Median <- apply(dp, 2, function(x) median(x, na.rm = TRUE))
# 
# ggplot(dp2)+
#   geom_col(aes(x=Sample,y=MeanDP))+
#   geom_point(aes(x=Sample,y=Median))+
#   #geom_pointrange(aes(x=Sample,y=MeanDP, ymin=Min,ymax=Max))+
#   ylab("DP")+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=6.5))
# 
# dp3<-data.frame(row.names(dp),rowMeans(dp,na.rm = T))
# colnames(dp3)<-c("SNP","MeanDP")
# dp3$MedianDP <- apply(dp, 1, function(x) median(x, na.rm = TRUE))
# dp3$Max <- apply(dp, 1, function(x) max(x, na.rm = TRUE))
# dp3$Min <- apply(dp, 1, function(x) min(x, na.rm = TRUE))
# 
# ggplot(dp3)+
#   geom_col(aes(x=SNP,y=MeanDP))+
#   geom_point(aes(x=SNP,y=MedianDP))+
#   #geom_pointrange(aes(x=Sample,y=MeanDP, ymin=Min,ymax=Max))+
#   ylab("DP")+
#   ylim(0,50)+
#   theme_classic()+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())
# 
# mean((dp3$MeanDP))
# sd((dp3$MeanDP))
# 
# ggplot(dp3[dp3$MeanDP>14 & dp3$MeanDP<22,])+
#   geom_col(aes(x=SNP,y=MeanDP))+
#   geom_point(aes(x=SNP,y=MedianDP))+
#   #geom_pointrange(aes(x=Sample,y=MeanDP, ymin=Min,ymax=Max))+
#   ylab("DP")+
#   ylim(0,50)+
#   theme_classic()+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())
# 
# hist((dp3[dp3$MeanDP>14 & dp3$MeanDP<22,]$MeanDP))
# 
# GQ <- extract.gt(x, element = "GQ", as.numeric=TRUE)
# 
# 
indmiss<-read.table("out.imiss",sep="\t",header = T)
indmiss$INDV

ind<-c("S2_10", "S2_11", "S2_13", "S2_17", "S2_18", "S2_19", "S2_20", "S2_21", "S2_22", "S2_23", "S2_24", "S2_3", "S2_6", "S2_7", "S2_9", "AR_26", "AR_29", "AR_30", "AR_31", "AR_32", "AR_33", "AR_34", "AR_36", "AR_39", "AR_41", "AR_42", "AR_43", "AR_44", "AR_45", "AR_47", "AR_48", "S1_49", "S1_51", "S1_52", "S1_53", "S1_56", "S1_57", "S1_58", "S1_59", "S1_60", "S1_66", "S1_67", "S1_68", "S1_70", "S1_71", "S1_72", "S4_74", "S4_75", "S4_76", "S4_78", "S4_79", "S4_80", "S4_81", "S4_84", "S4_86", "S4_87", "S4_88", "S4_89", "S4_90", "S4_91", "S4_92", "S4_94", "S5_1", "S5_10", "S5_11", "S5_12", "S5_2", "S5_3", "S5_4", "S5_5", "S5_6", "S5_7", "S5_8", "S5_9", "S3_66", "S3_67", "S3_68", "S3_69", "S3_70", "S3_71", "S3_72", "S3_73", "S3_74", "S3_75", "S3_76", "OS_141", "OS_142", "OS_143", "OS_144", "OS_145", "OS_146", "OS_147", "OS_148", "OS_149", "OS_150", "OS_151")
samp<-c("10_S2","11_S2","13_S2","17_S2","18_S2","19_S2","20_S2","21_S2","22_S2","23_S2","24_S2","3_S2","6_S2","7_S2","9_S2","26_AR","29_AR","30_AR","31_AR","32_AR","33_AR","34_AR","36_AR","39_AR","41_AR","42_AR","43_AR","44_AR","45_AR","47_AR","48_AR","49_S1","51_S1","52_S1","53_S1","56_S1","57_S1","58_S1","59_S1","60_S1","66_S1","67_S1","68_S1","70_S1","71_S1","72_S1","74_S4","75_S4","76_S4","78_S4","79_S4","80_S4","81_S4","84_S4","86_S4","87_S4","88_S4","89_S4","90_S4","91_S4","92_S4","94_S4","1_S5","10_S5","11_S5","12_S5","2_S5","3_S5","4_S5","5_S5","6_S5","7_S5","8_S5","9_S5","66_S3","67_S3","68_S3","69_S3","70_S3","71_S3","72_S3","73_S3","74_S3","75_S3","76_S3","141_OS","142_OS","143_OS","144_OS","145_OS","146_OS","147_OS","148_OS","149_OS","150_OS","151_OS")
pop<-c(rep("S2",15),rep("AR",16),rep("S1",15),rep("S4",16),rep("S5",12),rep("S3",11),rep("OS",11))

year<-c(rep("2018",62),rep("2021",34))
indmiss$Sample<-ind
indmiss$Sample2<-samp
indmiss$Pop<-pop
indmiss$Year<-year
indmiss<-indmiss[order(indmiss$F_MISS,decreasing = T),]
indmiss$Sample2<-factor(indmiss$Sample2,levels = indmiss$Sample2)
ggplot(data=indmiss)+
  geom_col(aes(x = Sample2, y = (F_MISS),fill=Year))+
  geom_abline(slope=0,intercept = 0.5,col="black",lwd=1)+
  xlab("")+
  ylab("Fraction of sites with missing data")+
  theme_classic()+
  theme(axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=14,face="bold"),
        axis.text.y = element_text(size=14,face="bold"),
        axis.text.x = element_text(size=14,angle = 270, vjust = 0.2,hjust = 1,face="bold"))

aggregate(F_MISS~Pop,indmiss,FUN=mean)

#New filtered file
x<-read.vcfR("../P.borealis_stacks.g5mac3dp10ind50g95mdp20p1RandOnePerContig.vcf")

queryMETA(x)
dp <- extract.gt(x, element = "DP", as.numeric=TRUE)
GQ <- extract.gt(x, element = "GQ", as.numeric=TRUE)

summary(dp)
min(dp,na.rm = T)
max(dp,na.rm = T)

mean(colMeans(dp,na.rm = T))
median(colMeans(dp,na.rm = T))
sd(colMeans(dp,na.rm = T))
hist(dp,xlab="Sequencing depth")
min(apply(dp, 2, function(x) min(x, na.rm = TRUE)))
max(apply(dp, 2, function(x) max(x, na.rm = TRUE)))

summary(GQ)
hist(GQ)
min(apply(GQ, 2, function(x) min(x, na.rm = TRUE)))
max(apply(GQ, 2, function(x) max(x, na.rm = TRUE)))

geno <- extract.gt(x) # Character matrix Containing the genotypes
position <- getPOS(x) # Positions in bp
chromosome <- getCHROM(x) # Chromosome information
pos_loc<-paste(chromosome,position,sep="_")

link<-table(chromosome)
mean(link)
sum(link>1)/length(link)

pop2<-c(rep("S2",9),rep("AR",10),rep("S1",10),rep("S4",7),rep("S5",11),rep("S3",10),rep("OS",11))
samp2<-c("11_S2","18_S2","19_S2","20_S2","21_S2","22_S2","24_S2","3_S2","9_S2","29_AR","30_AR","31_AR","36_AR","39_AR","42_AR","43_AR","44_AR","45_AR","48_AR","49_S1","51_S1","52_S1","57_S1","59_S1","66_S1","67_S1","68_S1","70_S1","72_S1","74_S4","75_S4","80_S4","81_S4","88_S4","89_S4","92_S4","1_S5","10_S5","11_S5","12_S5","2_S5","3_S5","4_S5","5_S5","6_S5","8_S5","9_S5","67_S3","68_S3","69_S3","70_S3","71_S3","72_S3","73_S3","74_S3","75_S3","76_S3","141_OS","142_OS","143_OS","144_OS","145_OS","146_OS","147_OS","148_OS","149_OS","150_OS","151_OS")
year2<-c(rep("2018",36),rep("2021",32))

#prepare sample for OutFLANK analysis
G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno),dimnames = list(pos_loc,colnames(geno)))

G[geno %in% c("0/0")] <- 0
G[geno %in% c("0/1", "1/0")] <- 1
G[geno %in% c("1/1")] <- 2

sum(is.na(geno))/(nrow(geno)*ncol(geno))

G[is.na(G)]<-9

#separate out the temporal samples
geno_18<-geno[,c(1:36)]
colnames(geno_18)

geno_21<-geno[,c(37:68)]
colnames(geno_21)

sum(is.na(geno_18))/(nrow(geno_18)*ncol(geno_18))

sum(is.na(geno_21))/(nrow(geno_21)*ncol(geno_21))

nll_18<-NULL
for(i in 1:nrow(geno_18)){
  if(length(unique(geno_18[i,]))==1){
    nll_18<-c(nll_18,i)
  }
}

nll_21<-NULL
for(i in 1:nrow(geno_21)){
  if(length(unique(geno_21[i,]))==1){
    nll_21<-c(nll_21,i)
  }
}

G18 <- matrix(NA, nrow = nrow(geno_18), ncol = ncol(geno_18),dimnames = list(pos_loc,colnames(geno_18)) )

G18[geno_18 %in% c("0/0")] <- 0
G18[geno_18 %in% c("0/1", "1/0")] <- 1
G18[geno_18 %in% c("1/1")] <- 2

G18[is.na(G18)]<-9

G21 <- matrix(NA, nrow = nrow(geno_21), ncol = ncol(geno_21),dimnames = list(pos_loc,colnames(geno_21)) )

G21[geno_21 %in% c("0/0")] <- 0
G21[geno_21 %in% c("0/1", "1/0")] <- 1
G21[geno_21 %in% c("1/1")] <- 2

G21[is.na(G21)]<-9

SNPmat18<-(t(G18))

SNPmat21<-(t(G21))

#separate out just the Skjálfandi samples
G_Skj<-G[,c(1:9,20:68)]

G18_Skj<-G[,c(1:9,20:36)]

G21_Skj<-G[,c(37:68)]

SNPmat_Skj<-(t(G_Skj))

SNPmat18_Skj<-(t(G18_Skj))

SNPmat21_Skj<-(t(G21_Skj))

#smartPCA
spca_tab<-read.table("../P.borealis_stacks.g5mac3dp10ind50g95mdp20p1RandOnePerContig.evec")
spca_Per_exp<-read.table("../P.borealis_stacks.g5mac3dp10ind50g95mdp20p1RandOnePerContig.eval")
spca_tab$V7<-pop2
colnames(spca_tab)<-c("sample.id","EV1","EV2","EV3", "EV4","EV5","Pop")
spca_tab$order<-c(rep(3,9),rep(1,10),rep(2,10),rep(5,7),rep(6,11),rep(4,10),rep(7,11))
spca_tab<-spca_tab[order(spca_tab$order),]
spca_tab$Pop<-factor(spca_tab$Pop,levels=unique(spca_tab$Pop))
  
plot(spca_Per_exp$V1[1:10],type="b",ylab=("Eigenvalue"),xlab="PC")
spca_Per_exp$V1<-(spca_Per_exp$V1/sum(spca_Per_exp$V1))*100

ggplot(spca_tab,aes(EV1,EV2,color=Pop,shape=Pop)) +
  xlab(paste("PC1 (",round(spca_Per_exp[1,],2),"%)",sep="")) + 
  ylab(paste("PC2 (",round(spca_Per_exp[2,],2),"%)",sep="")) + 
  stat_ellipse(data = spca_tab[spca_tab$Pop==c("OS"),], level=0.75,size=2)+
  geom_point(size=5) +
  #geom_point(size=3, subset = .(label == 'point'))+
  scale_shape_manual(name="Pop", labels=unique(spca_tab$Pop)[order(unique(spca_tab$Pop))], values=shp)+
  scale_color_manual(name="Pop", labels=unique(spca_tab$Pop)[order(unique(spca_tab$Pop))], values=cl)+
  labs(color="") + 
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))

#smartPCA redone w/ no MAF
spca_tab<-read.table("P.borealis_stacks.g5mac3dp10ind50g95mdp20p1.evec")
spca_Per_exp<-read.table("P.borealis_stacks.g5mac3dp10ind50g95mdp20p1.eval")
plot(spca_Per_exp$V1[1:10],type="b",ylab=("Eigenvalue"),xlab="PC")
spca_Per_exp$V1<-(spca_Per_exp$V1/sum(spca_Per_exp$V1))*100
spca_tab$V7<-pop2
colnames(spca_tab)<-c("sample.id","EV1","EV2","EV3", "EV4","EV5","Pop")
spca_tab$Pop<-factor(spca_tab$Pop, levels = c("AR","S1","S2","S3","S4","S5","OS"))
#spca_tab$order<- c(rep(1,25),rep(2,12),rep(3,5),rep(4,11),rep(9,3),rep(5,12),rep(6,8),rep(7,2),rep(8,13),10)
#spca_tab<-spca_tab[order(spca_tab$order),]

barplot(spca_Per_exp$V1)

plot(spca_Per_exp$V1[1:10],type="b",ylab=("Percent variance explained"),xlab="PC")

ggplot(spca_tab,aes(EV1,EV2,color=Pop,shape=Pop)) +
  xlab(paste("PC1 (",round(spca_Per_exp[1,],2),"%)",sep="")) + 
  ylab(paste("PC2 (",round(spca_Per_exp[2,],2),"%)",sep="")) + 
  stat_ellipse(level=0.75,size=2)+
  geom_point(size=5) +
  #geom_point(size=3, subset = .(label == 'point'))+
  scale_shape_manual(name="Pop", labels=unique(spca_tab$Pop)[order(unique(spca_tab$Pop))], values=shp)+
  scale_color_manual(name="Pop", labels=unique(spca_tab$Pop)[order(unique(spca_tab$Pop))], values=cl)+
  labs(color="") + 
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))

ggplot(spca_tab,aes(EV1,EV3,color=Pop,shape=Pop)) +
  xlab(paste("PC1 (",round(spca_Per_exp[1,],2),"%)",sep="")) + 
  ylab(paste("PC3 (",round(spca_Per_exp[3,],2),"%)",sep="")) + 
  stat_ellipse(level=0.75,size=2)+
  geom_point(size=5) +
  #geom_point(size=3, subset = .(label == 'point'))+
  scale_shape_manual(name="Pop", labels=unique(spca_tab$Pop)[order(unique(spca_tab$Pop))], values=shp)+
  scale_color_manual(name="Pop", labels=unique(spca_tab$Pop)[order(unique(spca_tab$Pop))], values=cl)+
  labs(color="") + 
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))

ggplot(spca_tab,aes(EV1,EV4,color=Pop,shape=Pop)) +
  xlab(paste("PC1 (",round(spca_Per_exp[1,],2),"%)",sep="")) + 
  ylab(paste("PC4 (",round(spca_Per_exp[4,],2),"%)",sep="")) + 
  stat_ellipse(level=0.75,size=2)+
  geom_point(size=5) +
  #geom_point(size=3, subset = .(label == 'point'))+
  scale_shape_manual(name="Pop", labels=unique(spca_tab$Pop)[order(unique(spca_tab$Pop))], values=shp)+
  scale_color_manual(name="Pop", labels=unique(spca_tab$Pop)[order(unique(spca_tab$Pop))], values=cl)+
  labs(color="") + 
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))

ggplot(spca_tab,aes(EV1,EV5,color=Pop,shape=Pop)) +
  xlab(paste("PC1 (",round(spca_Per_exp[1,],2),"%)",sep="")) + 
  ylab(paste("PC5 (",round(spca_Per_exp[5,],2),"%)",sep="")) + 
  stat_ellipse(level=0.75,size=2)+
  geom_point(size=5) +
  #geom_point(size=3, subset = .(label == 'point'))+
  scale_shape_manual(name="Pop", labels=unique(spca_tab$Pop)[order(unique(spca_tab$Pop))], values=shp)+
  scale_color_manual(name="Pop", labels=unique(spca_tab$Pop)[order(unique(spca_tab$Pop))], values=cl)+
  labs(color="") + 
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))


#ADMIXTURE visualization

cv<-read.csv("admixture_CV.csv")

ggplot(cv)+
  geom_line(aes(K,CV),size=1)+
  #geom_vline(xintercept = 4,linetype="dashed")+
  scale_x_continuous(breaks = c(1:15))+
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))

#ADMIXTURE visualization with no MAF filtration and random SNP per contig

cv<-read.csv("../ADMIXTURE_CV3.txt",sep="\t")

ggplot(cv)+
  geom_line(aes(K,CV),size=1)+
  #geom_vline(xintercept = 4,linetype="dashed")+
  scale_x_continuous(breaks = c(1:10))+
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))

tbl2=read.table("../P.borealis_stacks.g5mac3dp10ind50g95mdp20p1RandOnePerContig.2.Q")
tbl2$pop<-pop2
row.names(tbl2)<-samp2
tbl2<-tbl2[order(tbl2[,3],tbl2[,1],tbl2[,2]),]
tbl2<-tbl2[c(which(tbl2$pop=="AR"),which(tbl2$pop=="S1"),which(tbl2$pop=="S2"),which(tbl2$pop=="S3"),which(tbl2$pop=="S4"),which(tbl2$pop=="S5"),which(tbl2$pop=="OS")),]
barplot(t(as.matrix(tbl2)), col=c("red","blue"),xlab=NA, ylab="Ancestry", border=NA,las=2,cex.lab=1.5,cex.names = .75,cex.axis = 1.5)

tbl8=read.table("../P.borealis_stacks.g5mac3dp10ind50g95mdp20p1RandOnePerContig.8.Q")
tbl8$pop<-pop2
row.names(tbl8)<-samp2
tbl8<-tbl8[order(tbl8[,9],tbl8[,1],tbl8[,2],tbl8[,3],tbl8[,4],tbl8[,5],tbl8[,6],tbl8[,7],tbl8[,8]),]
tbl8<-tbl8[c(which(tbl8$pop=="AR"),which(tbl8$pop=="S1"),which(tbl8$pop=="S2"),which(tbl8$pop=="S3"),which(tbl8$pop=="S4"),which(tbl8$pop=="S5"),which(tbl8$pop=="OS")),]
barplot(t(as.matrix(tbl8)), col=c("red","orange","yellow","green","lightblue","darkblue","violet","black"),xlab=NA, ylab="Ancestry", border=NA,las=2,cex.lab=1.5,cex.names = .75,cex.axis = 1.5)

Qmean<-data.frame(aggregate(V1~pop,tbl2,FUN=mean),V2=aggregate(V2~pop,tbl2,FUN=mean)$V2)
#Qmean<-data.frame(V1=Qmean$V1,V2=Qmean$V2,pop=c("AR","S3","S4","S1","S2","OS","S5"))
row.names(Qmean)<-Qmean$pop
Qmean<-Qmean[c(which(Qmean$pop=="AR"),which(Qmean$pop=="S1"),which(Qmean$pop=="S2"),which(Qmean$pop=="S3"),which(Qmean$pop=="S4"),which(Qmean$pop=="S5"),which(Qmean$pop=="OS")),]
Qmean$pop<-factor(Qmean$pop,level=Qmean$pop)

ggplot(data = Qmean)+
  geom_col(aes(x=pop,y=V2),fill="orange")+
  #geom_line(inherit.aes = F,aes(x=Qmean$pop,y=c(NA,Qmean$V2[2],Qmean$V2[3],Qmean$V2[4],Qmean$V2[5],Qmean$V2[6],Qmean$V2[7])),group=1) +
  geom_line(inherit.aes = F,aes(x=Qmean$pop,y=c(NA,NA,NA,Qmean$V2[4],sum(Qmean$V2[4],(Qmean$V2[6]-Qmean$V2[4])/2),Qmean$V2[6],Qmean$V2[7])),group=1,size=1.5) +
  xlab("")+
  ylab("Mean offshore ancestry assignment")+
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))

##################################################################
#indmiss<-read.table("postfiltmiss.imiss",sep="\t",header = T)
indmiss<-read.table("../postfiltRandmiss.imiss",sep="\t",header = T)
year2<-c(rep("2018",36),rep("2021",32))
indmiss$Sample<-samp2
indmiss$Pop<-pop2
indmiss$Year<-year2
indmiss<-indmiss[order(indmiss$F_MISS,decreasing = T),]
indmiss$Sample<-factor(indmiss$Sample,levels = indmiss$Sample)

ggplot(data=indmiss)+
  geom_col(aes(x = Sample, y = (F_MISS),fill=Year))+
  geom_abline(slope=0,intercept = 0.5,col="black",lwd=1)+
  xlab("")+
  ylab("Fraction of sites with missing data")+
  ylim(0,1)+
  theme_classic()+
  theme(axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=14,face="bold"),
        axis.text.y = element_text(size=14,face="bold"),
        axis.text.x = element_text(size=14,angle = 270, vjust = 0.2,hjust = 1,face="bold"))

aggregate(N_MISS~Pop,indmiss,FUN=max)
aggregate(N_MISS~Pop,indmiss,FUN=min)
aggregate(F_MISS~Pop,indmiss,FUN=mean)

#Fjarlægð
library(geo)
pos1 <- list(lat = c(65, 66), lon = c(-19, -20))
pos2 <- list(lat = c(62, 65), lon = c(-19, -20))

(dists <- arcdist(pos1, pos2,scale="km"))

NNC <- list(lat = 50.3, lon = 54.1)
SLE <- list(lat = 48.58, lon = 68.58)
ESS <- list(lat = 45.38, lon = 61.03)

arcdist(NNC,SLE,scale="km")         # pos1 and pos2 are lists of coordinates
arcdist(NNC,ESS,scale="km")
arcdist(SLE, ESS,scale="km")

std<-read.csv("stodvar.csv")
std$stod[2]
AR <- list(lat = c(as.numeric(std$kastad_breidd)[2], as.numeric(std$hift_breidd)[2]), lon = c(as.numeric(std$kastad_lengd)[2], as.numeric(std$hift_lengd)[2]))
std$stod[4]
S1<- list(lat = c(as.numeric(std$kastad_breidd)[4], as.numeric(std$hift_breidd)[4]), lon = c(as.numeric(std$kastad_lengd)[4], as.numeric(std$hift_lengd)[4]))
std$stod[3]
S2<- list(lat = c(as.numeric(std$kastad_breidd)[3], as.numeric(std$hift_breidd)[3]), lon = c(as.numeric(std$kastad_lengd)[3], as.numeric(std$hift_lengd)[3]))
std$stod[7]
S3<- list(lat = c(as.numeric(std$kastad_breidd)[7], as.numeric(std$hift_breidd)[7]), lon = c(as.numeric(std$kastad_lengd)[7], as.numeric(std$hift_lengd)[7]))
std$stod[1]
S4<- list(lat = c(as.numeric(std$kastad_breidd)[1], as.numeric(std$hift_breidd)[1]), lon = c(as.numeric(std$kastad_lengd)[1], as.numeric(std$hift_lengd)[1]))
std$stod[6]
S5<- list(lat = c(as.numeric(std$kastad_breidd)[6], as.numeric(std$hift_breidd)[6]), lon = c(as.numeric(std$kastad_lengd)[6], as.numeric(std$hift_lengd)[6]))
std$stod[5]
OS <- list(lat = c(as.numeric(std$kastad_breidd)[5], as.numeric(std$hift_breidd)[5]), lon = c(as.numeric(std$kastad_lengd)[5], as.numeric(std$hift_lengd)[5]))

sta<-data.frame(sta=c("AR","OS","S1", "S2", "S3", "S4","S5"), lat=c(mean(unlist(AR$lat)),mean(unlist(OS$lat)),mean(unlist(S1$lat)),mean(unlist(S2$lat)),mean(unlist(S3$lat)),mean(unlist(S4$lat)),mean(unlist(S5$lat))),lon=c(mean(unlist(AR$lon)),mean(unlist(OS$lon)),mean(unlist(S1$lon)),mean(unlist(S2$lon)),mean(unlist(S3$lon)),mean(unlist(S4$lon)),mean(unlist(S5$lon))))

DIST<-matrix(data=c(
  0,
  mean(arcdist(S1,S2,scale="km")),
  mean(arcdist(S1,S3,scale="km")),
  mean(arcdist(S1,S4,scale="km")),
  mean(arcdist(S1,S5,scale="km")),
  mean(arcdist(S1,OS,scale="km")),
  mean(arcdist(S1,AR,scale="km")),
  #
  mean(arcdist(S2,S1,scale="km")),
  0,
  mean(arcdist(S2,S3,scale="km")),
  mean(arcdist(S2,S4,scale="km")),
  mean(arcdist(S2,S5,scale="km")),
  mean(arcdist(S2,OS,scale="km")),
  mean(arcdist(S2,AR,scale="km")),
  #
  mean(arcdist(S3,S1,scale="km")),
  mean(arcdist(S3,S2,scale="km")),
  0,
  mean(arcdist(S3,S4,scale="km")),
  mean(arcdist(S3,S5,scale="km")),
  mean(arcdist(S3,OS,scale="km")),
  mean(arcdist(S3,AR,scale="km")),
  #
  mean(arcdist(S4,S1,scale="km")),
  mean(arcdist(S4,S2,scale="km")),
  mean(arcdist(S4,S3,scale="km")),
  0,
  mean(arcdist(S4,S5,scale="km")),
  mean(arcdist(S4,OS,scale="km")),
  mean(arcdist(S4,AR,scale="km")),
  #
  mean(arcdist(S5,S1,scale="km")),
  mean(arcdist(S5,S2,scale="km")),
  mean(arcdist(S5,S4,scale="km")),
  mean(arcdist(S5,S3,scale="km")),
  0,
  mean(arcdist(S5,OS,scale="km")),
  mean(arcdist(S5,AR,scale="km")),
  #
  mean(arcdist(OS,S1,scale="km")),
  mean(arcdist(OS,S2,scale="km")),
  mean(arcdist(OS,S4,scale="km")),
  mean(arcdist(OS,S3,scale="km")),
  mean(arcdist(OS,S5,scale="km")),
  0,
  mean(arcdist(OS,AR,scale="km")),
  #
  mean(arcdist(AR,S1,scale="km")),
  mean(arcdist(AR,S2,scale="km")),
  mean(arcdist(AR,S4,scale="km")),
  mean(arcdist(AR,S3,scale="km")),
  mean(arcdist(AR,S5,scale="km")),
  mean(arcdist(AR,OS,scale="km")),
  0),nrow=7)
  
colnames(DIST)<-c("S1","S2","S3","S4","S5","OS","AR")
row.names(DIST)<-c("S1","S2","S3","S4","S5","OS","AR")

FST<-read.csv("FST_WC_dist_RandSNP.csv")
FST<-as.matrix(FST[,-1])
FST[1,-1]<-FST[-1,1]
FST[2,-c(1,2)]<-FST[-c(1,2),2]
FST[3,-c(1,2,3)]<-FST[-c(1,2,3),3]
FST[4,-c(1,2,3,4)]<-FST[-c(1,2,3,4),4]
FST[5,-c(1,2,3,4,5)]<-FST[-c(1,2,3,4,5),5]
FST[6,-c(1,2,3,4,5,6)]<-FST[-c(1,2,3,4,5,6),6]
FST[-7,7]<-FST[7,-7]

FST<-FST/(1-FST)
#Öll gögn
(Mant_tot<-vegan::mantel(FST,DIST,method="spearman",permutations=1000))
plot(density(Mant_tot$perm,),xlab="Mantel",ylab="Density",main="")
abline(v=Mant_tot$statistic, col="red")
ade4::mantel.rtest(as.dist(FST),as.dist(DIST),)
vegan::mantel(DIST,FST,permutations=999)
ade4::mantel.rtest(as.dist(DIST),as.dist(FST))

#2018 - without Arnarfjarðar
(Mant_18<-vegan::mantel(FST[c(1,2,4),c(1,2,4)],DIST[c(1,2,4),c(1,2,4)],method="spearman",permutations=999))
vegan::mantel(DIST[c(1,2,4),c(1,2,4)],FST[c(1,2,4),c(1,2,4)],permutations=999)
plot(density(Mant_18$perm,),xlab="Mantel",ylab="Density",main="")
abline(v=Mant_18$statistic, col="red")

#2021
(Mant_21<-vegan::mantel(FST[c(3,5,6),c(3,5,6)],DIST[c(3,5,6),c(3,5,6)],method="spearman",permutations=999))
vegan::mantel(DIST[c(3,5,6),c(3,5,6)],FST[c(3,5,6),c(3,5,6)],permutations=999)
plot(density(Mant_21$perm,),xlab="Mantel",ylab="Density",main="")
abline(v=Mant_21$statistic, col="red")

#2018 - all
(Mant_18T<-vegan::mantel(FST[c(1,2,4,7),c(1,2,4,7)],DIST[c(1,2,4,7),c(1,2,4,7)],method="spearman",permutations=999))
plot(density(Mant_18T$perm,),xlab="Mantel",ylab="Density",main="")
abline(v=Mant_18T$statistic, col="red")


#All data except Arnarfjarðar
(Mant_sAr<-vegan::mantel(FST[-7,-7],DIST[-7,-7],method="spearman",permutations=1000))
vegan::mantel(DIST[-7,-7],FST[-7,-7],method="spearman",permutations=1000)
ade4::mantel.rtest(as.dist(FST[-7,-7]),as.dist(DIST[-7,-7]))
ade4::mantel.rtest(as.dist(DIST[-7,-7]),as.dist(FST[-7,-7]))
plot(density(Mant_sAr$perm,),xlab="Mantel",ylab="Density",main="")
abline(v=Mant_sAr$statistic, col="red")

#Drawn distance from OS against FST
plot(as.dist(DIST),as.dist(FST))
plot(as.dist(DIST[c(1,2,4),c(1,2,4)]),as.dist(FST[c(1,2,4),c(1,2,4)]))
plot(as.dist(DIST[c(1,2,4,7),c(1,2,4,7)]),as.dist(FST[c(1,2,4,7),c(1,2,4,7)]))
plot(as.dist(DIST[c(4,5,6),c(4,5,6)]),as.dist(FST[c(4,5,6),c(4,5,6)]))

plot(DIST[6,-6],FST[6,-6])
plot(DIST[6,-c(1,2,4,6,7)],FST[6,-c(1,2,4,6,7)])
FST<-as.matrix(FST)

Dist_FST<-cbind(c(DIST),c(FST))
Dist_FST<-as.data.frame(Dist_FST)
colnames(Dist_FST)<-c("km","FST")
Dist_FST$Site<-rep(colnames(DIST),7)
Dist_FST$Site2<-c(rep(colnames(DIST)[1],7),rep(colnames(DIST)[2],7),rep(colnames(DIST)[3],7),rep(colnames(DIST)[4],7),rep(colnames(DIST)[5],7),rep(colnames(DIST)[6],7),rep(colnames(DIST)[7],7))
Dist_FST<-Dist_FST[-c(1,8,9,15,16,17,22,23,24,25,29,30,31,32,33,36:41,43:49),]

ggplot(Dist_FST,aes(km,FST,color=Site))+
  #geom_line(aes(km,FST),linewidth=0.75)+
  #geom_point(aes(km,FST,color=Year),size=4)+
  geom_point(size=4)+
  stat_smooth(method="lm",se=F)+
  #geom_text(aes(km,FST,label=row.names(OSdist[-6,])),nudge_x=-1.5,nudge_y=0.0005,size=5)+
  xlab("Distance (km)")+
  ylab(expression('F'[ST]/(1-'F'[ST])))+
  theme_classic()+
  theme(text = element_text(size=20))

ggplot(Dist_FST,aes(km,FST))+
  #geom_line(aes(km,FST),linewidth=0.75)+
  #geom_point(aes(km,FST,color=Year),size=4)+
  geom_point(size=4)+
  stat_smooth(method="lm",se=F)+
  #geom_text(aes(km,FST,label=row.names(OSdist[-6,])),nudge_x=-1.5,nudge_y=0.0005,size=5)+
  xlab("Distance (km)")+
  ylab(expression('F'[ST]/(1-'F'[ST])))+
  theme_classic()+
  theme(text = element_text(size=20))

Dist_FST<-cbind(c(DIST[-7,-7]),c(FST[-7,-7]))
Dist_FST<-as.data.frame(Dist_FST)
colnames(Dist_FST)<-c("km","FST")
Dist_FST$Site<-rep(colnames(DIST)[-7],6)
Dist_FST$Site2<-c(rep(colnames(DIST)[1],6),rep(colnames(DIST)[2],6),rep(colnames(DIST)[3],6),rep(colnames(DIST)[4],6),rep(colnames(DIST)[5],6),rep(colnames(DIST)[6],6))
Dist_FST<-Dist_FST[-c(1,7,8,13,14,15,19,20,21,22,25,26,27,28,29,31,32,33,34,35,36),]
ggplot(Dist_FST,aes(km,FST,color=Site))+
  #geom_line(aes(km,FST),linewidth=0.75)+
  #geom_point(aes(km,FST,color=Year),size=4)+
  geom_point(size=4)+
  stat_smooth(method="lm",se=F)+
  stat_ellipse(level=0.75,size=2)+
  #geom_text(aes(km,FST,label=row.names(OSdist[-6,])),nudge_x=-1.5,nudge_y=0.0005,size=5)+
  xlab("Distance (km)")+
  ylab(expression('F'[ST]/(1-'F'[ST])))+
  theme_classic()+
  theme(text = element_text(size=20))

ggplot(Dist_FST,aes(km,FST))+
  #geom_line(aes(km,FST),linewidth=0.75)+
  #geom_point(aes(km,FST,color=Year),size=4)+
  geom_point(size=4)+
  stat_smooth(method="lm",se=F)+
  stat_ellipse(inherit.aes = F,data=Dist_FST[Dist_FST$Site%in%c("AR","OS") | Dist_FST$Site2%in%c("AR","OS"),],aes(km,FST,color=Site),level=0.95,size=2)+
  #geom_text(aes(km,FST,label=row.names(OSdist[-6,])),nudge_x=-1.5,nudge_y=0.0005,size=5)+
  xlab("Distance (km)")+
  ylab(expression('F'[ST]/(1-'F'[ST])))+
  theme_classic()+
  theme(text = element_text(size=20))

OSdist<-data.frame(t(rbind(DIST[6,-6],FST[6,-6])))
colnames(OSdist)<-c("km","FST")
OSdist$Year<-factor(c(2018,2018,2021,2018,2021,2018),levels=c(2018,2021))

ggplot(OSdist[-6,])+
  geom_line(aes(km,FST),linewidth=0.75)+
  geom_point(aes(km,FST,color=Year),size=4)+
  geom_text(aes(km,FST,label=row.names(OSdist[-6,])),nudge_x=-1.25,nudge_y=0.00025,size=5)+
  xlab("Distance from offshore site OS (km)")+
  ylab(expression('F'[ST]/(1-'F'[ST])))+
  theme_classic()+
  theme(text = element_text(size=20))

Sdist<-data.frame(t(rbind(DIST[1,],FST[1,])))
colnames(Sdist)<-c("km","FST")
Sdist$Year<-factor(c(2018,2018,2021,2018,2021,2021,2018),levels=c(2018,2021))

ggplot(Sdist[-7,])+
  geom_line(aes(km,FST),linewidth=0.5)+
  geom_point(aes(km,FST,color=Year),size=4)+
  geom_text(aes(km,FST,label=row.names(Sdist[-7,])),nudge_y=0.0009,size=5)+
  xlab("Distance from nearshore")+
  ylab(expression('F'[ST]))+
  theme_classic()+
  theme(text = element_text(size=20))

#Visualize IBD 
plot(log(as.dist(DIST[-7,-7])),as.dist(FST[-7,-7]))

arcdist(lat=AR[[1]][1],lon=AR[[2]][1],lat1=AR[[1]][2],lon1=AR[[2]][2])

env<-std[,c(6:10,13:15,17)]
row.names(env)<-c("S4","AR","S2","S1","OS","S5","S3")
env<-env[c(4,3,7,1,6,5,2),]
FST[is.na(FST)]<-0
FST<-data.frame(FST)
vegan::adonis2(FST~kastad_breidd*kastad_lengd*hift_breidd*hift_lengd*botnhiti*yfirbordshiti*lofthiti*toglengd,data=env, permutations = 99)
vegan::adonis2(FST~botnhiti*yfirbordshiti*lofthiti,data=env)

require(adegenet)
require(StAMPP)
G2 <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno),dimnames = list(pos_loc,colnames(geno)))
G2[geno %in% c("0/0")] <- 0
G2[geno %in% c("0/1", "1/0")] <- 1
G2[geno %in% c("1/1")] <- 2

Gl<-new("genlight",gen=t(G2))
pop(Gl)<-pop2
FST_p<-stamppFst(Gl, nboots = 1000, percent = 95, nclusters = 1)
FST_p$Fsts
FST_p$Pvalues
btstrp<-FST_p$Bootstraps
Pop1<-btstrp[,1]
Pop2<-btstrp[,2]
btstrp<-as.matrix(btstrp[,-c(1,2)])
tmp<-as.numeric(btstrp)
btstrp<-matrix(tmp,nrow=21)
str(btstrp)

require(matrixStats)
bt_summary<-data.frame(Pop1=Pop1, Pop2=Pop2, BootMean=rowMeans(btstrp[,-c(1,2)]),BootMedian=rowMedians(btstrp[,-c(1,2)]),BootVars=rowVars(btstrp[,-c(1,2)]))

######################################################
# check env. data in Skjálfandi
######################################################
x<-read.csv("Skjalfandi_ps.csv")
summary(x$latitude)
summary(x$longitude)

#S1 
S1_env<-x[x$latitude>66.0654 & x$latitude<66.07,]
#S2
S2_env<-x[x$latitude>66.188 & x$latitude<66.25,]
#S3
S3_env<-x[x$latitude>66.265 & x$latitude<66.3,]
#S4
S4_env<-x[x$latitude>66.37 & x$latitude<66.4,]
#S5
S5_env<-x[x$latitude>66.44 & x$latitude<66.5,]

summary(S1_env$salinity)
summary(S2_env$salinity)
summary(S3_env$salinity)
summary(S4_env$salinity)
summary(S5_env$salinity)


################################################################
# Check Bayescan detection
################################################################
#Data prep (convert from VCF to GESTE)
# if (!require("devtools")) install.packages("devtools")
# devtools::install_github("thierrygosselin/radiator")
require(radiator)
require(ggplot2)
kampalampi <- genomic_converter(
  data = "../P.borealis_stacks.g5mac3dp10ind50g95mdp20p1RandOnePerContig.vcf",
  output = c("bayescan"), filename="P.borealis.g5mac3dp10ind50g95mdp20p1RandOnePerContig.geste",
  parallel.core = 1L)

kampalampiP <- genomic_converter(
  data = "../P.borealis_stacks.g5mac3dp10ind50g95mdp20p1RandOnePerContig.vcf",
  strata = "strata.tsv",
  output = c("bayescan"), filename="P.borealis.g5mac3dp10ind50g95mdp20p1RandOnePerContig_pop.geste",
  parallel.core = 1L)

require(coda)
chain<-read.table("../P.borealis.g5mac3dp10ind50g95mdp20p1RandOnePerContig.g_fst.txt",header=TRUE)
chain<-mcmc(chain,thin=10)
#plot(chain)
summary(chain)
autocorr.diag(chain)

BayesMCMC<-read.table("../P.borealis.g5mac3dp10ind50g95mdp20p1RandOnePerContig.g.sel",colClasses = "numeric")

#bayescan=read.table("P.borealis.g5mac3dp10ind50g95mdp20p1.PGD_geste_bayescan_fst.txt") 
bayescan=read.table("../P.borealis.g5mac3dp10ind50g95mdp20p1RandOnePerContig.g_fst.txt") 
#SNPb=read.table("../P.borealis_stacks.g5mac3dp10ind50g95mdp20p1_id_vcf.txt",header=FALSE) 
SNPb<-paste(chromosome,position-1,sep=":")
#SNPbs<-SNPb[-c(127,1449,2010)] #Remove filtered markers
bayescan=cbind(SNPb, bayescan) 
colnames(bayescan)=c("SNP","PROB","LOG_PO","Q_VALUE","ALPHA","FST")  
bayescan$SELECTION <- ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE<=0.05,"diversifying",ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE>0.05,"neutral","balancing")) 
bayescan$SELECTION<- factor(bayescan$SELECTION)
levels(bayescan$SELECTION) 
summary(bayescan$SELECTION)

bayescan$LOG10_Q <- -log10(bayescan$Q_VALUE)

ggplot(bayescan,aes(x=LOG10_Q,y=FST, label=bayescan$POS))+
  geom_point(aes(fill=bayescan$SELECTION), pch=21, size=2)+ 
  #geom_text()+ 
  scale_fill_manual(name="Selection",values=c("white","red","orange"))+ 
  labs(x="Log(q-value)")+ 
  labs(y="Fst")+ 
  theme(axis.title=element_text(size=12, family="Helvetica",face="bold"), legend.position="none")+ 
  theme(axis.text.x=element_text(colour="black"))+ 
  theme(axis.text.y=element_text(colour="black",size=12))+ 
  theme(axis.text.x=element_text(colour="black",size=12))+ 
  theme(panel.border = element_rect(colour="black", fill=NA, size=3),  
        axis.title=element_text(size=18,colour="black",family="Helvetica",face="bold")) +
  theme_classic()

require(pcadapt)
#path_to_file <- "../P.borealis_stacks.g5mac3dp10ind50g95mdp20p1.recode.vcf"
path_to_file <- "../P.borealis_stacks.g5mac3dp10ind50g95mdp20p1RandOnePerContig.vcf"
filename <- read.pcadapt(path_to_file, type = "vcf")
z <- pcadapt(input = filename, K = 20, min.maf = 0.001) 
plot(z, option = "screeplot")
plot(z, option="scores",pop=pop2)
plot(z, option="scores", i = 1, j = 3,pop=pop2)
plot(z, option="scores", i = 1, j = 4,pop=pop2)
z <- pcadapt(input = filename, K = 2, min.maf = 0.001)
summary(z)
plot(z,option="manhattan")
plot(z,option="qqplot")
hist(z$pvalues, zlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(z, option = "stat.distribution")

require(qvalue)
qval <- qvalue(z$pvalues)$qvalues
alpha <- 0.05
outliers_q <- which(qval < alpha)
length(outliers_q)
row.names(geno[outliers_q,])

padj <- p.adjust(z$pvalues,method="BH")
alpha <- 0.05
outliers_BH <- which(padj < alpha)
length(outliers_BH)

padj <- p.adjust(z$pvalues,method="bonferroni")
alpha <- 0.05
outliers_Bf <- which(padj < alpha)
length(outliers_Bf)

row.names(geno[outliers_q,])
row.names(geno[outliers_BH,])
row.names(geno[outliers_Bf,])
#row.names(geno[as.numeric(row.names(OutLoc)),])

bayescan2<-bayescan[bayescan$SNP%in%row.names(geno[outliers_q,]),]
summary(bayescan2$SELECTION)
ggplot(bayescan2,aes(x=LOG10_Q,y=FST, label=SNP))+
  geom_point(aes(fill=SELECTION), pch=21, size=2)+ 
  #geom_text()+ 
  scale_fill_manual(name="Selection",values=c("white","red","orange"))+ 
  labs(x="Log(q-value)")+ 
  ylab(expression('F'[ST]))+
  theme(axis.title=element_text(size=12, family="Helvetica",face="bold"), legend.position="none")+ 
  theme(axis.text.x=element_text(colour="black"))+ 
  theme(axis.text.y=element_text(colour="black",size=12))+ 
  theme(axis.text.x=element_text(colour="black",size=12))+ 
  theme(panel.border = element_rect(colour="black", fill=NA, size=3),  
        axis.title=element_text(size=18,colour="black",family="Helvetica",face="bold")) +
  theme_classic()
  

#Subset of outlier alleles
outliers<-bayescan2[bayescan2$SELECTION=="balancing",]$SNP #SNPAdapt og Bayescan

outliers<-gsub(":","_",outliers)

pos_loc<-paste(chromosome,position-1,sep="_")
require(hierfstat)
G2 <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno),dimnames = list(pos_loc,colnames(geno)))
G2[geno %in% c("0/0")] <- 11
G2[geno %in% c("0/1", "1/0")] <- 12
G2[geno %in% c("1/1")] <- 22

G3<-data.frame(t(G2))
G3$Pop<-pop2
G3$Ind<-samp2
G3<-G3[,c(1100,1101,1:1099)]
#G3<-G3[,c(2365,1:2364)]

AlStat<-basic.stats(G3)
str(AlStat)
FIS<-AlStat$Fis
FIS_out<-FIS[row.names(FIS)%in% paste("X",outliers,sep=""),]
weirdo<-data.frame(Pop=G3$Pop,G3[,colnames(G3)%in%paste("X",outliers,sep="")])
PopFreq<-AlStat$pop.freq
PopFreqOut<-PopFreq[names(PopFreq)%in% paste("X",outliers,sep="")]

#Look at all Bayescan alleles flagged as being under balancing selection

perLoc<-AlStat$perloc
perLocOut<-perLoc[row.names(perLoc)%in% paste("X",outliers,sep=""),]
barplot(perLocOut$Fis, ylab=expression("F"[IS]),names.arg = c(row.names(perLocOut)),xlab = "Outlier alleles")
barplot(perLocOut$Ho, ylab=expression("H"[O]),names.arg = c(row.names(perLocOut)),xlab = "Outlier alleles")

write.table(perLoc,file = "SummaryStatsAll.csv",sep = "\t",row.names = T)
row.names(perLocOut[which(perLocOut$Fis==max(perLocOut$Fis)),])
PopFreqOut[names(PopFreqOut)==row.names(perLocOut[which(perLocOut$Fis==max(perLocOut$Fis)),])]

require(HardyWeinberg)
require(dplyr)
Gout<-G3[,colnames(G3)%in%paste("X",outliers,sep="")]

Ghyb<-G3[G3$Ind%in%c("70_S1","6_S5","43_AR","10_S5","70_S3","76_S3","67_S3"),]
GhybOut<-Ghyb[,colnames(Ghyb)%in%paste("X",outliers,sep="")]
GhybOut
GhybOut$Pop<-Ghyb$Pop
GhybOut$Ind<-Ghyb$Ind
GhybOut<-GhybOut[,c(10,11,1:9)]

for(i in 1:ncol(Gout)){
  z<-c(unname(summary(as.factor(Gout[,i])))[1],unname(summary(as.factor(Gout[,i])))[2],unname(summary(as.factor(Gout[,i])))[3])
  names(z)<-c("AA","AB","BB")
  print(colnames(Gout)[i])
  print(HWExact(z)$pval)
}

require(reshape2)

#When looking at all Bayescan identified alleles, need to remove one fixed allele

PopFreqOutA<-NULL
for(i in 1:length(PopFreqOut)){
  PopFreqOutA<-cbind(PopFreqOutA,PopFreqOut[[i]][2,])
}
colnames(PopFreqOutA)<-names(PopFreqOut)
Pop_afreq<-data.frame(PopFreqOutA)

AF<-data.frame(Pop_afreq)
AF$Pop<-row.names(Pop_afreq)
AF$Pop<-factor(row.names(Pop_afreq),levels=c("AR","S1","S2","S3","S4","S5","OS"))

df <- melt(AF,  id.vars = 'Pop', variable.name = 'Alleles')
df_out<-df[df$Alleles%in%paste("X",outliers,sep=""),]

df_out$Alleles<-factor(df_out$Alleles,levels=df_out[df_out$Pop=="OS",][order(df_out[df_out$Pop=="OS",]$value),2])

ggplot(df_out)+
  #geom_col(aes(Pop,value,group=Alleles,fill=Alleles))+
  geom_bar(aes(Pop,value,group=Alleles,fill=Alleles),stat="identity", width=0.75, position = "dodge")+
  xlab("")+
  #ggtitle("PCAdapt outliers")+
  ylab("Outlier allele frequencies at each sample site")+
  theme_classic()+
  #theme(legend.position = "none",text = element_text(size=20))+
  theme(axis.title.y=element_text(size=15))

ggplot(df_out)+
  geom_line(aes(Pop,value,group=Alleles,colour=Alleles),linewidth=1)+
  xlab("")+
  #ggtitle("PCAdapt outliers")+
  ylab("Outlier allele frequencies at each sample site")+
  theme_classic()+
  #theme(legend.position = "none",text = element_text(size=20))+
  theme(text = element_text(size=20))+
  theme(axis.title.y=element_text(size=15))

######################################################
#Environmental association
######################################################
library(mar)
con <- connect_mar (dbname = "sjor")

llist <- c("B13-2018", "B4-2021")

st.info <- ##get station data
  les_stod(con) %>%
  left_join(les_syni(con), by ='stod_id') %>%
  filter(leidangur %in% c('B13-2018','B4-2021'),
         stod_nr %in% c(1199,1504,1514,153,1513,152,163)) %>%
  select(synis_id, ar, stod_nr, leidangur, toglengd,lon = kastad_lengd, lon1= hift_lengd,lat = kastad_breidd, lat1= hift_breidd,botndypi_kastad, botndypi_hift, togbyrjun, vir_uti, timi, hlerabil, veidarfaeri_nr, toghradi)

st.info<-data.frame(st.info)

st.info$meanDypt<-rowMeans(cbind(st.info$botndypi_kastad,st.info$botndypi_hift))

colnames(st.info)
data.frame(st.info$stod_nr,st.info$leidangur,st.info$ar,st.info$meanDypt)

umhv<-read.csv("umhverfi.csv")
env<-data.frame(sta=umhv$stod, btm_lat=umhv$breiddargrada, btm_tmp=umhv$botnhiti, srf_tmp=umhv$yfirbordshiti, sal=umhv$selta, dpth=umhv$dypt)
#env<-env[c(4,3,6),]  
rownames(env)<-env$sta
env<-env[,-c(1,2)]
#env<-env[c(1,3,4,5,6,7,2),]
barplot(env$btm_tmp,names.arg = row.names(env), main="Bottom temp. (C°)")
barplot(env$srf_tmp,names.arg = row.names(env), main="Surface temp. (C°)")
barplot(env$sal,names.arg = row.names(env), ylim=c(34,35),main="P.O.S.")
barplot(env$dpth,names.arg = row.names(env), main="Depth (m)")

EnvCor<-cor(env,method = "spearman")
diag(EnvCor)<-0
corrplot::corrplot(EnvCor,type="lower")

env<-env[,-4] #removed Depth due to high level of negative correlation with temperature and positive correlation with salinity

require(vegan)
Kampa.rda <- rda(Pop_afreq~., data=env)
Kampa.rda
ordiplot(Kampa.rda,type="text") #species points are shown as crosses, environmental points are shown as circles
#Constrained Proportion: variance of Y explained by X (68%)
#Unconstrained Proportion: unexplained variance in Y (32%)
#"The included environmental variables explain 68% of the variation in fish community composition across sites.”

# "inertia" = general term for variation in the data
# "species" = term for the variables
# "sites" = term for the observations

KampaLim.rda <- rda(Pop_afreq~btm_tmp + srf_tmp + sal, data=env)
KampaLim.rda

KampaLimT.rda <- rda(Pop_afreq~btm_tmp + srf_tmp, data=env)
KampaLimT.rda
ordiplot(KampaLimT.rda,type="text")

KampaSrf.rda <- rda(Pop_afreq~srf_tmp, data=env)
KampaSrf.rda
ordiplot(KampaSrf.rda,cex=3,type="text")

KampaBtm.rda <- rda(Pop_afreq~btm_tmp, data=env)
KampaBtm.rda
ordiplot(KampaBtm.rda,type="text")

KampaDpth.rda <- rda(Pop_afreq~dpth, data=env)
KampaDpth.rda
ordiplot(KampaDpth.rda,type="text")

KampaSal.rda <- rda(Pop_afreq~sal, data=env)
KampaSal.rda
ordiplot(KampaSal.rda)

AF<-data.frame(Pop_afreq)
AF$Pop<-row.names(Pop_afreq)
AF$Pop<-factor(row.names(Pop_afreq),levels=c("AR","S1","S2","S3","S4","S5","OS"))

df <- melt(AF,  id.vars = 'Pop', variable.name = 'Alleles')
df$Pop <- factor(df$Pop, levels=c("AR","S1","S2","S3","S4","S5","OS"))
dfOut<-df[df$Alleles%in%row.names(R2)[row.names(R2)%in%paste("X",outliers,sep="")],]

ggplot(dfOut)+
  geom_line(aes(Pop,value,group=Alleles,colour=Alleles),linewidth=1)+
  xlab("")+
  ylab("Environment associated allele frequencies")+
  theme_classic()+
  theme(text = element_text(size=20))+
  theme(axis.title.y=element_text(size=15))

ggplot(dfOut[dfOut$Pop%in%c("AR","S1","S2","S4"),])+
  geom_line(aes(Pop,value,group=Alleles,colour=Alleles),linewidth=1)+
  xlab("")+
  ylab("Temperature associated allele frequencies")+
  theme_classic()+
  theme(legend.position = "none",text = element_text(size=20))+
  theme(axis.title.y=element_text(size=15))

ggplot(dfOut[dfOut$Pop%in%c("S3","S5","OS"),])+
  geom_line(aes(Pop,value,group=Alleles,colour=Alleles),linewidth=1)+
  xlab("")+
  ylab("Temperature associated allele frequencies")+
  theme_classic()+
  theme(legend.position = "none",text = element_text(size=20))+
  theme(axis.title.y=element_text(size=15))

perLocAdapt<-perLoc[row.names(perLoc)%in%row.names(R2)[(row.names(R2)%in%paste("X",outliers,sep=""))],]
barplot(perLocAdapt$Fis, ylab=expression("F"[IS]),names.arg = c(row.names(perLocAdapt)),xlab = "Outlier alleles")
barplot(perLocAdapt$Ho, ylab=expression("H"[O]),ylim=c(0,0.5), names.arg = c(row.names(perLocAdapt)),xlab = "Outlier alleles")
barplot(perLocAdapt$Hs, ylab=expression("H"[S]),ylim=c(0,0.5),names.arg = c(row.names(perLocAdapt)),xlab = "Outlier alleles")

###################################################################################################
#L50 gögn frá OS, S1, og S4
###################################################################################################
setwd("C:/Users/aki/Desktop/P.borealis/")
l50<-read.table("l50.skjalfandi.csv",sep=";",header=T)
l50_S1<-l50[l50$reg=="S1",]
l50_S4<-l50[l50$reg=="S4",]
l50_OS<-l50[l50$reg=="OS",]
??kolmogorov
ks.test(l50_S1$L50,l50_S4$L50)
ks.test(l50_S1$L50,l50_OS$L50)
ks.test(l50_S4$L50,l50_OS$L50)
ks.test(l50_S1$L50,l50_S4$L50,alternative="greater")
ks.test(l50_S4$L50,l50_S1$L50,alternative="greater")
ks.test(l50_S1$L50,l50_OS$L50,alternative="greater")
ks.test(l50_OS$L50,l50_S1$L50,alternative="greater")
ks.test(l50_S4$L50,l50_OS$L50,alternative="greater")
ks.test(l50_OS$L50,l50_S4$L50,alternative="greater")

t.test(l50_S4$L50,l50_S1$L50,paired = T)
t.test(l50_S1$L50,l50_S4$L50,paired = T,alternative="greater")
t.test(l50_S1$L50,l50_S4$L50,paired = T,alternative="less")

t.test(l50_OS$L50,l50_S1$L50,paired = T)
t.test(l50_OS$L50,l50_S4$L50,paired = T)
t.test(l50_OS$L50,l50_S1$L50,alternate="g",paired = T)
t.test(l50_OS$L50,l50_S4$L50,alternate="g",paired = T)
