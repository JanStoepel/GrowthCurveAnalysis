#' Inputs: csv
#' -------
#' output: File supplement/MainEffects.txt and supplement/InterActions.txt
#' output: Graphic file fig/2014-11-24/Interaction.TM35.TM05.pdf
library(gdata)
library(ggplot2)
library(plyr)
library(dplyr)
library(stringr)
library(tidyr)
FP <- file.path("YOUR_DESTINATION_FOLDER_HERE","YOUR_DESTINATION_NAME_HERE")
AUC <- read.csv("YOUR_SOURCE_DATA_HERE.csv",as.is=TRUE)
Layout <- read.csv("YOUR_LAYOUT_FILE_HERE.csv",as.is=TRUE)
X <- merge(Layout,AUC)
# Remove blanks and set up factors to indicate deletions
X <- filter(X,ORF!="Blank")
X <- mutate(X,CTF8Deleted=grepl("ctf8",word(ORF,-1)),OrfDeleted=word(ORF))
X$OrfDeleted[!grepl("Y",X$OrfDeleted)] <- "None"
X <- filter(X,!(OrfDeleted=="YGL066W∆" & Plate=="GC_07"))
X$CTF8Deleted <- factor(X$CTF8Deleted)
X$OrfDeleted <- factor(X$OrfDeleted)
# fit the model and extract summaries
library(lmerTest)
m1 <- lmer(TM35-TM05~CTF8Deleted*OrfDeleted+(1|Plate),data=X)
#m2 <- lmer(log10(TM35-TM05)~CTF8Deleted*OrfDeleted+(1|Plate),data=X)
z <- summary(m1)
z0 <- coef(z)
z1 <- confint(m1,method="Wald")
zz <- cbind(z0[,c(1,5)],z1)[,c(1,3,4,2)]
# write the main effects to a file
File <- "supplement/MainEffects.txt"
XX <- zz[1:(nlevels(X$OrfDeleted)+1),]
rownames(XX)[1:2] <- c("WT Level","CTF8")
rownames(XX) <- gsub("OrfDeleted","",rownames(XX))
rownames(XX)[-1] <- paste(rownames(XX)[-1],"- WT")
XX <- data.frame(Effect=rownames(XX),XX,check.names=FALSE)
names(XX)[5] <- "Pvalue"
XX$Qvalue <- p.adjust(XX[[5]],method="fdr")
write.table(XX,File,quote=F,row.names=FALSE,sep="\t")
# write the interactions to a file
File <- "supplement/InterActions.txt"
XX <- zz[nlevels(X$OrfDeleted)+2:nlevels(X$OrfDeleted),]
rownames(XX) <- gsub("CTF8DeletedTRUE:OrfDeleted","",rownames(XX))
XX <- data.frame("CTF8 and"=rownames(XX),XX,check.names=FALSE)
names(XX)[5] <- "Pvalue"
XX$Qvalue <- p.adjust(XX[[5]],method="fdr")
write.table(XX,File,quote=F,row.names=FALSE,sep="\t")
# create graphic of the interactions
XX <- zz[nlevels(X$OrfDeleted)+2:nlevels(X$OrfDeleted),1:3]
XX <- cbind(XX,zz[3:(nlevels(X$OrfDeleted)+1),1:3])
rownames(XX) <- gsub("CTF8DeletedTRUE:OrfDeleted","",rownames(XX))
XX <- as.data.frame(XX[order(XX[,1]),])
names(XX) <- c("Interaction","LB","UB","MainEffect","LBm","UBm")
XX <- cbind(Orf=factor(rownames(XX),levels=rownames(XX)),XX)
XX$WT <- 0
XX$CTF8 <- zz[2,1]
XX$LBc <- zz[2,2]
XX$UBc <- zz[2,3]
X2 <- gather(select(XX,-contains("B")),What,Value,Interaction,MainEffect,CTF8,WT)
X1 <- read.xls("cleanData/2014-11-03/GC_JS_2015_01.xlsx",sheet=1,as.is=TRUE,skip=1)[[1]]
X0 <- as.character(X2$Orf)
X0 <- substring(X0,1,nchar(X0)-1)
ii <- X0 %in% X1
X2 <- droplevels(X2[ii,])
X0 <- as.character(XX$Orf)
X0 <- substring(X0,1,nchar(X0)-1)
ii <- X0 %in% X1
XX <- droplevels(XX[ii,])
P1 <- ggplot(X2,aes(Value,Orf)) + geom_segment(aes(x=LB,xend=UB,y=Orf,yend=Orf),data=XX,size=1.5) +
geom_point(aes(colour=What,shape=What),size=3) + xlab("Effect on AUC TM35-TM05") +
ggtitle("Interaction between Orf and CTF8 (with 95%CI)")
ggsave("Interaction.TM35.TM05.pdf",P1,path=FP,height=20,width=10)
ggsave("Interaction.TM35.TM05.png",P1,path=FP,height=20,width=10)
P2 <- ggplot(X2,aes(Value,Orf)) + geom_segment(aes(x=LBc,xend=UBc,y=Orf,yend=Orf),colour="blue",data=XX,size=0.5) +
geom_segment(aes(x=LBm,xend=UBm,y=Orf,yend=Orf),colour="green",data=XX,size=0.5) +
geom_segment(aes(x=LB,xend=UB,y=Orf,yend=Orf),data=XX,size=1,colour="red") +
geom_point(aes(colour=What,shape=What),size=4) + xlab("Effect on AUC TM35-TM05") +
ggtitle("Interaction between Orf and CTF8 (with 95%CI)")
ggsave("Interaction.TM35.TM05.alt.pdf",P2,path=FP,height=20,width=10)
ggsave("Interaction.TM35.TM05.alt.png",P2,path=FP,height=20,width=10)
