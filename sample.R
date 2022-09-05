#the first data set 
load("x") # x_{ij} in eq. (19)
load("x1") # x_{kj} in eq. (19)
index <- match(substring(colnames(x1),1,20), substring(colnames(x),1,20)) #select common samples
class <- c(rep("T",255),rep("N",71)) #define class label, T:Tumor, N:normal tissue
x1 <- x1[,!is.na(index)] #select common samples in x_{ij} with x_{kj}
class <- class[!is.na(index)] #ignore missing class labels
index <- index[!is.na(index)] #ignore missing samples
class0 <- rep(-1,length(class)) #numerize labels
x1[,-1] <- x1[,-1] - min(x1[,-1]) #make values positive 
class0[class=="T"] <- 1 #numerize labels for tumor 
class0[class0>0] <- class0[class0>0]/sum(class0>0) #y_j in eq. (29) for normal tissues
class0[class0<0] <- class0[class0<0]/sum(class0<0) #y_j in eq. (29) for tumors
#---- compute SVD ---
Z <-matrix(0,dim(x1)[1],dim(x)[1])
for (i in c(2:length(index)))
{
    cat(i, " ")
    Z <- Z + outer(as.vector(x1[,i]),as.vector(x[,index[i]]),"*")
} 
SVD <- svd(Z)
#----- compute SVD ----- 
U <-t(x1[,-1]) %*%  SVD$u 
V <-t(x[,index[-1]]) %*%  SVD$v
P1 <- pchisq(scale(data.matrix(x[,index][,-1]) %*% class0[-1])^2,1,lower.tail=F) 
P <- pchisq(scale(SVD$v[,2])^2,1,lower.tail=F)
pdf(file="qqplot_TCGA_mRNA.pdf") #Fig. 4(A)
qqplot(P,P1,pch=16,cex=0.5,xlab="TD",ylab="projection");abline(0,1,col=2)
dev.off()
table(p.adjust(P,"BH")<0.01,p.adjust(P1,"BH")<0.01) #Table 1
#      FALSE  TRUE
#FALSE 19447    17
#TRUE     11    61



P1 <- pchisq(scale(scale(data.matrix(x1[,-1])) %*% class0[-1])^2,1,lower.tail=F) 
P <- pchisq(scale(SVD$u[,2])^2,1,lower.tail=F)
pdf(file="qqplot_TCGA_miRNA.pdf") #Fig. 4(B)
qqplot(P,P1,pch=16,cex=0.5,xlab="TD",ylab="projection");abline(0,1,col=2)
dev.off()
table(p.adjust(P,"BH")<0.01,p.adjust(P1,"BH")<0.01) #Table 2

#       FALSE TRUE
#FALSE   812    2
#TRUE      0   11

#shffuled set 
U_all <- NULL
V_all <- NULL
for (i in c(1:100)) #it will take several hours 
{
    cat("\n",i," ")
    Z0 <-matrix(0,dim(x1)[1],dim(x)[1])
    x10 <- apply(x1,2,sample)
    x0 <- apply(x[,index],2,sample)
    for (ii in c(2:length(index)))
    {
        cat(ii, " ")
        Z0 <- Z0 + outer(as.numeric(as.vector(x10[,ii])),as.numeric(as.vector(x0[,ii])),"*")
    } 
    SVD <-svd(Z0)
    U_all <- cbind(U_all,SVD$u[,2])
    V_all <- cbind(V_all,SVD$v[,2])
  
}
SVD <- svd(Z)
SUM <- apply(abs(x1[,-1]),1,sum)
#RANK <- rank(-(scale(c(SVD$u[,2],unlist(U_all))))^2)[1:length(SVD$u[,2])] #for all
RANK <- rank(-(scale(c(SVD$u[rank(-SUM)<=500,2],unlist(U_all))))^2)[1:length(SVD$u[rank(-SUM)<=500,2])] #for top expressed 500 
P0 <- RANK/(length(U_all)*1.01)
#pdf(file="KA-Lok-2_hist_miRNA.pdf") #for all, Fig 1(A)
pdf(file="KA-Lok-2_hist_miRNA_top_500.pdf") #for top expressed 500, Fig. 1(B)
hist(1-P0,breaks=100,xlab="1-P",main="")
dev.off()
P <- pchisq(scale(SVD$u[,2])^2,1,lower.tail=F)
table(p.adjust(P0,"BH")<0.1,(p.adjust(P,"BH")<0.01)[order(-SUM)<=500]) #Table 6 
#       FALSE TRUE
#FALSE   488    0
#TRUE      1   11

SUM <- apply(abs(x[,index][,-1]),1,sum)
RANK <- rank(-(scale(c(SVD$v[rank(-SUM)<=3000,2],unlist(V_all))))^2)[1:length(SVD$v[rank(-SUM)<=3000,2])] #for  top expressed 3000
#RANK <- rank(-(scale(c(SVD$v[,2],unlist(V_all))))^2)[1:length(SVD$v[,2])] #for all
P0 <- RANK/(length(V_all)*1.01)
pdf(file="KA-Lok-2_hist_mRNA.pdf") #for all, Fig 2(A)
#pdf(file="KA-Lok-2_hist_mRNA_top_3000.pdf") #for top expressed 3000, Fig. 2(B)
hist(1-P0,breaks=100,xlab="1-P",main="")
dev.off()


P <- pchisq(scale(SVD$v[,2])^2,1,lower.tail=F)
table(p.adjust(P0,"BH")<0.1,(p.adjust(P,"BH")<0.01)[rank(-SUM)<=3000]) #Table 7
#     FALSE TRUE
#FALSE  2928    3
#TRUE      0   69

#second data set 
#download two files from GSD16441 to the current directory
x <- read.csv("GSE16441-GPL6480_series_matrix.txt.gz",sep="\t",comment.char="!")
x1 <- read.csv("GSE16441-GPL8659_series_matrix.txt.gz",sep="\t",comment.char="!")
#----- compute SVD ------
Z <-matrix(0,dim(x1)[1],dim(x)[1])
for (i in c(2:35))
{
    cat(i, " ")
    Z <- Z + outer(as.vector(x1[,i]),as.vector(x[,i]),"*")
} 
SVD <- svd(Z)
#---- compute SVD ----
class <- rep(c(-1,1),each=17) #class labels 
P1 <- pchisq(scale(data.matrix(x1[,-1]) %*% class)^2,1,lower.tail=F)
P<- pchisq(scale(SVD$u[,2])^2,1,lower.tail=F)
pdf(file="qqplot_GEO_miRNA.pdf") #Fig 4(D)
qqplot(P,P1,pch=16,cex=0.5,xlab="TD",ylab="projection");abline(0,1,col=2)
dev.off()
table(p.adjust(P,"BH")<0.01,p.adjust(P1,"BH")<0.01) #Table 4
#       FALSE TRUE
#FALSE   316    0
#TRUE      0    3

P1 <- pchisq(scale(data.matrix(x[,-1]) %*% class)^2,1,lower.tail=F)
P<- pchisq(scale(SVD$v[,2])^2,1,lower.tail=F)
pdf(file="qqplot_GEO_mRNA.pdf") #Fig 4(C)
qqplot(P,P1,pch=16,cex=0.5,xlab="TD",ylab="projection");abline(0,1,col=2)
dev.off()
table(p.adjust(P,"BH")<0.01,p.adjust(P1,"BH")<0.01) #Table 3
#      FALSE  TRUE
#FALSE 33781     8
#TRUE     23   186

#shuffled data set 
U_all <- NULL
V_all <- NULL
for (i in c(1:100)) #it will take one to two hours 
{
    cat("\n",i," ")
    Z0 <-matrix(0,dim(x1)[1],dim(x)[1])
    x10 <- apply(x1,2,sample)
    x0 <- apply(x,2,sample)
    for (ii in c(2:35))
    {
        cat(ii, " ")
        Z0 <- Z0 + outer(as.numeric(as.vector(x10[,ii])),as.numeric(as.vector(x0[,ii])),"*")
    } 
    SVD <-svd(Z0)
    U_all <- cbind(U_all,SVD$u[,2])
    V_all <- cbind(V_all,SVD$v[,2])
}
SVD <- svd(Z)
RANK <- rank(-(scale(c(SVD$u[,2],unlist(U_all))))^2)[1:length(SVD$u[,2])]
P0 <- RANK/(length(U_all)*1.01)
pdf(file="Ka-Lok-2_hist_miRNA_GEO.pdf") #Fig. 6(A)
hist(1-P0,breaks=100,xlab="1-P",main="")
dev.off()

RANK <- rank(-(scale(c(SVD$v[,2],unlist(V_all))))^2)[1:length(SVD$v[,2])]
P0 <- RANK/(length(V_all)*1.01)
pdf(file="Ka-Lok-2_hist_mRNA_GEO.pdf") #Fig. 6(B)
hist(1-P0,breaks=100,xlab="1-P",main="")
dev.off()

P<- pchisq(scale(SVD$v[,2])^2,1,lower.tail=F)
table(p.adjust(P0,"BH")<0.1,p.adjust(P,"BH")<0.01) #Table 8
#    FALSE  TRUE
#FALSE 33736     0
#TRUE     53   209

#the third data set

#download a file to the current directory from GSE147507 
x <- read.csv("GSE147507_RawReadCounts_Human.tsv.gz",sep="\t",comment.char="!")

load("Z") #x_{ijkm} in eq. (15)
Z <- apply(Z,2:4,scale) #standardize x_{ijkm} 
require(rTensor)
HOSVD <- hosvd(as.tensor(Z),c(10,5,2,3)) #apply HOSVD as in eq. (15)
P <- pchisq(scale(HOSVD$U[[1]][,5])^2,1,lower.tail=F) 
A1 <- outer(rep(1,5),outer(c(1,-1),rep(1,3)))
dim(A1)<- 5*2*3
dim(Z) <- c(21797,30)
PP <- pchisq(scale(Z %*% A1)^2,1,lower.tail=F)
table(p.adjust(P,"BH")<0.01,p.adjust(PP,"BH")<0.01) #Table 5
#     FALSE  TRUE
#FALSE 21582    52
#TRUE     60   103
pdf(file="qqplot_SARS-CoV-2.pdf") #Fig. 5
qqplot(P,PP,pch=16,cex=0.5,xlab="TD",ylab="projection");abline(0,1,col=2)
dev.off()

#shuffled data set 
dim(Z) <- c(21797,5,2,3)
set.seed(0)
U_all <- NULL
for (i in c(1:100))
{
    cat(i," ")
    Z0 <- apply(Z,2:4,sample)
    HOSVD0 <- hosvd(as.tensor(Z0),c(5,1,1,1))
    U_all <- cbind(U_all,HOSVD0$U[[1]][,5])
}
SUM <- apply(abs(Z),1,sum)
RANK <- rank(-(scale(c(HOSVD$U[[1]][rank(-SUM)<=2780,5],unlist(U_all))))^2)[1:length(HOSVD$U[[1]][rank(-SUM)<=2780,5])] #for top expressed 2780 
#RANK <- rank(-(scale(c(HOSVD$U[[1]][,5],unlist(U_all))))^2)[1:length(HOSVD$U[[1]][,5])] #for all
P0 <- RANK/(length(U_all)*1.01)
#pdf(file="SARS-CoV-2_hist_top2780.pdf") #for top expressed 2780, Fig. 3(A)
pdf(file="SARS-CoV-2_hist.pdf") #for all, Fig. 3(B)
hist(1-P0,breaks=100,xlab="1-P",main="")
dev.off()

table(p.adjust(P0,"BH")<0.1,(p.adjust(P,"BH")<0.01)[rank(-SUM)<=2780]) #Table 9
#    ALSE TRUE
#FALSE  2617  115
#TRUE      0   48
