#!/usr/bin/env Rscript
C <- read.csv("data/Compound.csv", header=TRUE)
P <- read.csv("data/Protein.csv", header=TRUE)
Y <- read.csv("data/Activity.csv", header=TRUE)
SMOTE_Y <- read.csv("data/Balanced_Activity.csv", header=TRUE)
library(C50)
library(unbalanced)
library(corrplot)
library(caret)
library(prospectr)
library(gridExtra)
library(Rcpi)

# This is where you should use for loop and please create a vector to
# store results instead of a,b,c,d...  Writing code this way is
# error-prone and very hard to debug because the code is redundant.

# At first glance, the code writing style is totally inconsistent.
# Please consult
# https://google-styleguide.googlecode.com/svn/trunk/Rguide.xml#filenames
# and try to do your best to follow the guideline.

# Balancing Data with SMOTE Algorithm
set.seed(1)
results <- list(10)
for ( i in 1:10) {
  results[[i]] <- ubBalance(X=cbind(C,P), Y=factor(Y$Before_SMOTE),
                           type="ubSMOTE", perc.over = 100,
                           perc.under=200, verbose=TRUE)  
}
set.seed(4)
smote_results <- lapply(results, function(x) {
  a <- cbind(x$Y, x$X)
  a <- a[order(x$Y),]
  x <- a[-1]
  return(x)

})
average_SMOTE <- Reduce("+", smote_results) / length(smote_results)

## Removeing Interacorrelation with the cutoff of 0.7
set.seed(200)
compoundSMOTE <- average_SMOTE[, 1:44]
proteinSMOTE <- average_SMOTE[, 45:443]
<<<<<<< HEAD
=======
set.seed(400)
>>>>>>> 2eba7529991c26e719b83ffb0fbd808e871e6636
compound_corr <- cor(compoundSMOTE)
dim(compound_corr)
compound_corr <- compound_corr[1:44, 1:44]
compound_highCorr <- findCorrelation(compound_corr, cutoff = .70)
length(compound_highCorr)
highCorrRemoveCompound <- compoundSMOTE[, -compound_highCorr]
Balanced_Activity <- SMOTE_Y[,1]
InputKennardCompound <- cbind(Balanced_Activity, highCorrRemoveCompound)
proteinSMOTE = proteinSMOTE[, -nearZeroVar(proteinSMOTE)]
protein_corr <- cor(proteinSMOTE)
dim(protein_corr)
protein_corr <- protein_corr[1:390, 1:390]
protein_highCorr <- findCorrelation(protein_corr, cutoff = .70)
length(protein_highCorr)
highCorrRemoveProtein <- proteinSMOTE[, -protein_highCorr]
highCorrRemoveProtein <- data.frame(highCorrRemoveProtein)
InputKennardProtein <- cbind(Balanced_Activity, highCorrRemoveProtein)
C_Inactive <- subset(InputKennardCompound, SMOTE_Y$Balanced_Activity == "Inactive")
C_Active <- subset(InputKennardCompound, SMOTE_Y$Balanced_Activity == "Active")
P_Inactive <- subset(InputKennardProtein, SMOTE_Y$Balanced_Activity == "Inactive")
P_Active <- subset(InputKennardProtein, SMOTE_Y$Balanced_Activity == "Active")
# Spliting the data set into training set and testing test using Kennard Stone Algorithm
<<<<<<< HEAD
set.seed(3422)
=======
>>>>>>> 2eba7529991c26e719b83ffb0fbd808e871e6636
data <- list(C_Inactive = C_Inactive,
             C_Active = C_Active,
             P_Inactive = P_Inactive,
             P_Active = P_Active)
kenStone <- lapply(data, function(x) {
  sel <- kenStone(x[-1], k = 70, metric="maha", pc=2)
  train <- x[sel$model, ]
  test <- x[sel$test, ]
  return(list(train=train,test=test))  
})
<<<<<<< HEAD

C_Train <- rbind(kenStone$C_Inactive$train, kenStone$C_Active$train)
C_Test <- rbind(kenStone$C_Inactive$test, kenStone$C_Active$test)
P_Train <- rbind(kenStone$P_Inactive$train, kenStone$P_Active$train)
P_Test <- rbind(kenStone$P_Inactive$test, kenStone$P_Active$test)
# Preparing input data for Proteocheometrics modeling
=======
# Preparing input data for Proteocheometrics modeling
C_Train <- rbind(kenStone$C_Inactive$train, kenStone$C_Active$train)
C_Test <- rbind(kenStone$C_Inactive$test, kenStone$C_Active$test)
P_Train <- rbind(kenStone$P_Inactive$train, kenStone$P_Active$train)
P_Test <- rbind(kenStone$P_Inactive$test, kenStone$P_Active$test)
>>>>>>> 2eba7529991c26e719b83ffb0fbd808e871e6636
set.seed(22)
c <- C_Train
p <- P_Train
train_c <- c[-1]
train_p <- p[-1]
crossTrain <- getCPI(train_c, train_p, type = "tensorprod")
dfcompound <- names(data.frame(C_Train[,2:24]))
dfprotein <- names(data.frame(P_Train[,2:86]))
CompoundNamecross <- rep(dfcompound, each=85)
ProteinNamecross <- rep(dfprotein,times=23)
label <- paste(CompoundNamecross, ProteinNamecross, sep="-")
colnames(crossTrain) <- label
Balanced_Activity_Train <- C_Train[1]
Balanced_Activity_Test <- C_Test[1]
CxP_Train <- cbind(Balanced_Activity_Train, data.frame(crossTrain))
protein_selfcrosstrain <- getCPI(train_p, train_p, type = "tensorprod")
df <- names(data.frame(P_Train[,2:86]))
ProteinName2 <- rep(df, times= 85)
ProteinName1 <- rep(df, each=85)
Label <- paste(ProteinName1, ProteinName2, sep="-")
colnames(protein_selfcrosstrain) <- Label
protein_selfcrosstrain <- data.frame(protein_selfcrosstrain)
index <- c(1, 87, 173, 259, 235, 431, 517, 603, 689, 775, 861, 947,
           1033, 1119, 1205, 1291, 1377, 1463, 1549, 1635, 1721, 1807,
           1893, 1979, 2065, 2151, 2237, 2323, 2409, 2495, 2581, 2667,
           2753, 2839, 2925, 3011, 3097, 3183, 3269, 3355, 3441, 3527, 
           3613, 3699, 3785, 3871, 3957, 4043, 4129, 4215, 4301, 4387,
           4473, 4559, 4645, 4731, 4817, 4903, 4989, 5075, 5161, 5247,
           5333, 5419, 5505, 5591, 5677, 5763, 5849, 5935, 6021, 6107,
           6193, 6279, 6365, 6451, 6537, 6623, 6709, 6795, 6881, 6967,
           7053, 7139, 7225)
protein_selfcrosstrain <- protein_selfcrosstrain[, -index]
TranposedIndexed_protein <- t(protein_selfcrosstrain)
index1 <- which(duplicated(TranposedIndexed_protein))
removed_protein_train <- TranposedIndexed_protein[-index1, ]
protein_finalselfcrosstrain <- t(removed_protein_train)
PxP_Train <- cbind(Balanced_Activity_Train, data.frame(protein_finalselfcrosstrain))
compound_selfcrosstrain <- getCPI(train_c, train_c, type = "tensorprod")
df <- names(data.frame(C_Train[,2:24]))
CompoundName2 <- rep(df, times=23)
CompoundName1 <- rep(df, each=23)
label <- paste(CompoundName1, CompoundName2, sep="-")
colnames(compound_selfcrosstrain) <- label
index <- c(1, 25, 49, 73, 97, 121, 145, 169, 193,
          217, 241, 265, 289, 313, 337, 361, 385,
          409, 433, 457, 481, 505, 529)
compound_selfcrosstrain <- compound_selfcrosstrain[, -index]
TranposedIndexed_compound <- t(compound_selfcrosstrain)
index3 <- which(duplicated(TranposedIndexed_compound))
removed_compound_train <- TranposedIndexed_compound[-index3, ]
compound_finalselfcrosstrain <- t(removed_compound_train)
CxC_Train <- cbind(Balanced_Activity_Train, data.frame(compound_finalselfcrosstrain))

## Feature Importance
input <- list(C_Train = C_Train,
              P_Train = P_Train,
              CxP_Train = CxP_Train, 
              PxP_Train = PxP_Train,
              CxC_Train = CxC_Train)
# C5.0 Algorithm is used to obtain feature importance of BT chemicals and proteins
set.seed(34)
Importance <- lapply(input, function(x) {
  data <- data.frame(x)
  Model <- C5.0(Balanced_Activity~., data = data, rules=TRUE)
  Importance <- C5imp(Model)
  return(Importance)
})
set.seed(333)
compound <- data.frame(Importance$C_Train)
protein <- data.frame(Importance$P_Train)
CxP <- data.frame(Importance$CxP_Train)
PxP <- data.frame(Importance$PxP_Train)
CxC <- data.frame(Importance$CxC_Train)
set.seed(1)
compoundtop10 <- head(compound,10) # top 10 compound
proteintop10 <- head(protein,10) # top 10 compound
CxPtop10 <- head(CxP,10) # top 10 Cross-terms
PxPtop10 <- head(PxP,10) # top 10 Protein-Protein Cross-terms
CxCtop10 <- head(CxC,10) # top 10 Compound-Compound Cross-terms
set.seed(2)
compound_labels <- c("SubFP5",
                     "SubFP133",
                     "SubFP8",
                     "SubFP151",
                     "SubFP1",
                     "SubFP295",
                     "SubFP12",
                     "SubFP16",
                     "SubFP28",
                     "SubFP84")
top10compound <- cbind(compound_labels, compoundtop10)
set.seed(3)
myDF2 <- cbind(Protein = rownames(proteintop10, proteintop10))
top10protein <- cbind(myDF2, proteintop10)
crossterms_labels <- c("SubFP133-X5z3scl2.lag18",
                       "SubFP84-z3scl2.lag4",
                       "SubFP84-X5z3scl3.lag14",
                       "SubFP96-X3z3scl3.lag14",
                       "SubFP1-z3scl1.lag3",
                       "SubFP1-z3scl1.lag12",
                       "SubFP1-z3scl1.lag13",
                       "SubFP1-z3scl1.lag15",
                       "SubFP1-z3scl1.lag19",
                       "SubFP1-z3scl1.lag20")
top10cross_terms <- cbind(crossterms_labels, CxPtop10)
PxP_labels <- c("X2z3scl2.lag6-X4z3scl2.lag6",
                "z3scl2.lag16-X5z3scl1.lag21",
                "z3scl1.lag3-z3scl1.lag12",
                "z3scl1.lag3-z3scl1.lag13",
                "z3scl1.lag3-z3scl1.lag15",
                "z3scl1.lag3-z3scl1.lag19",
                "z3scl1.lag3-z3scl2.lag20",
                "z3scl1.lag3-z3scl2.lag22",
                "z3scl1.lag3-z3scl2.lag23",
                "z3scl1.lag3-z3scl2.lag25")
top10PxP <- cbind(PxP_labels, PxPtop10)
CxC_labels <- c("SubFP16-SubFP143",
                "SubFP85-SubFP303",
                "SubFP296-SubFP300",
                "SubFP143-SubFP224",
                "SubFP12-SubFP302",
                "SubFP5-SubFP296",
                "SubFP12-SubFP28",
                "SubFP1-SubFP8",
                "SubFP86-SubFP151",
                "SubFP1-SubFP5")
top10CxC <- cbind(CxC_labels, CxCtop10)
set.seed(4323)
a <- data.frame(top10compound)
b <- data.frame(top10protein)
c <- data.frame(top10cross_terms)
d <- data.frame(top10PxP)
e <- data.frame(top10CxC)
set.seed(7) # Reordeing Feature usage from highest to lowest
a$compound_labels <- factor(a$compound_labels, levels =a[order(a$Overall), "compound_labels"])
b$Protein <- factor(b$Protein, levels =b[order(b$Overall), "Protein"])
c$crossterms_labels <- factor(c$crossterms_labels, level=c[order(c$Overall), "crossterms_labels"])
d$PxP_labels <- factor(d$PxP_labels, level=d[order(d$Overall), "PxP_labels"])
e$CxC_labels <- factor(e$CxC_labels, level=e[order(e$Overall), "CxC_labels"])
set.seed(9) # Creating the Feature Importance plot
pdf("Fig_Feature_Importance_1.pdf", width = 12, height = 6)
set.seed(10)
z <- ggplot(a, aes(x= Overall, y = compound_labels)) +  theme(axis.title.x=element_text(size=20,face="bold"),
                                                  plot.title = element_text(size=20, face="bold")) +
  geom_point(size=4, colour = "black", fill = "red", pch=21) + coord_fixed(ratio=30) +
  ggtitle("Compound") + xlab("Feature Usage") + ylab("") + theme_bw() + 
  scale_x_continuous(breaks = round(seq(min(0), max(100), by = 25),1), limits= c(0, 100))

y <- ggplot(b, aes(x=Overall, y=Protein)) + theme(axis.title.x=element_text(size=20,face="bold"),
                                              plot.title = element_text(size=20,face="bold")) +
  geom_point(size=4, colour = "black", fill = "blue", pch=21) + coord_fixed(ratio=30) +
  ggtitle("Protein") + xlab("Feature Usage") + ylab("") + theme_bw() + 
  scale_x_continuous(breaks = round(seq(min(0), max(100), by = 25),1), limits= c(0, 100))

x <- ggplot(c, aes(x=Overall, y=crossterms_labels)) + theme(axis.title.x=element_text(size=20,face="bold"),
                                                  plot.title = element_text(size=20,face="bold")) +
  geom_point(size=4, colour = "black", fill="green", pch=21) +coord_fixed(ratio=34) +
  ggtitle("Cross-Terms") + xlab("Feature Usage") + ylab("") + theme_bw() + 
  scale_x_continuous(breaks = round(seq(min(0), max(100), by = 25),1), limits= c(0, 100))

w <- ggplot(d, aes(x=Overall, y=PxP_labels)) + theme(axis.title.x=element_text(size=20, face="bold"),
                                                     plot.title = element_text(size=20,face="bold")) +
  geom_point(size=4, colour = "black", fill= "orange", pch=21) + coord_fixed(ratio=30) +
  ggtitle("Protein-Protein") + xlab("Feature Usage") + ylab("") + theme_bw() +
  scale_x_continuous(breaks = round(seq(min(0), max(100), by = 25), 1), limits=c(0, 100))
  
v <- ggplot(e, aes(x=Overall, y=CxC_labels)) + theme(axis.title.x=element_text(size=20, face="bold"),
                                                     plot.title = element_text(size=20, face="bold")) +
  geom_point(size=4, colour = "black", fill = "black",  pch=21) + coord_fixed(ratio=30) +
  ggtitle("Compound-Compound") + xlab("Feature Usage") + ylab("") + theme_bw() +
  scale_x_continuous(breaks = round(seq(min(0), max(100), by = 25), 1), limits=c(0, 100))


plot.list <- list(z,y,x)
args.list <- c(plot.list,list(ncol=3, nrow=1))
do.call(grid.arrange, args.list)

dev.off()
set.seed(34550)
pdf("Fig_Feature_Importance_2.pdf", width = 12, height = 6)
grid.arrange(w, v, ncol=2)
dev.off()
