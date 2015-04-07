#!/usr/bin/env Rscript
C <- read.csv("data/Compound.csv", header=TRUE)
P <- read.csv("data/Protein.csv", header=TRUE)
Y <- read.csv("data/Activity.csv", header=TRUE)
SMOTE_Y <- read.csv("data/Balanced_Activity.csv", header=TRUE)
library(unbalanced)
library(corrplot)
library(caret)
library(prospectr)
library(gridExtra)
library(Rcpi)
# Balancing Data with SMOTE Algorithm
set.seed(1)

# This is where you should use for loop and please create a vector to
# store results instead of a,b,c,d...  Writing code this way is
# error-prone and very hard to debug because the code is redundant.

# At first glance, the code writing style is totally inconsistent.
# Please consult
# https://google-styleguide.googlecode.com/svn/trunk/Rguide.xml#filenames
# and try to do your best to follow the guideline.

SMOTE <- ubBalance(X=cbind(C,P), Y=factor(Y$Before_SMOTE), 
                   type= "ubSMOTE", perc.over = 100,
                   perc.under=200, verbose=TRUE)
a <- cbind(SMOTE$Y, SMOTE$X)
a <- a[order(SMOTE$Y),]
a <- a[-1]
set.seed(2)
SMOTE <- ubBalance(X=cbind(C,P), Y=factor(Y$Before_SMOTE), 
                   type= "ubSMOTE", perc.over = 100,
                   perc.under=200, verbose=TRUE)
b <- cbind(SMOTE$Y, SMOTE$X, row.names=NULL)
b <- b[order(SMOTE$Y),]
b <- b[-1]
set.seed(3)
SMOTE <- ubBalance(X=cbind(C,P), Y=factor(Y$Before_SMOTE), 
                   type= "ubSMOTE", perc.over = 100,
                   perc.under=200, verbose=TRUE)
c <- cbind(SMOTE$Y, SMOTE$X, row.names=NULL)
c <- c[order(SMOTE$Y),]
c <- c[-1]
set.seed(4)
SMOTE <- ubBalance(X=cbind(C,P), Y=factor(Y$Before_SMOTE), 
                   type= "ubSMOTE", perc.over = 100,
                   perc.under=200, verbose=TRUE)
d <- cbind(SMOTE$Y, SMOTE$X, row.names=NULL)
d <- d[order(SMOTE$Y),]
d <- d[-1]
set.seed(5)
SMOTE <- ubBalance(X=cbind(C,P), Y=factor(Y$Before_SMOTE), 
                   type= "ubSMOTE", perc.over = 100,
                   perc.under=200, verbose=TRUE)
e <- cbind(SMOTE$Y, SMOTE$X, row.names=NULL)
e <- e[order(SMOTE$Y),]
e <- e[-1]
set.seed(6)
SMOTE <- ubBalance(X=cbind(C,P), Y=factor(Y$Before_SMOTE), 
                   type= "ubSMOTE", perc.over = 100,
                   perc.under=200, verbose=TRUE)
f <- cbind(SMOTE$Y, SMOTE$X, row.names=NULL)
f <- f[order(SMOTE$Y),]
f <- f[-1]
set.seed(7)
SMOTE <- ubBalance(X=cbind(C,P), Y=factor(Y$Before_SMOTE), 
                   type= "ubSMOTE", perc.over = 100,
                   perc.under=200, verbose=TRUE)
g <- cbind(SMOTE$Y, SMOTE$X, row.names=NULL)
g <- g[order(SMOTE$Y),]
g <- g[-1]
set.seed(8)
SMOTE <- ubBalance(X=cbind(C,P), Y=factor(Y$Before_SMOTE), 
                   type= "ubSMOTE", perc.over = 100,
                   perc.under=200, verbose=TRUE)
h <- cbind(SMOTE$Y, SMOTE$X, row.names=NULL)
h <- h[order(SMOTE$Y),]
h <- h[-1]
set.seed(9)
SMOTE <- ubBalance(X=cbind(C,P), Y=factor(Y$Before_SMOTE), 
                   type= "ubSMOTE", perc.over = 100,
                   perc.under=200, verbose=TRUE)
i <- cbind(SMOTE$Y, SMOTE$X, row.names=NULL)
i <- i[order(SMOTE$Y),]
i <- i[-1]
set.seed(10)
SMOTE <- ubBalance(X=cbind(C,P), Y=factor(Y$Before_SMOTE), 
                   type= "ubSMOTE", perc.over = 100,
                   perc.under=200, verbose=TRUE)
j <- cbind(SMOTE$Y, SMOTE$X, row.names=NULL)
j <- j[order(SMOTE$Y),]
j <- j[-1]
set.seed(100)
Total <- a + b + c + d + e + f + g + h + i + j
Average_SMOTE <- Total/10
## Removeing Interacorrelation with the cutoff of 0.7
set.seed(200)
compoundSMOTE <- Average_SMOTE[, 1:44]
proteinSMOTE <- Average_SMOTE[, 45:443]
set.seed(400)
compound_corr <- cor(compoundSMOTE)
dim(compound_corr)
compound_corr <- compound_corr[1:44, 1:44]
compound_highCorr <- findCorrelation(compound_corr, cutoff = .70)
length(compound_highCorr)
highCorrRemoveCompound <- compoundSMOTE[, -compound_highCorr]
Balanced_Activity <- SMOTE_Y[,1]
InputKennardCompound <- cbind(Balanced_Activity, highCorrRemoveCompound)
set.seed(400)
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
Inactive <- C_Inactive
x <- data.frame(Inactive)
sel <- kenStone(x[-1], k = 70, metric = "mahal", pc=2)
trainInactiveCompound <- x[sel$model, ]
testInactiveCompound <- x[sel$test, ]
Active <- C_Active
x <- data.frame(Active)
sel <- kenStone(x[-1], k = 70, metric = "mahal", pc=2)
trainActiveCompound <- x[sel$model, ]
testActiveCompound <- x[sel$test, ]
Inactive <- P_Inactive
x <- data.frame(Inactive)
sel <- kenStone(x[-1], k = 70, metric = "mahal", pc=2)
trainInactiveProtein <- x[sel$model, ]
testInactiveProtein <- x[sel$test, ]
Active <- P_Active
x <- data.frame(Active)
sel <- kenStone(x[-1], k = 70, metric = "mahal", pc=2)
trainActiveProtein <- x[sel$model, ]
testActiveProtein <- x[sel$test, ]
# Preparing input data for Proteocheometrics modeling
set.seed(567)
C_Train <- rbind(trainInactiveCompound, trainActiveCompound)
C_Test <- rbind(testInactiveCompound, testActiveCompound)
P_Train <- rbind(trainInactiveProtein, trainActiveProtein)
P_Test <- rbind(testInactiveProtein, testActiveProtein)
set.seed(22)
c <- C_Train
p <- P_Train
train_c <- c[-1]
train_p <- p[-1]
crossTrain <- getCPI(train_c, train_p, type = "tensorprod")
c <- C_Test
p <- P_Test
test_c <- c[-1]
test_p <- p[-1]
crossTest <- getCPI(test_c, test_p, type = "tensorprod")
Balanced_Activity_Train <- C_Train[1]
Balanced_Activity_Test <- C_Test[1]
CxP_Train <- cbind(Balanced_Activity_Train, data.frame(crossTrain))
CxP_Test <- cbind(Balanced_Activity_Test, data.frame(crossTest))
## Feature Importance
input <- list(C_Train = C_Train,
              P_Train = P_Train,
              CxP_Train = CxP_Train)
# C5.0 Algorithm is used to obtain feature importance of BT chemicals and proteins
set.seed(34440)
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
set.seed(1)
compoundtop10 <- head(compound,10) # top 10 compound
proteintop10 <- head(protein,10) # top 10 compound
CxPtop10 <- head(CxP,10) # top 10 Cross-terms
set.seed(2)
myDF1 <- cbind(Compound = rownames(compoundtop10, compoundtop10))
top10compound <- cbind(myDF1, compoundtop10)
set.seed(3)
myDF2 <- cbind(Protein = rownames(proteintop10, proteintop10))
top10protein <- cbind(myDF2, proteintop10)
crossterms_labels <- c("SubFPC96-X5z3scl3.lag19",
                       "SubFPC12-z3scl2.lag2",
                       "SubFPC28-X5z3scl3.lag3",
                       "SubFPC28-X3z3scl3.lag3",
                       "SubFPC1-z3scl1.lag15",
                       "SubFPC303-X5z3scl1.lag21",
                       "SubFPC2-z3scl1.lag13",
                       "SubFPC2-z3scl1.lag19",
                       "SubFPC2-z3scl1.lag20",
                       "SubFPC2-z3scl1.lag23")
top10cross_terms <- cbind(crossterms_labels, CxPtop10)
set.seed(4323)

set.seed(6)
a <- data.frame(top10compound)
b <- data.frame(top10protein)
c <- data.frame(top10cross_terms)
set.seed(7) # Reordeing Feature usage from highest to lowest
a$Compound <- factor(a$Compound, levels =a[order(a$Overall), "Compound"])
b$Protein <- factor(b$Protein, levels =b[order(b$Overall), "Protein"])
c$crossterms_labels <- factor(c$crossterms_labels, level=c[order(c$Overall), "crossterms_labels"])
set.seed(9) # Creating the Feature Importance plot
pdf("Feature Importance for BT.pdf", width = 12, height = 6)
set.seed(10)
x <- ggplot(a, aes(x= Overall, y = Compound)) +  theme(axis.title.x=element_text(size=20,face="bold"),
                                                  plot.title = element_text(size=20, face="bold")) +
  geom_point(size=4, colour = "red") + coord_fixed(ratio=30) +
  ggtitle("Compound") + xlab("Feature Usage") + ylab("") + theme_bw() + 
  scale_x_continuous(breaks = round(seq(min(0), max(100), by = 25),1), limits= c(0, 100))

y <- ggplot(b, aes(x=Overall, y=Protein)) + theme(axis.title.x=element_text(size=20,face="bold"),
                                              plot.title = element_text(size=20,face="bold")) +
  geom_point(size=4, colour = "blue") + coord_fixed(ratio=30) +
  ggtitle("Protein") + xlab("Feature Usage") + ylab("") + theme_bw() + 
  scale_x_continuous(breaks = round(seq(min(0), max(100), by = 25),1), limits= c(0, 100))

z <- ggplot(c, aes(x=Overall, y=crossterms_labels)) + theme(axis.title.x=element_text(size=20,face="bold"),
                                                  plot.title = element_text(size=20,face="bold")) +
  geom_point(size=4, colour = "green") +coord_fixed(ratio=34) +
  ggtitle("Cross-Terms") + xlab("Feature Usage") + ylab("") + theme_bw() + 
  scale_x_continuous(breaks = round(seq(min(0), max(100), by = 25),1), limits= c(0, 100))

grid.arrange(x, y, z,  ncol = 3)
dev.off()

