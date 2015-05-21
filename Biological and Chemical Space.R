#!/usr/bin/env Rscript
pdf("Fig_PCA.pdf", width = 20, height = 10)
Biological_Space <- read.csv("data/Biological_Space.csv", header=TRUE)
df <- data.frame(Biological_Space[,2:400])
df = df[, -nearZeroVar(df)]
set.seed(44321)
pca <- prcomp(df, retx=TRUE,scale.=TRUE)
scores <- pca$x[,1:3]
km <- kmeans(scores, center=3, nstart=5)
Protein <- Biological_Space[,1]
ProteinNumber <- c("10", "11", "12", "13", "14", "15", "16", "9", "8", "7", "6", "5", "4", "3", "2", "1")
ggdata <- data.frame(scores, Cluster=km$cluster)
ggdata <- cbind(ProteinNumber, ggdata)
set.seed(23)
x <- ggplot(ggdata) +
  geom_point(aes(x=PC1, y=PC2,
                 color=factor(Cluster)), size=5, shape=20) +
  ggtitle("A") +
  stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Cluster)),
               geom="polygon", level=0.95, alpha=0.2) +
  guides(color=guide_legend("Cluster"),fill=guide_legend("Cluster")) +
  geom_text(aes(x=PC1, y=PC2, label=ProteinNumber), size=8, hjust=0, vjust=-0.5, alpha=0.45) +
  theme(
    legend.position=("none"),
    plot.title = element_text(size=20, face="bold", colour="black", vjust = 2, hjust=-0.07),
    axis.text.y = element_text(size = 15),
    axis.ticks.length = unit(0.3, "cm"),
    axis.text.x = element_text(size = 15),
    legend.title=element_blank(),
    axis.title.x = element_text(color="black", size=20),
    axis.title.y = element_text(color="black", size=20)) + 
  scale_y_continuous(limits = c(-15,20), expand = c(0,0)) +
  scale_x_continuous(limits = c(-15,20), expand = c(0,0)) 
## Compounds
set.seed(42342)
Chemical_Space <- read.csv("data/Chemical_Space.csv", header=TRUE)
df <- data.frame(Chemical_Space[,2:45])
pca <- prcomp(df, retx=TRUE,scale.=TRUE)
set.seed(34233)
scores <- pca$x[,1:3]
km <- kmeans(scores, center=2, nstart=5)
Compound <- c("1", "2", "3", "4", "5", "6", "7",
              "8", "9", "10", "11", "12", "13",
              "14", "15", "16", "17", "18")
ggdata <- data.frame(scores, Cluster=km$cluster)
ggdata <- cbind(Compound, ggdata)
set.seed(10093)
y <- ggplot(ggdata) +
  geom_point(aes(x=PC1, y=PC2,
                 color=factor(Cluster)), size=5, shape=20) +
  ggtitle("B") +
  stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Cluster)),
               geom="polygon", level=0.95, alpha=0.2) +
  guides(color=guide_legend("Cluster"),fill=guide_legend("Cluster")) +
  geom_text(aes(x=PC1, y=PC2, label=Compound), size=8, hjust=0, vjust=-0.5, alpha=0.45) +
  theme(
    legend.position=("none"),
    plot.title = element_text(size=20, face="bold", colour="black", vjust = 2, hjust=-0.07),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.ticks.length = unit(0.3, "cm"),
    legend.title=element_blank(),
    axis.title.x = element_text(color="black", size=20),
    axis.title.y = element_blank()) + 
  scale_y_continuous(limits = c(-15,20), expand = c(0,0)) +
  scale_x_continuous(limits = c(-15,20), expand = c(0,0)) 

grid.arrange(x, y, nrow=1)
dev.off()