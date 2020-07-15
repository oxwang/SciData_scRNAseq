library(Seurat)
library(dplyr)

input <- '/.../'
Cells <- readRDS(paste0(input,'/data/cells_filtered.rds'))

# find HVG for sample A and sample B
Cells_A <- subset(x=Cells, subset = celltype == "HCC1395")
Cells_A <- FindVariableFeatures(Cells_A, selection.method = "vst", nfeatures = 2000)
hvg_A <- Cells_A@assays$RNA@var.features

Cells_B <- subset(x=Cells, subset = celltype == "HCC1395BL")
Cells_B <- FindVariableFeatures(Cells_B, selection.method = "vst", nfeatures = 2000)
hvg_B <- Cells_B@assays$RNA@var.features

# calculate mean expression of HVG for each data sets
library(dplyr)
mean_hvg_A <- matrix(nrow = length(hvg_A), ncol = length(unique(Cells_A$sample)),
                     dimnames = list(hvg_A, unique(Cells_A$sample)))
var_hvg_A <- matrix(nrow = length(hvg_A), ncol = length(unique(Cells_A$sample)),
                     dimnames = list(hvg_A, unique(Cells_A$sample)))

for (i in 1:length(unique(Cells_A$sample))) {
  t <- unique(Cells_A$sample)[i]
  mean_hvg_A[,i] <- Matrix::rowMeans(Cells_A@assays$RNA@data[hvg_A,which(Cells_A$sample == t)], na.rm = T)
  var_hvg_A[,i] <- matrixStats::rowVars(as.matrix(Cells_A@assays$RNA@data[hvg_A,which(Cells_A$sample == t)]), na.rm = T)
}

mean_hvg_B <- matrix(nrow = length(hvg_B), ncol = length(unique(Cells_B$sample)),
                     dimnames = list(hvg_B, unique(Cells_B$sample)))
var_hvg_B <- matrix(nrow = length(hvg_B), ncol = length(unique(Cells_B$sample)),
                     dimnames = list(hvg_B, unique(Cells_B$sample)))

for (i in 1:length(unique(Cells_B$sample))) {
  t <- unique(Cells_B$sample)[i]
  mean_hvg_B[,i] <- Matrix::rowMeans(Cells_B@assays$RNA@data[hvg_B,which(Cells_B$sample == t)], na.rm = T)
  var_hvg_B[,i] <- matrixStats::rowVars(as.matrix(Cells_B@assays$RNA@data[hvg_B,which(Cells_B$sample == t)]), na.rm = T)
}

## correlation ----
cor.A.mean <- cor(mean_hvg_A)
cor.A.var <- cor(var_hvg_A)
cor.B.mean <- cor(mean_hvg_B)
cor.B.var <- cor(var_hvg_B)

par(ps=8, lwd=0.75)
par(mfrow=c(2,2))
corrplot::corrplot(cor.A.mean, method = "number", order="hclust", rect.lwd = 1,
         addrect=2, number.cex = 0.8, tl.cex = 1, cl.cex = 0.9, cl.lim = c(0, 1),
         mar = c(0, 0, 0, 0), tl.col="black")
corrplot::corrplot(cor.A.var, method = "number", order="hclust", rect.lwd = 1,
         addrect=2, number.cex = 0.8, tl.cex = 1, cl.cex = 0.9, cl.lim = c(0, 1),
         mar = c(0, 0, 0, 0), tl.col="black")
corrplot::corrplot(cor.B.mean, method = "number", order="hclust", rect.lwd = 1,
                   addrect=2, number.cex = 0.8, tl.cex = 1, cl.cex = 0.9, cl.lim = c(0, 1),
                   mar = c(0, 0, 0, 0), tl.col="black")
corrplot::corrplot(cor.B.var, method = "number", order="hclust", rect.lwd = 1,
                   addrect=2, number.cex = 0.8, tl.cex = 1, cl.cex = 0.9, cl.lim = c(0, 1),
                   mar = c(0, 0, 0, 0), tl.col="black")

# cell percentage in each data sets
df1 <- Cells_A[[]] %>%
  as_tibble() %>%
  select(sample, color) %>%
  group_by(sample, color) %>%
  summarise(Freq = n()) %>%
  mutate(cells = "cells")
df2 <- Cells_B[[]] %>%
  as_tibble() %>%
  select(sample, color) %>%
  group_by(sample, color) %>%
  summarise(Freq = n()) %>%
  mutate(cells = "cells")

library(scales)
ggplot(df1, aes(cells, Freq, fill=sample)) +
  geom_bar(position="fill", stat="identity", width = 100) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = df1$color) +
  ylab("Percentage of cells") + labs(title = paste0(sum(df1$Freq), " cells")) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = "black", size = 9),
        axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 9),
        legend.key.size = unit(0.4, "cm"),legend.spacing.y = unit(0.4, 'cm'),
        legend.justification = c(0,1),
        legend.text = element_text(size = 9), legend.title = element_blank())

ggplot(df2, aes(cells, Freq, fill=sample)) +
  geom_bar(position="fill", stat="identity", width = 100) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = df2$color) +
  ylab("Percentage of cells") + labs(title = paste0(sum(df2$Freq), " cells")) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = "black", size = 9),
        axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 9),
        legend.key.size = unit(0.4, "cm"),legend.spacing.y = unit(0.4, 'cm'),
        legend.justification = c(0,1),
        legend.text = element_text(size = 9), legend.title = element_blank())
