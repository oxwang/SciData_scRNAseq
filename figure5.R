library(tidyverse)
library(ggrepel)
library(scales)

input <- '/.../'
Cells <- readRDS(paste0(input,'/data/cells_filtered.rds'))
Cells.unfiltered <- readRDS(paste0(input,'/data/cells_unfiltered.rds'))

sample_order <- unique(Cells$sample)

keeped_cells <- ifelse(colnames(Cells.unfiltered) %in% colnames(Cells),"high_quality","low_quality")
Cells.unfiltered$keeped_cells <- keeped_cells

# Plot (number of genes)
df.p4 <- Cells.unfiltered[[]]
df.p4$sample <- factor(df.p4$sample, levels = unique(df.p4$sample))
df.p4 <- df.p4 %>%
  group_by(sample, keeped_cells, color) %>%
  summarise(subtotal=mean(nFeature_RNA)) %>%
  spread(keeped_cells, subtotal)

p1 <- ggplot(df.p4, aes(high_quality, low_quality, color=sample)) +
  geom_vline(xintercept = 4500, color = "grey", linetype=2) +
  geom_hline(yintercept = 2500, color = "grey", linetype=2) +
  geom_label_repel(aes(label=ifelse(high_quality>4500|low_quality>2500, as.character(sample),'')), size=9/.pt) +
  geom_point()+
  labs(title = "Mean number of genes") + xlab("High quality cells") + ylab("Low quality cells (filtered)") +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = df.p4$color) +
  theme_bw(base_size = 9) +
  guides(color=guide_legend(ncol=2)) +
  theme(plot.title = element_text(hjust = 0.5, size = 9,face = "bold"),
        axis.text = element_text(size = 9),legend.title = element_blank(),
        panel.grid = element_blank())

# Plot (read depth)
df.p3 <- Cells.unfiltered[[]]
df.p3$sample <- factor(df.p3$sample, levels = unique(df.p3$sample))
df.p3 <- df.p3 %>% group_by(sample, keeped_cells, color) %>% summarise(subtotal=mean(nCount_RNA)) %>%
  spread(keeped_cells, subtotal)

p2 <- ggplot(df.p3, aes(high_quality, low_quality, color=sample)) +
  geom_vline(xintercept = 1e+05, color = "grey", linetype=2) +
  geom_hline(yintercept = 1e+05, color = "grey", linetype=2) +
  geom_label_repel(aes(label=ifelse(high_quality>1e+05 , as.character(sample),'')), size=9/.pt) +
  geom_point() +
  scale_colour_manual(values = df.p3$color) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Mean number of read counts") + xlab("High quality cells") + ylab("Low quality cells (filtered)") +
  theme_bw(base_size = 9) +
  guides(color=guide_legend(ncol=2)) +
  theme(plot.title = element_text(hjust = 0.5, size = 9,face = "bold"),
        axis.text = element_text(size = 9), legend.title = element_blank(),
        panel.grid = element_blank())



# Plot (percentage of MT)
df.p5 <- Cells.unfiltered[[]]
df.p5$sample <- factor(df.p5$sample, levels = unique(df.p5$sample))
df.p5 <- df.p5 %>% group_by(sample, keeped_cells, color) %>% summarise(subtotal=mean(percent.mt)) %>%
  spread(keeped_cells, subtotal)

p3 <- ggplot(df.p5, aes(high_quality, low_quality, color=sample)) +
  geom_hline(yintercept = 10, color = "grey", linetype=2) +
  geom_vline(xintercept = 5, color = "grey", linetype=2) +
  geom_label_repel(aes(label=ifelse(high_quality<5 & low_quality<10, as.character(sample),'')), size=9/.pt) +
  geom_point()+
  scale_colour_manual(values = df.p5$color) +
  labs(title = "Mean percentage of mitocondrial genes") + xlab("High quality cells") + ylab("Low quality cells (filtered)") +
  theme_bw(base_size = 9) +
  guides(color=guide_legend(ncol=2)) +
  theme(plot.title = element_text(hjust = 0.5, size = 9,face = "bold"),
        axis.text = element_text(size = 9),legend.title = element_blank(),
        panel.grid = element_blank())

# plot cell percentage
df1 <- as.data.frame(table(Cells.unfiltered$keeped_cells))
df1$cells <- "Cells"
df1$Var1 <- factor(df1$Var1, levels = rev(levels(df1$Var1)))

p4 <- ggplot(df1, aes(cells, Freq, fill=Var1)) +
  geom_bar(position="fill", stat="identity", width = 100) +
  geom_text(label=df1$Freq[1], y=df1$Freq[1]/(sum(df1$Freq)*2), color="white", size=9/.pt) +
  geom_text(label=df1$Freq[2], y=(1-df1$Freq[2]/(sum(df1$Freq)*2)), color="white", size=9/.pt) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  scale_fill_manual(name = '', labels = c('low quality', "high quality"), values = c('darkred',"darkgreen")) +
  ylab("Percentage of cells") + labs(title = paste0(sum(df1$Freq), " cells")) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = "black", size = 9),
        legend.position = "bottom", legend.key.size = unit(0.4, "cm"),legend.spacing.y = unit(0.4, 'cm'),
        axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 9,face = "bold"),
        legend.text = element_text(size = 9)) +
  guides(fill=guide_legend(nrow=2))


