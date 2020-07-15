library(Seurat)
library(infercnv)
library(tidyverse)
library(harmony)
library(ggdendro)

input <- '/.../'
cells.filtered <- readRDS(paste0(input,'/data/cells_filtered.rds'))

#############################################
# inferCNV analysis
#############################################
# downsample cells in each datasets
table(cells.filtered$sample)
large <- subset(cells.filtered, subset = sample %in% unique(cells.filtered$sample)[c(1:4,6:8,10:12)])
df <- large[[]] %>% rownames_to_column('id') %>% group_by(sample) %>% sample_n(1000) %>% column_to_rownames('id')
large2 <- subset(cells.filtered, cells = rownames(df))
small <- subset(cells.filtered, subset = sample %in% unique(cells.filtered$sample)[c(13,15,17,19)])
cells.subsample <- merge(large2, small)

counts_matrix = as.matrix(cells.subsample@assays$RNA@counts)
write.table(counts_matrix, file='data/fig3_all.matrix', quote=F, sep="\t")
metadata <- as.data.frame(cells.subsample[[]][, 4, drop=FALSE])
write.table(metadata, file='data/fig3_all.matrix.Annotations.txt', quote=F, sep="\t", col.names = FALSE)

# run infercnv
library(infercnv)
# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="./data/fig3_all.matrix",
                                    annotations_file="./data/fig3_all.matrix.Annotations.txt",
                                    delim="\t",
                                    gene_order_file="./data/fig3_ensemble_id_gen_pos.txt",
                                    ref_group_names=c("10X_LLU_B")
)

output_dir <- "result"
infercnv_obj = infercnv::run(infercnv_obj,
                             analysis_mode='subclusters',
                             hclust_method='ward.D2',
                             tumor_subcluster_pval=0.05,
                             tumor_subcluster_partition_method='qnorm',
                             cutoff=0.1,
                             out_dir=output_dir,
                             plot_steps=T,
                             denoise=T,
                             HMM=T,
                             png_res=300,
                             num_threads=100
)

# Extracting HMM features - map_metadata_from_infercnv.txt
seurat_obj = add_to_seurat(infercnv_output_path=output_dir,
                           seurat_obj=cells.subsample,
                           top_n=10
)

# Extracting the dendrogram plot data
dend <- read.newick("result/infercnv.observations_dendrogram.txt") %>%
        as.hclust.phylo()
ddata <- dendro_data(dend, type = "rectangle")
ggplot(ddata$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), size=0.5/.pt) +
  coord_flip() +
  scale_y_reverse(expand = c(0.2, 0)) +
  theme_dendro()

#############################################
# Seurat analysis
#############################################
cells <- cells.subsample %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = pbmc@var.genes, npcs = 20, verbose = FALSE)
cells <- cells %>% RunHarmony(group.by.vars = "sample")
cells <- RunUMAP(cells, reduction = "harmony", dims = 1:30)
cells <- cells %>% FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution=0.001)
cells$clusters <- ifelse(cells$seurat_clusters == 0, "HCC1395BL", "HCC1395")

DimPlot(object = cells, reduction = "harmony", pt.size = .1, group.by = "clusters",
              do.return = TRUE) +
  scale_color_manual(values = c("#85144b","#0074D9")) +
  border() +
  labs(x="UMAP_1", y="UMAP_2") +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 9, colour = "black"),
        legend.position = "bottom")
