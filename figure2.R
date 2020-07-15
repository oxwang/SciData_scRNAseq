

###################################################################
## low quality cell filtering and logCPM normalization by Seurat ##
###################################################################


library(Seurat)
library(Matrix)

input <- '/.../'
output <- '/.../'

load(paste0(input,'/data/fig2_da_counts.rdata'))
load(paste0(input,'/data/fig2_da_cls.rdata'))
gene_anno <- cls[[21]]

batch <- c("10X_LLU_A", "10X_NCI_A", "10X_NCI_M_A", "10X_LLU_B", "10X_NCI_B", "10X_NCI_M_B", "10X_LLU_Mix10", "10X_NCI_Mix5", "10X_NCI_Mix5_F", "10X_NCI_M_Mix5", "10X_NCI_M_Mix5_F", "10X_NCI_M_Mix5_F2","C1_FDA_HT_A", "C1_LLU_A", "ICELL8_SE_A", "ICELL8_PE_A", "C1_FDA_HT_B", "C1_LLU_B", "ICELL8_SE_B", "ICELL8_PE_B")
col_type <- c('darkturquoise','darkorchid','hotpink','firebrick','seagreen4','burlywood','darkorange','deeppink','#C7D591','#DAD9CE',
'#DC8857','#84C1E1','#9D83DC','#B049E2','#D9AACB','#E566A9','#80E570','#6C798A','#E1E058','#76E0C1')


raw_cell <- c()
mito_cell <- c()
process_cell <- c()
cca_process <- counts_norm <- vector('list',length=20)

for(i in 1:length(counts))
{
	celltype <- cls[[i]]
	celltype <- celltype[colnames(counts[[i]]),]
	colnames(counts[[i]]) <- paste(colnames(counts[[i]]),i,sep="_")
	rownames(counts[[i]]) <- gene_anno[rownames(counts[[i]]),2]
	raw_cell <- c(raw_cell,ncol(counts[[i]]))
	meta_da <- cbind.data.frame(batch[i],celltype[,2])
	rownames(meta_da) <- colnames(counts[[i]])
	colnames(meta_da) <- c('batch','celltype')
	tmpCells <- CreateSeuratObject(counts[[i]],min.cells = 3,min.features = 200, meta.data=meta_da)
	mito_features <- grep(pattern = "^MT-", x = rownames(x = tmpCells), value = TRUE)
	mito_perct <- colSums(counts[[i]][mito_features,])/colSums(counts[[i]])
	mito_cell <- c(mito_cell,sum(mito_perct >= 0.1))
	percent_mito <- colSums(x = GetAssayData(object = tmpCells, slot = 'counts')[mito_features, ]) / colSums(x = GetAssayData(object = tmpCells, slot = 'counts'))
	tmpCells[['percent_mito']] <- percent_mito
	Total_mRNAs <- tmpCells[["nCount_RNA"]]$nCount_RNA
	mupper_bound <- 10^(mean(log10(Total_mRNAs)) + 2*sd(log10(Total_mRNAs)))
	mlower_bound <- 10^(mean(log10(Total_mRNAs)) - 2*sd(log10(Total_mRNAs)))
	Total_Genes <- tmpCells[["nFeature_RNA"]]$nFeature_RNA
	gupper_bound <- 10^(mean(log10(Total_Genes)) + 2*sd(log10(Total_Genes)))
	glower_bound <- 10^(mean(log10(Total_Genes)) - 2*sd(log10(Total_Genes)))

	tmpCells <- subset(x = tmpCells, subset = nFeature_RNA > glower_bound & nFeature_RNA < gupper_bound & nCount_RNA > mlower_bound & nCount_RNA < mupper_bound & percent_mito < 0.1)
	tmpCells <- NormalizeData(tmpCells,verbose=F)

## generate raw counts data for bbknn
	counts[[i]] <- GetAssayData(object=tmpCells, slot='counts')

## generate logCPM normalized data for uncorrected
	counts_norm[[i]] <- GetAssayData(object=tmpCells, slot='data')
	process_cell <- c(process_cell,ncol(counts_norm[[i]]))

##	generate cca_process data for seurat integration
	cca_process[[i]] <- tmpCells

## generate merged Cells data for Harmony
	if(i == 1) {Cells <- tmpCells}
	else {Cells <- merge(Cells,tmpCells)}	
}


## generate filtered gene counts matrix for bbknn

for(i in 1:length(counts))
{
	gene_counts <- counts[[i]]
	gene_counts <- cbind.data.frame(rownames(gene_counts),as.matrix(gene_counts))
	colnames(gene_counts) <- c("gene",colnames(gene_counts)[-1])
	write.table(gene_counts,file=paste0(output,'/batch/',batch[i],'.txt'),sep="\t",row.names=F,quote=F)
}


hvg_num <- 2000

######################
## seruat v3 method ##
######################


cca_da <- vector('list',length(cca_process))
for(i in 1:length(cca_da))
{
	cca_da[[i]] <- FindVariableFeatures(cca_process[[i]],nfeatures=hvg_num)
}
cca_anchors <- FindIntegrationAnchors(object.list=cca_da,anchor.features=hvg_num,k.filter=30)
cca_integrated <- IntegrateData(cca_anchors)
cca_integrated <- ScaleData(cca_integrated,verbose=F)
cca_integrated <- RunPCA(cca_integrated, npcs=30, verbose=F)
cca_integrated <- RunTSNE(cca_integrated,reduction='pca',dims=1:30)
cca_integrated <- RunUMAP(cca_integrated,reduction='pca',dims=1:30)
save(cca_integrated,file=paste0(output,'fig2_seurat.Rdata'))


####################
## Hormony method ##
####################


library(harmony)
library(Rtsne)

harmony_data <- Cells
harmony_data <- FindVariableFeatures(harmony_data,nfeatures=hvg_num)
harmony_data <- ScaleData(harmony_data,verbose=FALSE)
harmony_data <- RunPCA(harmony_data,pc.genes=harmony_data@var.genes,npcs=30,verbose=FALSE)
harmony_data <- RunHarmony(harmony_data,"batch")
harmony_data <- RunTSNE(harmony_data,reduction='harmony',dims=1:30)
harmony_data <- RunUMAP(harmony_data,reduction='harmony',dims=1:30)
save(harmony_data,file=paste0(output,'fig2_harmony.Rdata'))


#############
## fastMNN ##
#############

library(SeuratWrappers)
library(batchelor)

mnn_data <- Cells
mnn_data <- FindVariableFeatures(mnn_data,nfeatures=hvg_num)
batch_id <- rep(batch,process_cell)

mnn_data <- RunFastMNN(object.list=SplitObject(mnn_data, split.by = 'batch'),features=hvg_num)
mnn_data <- RunTSNE(mnn_data,reduction='mnn',dims=1:30)
mnn_data <- RunUMAP(mnn_data,reduction='mnn',dims=1:30)

temp_hvgs <- VariableFeatures(mnn_data)
mnn_process <- GetAssayData(mnn_data,'data')
mnn_data <- batchelor::fastMNN(as.matrix(mnn_process), batch=factor(batch_id), subset.row=temp_hvgs)

mnn_corrected <- assay(mnn_data,'reconstructed')
meta_da <- data.frame(batch_id)
rownames(meta_da) <- colnames(mnn_corrected)
colnames(meta_da) <- 'batch'
mnn_data <- CreateSeuratObject(mnn_corrected,min.cells = 0,min.features = 0, meta.data=meta_da)
mnn_data <- ScaleData(mnn_data,verbose=F)
mnn_data <- RunPCA(mnn_data, npcs=30, verbose=F,features=rownames(mnn_corrected))
mnn_embeddings <- Embeddings(mnn_data[['pca']])
mnn_data <- RunTSNE(mnn_data,reduction='pca',dims=1:30)
mnn_data <- RunUMAP(mnn_data,reduction='pca',dims=1:30)
save(mnn_data,file=paste0(output,'fig2_mnn.Rdata'))


###########
## BBKNN ##
###########

library(ggplot2)

temp_output <- '/genomics/1_Projects/FDA_QC/Manuscript/sc_Scientific_Data/Figure/batch/bbknn/'
umap <- read.csv(paste0(temp_output,'umap_coord_bbknn.csv'),row.names=1)
temp_batch <- table(umap[,3])[batch]
colnames(umap) <- c('UMAP_1','UMAP_2','sample')

p7 <- ggplot(umap,aes(x=UMAP_1, y=UMAP_2, color=sample)) + geom_point(size=0.2) + labs(title='BBKNN') + 
scale_color_manual(values=col_type) + theme_classic() + guides(colour = guide_legend(override.aes = list(size=3))) +
theme(plot.title = element_text(size = 9,hjust=0.5, face='bold'), axis.title=element_text(size=9), axis.text=element_text(size=9), legend.position='none',panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))		


#################
## uncorrected ##
#################

common_genes <- rownames(counts_norm[[1]])
for(i in 2:length(counts))
{
	common_genes <- intersect(common_genes,rownames(counts_norm[[i]]))
}

common_counts_norm <- vector('list',length(counts_norm))
for(i in 1:length(counts_norm))
{	
	common_counts_norm[[i]] <- counts_norm[[i]][common_genes,]
}


unc_process <- do.call(cbind,common_counts_norm)
meta_da <- c()

for(i in 1:length(common_counts_norm))
{
	celltype <- cls[[i]]
	rownames(celltype) <- paste0(rownames(celltype),'_',i)
	celltype <- celltype[colnames(common_counts_norm[[i]]),]
	temp <- cbind.data.frame(batch[i],celltype[,2])
	meta_da <- rbind(meta_da,temp)
}

rownames(meta_da) <- colnames(unc_process)
colnames(meta_da) <- c('batch','celltype')
unc_da <- CreateSeuratObject(unc_process,min.cells = 3,min.features = 200, meta.data=meta_da)

unc_out <- FindVariableFeatures(unc_da,nfeatures=hvg_num)
unc_out <- ScaleData(unc_out,verbose=F)
unc_out <- RunPCA(unc_out, npcs=30, verbose=F)
unc_out <- RunTSNE(unc_out,reduction='pca',dims=1:30)
unc_out <- RunUMAP(unc_out,reduction='pca',dims=1:30)
save(unc_out,file=paste0(output,'fig2_unc.Rdata'))




library(ggpubr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

p1 <- DimPlot(harmony_data, reduction = "tsne", group.by = "batch",cols=col_type)
p2 <- DimPlot(harmony_data, reduction = "umap", group.by = "batch",cols=col_type)
p3 <- DimPlot(unc_out, reduction = "tsne", group.by = "batch",cols=col_type)
p4 <- DimPlot(unc_out, reduction = "umap", group.by = "batch",cols=col_type)
p5 <- DimPlot(mnn_data, reduction = "tsne", group.by = "batch",cols=col_type)
p6 <- DimPlot(mnn_data, reduction = "umap", group.by = "batch",cols=col_type)
p8 <- DimPlot(cca_integrated, reduction = "tsne", group.by = "batch",cols=col_type)
p9 <- DimPlot(cca_integrated, reduction = "umap", group.by = "batch",cols=col_type)

p1 <- p1 + labs(title='Harmony') + theme(axis.text = element_text(size = 9), 
plot.title = element_text(size = 9,hjust=0.5), axis.title.x = element_text(size = 9), 
axis.title.y = element_text(size = 9), legend.title = element_text(size = 9), 
legend.text = element_text(size = 9), panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))

p2 <- p2 + labs(title='Harmony') + theme(axis.text = element_text(size = 9), 
plot.title = element_text(size = 9,hjust=0.5), axis.title.x = element_text(size = 9), 
axis.title.y = element_text(size = 9), legend.position = 'none', panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))

p3 <- p3 + labs(title='Uncorrected') + theme(axis.text = element_text(size = 9), 
plot.title = element_text(size = 9,hjust=0.5), axis.title.x = element_text(size = 9), 
axis.title.y = element_text(size = 9), legend.position = 'none', panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))

p4 <- p4 + labs(title='Uncorrected') + theme(axis.text = element_text(size = 9), 
plot.title = element_text(size = 9,hjust=0.5), axis.title.x = element_text(size = 9), 
axis.title.y = element_text(size = 9), legend.position = 'none', panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))

p5 <- p5 + labs(title='fastMNN') + theme(axis.text = element_text(size = 9), 
plot.title = element_text(size = 9,hjust=0.5), axis.title.x = element_text(size = 9), 
axis.title.y = element_text(size = 9), legend.position = 'none', panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))

p6 <- p6 + labs(title='fastMNN') + theme(axis.text = element_text(size = 9), 
plot.title = element_text(size = 9,hjust=0.5), axis.title.x = element_text(size = 9), 
axis.title.y = element_text(size = 9), legend.position = 'none', panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))

p8 <- p8 + labs(title='Seurat V3') + theme(axis.text = element_text(size = 9), 
plot.title = element_text(size = 9,hjust=0.5), axis.title.x = element_text(size = 9), 
axis.title.y = element_text(size = 9), legend.position = 'none', panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))

p9 <- p9 + labs(title='Seurat V3') + theme(axis.text = element_text(size = 9), 
plot.title = element_text(size = 9,hjust=0.5), axis.title.x = element_text(size = 9), 
axis.title.y = element_text(size = 9), legend.position = 'none', panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))

p_legend <- get_legend(p1)
p1 <- p1 + theme(legend.position = 'none')

p10 <- plot_grid(p4,p2,p6,p7,p9, ncol=2, labels=c('a','b','c','d','e'))
p11 <- plot_grid(p10,p_legend,ncol=2, rel_widths=c(4,1))

jpeg(paste0(output,'/figure2.jpeg'),width=10,height=10,res=300,units='in')
plot(p11)
dev.off()


pdf(paste0(output,'/figure2.pdf'),width=10,height=10)
plot(p11)
dev.off()





