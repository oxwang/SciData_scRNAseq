
library(ggpubr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

input <- '/.../'
output <- '/.../'

load(paste0(input,'/data/fig4_data.rdata'))
cell_num <- data[[1]]
gene_exprs_num_10x <- data[[2]]
mix_cell_perct <- data[[3]]
gene_exprs_num <- data[[4]]


p1 <- ggplot(cell_num, aes(x=Platform,y=Cell_num,fill=Pipeline)) + geom_bar(stat="identity",position=position_dodge()) + 
labs(title="Number of selected cells",x="", y="") + scale_y_continuous(limits= c(0, 7000), breaks = seq(0, 7000, by = 2000)) +
theme(plot.title = element_text(hjust=0.5, size = 9), legend.position=c(0.84,0.83), legend.title = element_blank(), 
legend.text=element_text(size = 9), axis.text.x = element_text(angle=30, size = 9,hjust=1), axis.text.y = element_text(size = 9), 
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))


p2 <- ggplot(gene_exprs_num_10x, aes(x=Platform,y=Expressed_gene,fill=Pipeline)) + geom_boxplot(position=position_dodge(),outlier.size=0.2) + 
labs(title="Number of expressed genes per cell",x="", y="") + scale_y_continuous(limits=c(0,12000), breaks = seq(0, 12000, by = 2000)) +
theme(plot.title = element_text(hjust=0.5, size = 9), legend.position=c(0.84,0.83), legend.title = element_blank(), 
legend.text=element_text(size = 9), axis.text.x = element_text(angle=30, size = 9,hjust=1), axis.text.y = element_text(size = 9), 
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))


p3 <- ggplot(mix_cell_perct, aes(x=Platform,y=Perct,fill=Pipeline)) + geom_bar(stat="identity",position=position_dodge()) +
labs(title="Perctentage of HCC1395 cells",x="", y="") + scale_y_continuous(limits=c(0,0.12), breaks = seq(0, 0.12, by = 0.02)) +
theme(plot.title = element_text(hjust=0.5, size = 9), legend.position=c(0.84,0.83), legend.title = element_blank(), 
legend.text=element_text(size = 9), axis.text.x = element_text(angle=30, size = 9,hjust=1), axis.text.y = element_text(size = 9), 
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))


p4 <- ggplot(gene_exprs_num, aes(x=Platform,y=Expressed_gene,fill=Pipeline)) + geom_boxplot(position=position_dodge(),outlier.size=0.2) + 
labs(title="Number of expressed genes per cell",x="", y="") + scale_y_continuous(limits=c(0,14000), breaks = seq(0, 14000, by = 2000)) +
theme(plot.title = element_text(hjust=0.5, size = 9), legend.position=c(0.85,0.87), legend.title = element_blank(), 
legend.text=element_text(size = 9), axis.text.x = element_text(angle=30, size = 9,hjust=1), axis.text.y = element_text(size = 9), 
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))


jpeg(paste0(output,'figure4.jpeg'),res=300,width=10,height=10,unit='in')
plot_grid(p1,p2,p3,p4,labels=c('a','b','c','d'))
dev.off()


pdf(paste0(output,'figure4.pdf'),width=10,height=10)
plot_grid(p1,p2,p3,p4,labels=c('a','b','c','d'))
dev.off()














