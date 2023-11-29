.libPaths("/home/zhouxin/R/x86_64-pc-linux-gnu-library/4.2")


library(Seurat)
library(ggplot2)

seu_POD <- Read10X("~/liumin/POD/outs/filtered_feature_bc_matrix")
seu_CARM <- Read10X("~/liumin/CAR-Ms/outs/filtered_feature_bc_matrix")
seu_CARTM <- Read10X("~/liumin/CAR-Trem2-Ms/outs/filtered_feature_bc_matrix")

POD <- CreateSeuratObject(seu_POD, project = "POD", min.cells = 3, min.features = 200)
CARM <- CreateSeuratObject(seu_CARM, project = "CAR-Ms", min.cells = 3, min.features = 200)
CARTM <- CreateSeuratObject(seu_CARTM, project = "CAR-Trem2-Ms", min.cells = 3, min.features = 200)


seu <- merge(POD, c(CARM, CARTM), add.cell.ids = c("POD", "CARM", "CARTM"))

table(seu$orig.ident)

seu$orig.ident <- factor(seu$orig.ident, levels = c("POD", "CAR-Ms", "CAR-Trem2-Ms"))


seu[["percent.mt"]] <- PercentageFeatureSet(object = seu, pattern = "^mt-")

pdf(file = "1.pdf", width = 10, height = 5)
VlnPlot(object = seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()

seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

table(seu$orig.ident)

rownames(seu)

seu$percent.mt <- NULL

colnames(seu)



DoubletFinder_multisamples <- function(seu_list){
  library(tidyverse)
  library(Seurat)
  library(DoubletFinder)
  library(conflicted)
  conflict_prefer("select", "dplyr")
  for (idx in 1:length(seu_list)) {
    data <- NormalizeData(seu_list[[idx]])
    data <- ScaleData(data, verbose = FALSE)
    data <- FindVariableFeatures(data, verbose = FALSE)
    data <- RunPCA(data, npcs = 40, verbose = FALSE)
    data <- RunUMAP(data, reduction = "pca", dims = 1:30)
    data <- FindNeighbors(data, reduction = "pca", dims = 1:30)
    data <- FindClusters(data, resolution = 0.5)
    sweep.res.list <- paramSweep_v3(data, PCs = 1:30, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    ## Homotypic Doublet Proportion Estimate 
    annotations <- data@meta.data$ClusteringResults
    homotypic.prop <- modelHomotypic(annotations)   
    #10X平台按每增加1000个细胞，双细胞比率增加千分之8来计算。
    nExp_poi <- round((ncol(data)*8*1e-6)*nrow(data@meta.data)) 
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    ## Run DoubletFinder with varying classification stringencies 
    data <- doubletFinder_v3(data, PCs = 1:30, pN = 0.25, pK = 0.09, 
                             nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    ## save results
    seu_list[[idx]]$doubFind_res = data@meta.data %>% select(contains('DF.classifications'))
    seu_list[[idx]]$doubFind_score = data@meta.data %>% select(contains('pANN'))
  }
  seu_list
}

table(seu$orig.ident)

seu_list <- SplitObject(seu, split.by = "orig.ident")
seu_list <- DoubletFinder_multisamples(seu_list = seu_list)
seu <- merge(seu_list$POD, c(seu_list$`CAR-Ms`, seu_list$`CAR-Trem2-Ms`))


table(seu$orig.ident, seu$doubFind_res)

seu <- subset(seu, doubFind_res == "Singlet")

table(seu$orig.ident)

seu$orig.ident <- factor(seu$orig.ident, levels = c("POD", "CAR-Ms", "CAR-Trem2-Ms"))


saveRDS(seu, file = "seu.rds")

##########################################################################################

seu$Group <- seu$orig.ident
seu$RNA_snn_res.0.05 <- NULL
seu$RNA_snn_res.0.1 <- NULL
seu$RNA_snn_res.0.3 <- NULL
seu$RNA_snn_res.0.5 <- NULL
seu$seurat_clusters <- NULL
seu@reductions$pca <- NULL
seu@reductions$umap <- NULL
seu@reductions$harmony <- NULL


seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)

#library(harmony)
#seu <- RunHarmony(seu, group.by.vars = "Group")


seu <- RunUMAP(seu, dims = 1:30, reduction = "pca")
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.05)

DimPlot(seu, split.by = "orig.ident", label = T, reduction = "umap")

degs <- FindAllMarkers(seu, logfc.threshold = 0.25, min.pct = 0.5, only.pos = T)

write_csv(degs, file = "degs.csv")


marker <- c("Col1a1", "Dcn", "Itgbl1", "Lyz2", "Adgre1", "Cd68", "Cd3d", "Cd28", "Il2rb", "Krt14", "Krt5", "Acta2", "Pecam1", "Cdh5", "Lyve1", "S100a9", "S100a8")

DotPlot(seu, features = c("Adgre1", "Cd68", "Trem2"), assay = "RNA", group.by = "seurat_clusters", cols = c("gray90", "red")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "italic"))

DotPlot(seu, features = c("Col1a1", "Dcn", "Col3a1"), assay = "RNA", group.by = "seurat_clusters", cols = c("gray90", "red")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "italic"))

DotPlot(seu, features = c("Cd3d", "Il2rb"), assay = "RNA", group.by = "seurat_clusters", cols = c("gray90", "red")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "italic"))

DotPlot(seu, features = c("Krt14", "Krt5"), assay = "RNA", group.by = "seurat_clusters", cols = c("gray90", "red")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "italic"))

DotPlot(seu, features = c("Acta2"), assay = "RNA", group.by = "seurat_clusters", cols = c("gray90", "red")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "italic"))



DotPlot(seu, features = c("Pecam1", "Cdh5", "Lyve1"), assay = "RNA", group.by = "seurat_clusters", cols = c("gray90", "red")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "italic"))

DotPlot(seu, features = c("Ccl3", "Hdc"), assay = "RNA", group.by = "seurat_clusters", cols = c("gray90", "red")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "italic"))


FeaturePlot(seu, features = "Dpp4", split.by = "orig.ident")

VlnPlot(seu, features = "Cd3e", pt.size = 0)

VlnPlot(seu, features = "Dpp4", pt.size = 0.1, split.by = "orig.ident", group.by = "seurat_clusters")

cols_clusters <- RColorBrewer::brewer.pal(n = 12, name = "Paired")

DotPlot(seu, features = "Dpp4", split.by = "orig.ident", cols = cols_clusters)

VlnPlot(cluster0, features = "Dpp4", split.by = "orig.ident", cols = cols_clusters)

cluster0 <- subset(seu, seurat_clusters == "0")


cluster0$seurat_clusters <- factor(cluster0$seurat_clusters, levels = "0")

VlnPlot(seu, features = "Dpp4", pt.size = 0.1, split.by = "orig.ident")


DotPlot(seu, features = c("Col1a1", "Col3a1","Dcn"))
VlnPlot(seu, features = c("Col1a1", "Col3a1","Dcn"), split.by = "orig.ident", pt.size = 0)



expr <- cluster0@assays$RNA$data
library(dplyr)
gene_expression <- expr %>% .["Dpp4",] %>% as.data.frame()
colnames(gene_expression) <- "Dpp4"
gene_expression$cell <- rownames(gene_expression)


gene_expression_sel <- gene_expression[which(gene_expression$Dpp4>0),]
qqqq <- cluster0[,rownames(gene_expression_sel)]

table(qqqq$orig.ident)
table(cluster0$orig.ident)


proportion_barplot <- function(pbmc, cluster1, cluster2, col.cluster1){
  proportion <- as.data.frame(prop.table(x = table(pbmc@meta.data[[cluster1]], pbmc@meta.data[[cluster2]]), margin = 2))
  colnames(proportion) <- c(cluster1, cluster2, "proportion")
  ggplot(proportion, aes_string(x = cluster2, group = cluster1)) + 
    geom_bar(aes_string(y="proportion", fill=cluster1), stat = "identity", width = 0.8)+ 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_y_continuous(limits=c(0, 1)) + 
    scale_x_discrete(limits = levels(pbmc@meta.data[[cluster2]])) + 
    scale_fill_manual(values = col.cluster1) +  
    cowplot::theme_cowplot() +
    ggtitle(label = "Proportion of cluster") # + coord_flip()
}


proportion_barplot(seu, cluster1 = "seurat_clusters", cluster2 = "orig.ident", col.cluster1 = cluster_31cols)

cluster_31cols <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00", 
                    '#B53E2B')






######################

seu@meta.data$cell.type <- seu@meta.data$seurat_clusters
Idents(seu) <- "cell.type"
seu <- RenameIdents(seu, '0' = "Fib", '1' = "Mac", '2' = "T", '3' = "Ep", '4' = "Myfib", '5' = "Ec", 
                    '6' = "Ly-Ec", '7' = "Neu")
seu$cell.type <- Idents(seu)
table(seu$cell.type)


#####

pdf(file = "2.pdf", height = 6, width = 7)
DimPlot(seu, reduction = "umap", group.by = "cell.type", label = T, repel = T, label.size = 5)
dev.off()

pdf(file = "2-1.pdf", height = 6, width = 7)
DimPlot(seu, reduction = "umap", group.by = "Group", label = T, repel = T, label.size = 5)
dev.off()

#####

marker <- c("Col1a1", "Dcn", "Itgbl1", "Lyz2", "Adgre1", "Cd68", "Cd3d", "Cd28", "Il2rb", "Krt14", "Krt5", "Acta2", "Pecam1", "Cdh5", "Lyve1", "S100a9", "S100a8")

pdf(file = "3.pdf", width = 9, height = 5)
DotPlot(seu, features = marker, assay = "RNA", group.by = "cell.type", cols = c("gray90", "red")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "italic"))
dev.off()


#####

pdf(file = "4.pdf", width = 5, height = 10)
proportion_barplot(seu, cluster1 = "cell.type", cluster2 = "Group", col.cluster1 = cols_clusters)
dev.off()

###



seu_fib <- subset(seu, cell.type == "Fib")

seu_fib <- NormalizeData(seu_fib)
seu_fib <- FindVariableFeatures(seu_fib)
seu_fib <- ScaleData(seu_fib)
seu_fib <- RunPCA(seu_fib)

#library(harmony)
#seu <- RunHarmony(seu, group.by.vars = "Group")


seu_fib <- RunUMAP(seu_fib, dims = 1:30, reduction = "pca")
seu_fib <- FindNeighbors(seu_fib, dims = 1:30)
seu_fib <- FindClusters(seu_fib, resolution = 0.1)


library(Seurat)
pdf(file = "5.pdf", width = 6, height = 5)
DimPlot(seu_fib, reduction = "umap", label = T)
dev.off()

proportion_barplot(seu_fib, cluster1 = "seurat_clusters", cluster2 = "Group", col.cluster1 = cols_clusters)
VlnPlot(seu_fib, features = "Dpp4", split.by = "Group")
FeaturePlot(seu_fib, features = "Dpp4", split.by = "Group", label = T)

VlnPlot(seu_fib, features = "Angptl1", split.by = "Group")

VlnPlot(seu_fib, features = "Apcdd1", split.by = "Group")
VlnPlot(seu_fib, features = "Col18a1", split.by = "Group")
VlnPlot(seu_fib, features = "Col13a1", split.by = "Group")



seu_fib@meta.data$cell.type2 <- seu_fib@meta.data$seurat_clusters
Idents(seu_fib) <- "cell.type2"
seu_fib <- RenameIdents(seu_fib, '0' = "Fib0", '1' = "Fib1", '2' = "Fib2", '3' = "Fib3")
seu_fib$cell.type2 <- Idents(seu_fib)
table(seu_fib$cell.type2)



seu2 <- seu


seu$cell.type[match(colnames(seu_fib),colnames(seu))] =  seu_fib$cell.type2 


seu$cell.type3 <- "1"

Idents(seu_fib)



expr <- seu_fib@assays$RNA$data
library(dplyr)
gene_expression <- expr %>% .["Dpp4",] %>% as.data.frame()
colnames(gene_expression) <- "Dpp4"
gene_expression$cell <- rownames(gene_expression)


gene_expression_sel <- gene_expression[which(gene_expression$Dpp4>0),]
qqqq <- seu_fib[,rownames(gene_expression_sel)]

table(qqqq$Group)

####

degs.fibs <- FindAllMarkers(seu_fib, only.pos = T, logfc.threshold = 0.25, min.pct = 0.5)





Degs_sig_subset <- subset(degs.fibs, p_val_adj < 0.01)
enrich_data <- Degs_sig_subset[, c("cluster","gene")]
library(org.Mm.eg.db)
library(clusterProfiler)
go_subset_id <- bitr(enrich_data[,2],
                     fromType="SYMBOL",
                     toType=c("ENTREZID","ENSEMBL","SYMBOL"),
                     OrgDb="org.Mm.eg.db")
library(dplyr)
go_subset <- inner_join(enrich_data,go_subset_id,by=c("gene"="SYMBOL"))
GO_CompareCluster_Reslut <- compareCluster(gene~cluster, 
                                           data=go_subset, 
                                           fun="enrichGO", 
                                           OrgDb = org.Mm.eg.db, 
                                           keyType = 'SYMBOL', 
                                           ont = "ALL", 
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05, 
                                           qvalueCutoff  = 0.2)
conflicts_prefer(clusterProfiler::simplify)
GO_filtered <- simplify(GO_CompareCluster_Reslut, 
                        cutoff = 0.5, 
                        by = "p.adjust")

pdf(file = "6.pdf", height = 15, width = 10)
dotplot(GO_filtered, x = "cluster", color = "p.adjust", 
        showCategory = 6, split = "ONTOLOGY", label_format=150) + 
  facet_grid(ONTOLOGY~., scales = "free")
dev.off()

FeaturePlot(seu_fib, features = c("Apcdd1", "Col18a1", "Col13a1"))

pdf(file = "7.pdf", width = 10, height = 5)
VlnPlot(seu_fib, features = c("Apcdd1", "Col18a1"), pt.size = 0)
dev.off()


seu_fib1_cluster1 <- subset(seu_fib, seurat_clusters == "1")
table(seu_fib1_cluster1$Group)

#####
mypie_plot <- function(pbmc, group.by = NULL, split.by = NULL, col.group = NULL, ncol = NULL){
  require(ggrepel)
  pie_data <- as.data.frame(prop.table(x = table(pbmc@meta.data[[group.by]], pbmc@meta.data[[split.by]]), margin = 2))
  colnames(pie_data) <- c("cluster", "group",  "proportion")
  pie_data$cluster <- factor(pie_data$cluster, levels = levels(pbmc@meta.data[[group.by]]))
  pie_data$group <- factor(pie_data$group, levels = levels(pbmc@meta.data[[split.by]]))
  ggplot(data=pie_data, mapping=aes_string(x="1", y = "proportion", fill="cluster"))+
    geom_bar(stat="identity", width=1,position='stack', size=1)+
    coord_polar("y", start=0)+
    facet_wrap(~group, ncol = ncol)+
    scale_fill_manual(values=col.group)+
    geom_text_repel(stat="identity",aes_string(y="proportion", label = quote(paste0(round(proportion*100, 1), "%"))), size=3,
                    position=position_stack(vjust = 0.5)) +
    theme_minimal()+ 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank(),
          axis.ticks = element_blank(),
          strip.text = element_text(size = 12))
}


#####

mypie_plot(pbmc = seu_fib1_cluster1, group.by = "Group", col.group = cols_clusters[1:6])

proportion_barplot(seu_fib1_cluster1, cluster1 = "Group", cluster2 = "seurat_clusters", col.cluster1 = cols_clusters)

library(plotrix)

pie3D(table(seu_fib1_cluster1$Group),labels = c("POD", "CAR-Ms", "CAR-Trem2-Ms"),explode = 0.1)

hee <- seu_fib1_cluster1$Group
hee <- as.numeric(hee)
piepercent = paste(round(100*(hee)/sum(hee)), "%")
pie(table(hee), labels=piepercent, col=cols_Set1, family='GB1')

aaa <- table(seu_fib1_cluster1$Group)
names <- c("POD", "CAR-Ms", "CAR-Trem2-Ms")
cols = c("#ED1C24","#22B14C","#FFC90E")

pdf(file = "8.pdf", width = 5, height = 5)
pie(aaa, labels=names, col=cols)
dev.off()

###


VlnPlot(seu_fib1_cluster1, features = "Crabp1", group.by = "Group")

pdf(file = "9.pdf", width = 13, height = 10)
VlnPlot(seu_fib1_cluster1, features = c("Dpp4", "Col1a1", "Col3a1", "Tgfb1", "En1", "Myc", "Twist1", "Bmp4", "Trps1"), split.by = "Group", cols = cluster_31cols, pt.size = 0)
dev.off()






FeaturePlot(seu_fib1_cluster1, features ="Tgfb1", split.by = "Group")


VlnPlot(seu_fib1_cluster1, features ="Dkk3", split.by = "Group")



Idents(seu_fib1_cluster1) <- "Group"
degs.fibs.cluster1 <- FindAllMarkers(seu_fib1_cluster1, logfc.threshold = 0.25, min.pct = 0.5, only.pos = F)

write_csv(degs.fibs.cluster1, "degs.cluster1.csv")

Degs_sig <- subset(degs.fibs.cluster1, p_val_adj < 0.01)
library(dplyr)
Degs_sig %>% group_by(gene) %>% top_n(1, avg_log2FC) -> Degs_sig_no_dup
Degs_sig_no_dup  %>% group_by(cluster) %>% top_n(5, -p_val_adj) -> top5
top5 %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

library(scRNAtoolVis)
cols_clusters <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
pdf(file = "10.pdf", width = 8, height = 8)
AverageHeatmap(object = seu_fib1_cluster1, group.by = "Group", annoCol = T, myanCol = cols_clusters[2:4],
               markerGene = Degs_sig_no_dup$gene, row_title = "", 
               showRowNames = F, markGenes = top5$gene)

dev.off()


Degs_sig %>% group_by(gene) %>% top_n(1, avg_log2FC) -> Degs_sig_no_dup
Degs_sig_no_dup  %>% group_by(cluster) %>% top_n(10, -p_val_adj) -> top10
top10 %>% group_by(cluster) %>% top_n(10, avg_log2FC) -> top10

library(scRNAtoolVis)
cols_clusters <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
pdf(file = "10-1.pdf", width = 8, height = 8)
AverageHeatmap(object = seu_fib1_cluster1, group.by = "Group", annoCol = T, myanCol = cols_clusters[2:4],
               markerGene = Degs_sig_no_dup$gene, row_title = "", 
               showRowNames = F, markGenes = top10$gene)

dev.off()


#####

seu_ec <- subset(seu, cell.type == "Ec")
table(seu_ec$orig.ident)

table(seu$cell.type)

seu_ec <- NormalizeData(seu_ec)
seu_ec <- FindVariableFeatures(seu_ec)
seu_ec <- ScaleData(seu_ec)
seu_ec <- RunPCA(seu_ec)

#library(harmony)
#seu <- RunHarmony(seu, group.by.vars = "Group")


seu_ec <- RunUMAP(seu_ec, dims = 1:30, reduction = "pca")
seu_ec <- FindNeighbors(seu_ec, dims = 1:30)
seu_ec <- FindClusters(seu_ec, resolution = 0.3)



pdf(file = "11.pdf", width = 5, height = 4)
DimPlot(seu_ec, reduction = "umap", label = T)
dev.off()

########
Idents(seu_ec)
degs.ec <- FindAllMarkers(seu_ec, only.pos = T, logfc.threshold = 0.25, min.pct = 0.3)


Degs_sig_subset <- subset(degs.ec, p_val_adj < 0.01)
enrich_data <- Degs_sig_subset[, c("cluster","gene")]
library(org.Mm.eg.db)
library(clusterProfiler)
go_subset_id <- bitr(enrich_data[,2],
                     fromType="SYMBOL",
                     toType=c("ENTREZID","ENSEMBL","SYMBOL"),
                     OrgDb="org.Mm.eg.db")
library(dplyr)
go_subset <- inner_join(enrich_data,go_subset_id,by=c("gene"="SYMBOL"))
GO_CompareCluster_Reslut <- compareCluster(gene~cluster, 
                                           data=go_subset, 
                                           fun="enrichGO", 
                                           OrgDb = org.Mm.eg.db, 
                                           keyType = 'SYMBOL', 
                                           ont = "ALL", 
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05, 
                                           qvalueCutoff  = 0.2)
conflicts_prefer(clusterProfiler::simplify)
GO_filtered <- simplify(GO_CompareCluster_Reslut, 
                        cutoff = 0.5, 
                        by = "p.adjust")

pdf(file = "12.pdf", height = 15, width = 8)
dotplot(GO_filtered, x = "cluster", color = "p.adjust", 
        showCategory = 6, split = "ONTOLOGY", label_format=150) + 
  facet_grid(ONTOLOGY~., scales = "free")
dev.off()



####
Degs_sig_subset %>% group_by(gene) %>% top_n(1, avg_log2FC) -> Degs_sig_no_dup
Degs_sig_no_dup  %>% group_by(cluster) %>% top_n(5, -p_val_adj) -> top5
top5 %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

library(scRNAtoolVis)
cols_clusters <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
pdf(file = "13.pdf", width = 8, height = 8)
AverageHeatmap(object = seu_ec, group.by = "seurat_clusters", annoCol = T, myanCol = cols_clusters[2:4],
               markerGene = Degs_sig_no_dup$gene, row_title = "", 
               showRowNames = F, markGenes = top5$gene)

dev.off()

write_csv(Degs_sig_subset, file = "degs.ec.csv")


VlnPlot(seu_ec, features = c("Adamt51"), pt.size = 0, split.by = "Group")

seu_ec_cluster2 <- subset(seu_ec, seurat_clusters == "2")
Idents(seu_ec_cluster2) <- "Group"
table(seu_ec_cluster2$orig.ident)
deg.ec.cluster2 <- FindAllMarkers(seu_ec_cluster2, logfc.threshold = 0.25, min.pct = 0.3, only.pos = F)
write_csv(deg.ec.cluster2, file = "degs.ec.cluster2.csv")



#####
deg.ec.cluster2 %>% group_by(gene) %>% top_n(1, avg_log2FC) -> Degs_sig_no_dup
Degs_sig_no_dup  %>% group_by(cluster) %>% top_n(5, -p_val_adj) -> top5
top5 %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

library(scRNAtoolVis)
cols_clusters <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
pdf(file = "14.pdf", width = 8, height = 8)
AverageHeatmap(object = seu_ec_cluster2, group.by = "Group", annoCol = T, myanCol = cols_clusters[2:4],
               markerGene = Degs_sig_no_dup$gene, row_title = "", 
               showRowNames = F, markGenes = top5$gene)
dev.off()

###
deg.ec.cluster2 <- FindAllMarkers(seu_ec_cluster2, logfc.threshold = 0, min.pct = 0.3, only.pos = T)

Degs_sig_subset <- subset(deg.ec.cluster2, p_val_adj < 0.01)
enrich_data <- Degs_sig_subset[, c("cluster","gene")]
library(org.Mm.eg.db)
library(clusterProfiler)
go_subset_id <- bitr(enrich_data[,2],
                     fromType="SYMBOL",
                     toType=c("ENTREZID","ENSEMBL","SYMBOL"),
                     OrgDb="org.Mm.eg.db")
library(dplyr)
go_subset <- inner_join(enrich_data,go_subset_id,by=c("gene"="SYMBOL"))
GO_CompareCluster_Reslut <- compareCluster(gene~cluster, 
                                           data=go_subset, 
                                           fun="enrichGO", 
                                           OrgDb = org.Mm.eg.db, 
                                           keyType = 'SYMBOL', 
                                           ont = "ALL", 
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05, 
                                           qvalueCutoff  = 0.2)
conflicts_prefer(clusterProfiler::simplify)
GO_filtered <- simplify(GO_CompareCluster_Reslut, 
                        cutoff = 0.5, 
                        by = "p.adjust")

pdf(file = "12.pdf", height = 15, width = 8)
dotplot(GO_filtered, x = "cluster", color = "p.adjust", 
        showCategory = 6, split = "ONTOLOGY", label_format=150) + 
  facet_grid(ONTOLOGY~., scales = "free")
dev.off()


####
seu_ec_cluster0 <- subset(seu_ec, seurat_clusters == "0")
Idents(seu_ec_cluster0) <- "Group"
table(seu_ec_cluster0$orig.ident)
deg.ec.cluster0 <- FindAllMarkers(seu_ec_cluster0, logfc.threshold = 0.25, min.pct = 0.3, only.pos = F)
write_csv(deg.ec.cluster0, file = "degs.ec.cluster0.csv")

seu_ec_cluster1 <- subset(seu_ec, seurat_clusters == "1")
Idents(seu_ec_cluster1) <- "Group"
table(seu_ec_cluster1$orig.ident)
deg.ec.cluster1 <- FindAllMarkers(seu_ec_cluster1, logfc.threshold = 0.25, min.pct = 0.3, only.pos = F)
write_csv(deg.ec.cluster1, file = "degs.ec.cluster1.csv")

#####
deg.ec.cluster2 %>% group_by(gene) %>% top_n(1, avg_log2FC) -> Degs_sig_no_dup
Degs_sig_no_dup  %>% group_by(cluster) %>% top_n(5, -p_val_adj) -> top5
top5 %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

library(scRNAtoolVis)
cols_clusters <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
pdf(file = "14.pdf", width = 8, height = 8)
AverageHeatmap(object = seu_ec_cluster2, group.by = "Group", annoCol = T, myanCol = cols_clusters[2:4],
               markerGene = Degs_sig_no_dup$gene, row_title = "", 
               showRowNames = F, markGenes = top5$gene)
dev.off()

###
deg.ec.cluster2 <- FindAllMarkers(seu_ec_cluster2, logfc.threshold = 0, min.pct = 0.3, only.pos = T)

Degs_sig_subset <- subset(deg.ec.cluster2, p_val_adj < 0.01)
enrich_data <- Degs_sig_subset[, c("cluster","gene")]
library(org.Mm.eg.db)
library(clusterProfiler)
go_subset_id <- bitr(enrich_data[,2],
                     fromType="SYMBOL",
                     toType=c("ENTREZID","ENSEMBL","SYMBOL"),
                     OrgDb="org.Mm.eg.db")
library(dplyr)
go_subset <- inner_join(enrich_data,go_subset_id,by=c("gene"="SYMBOL"))
GO_CompareCluster_Reslut <- compareCluster(gene~cluster, 
                                           data=go_subset, 
                                           fun="enrichGO", 
                                           OrgDb = org.Mm.eg.db, 
                                           keyType = 'SYMBOL', 
                                           ont = "ALL", 
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05, 
                                           qvalueCutoff  = 0.2)
conflicts_prefer(clusterProfiler::simplify)
GO_filtered <- simplify(GO_CompareCluster_Reslut, 
                        cutoff = 0.5, 
                        by = "p.adjust")

pdf(file = "12.pdf", height = 15, width = 8)
dotplot(GO_filtered, x = "cluster", color = "p.adjust", 
        showCategory = 6, split = "ONTOLOGY", label_format=150) + 
  facet_grid(ONTOLOGY~., scales = "free")
dev.off()





seu_ec@meta.data$cell.type2 <- seu_ec@meta.data$seurat_clusters
Idents(seu_ec) <- "cell.type2"
seu_ec <- RenameIdents(seu_ec, '0' = "Ec0", '1' = "Ec1", '2' = "Ec2")
seu_ec$cell.type2 <- Idents(seu_ec)
table(seu_ec$cell.type2)




seu_cellchat <- merge(seu_fib, seu_ec)





table(seu_cellchat$cell.type2)


seu_cellchat$cell.type3 <- NULL
seu_cellchat$cell.type





table(seu_cellchat$orig.ident)











################################


POD <- subset(seu_cellchat, Group == "POD")
CARM <- subset(seu_cellchat, Group == "CAR-Ms")
CARTM <- subset(seu_cellchat, Group == "CAR-Trem2-Ms")


#####POD
library(CellChat)
data.input.POD <- GetAssayData(POD, assay = "RNA", slot = "data")
Idents(POD) <- POD$cell.type2
labels <- Idents(POD)
identity.POD <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.POD <- createCellChat(object = data.input.POD)
cellchat.seu.POD <- addMeta(cellchat.seu.POD, meta = identity.POD, 
                             meta.name = "labels")
cellchat.seu.POD <- setIdent(cellchat.seu.POD, ident.use = "labels") 

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.POD@DB <- CellChatDB.use

cellchat.seu.POD <- subsetData(cellchat.seu.POD)
cellchat.seu.POD <- identifyOverExpressedGenes(cellchat.seu.POD)
cellchat.seu.POD <- identifyOverExpressedInteractions(cellchat.seu.POD)
cellchat.seu.POD <- projectData(cellchat.seu.POD, PPI.mouse)
cellchat.seu.POD <- computeCommunProb(cellchat.seu.POD)
cellchat.seu.POD <- filterCommunication(cellchat.seu.POD, min.cells = 10)
df.net.POD <- subsetCommunication(cellchat.seu.POD)
write.csv(df.net.POD, "net.LR.POD.csv")

cellchat.seu.POD <- computeCommunProbPathway(cellchat.seu.POD)
df.netp.POD <- subsetCommunication(cellchat.seu.POD, slot.name = "netP")
write.csv(df.netp.POD, "net.patway.POD.csv")

cellchat.seu.POD <- aggregateNet(cellchat.seu.POD)
cellchat.seu.POD <- netAnalysis_computeCentrality(cellchat.seu.POD)

saveRDS(cellchat.seu.POD, file = "cellchat.POD.rds")



######


data.input.CARM <- GetAssayData(CARM, assay = "RNA", slot = "data")
Idents(CARM) <- CARM$cell.type2
labels <- Idents(CARM)
identity.CARM <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.CARM <- createCellChat(object = data.input.CARM)
cellchat.seu.CARM <- addMeta(cellchat.seu.CARM, meta = identity.CARM, 
                            meta.name = "labels")
cellchat.seu.CARM <- setIdent(cellchat.seu.CARM, ident.use = "labels") 

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.CARM@DB <- CellChatDB.use

cellchat.seu.CARM <- subsetData(cellchat.seu.CARM)
cellchat.seu.CARM <- identifyOverExpressedGenes(cellchat.seu.CARM)
cellchat.seu.CARM <- identifyOverExpressedInteractions(cellchat.seu.CARM)
cellchat.seu.CARM <- projectData(cellchat.seu.CARM, PPI.mouse)
cellchat.seu.CARM <- computeCommunProb(cellchat.seu.CARM)
cellchat.seu.CARM <- filterCommunication(cellchat.seu.CARM, min.cells = 10)
df.net.CARM <- subsetCommunication(cellchat.seu.CARM)
write.csv(df.net.CARM, "net.LR.CARM.csv")

cellchat.seu.CARM <- computeCommunProbPathway(cellchat.seu.CARM)
df.netp.CARM <- subsetCommunication(cellchat.seu.CARM, slot.name = "netP")
write.csv(df.netp.CARM, "net.patway.CARM.csv")

cellchat.seu.CARM <- aggregateNet(cellchat.seu.CARM)
cellchat.seu.CARM <- netAnalysis_computeCentrality(cellchat.seu.CARM)

saveRDS(cellchat.seu.CARM, file = "cellchat.CARM.rds")



#####
showDatabaseCategory(CellChatDB)
data.input.CARTM <- GetAssayData(CARTM, assay = "RNA", slot = "data")
Idents(CARTM) <- CARTM$cell.type2
labels <- Idents(CARTM)
identity.CARTM <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.CARTM <- createCellChat(object = data.input.CARTM)
cellchat.seu.CARTM <- addMeta(cellchat.seu.CARTM, meta = identity.CARTM, 
                             meta.name = "labels")
cellchat.seu.CARTM <- setIdent(cellchat.seu.CARTM, ident.use = "labels") 

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.CARTM@DB <- CellChatDB.use

cellchat.seu.CARTM <- subsetData(cellchat.seu.CARTM)
cellchat.seu.CARTM <- identifyOverExpressedGenes(cellchat.seu.CARTM)
cellchat.seu.CARTM <- identifyOverExpressedInteractions(cellchat.seu.CARTM)
cellchat.seu.CARTM <- projectData(cellchat.seu.CARTM, PPI.mouse)
cellchat.seu.CARTM <- computeCommunProb(cellchat.seu.CARTM)
cellchat.seu.CARTM <- filterCommunication(cellchat.seu.CARTM, min.cells = 10)
df.net.CARTM <- subsetCommunication(cellchat.seu.CARTM)
write.csv(df.net.CARTM, "net.LR.CARTM.csv")

cellchat.seu.CARTM <- computeCommunProbPathway(cellchat.seu.CARTM)
df.netp.CARTM <- subsetCommunication(cellchat.seu.CARTM, slot.name = "netP")
write.csv(df.netp.CARTM, "net.patway.CARTM.csv")

cellchat.seu.CARTM <- aggregateNet(cellchat.seu.CARTM)
cellchat.seu.CARTM <- netAnalysis_computeCentrality(cellchat.seu.CARTM)

saveRDS(cellchat.seu.CARTM, file = "cellchat.CARTM.rds")


#####
object.list <- list(POD = cellchat.seu.POD, 
                    CARM = cellchat.seu.CARM, 
                    CARTM = cellchat.seu.CARTM)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(cellchat, file = "cellchat.rds")


####
cols_Set1 <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
p1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), 
                          color.use = cols_Set1[2:4])
p2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), 
                          measure = "weight", color.use = cols_Set1[2:4])
gg <- p1+p2
pdf(file = "30.pdf", width = 10, height = 10)
gg
dev.off()

######
pdf(file = "31.pdf", height = 10, width = 10)
groupSize1 <- as.numeric(table(cellchat.seu.POD@idents))
netVisual_circle(cellchat.seu.POD@net$weight, vertex.weight = groupSize1, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength - POD", 
                 arrow.width = 0.7, arrow.size = 0.7)
dev.off()

pdf(file = "32.pdf", height = 10, width = 10)
groupSize2 <- as.numeric(table(cellchat.seu.CARM@idents))
netVisual_circle(cellchat.seu.CARM@net$weight, vertex.weight = groupSize2, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength - CARM", 
                 arrow.width = 0.7, arrow.size = 0.7)
dev.off()

pdf(file = "33.pdf", height = 10, width = 10)
groupSize3 <- as.numeric(table(cellchat.seu.CARTM@idents))
netVisual_circle(cellchat.seu.CARTM@net$weight, vertex.weight = groupSize3, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength - CARTM", 
                 arrow.width = 0.7, arrow.size = 0.7)
dev.off()

cellchat@meta$labels
cellchat@meta$datasets

###
pdf(file = "34CARMvs.POD.pdf", height = 10, width = 10)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",
                          comparison = c(1,2), 
                          title.name = "CARM vs. POD", 
                          margin = 0.2, arrow.width = 0.5, arrow.size = 0.7)
dev.off()

pdf(file = "35CARTMvs.CARM.pdf", height = 10, width = 10)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",
                          comparison = c(2,3), 
                          title.name = "CARTM vs. CARM", 
                          margin = 0.2, arrow.width = 0.5, arrow.size = 0.7)
dev.off()

pdf(file = "36CARTMvs.POD.pdf", height = 10, width = 10)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",
                          comparison = c(1,3), 
                          title.name = "CARTM vs. POD", 
                          margin = 0.2, arrow.width = 0.5, arrow.size = 0.7)
dev.off()

######


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + 
    colSums(x@net$count)-diag(x@net$count)})

# control the dot size in the different datasets
weight.MinMax <- c(min(num.link), max(num.link)) 
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], 
                                               title = names(object.list)[i], 
                                               weight.MinMax = weight.MinMax)
}

pdf(file = "37.pdf", width = 15, height = 5)
patchwork::wrap_plots(plots = gg, nrow = 1)
dev.off()
######


pairLR.use <- extractEnrichedLR(cellchat, signaling = c("TGFb"))
netVisual_bubble(cellchat, sources.use = c(6), targets.use = c(1), remove.isolate = F)
cellchat@netP$POD$pathways




####

pdf(file = "38.pdf", width = 7, height = 20)
netVisual_bubble(cellchat, sources.use = 6, targets.use = 1,  comparison = c(1, 2, 3), angle.x = 45)
dev.off()



####

pdf(file = "39-1.pdf", height = 10, width = 5)
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(1, 3))
dev.off()

pdf(file = "39-2.pdf", height = 10, width = 5)
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(2, 3))
dev.off()

pdf(file = "39-3.pdf", height = 10, width = 5)
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(1, 2))
dev.off()




####

pdf(file = "25ec.prop.pdf", width = 10, height = 10)
proportion_barplot(seu_ec, cluster1 = "seurat_clusters", cluster2 = "Group", col.cluster1 = cols_clusters)
dev.off()

####


seu_fib_myfib <- subset(seu, cell.type %in% c("Fib", "Myfib"))

seu_fib_myfib$cell.type <- factor(seu_fib_myfib$cell.type, levels = c("Fib", "Myfib"))
#############
Idents(seu_fib_myfib)


data <- as(as.matrix(seu_fib_myfib@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = seu_fib_myfib@meta.data)
gene_annotation=data.frame(gene_short_name = rownames(seu_fib_myfib[["RNA"]]),
                           stringsAsFactors=F)
rownames(gene_annotation)<-gene_annotation$gene_short_name
fd <- new("AnnotatedDataFrame", data = gene_annotation)
HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       expressionFamily=negbinomial.size())
library(dplyr)
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)




Degs_subset <- FindAllMarkers(seu_fib_myfib, only.pos = T, logfc.threshold = 0.25, min.pct = 0.3)

HSMM_ordering_genes=unique(Degs_subset$gene)
HSMM <- setOrderingFilter(HSMM,ordering_genes = HSMM_ordering_genes) 

HSMM <- reduceDimension(HSMM,
                        max_components = 3,
                        reduction_method = "DDRTree") 
HSMM <- orderCells(HSMM)


cols_Set1 <- RColorBrewer::brewer.pal(n = 9, name = "Set1")


pdf(file = "26-1.pdf", width = 10, height = 10)
plot_cell_trajectory(HSMM, color_by = "cell.type", show_branch_points = T) + scale_color_manual(values = cols_Set1)
dev.off()


pdf(file = "26-2.pdf", width = 10, height = 10)
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
dev.off()

pdf(file = "26-3.pdf", width = 10, height = 10)
plot_cell_trajectory(HSMM, color_by = "cell.type")  + scale_color_manual(values = cols_Set1)+
  facet_wrap(~Group, nrow = 2)
dev.off()

saveRDS(HSMM, file = "HSMM.max3.rds")

#####
HSMM <- reduceDimension(HSMM,
                        max_components = 2,
                        reduction_method = "DDRTree") 
HSMM <- orderCells(HSMM)


cols_Set1 <- RColorBrewer::brewer.pal(n = 9, name = "Set1")


pdf(file = "27-1.pdf", width = 10, height = 10)
plot_cell_trajectory(HSMM, color_by = "cell.type", show_branch_points = T) + scale_color_manual(values = cols_Set1)
dev.off()


pdf(file = "27-2.pdf", width = 10, height = 10)
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
dev.off()

pdf(file = "27-3.pdf", width = 10, height = 10)
plot_cell_trajectory(HSMM, color_by = "cell.type")  + scale_color_manual(values = cols_Set1)+
  facet_wrap(~Group, nrow = 1)
dev.off()

saveRDS(HSMM, file = "HSMM.max2.rds")



###
cds_DGT_pseudotimegenes <- differentialGeneTest(HSMM,fullModelFormulaStr = "~sm.ns(Pseudotime)")
cds_DGT_pseudotimegenes_sig <- subset(cds_DGT_pseudotimegenes, qval < 0.01)
saveRDS(cds_DGT_pseudotimegenes, file = "cds_DGT_pseudotimegenes.rds")

Time_genes <- cds_DGT_pseudotimegenes_sig %>% pull(gene_short_name) %>% as.character()
p3I <- plot_pseudotime_heatmap(HSMM[Time_genes,], 
                               num_cluster = 4, 
                               show_rownames = T, 
                               return_heatmap = T)
ggsave(p3I, filename = "Figure3I.pdf", height = 20, width = 7)



Idents(seu)
VlnPlot(seu, features = "", split.by = "Group", pt.size = 0)

pdf(file = "Lrg1.pdf", width = 10, height = 5)
VlnPlot(seu_ec, features = "Lrg1", split.by = "Group", pt.size = 0)
dev.off()


pdf(file = "Gata4.pdf", width = 10, height = 5)
VlnPlot(seu_ec, features = "Gata4", split.by = "Group", pt.size = 0)
dev.off()


seu_ec <- FindClusters(seu_ec, resolution = 0.2)
DimPlot(seu_ec)




degs.ec <- FindAllMarkers(seu_ec, only.pos = T, logfc.threshold = 0.25, min.pct = 0.3)


Degs_sig_subset <- subset(degs.ec, p_val_adj < 0.01)
enrich_data <- Degs_sig_subset[, c("cluster","gene")]
library(org.Mm.eg.db)
library(clusterProfiler)
go_subset_id <- bitr(enrich_data[,2],
                     fromType="SYMBOL",
                     toType=c("ENTREZID","ENSEMBL","SYMBOL"),
                     OrgDb="org.Mm.eg.db")
library(dplyr)
go_subset <- inner_join(enrich_data,go_subset_id,by=c("gene"="SYMBOL"))
GO_CompareCluster_Reslut <- compareCluster(gene~cluster, 
                                           data=go_subset, 
                                           fun="enrichGO", 
                                           OrgDb = org.Mm.eg.db, 
                                           keyType = 'SYMBOL', 
                                           ont = "ALL", 
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05, 
                                           qvalueCutoff  = 0.2)
conflicts_prefer(clusterProfiler::simplify)
GO_filtered <- simplify(GO_CompareCluster_Reslut, 
                        cutoff = 0.5, 
                        by = "p.adjust")

pdf(file = "12-2.pdf", height = 15, width = 8)
dotplot(GO_filtered, x = "cluster", color = "p.adjust", 
        showCategory = 6, split = "ONTOLOGY", label_format=150) + 
  facet_grid(ONTOLOGY~., scales = "free")
dev.off()


seu_ec_1 <- subset(seu_ec, seurat_clusters == "1")
Idents(seu_ec_1) <- "Group"

degs111 <- FindAllMarkers(seu_ec_1, only.pos = T, logfc.threshold = 0.25, min.pct = 0.5)


Idents(seu_ec) <- "Group"

deg2 <- FindAllMarkers(seu_ec, min.pct = 0.3, logfc.threshold = 0.25, only.pos = T)

write_csv(deg2, file = "内皮不分亚群组间差异.csv")




netAnalysis_signalingRole_network(cellchat.seu.POD, signaling = c("Ccl"))


Idents(seu_ec) <- "seurat_clusters"
pdf(file = "28.pdf", width = 5, height = 8)
DotPlot(seu_ec, features = "Cox4i2", split.by = "Group", cols = cluster_31cols)
dev.off()


VlnPlot(seu, features = "Cox4i2", group.by = "cell.type")
pdf(file = "29.pdf", width = 5, height = 8)
DotPlot(seu, features = "Lrg1", group.by = "cell.type", cols = cluster_31cols)
dev.off()


VlnPlot(seu_ec, features = "Cox4i2", pt.size = 0)

VlnPlot(seu, features = "Trem2", split.by = "Group")
DotPlot(seu, features = "Trem2", split.by = "Group", cols = cluster_31cols)


seu_mac <- subset(seu, cell.type == "Mac")
seu_mac$cell.type <- factor(seu_mac$cell.type, levels = c("Mac"))
VlnPlot(seu_mac, features = "Trem2", split.by = "Group")


seu_mac <- NormalizeData(seu_mac)
seu_mac <- FindVariableFeatures(seu_mac)
seu_mac <- ScaleData(seu_mac)
seu_mac <- RunPCA(seu_mac)



seu_mac <- RunUMAP(seu_mac, dims = 1:30, reduction = "pca")
seu_mac <- FindNeighbors(seu_mac, dims = 1:30)
seu_mac <- FindClusters(seu_mac, resolution = 0.05)

pdf(file = "5.pdf", width = 6, height = 5)
DimPlot(seu_mac, reduction = "umap", label = T)
dev.off()


proportion_barplot(seu_mac, cluster1 = "seurat_clusters", cluster2 = "Group", col.cluster1 = cols_clusters)

VlnPlot(seu_mac, features = "Trem2", split.by = "Group")
DotPlot(seu_mac, features = "Trem2", split.by = "Group", cols = cols_clusters)



VlnPlot(seu_fib, features = "Dpp4", split.by = "Group")


pathways.show <- c("TGFb") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) 
# 控制不同数据集的边的权重
pdf(file = "40TGFb.pdf", width = 10, height = 10)
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

cellchat@netP$POD$pathways
cellchat@netP$CARM$pathways
cellchat@netP$CARTM$pathways


cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NL", "LS")) # set factor level
plotGeneExpression(cellchat, signaling = "TGFb", split.by = "datasets", colors.ggplot = T)


VlnPlot(seu_ec, features = "Lrg1", split.by = "Group")
pdf(file = "lrg1_dotplot.pdf", width = 5, height = 9)
DotPlot(seu_ec, features = "Lrg1", split.by = "Group", cols = cols_Set1)
dev.off()


Cellratio <- prop.table(table(Idents(seu_fib), seu_fib$orig.ident), margin = 2)
Cellratio <- data.frame(Cellratio)
library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]


seu_fib$seurat_clusters



###########################################################2023-9-21
HSMM <- readRDS("~/liumin/Results/HSMM.max2.rds")
library(monocle)
pdf(file = "成纤维细胞拟时序.pdf", width = 13, height = 5)
plot_cell_trajectory(HSMM, color_by = "cell.type")  + scale_color_manual(values = cols_Set1)+
  facet_wrap(~Group, nrow = 1)
dev.off()


rm(HSMM);gc()


seu <- readRDS("~/liumin/Results/seu.rds")
library(ggplot2)
library(Seurat)




