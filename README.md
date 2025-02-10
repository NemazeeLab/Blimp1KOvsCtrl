# Blimp1KOvsCtrl
# ScRNAseq 
dataset_loc <- '/Users/Li/Desktop/R_scrna_seq' 
ids <- c("ctl_sample_feature_bc_matrix","blimp_sample_feature_bc_matrix") 
ctl.data <- Read10X_h5(file.path(dataset_loc, ids[1], "ctl_sample_feature_bc_matrix.h5"), use.names = T)
blimp1.data <- Read10X_h5(file.path(dataset_loc, ids[2], "blimp_sample_feature_bc_matrix.h5"), use.names = T)
sdata.ctl <- CreateSeuratObject(counts = ctl.data, project = "Til_B_CTRL", min.cells = 3)
sdata.blimp1 <- CreateSeuratObject(counts = blimp1.data, project = "Til_B_BLIMP1", min.cells = 3)
sdata.ctl$type = "CTL"
sdata.blimp1$type = "BLIMP1"
alldata <- merge(sdata.ctl, sdata.blimp1 , add.cell.ids = c("CTL" , "BLIMP1"))
gc()
alldata <- PercentageFeatureSet(alldata, "^MT-", col.name = "percent_mito")
alldata <- PercentageFeatureSet(alldata, "^RP[SL]", col.name = "percent_ribo")
alldata <- PercentageFeatureSet(alldata, "^HB[^(P)]", col.name = "percent_hb")
alldata <- PercentageFeatureSet(alldata, "PECAM1|PF4", col.name = "percent_plat")
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
VlnPlot(alldata, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + NoLegend()
FeatureScatter(alldata, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
selected_c <- WhichCells(alldata, expression = nFeature_RNA > 200)
selected_f <- rownames(alldata)[Matrix::rowSums(alldata) > 3]
data.filt <- subset(alldata, features = selected_f, cells = selected_c)
dim(data.filt)
selected_mito <- WhichCells(data.filt, expression = percent_mito < 0.2)
data.filt <- subset(data.filt, cells = selected_mito)

dim(data.filt)

table(data.filt$orig.ident)
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")

VlnPlot(data.filt, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
  NoLegend()

data.filt = NormalizeData(data.filt)
data.filt <- CellCycleScoring(object = data.filt, g2m.features = cc.genes$g2m.genes,
                              s.features = cc.genes$s.genes)
VlnPlot(data.filt, features = c("S.Score", "G2M.Score"), group.by = "orig.ident",
        ncol = 4, pt.size = 0.1)
suppressMessages(require(DoubletFinder))

data.filt = FindVariableFeatures(data.filt, verbose = F)
data.filt = ScaleData(data.filt, vars.to.regress = c("nFeature_RNA", "percent_mito"),
                      verbose = F)
data.filt = RunPCA(data.filt, verbose = F, npcs = 20)
data.filt = RunUMAP(data.filt, dims = 1:10, verbose = F)
nExp <- round(ncol(data.filt) * 0.04)  # expect 4% doublets
data.filt <- doubletFinder_v3(data.filt, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

DF.name = colnames(data.filt@meta.data)[grepl("DF.classification", colnames(data.filt@meta.data))]
cowplot::plot_grid(ncol = 2, DimPlot(data.filt, group.by = "orig.ident") + NoAxes(),
                   DimPlot(data.filt, group.by = DF.name) + NoAxes())
VlnPlot(data.filt, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
data.filt = data.filt[, data.filt@meta.data[, DF.name] == "Singlet"]
dim(data.filt)

umap = cbind("Barcode" = rownames(Embeddings(object = data.filt, reduction = "umap")), Embeddings(object = data.filt, reduction = "umap"))
write.table(umap, file="/Users/Li/Desktop/R_scrna_seq/results/seurat_BlimpVSctl.csv", sep = ",", quote = F, row.names = F, col.names = T)
save(data.filt, umap, file = "seurat_BlimpVSctl.Rdata")
### to make Figure Figure S3a-c & Figure 3a-f
load(file = "seurat_BlimpVSctl.Rdata")
batch_ids <- data.frame(barcode = rownames(data.filt@meta.data), 
                        batch_id = sample(0:2, NROW(data.filt@meta.data), replace = TRUE),
                        stringsAsFactors = FALSE)
row.names(batch_ids) <- row.names(data.filt@meta.data)
data.filt <- AddMetaData(object = data.filt, metadata = batch_ids, col.name = NULL)
data.filt <- ScaleData(object = data.filt, vars.to.regress = 'batch_id')
data.filt <- RunPCA(object = data.filt, 
                    # features = seurat@assays$RNA@var.features, 
                    ndims.print = 1:12, 
                    nfeatures.print = 12)

set.seed(2020)
data.filt <- FindNeighbors(object = data.filt, dims = 1:12)

data.filt <- FindClusters(object = data.filt, 
                          reduction = "pca", 
                          dims = 1:12, 
                          resolution = 0.5,
                          random.seed = 2020)
DimPlot(data.filt, reduction = 'umap', 
        label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(data.filt, reduction = 'umap', pt.size = 0.5)
DimPlot(data.filt, group.by = "orig.ident", reduction = 'umap', 
        label = TRUE, pt.size = 0.5)

genes_to_check_Top = c('Cd79a','Ighd','Ms4a1','Mzb1','Jchain','Tnfrsf17','Fcgr1','Cd14','Treml4','Ace',
                       'Cd8a','Cd3d','Siglech','Ccr9','Cd200r3','Fcer1a',
                       'Tcrg-V6','Cd163l1','Tnfrsf9','Nkg7','Ccl5','Styk1','Kif13b','Mki67','Top2a')
DotPlot(data.filt, group.by = 'seurat_clusters',split.by = "orig.ident",cols = c("red","blue"),
        features = unique(genes_to_check_Top)) + RotatedAxis()
#to make Figure S3a-c

#Subclustering of B cells 
B_sub_m = data.filt[,data.filt@meta.data$seurat_clusters %in% c(0,9)]
sce.m = B_sub_m
sce.m <- NormalizeData(sce.m, normalization.method = "LogNormalize", scale.factor = 1e4) 
sce.m <- FindVariableFeatures(sce.m, selection.method = 'vst', nfeatures = 2000)
sce.m <- RunPCA(sce.m, features = VariableFeatures(object = sce.m)) 

sce.m <- FindNeighbors(sce.m, dims = 1:6)
sce.m <- FindClusters(sce.m, resolution = 1 )

head(Idents(sce.m), 2)
table(sce.m$seurat_clusters) 
sce.m <- RunUMAP(sce.m, dims = 1:6)
DimPlot(sce.m, reduction = 'umap')
DimPlot(sce. m, group.by = "orig.ident") + NoAxes()
> table(sce.m$orig.ident)
table(Idents(sce.m))
> table(Idents(sce.m),sce.m$orig.ident)
genes_to_check_fomzb = c('Ighd', 'Ccr7', 'Ighm', 'Mzb1')
DotPlot(sce.m, group.by = 'seurat_clusters',
        features = unique(genes_to_check_fomzb)) + RotatedAxis()
FeaturePlot(sce.m, features = c('Ighd', 'Ccr7', 'Ighm', 'Mzb1'))

df.sub.m.markers <- FindAllMarkers(object = sce.m, min.pct = 0.25, thresh.use = 0.25)
write.table(df.sub.m.markers, file="/Users/Li/Desktop/R_scrna_seq/results/seurat_df_sub_m_markers_BlimpVSctl.csv", sep = ",", quote = F, row.names = F, col.names = T)
df <- read.table("/Users/Li/Desktop/R_scrna_seq/results/sub.m/seurat_df_sub_m_markers_BlimpVSctl.csv", sep = ",", fill = TRUE, col.names=c("p_val","avg_log2FC","pct.1","pct.2","p_val_adj","cluster","gene"), header=T)
head(df)
df$label <- ifelse(df$p_val_adj<0.01,"adjust P-val<0.01","adjust P-val>=0.01")
head(df)
top10sig0 <- filter(df,cluster=="0") %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))
head(top10sig0)
top10sig1 <- filter(df,cluster=="1") %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))
head(top10sig1)
top10sig2 <- filter(df,cluster=="2") %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))
head(top10sig2)
top10sig3 <- filter(df,cluster=="3") %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))
head(top10sig3)
top10sig4 <- filter(df,cluster=="4") %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))
head(top10sig4)
top10sig5 <- filter(df,cluster=="5") %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))
head(top10sig5)
top10sig <- rbind(top10sig0,top10sig1,top10sig2,top10sig3,top10sig4,top10sig5)
df$size <- case_when(!(df$gene %in% top10sig$gene)~ 1,
                     df$gene %in% top10sig$gene ~ 2)
dt <- filter(df,size==1)
head(dt)
p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = cluster, y = avg_log2FC, color = label))
p
p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top10sig,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)
p
dfbar<-data.frame(x=c(0,1,2,3,4,5),
                  y=c(1.0,0.5,1.2,2.0,8.0,2.0))
dfbar1<-data.frame(x=c(0,1,2,3,4,5),
                   y=c(-2.0,-1.5,-1.5,-1.8,-5.0,-1.5))
p1 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)
p1
p2 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = dt,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top10sig,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)
p2
dfcol<-data.frame(x=c(0:5),
                  y=0,
                  label=c(0:5))
mycol <- c("#B2182B", "#D6604D", "#F4A582", "#D1E5F0", "#92C5DE","#4393C3")
p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=0.4,
                     color = "black",
                     fill = mycol,
                     alpha = 0.6,
                     show.legend = F)
p3
p4 <- p3+
  geom_text_repel(
    data=top10sig,
    aes(x=cluster,y=avg_log2FC,label=gene),
    force = 1.2,
    size = 2,
    arrow = arrow(length = unit(0.008, "npc"),
                  type = "open", ends = "last")
  )
p4
p5 <- p4 +
  scale_color_manual(name=NULL,
                     values = c("red","black"))
p5
p6 <- p5+
  labs(x="Cluster",y="average logFC")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =3,
            color ="white")
p6
p7 <- p6+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 15)
  )
p7
genes_to_check_Top = c('Fcer2a','Vpreb3','H3f3a','Tmem108','Ighd','Il4i1','Junb','Rpl35a','Rps29','Limd2',
                       'Rps15a', 'Rps16','Rps7','Rpl9','Ebf1','Rpl37a','Rps27','Fth1','Dusp2',
                       'S100a6','Apoe','Psap','Napsa','Nfatc1','Ptprj','Ms4a1','Cd9','Ighm','Cyp4f18',
                       'Ighg2c','Cxcr3','Ighg2b','Tbx21','Ass1','Fah','Cd52','Anxa2','Bhlhe41','Scimp','Cd86',
                       'Tnfrsf17','Eaf2','Jchain','Derl3','Cd28','Plpp5','Sdc1','Tent5c','Ccr10','Tigit',
                       'Ifi214','Ifit3','Ifit2','Usp18','Irf7','Ifi208','Oasl2','Irgm1','Oasl1','Ifi47','Tor3a')
DotPlot(sce.m, group.by = 'seurat_clusters',
        features = unique(genes_to_check_Top)) + RotatedAxis() + theme(axis.text.x = setting) + scale_color_gradient2(low = "blue", mid = "white", high = "red")

### VDJ combined with scRNAseq
CTL <- read.csv("/Users/Li/Desktop/VDJ/CTRL_03idCTRLBCELLtyVDJdt111522/CTL_filtered_contig_annotations.csv")
BLIMP1 <- read.csv("/Users/Li/Desktop/VDJ/BLIMP_01idBLIMPBCELLtyVDJdt111122/BLIMP_filtered_contig_annotations.csv")

contig_list_2 <- list(CTL, BLIMP1)
data("contig_list_2") #the data built into scRepertoire
head(contig_list_2[[1]])
combined.2 <- combineBCR(contig_list_2, samples = c( "CTL","BLIMP1"))
head(combined.2$CTL)
head(combined.2$BLIMP1)

subset.2 <- subsetContig(combined.2, name = "sample", 
                       variables = c("CTL","BLIMP1"))
quantContig(combined.2, cloneCall="strict", scale = T)
abundanceContig(combined.2, cloneCall = "gene", scale = F)
lengthContig(combined.2, cloneCall="aa", chain = "both") 
lengthContig(combined.2, cloneCall="aa", chain = "IGH") 
lengthContig(combined.2, cloneCall="aa", chain = "IGL")
compareClonotypes(combined.2, 
                  numbers = 10, 
                  samples = c("CTL","BLIMP1"), 
                  cloneCall="aa", 
                  graph = "alluvial")
scatterClonotype(combined.2, 
                 cloneCall ="aa", 
                 x.axis = "CTL", 
                 y.axis = "BLIMP1",
                 chain = "both",
                 dot.size = "total",
                 graph = "proportion",
                 exportTable = T)
options(max.print=999999)

scatterClonotype(combined.2, 
                 cloneCall ="aa", 
                 x.axis = "CTL", 
                 y.axis = "BLIMP1",
                 chain = "both",
                 dot.size = "total",
                 graph = "proportion")

sce.m <- get(load(file = 'sce.m.B.subset.Rdata'))
DimPlot(sce.m, label = T) + NoLegend()

sce.m <- combineExpression(combined.2, sce.m, 
                            cloneCall="gene", 
                            group.by = "sample", 
                            proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
head(sce.m@meta.data)
colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                                       "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                                       "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))
DimPlot(sce.m, group.by = "cloneType") +
  scale_color_manual(values = colorblind_vector(4), na.value="grey")

slot(sce.m, "meta.data")$cloneType <- factor(slot(sce.m, "meta.data")$cloneType, 
                                              levels = c("Hyperexpanded (100 < X <= 500)", 
                                                         "Large (20 < X <= 100)", 
                                                         "Medium (5 < X <= 20)", 
                                                         "Small (1 < X <= 5)", 
                                                         "Single (0 < X <= 1)", NA))
DimPlot(sce.m, group.by = "cloneType", split.by = "orig.ident") +
  scale_color_manual(values = colorblind_vector(3), na.value="grey") + 
  theme(plot.title = element_blank())
###to make Figure 4a,b
###
data <- GetAssayData(sce.m, assay = 'RNA', layer ='counts')
cell_metadata <- sce.m@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)

cds<- preprocess_cds(cds, method="PCA", num_dim = 20)
plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds,preprocess_method = "PCA", reduction_method ='UMAP')

cds.embed<-cds@int_colData$reducedDims$UMAP
int.embed <-Embeddings(sce.m,reduction="umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
cds<-cluster_cells(cds,reduction_method = "UMAP",k=8)

cds <- learn_graph(cds, learn_graph_control=list(ncenter=197), use_partition = FALSE)
cds<- order_cells(cds)
plot_cells(cds)

plot_cells(cds,color_cells_by = "seurat_clusters")
#to make Figure S4a
###
df <- read.table("/Users/Li/Desktop/R_scrna_seq/results/sub.m/df_sub.m_all.csv", sep = ",", fill = TRUE, col.names=c("p_val","avg_log2FC","pct.1","pct.2","p_val_adj","cluster","gene"), header=T)
dfsample <- split(df$gene,df$cluster)
length(dfsample)
dfsample$`0` = bitr(dfsample$`0`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`1` = bitr(dfsample$`1`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`2` = bitr(dfsample$`2`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`3` = bitr(dfsample$`3`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`4` = bitr(dfsample$`4`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`5` = bitr(dfsample$`5`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
genelist <- list("0" = dfsample$`0`$ENTREZID, 
                 "1" = dfsample$`1`$ENTREZID,
                 "2" = dfsample$`2`$ENTREZID,
                 "3" = dfsample$`3`$ENTREZID,
                 "4" = dfsample$`4`$ENTREZID,
                 "5" = dfsample$`5`$ENTREZID)
GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO",ont = "ALL", qvalueCutoff = 0.05, pvalueCutoff = 0.05, OrgDb = "org.Mm.eg.db")
dotplot(GOclusterplot)
KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG", organism="mmu")
dotplot(KEGGclusterplot)
###to make Figure S4b,d

###
genes_to_check_GC = c("Irf4",  "Myc","Dkc1", "Cd52","Mif", "Hspd1","Atp5b","Ran", "Cd83","Cd72",  
                      "Bcl6","Son","Tcf3","Pax5","Bach2","Icosl", "Gadd45a",    "Fam193a" )
DotPlot(sce.m, group.by = 'seurat_clusters',
        features = unique(genes_to_check_GC)) + RotatedAxis()+ scale_color_gradient2(low = "blue", mid = "white",high = "red") 
#to make Figure S4c


