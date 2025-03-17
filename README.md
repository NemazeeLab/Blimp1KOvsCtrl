# Blimp1KOvsCtrl<br>
# ScRNAseq set-up<br>
dataset_loc <- '/Users/Li/Desktop/R_scrna_seq' <br>
ids <- c("ctl_sample_feature_bc_matrix","blimp_sample_feature_bc_matrix") <br>
ctl.data <- Read10X_h5(file.path(dataset_loc, ids[1], "ctl_sample_feature_bc_matrix.h5"), use.names = T)<br>
blimp1.data <- Read10X_h5(file.path(dataset_loc, ids[2], "blimp_sample_feature_bc_matrix.h5"), use.names = T)<br>
sdata.ctl <- CreateSeuratObject(counts = ctl.data, project = "Til_B_CTRL", min.cells = 3)<br>
sdata.blimp1 <- CreateSeuratObject(counts = blimp1.data, project = "Til_B_BLIMP1", min.cells = 3)<br>
sdata.ctl$type = "CTL"<br>
sdata.blimp1$type = "BLIMP1"<br>
alldata <- merge(sdata.ctl, sdata.blimp1 , add.cell.ids = c("CTL" , "BLIMP1"))<br>
gc()<br>
alldata <- PercentageFeatureSet(alldata, "^MT-", col.name = "percent_mito")<br>
alldata <- PercentageFeatureSet(alldata, "^RP[SL]", col.name = "percent_ribo")<br>
alldata <- PercentageFeatureSet(alldata, "^HB[^(P)]", col.name = "percent_hb")<br>
alldata <- PercentageFeatureSet(alldata, "PECAM1|PF4", col.name = "percent_plat")<br>
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")<br>
VlnPlot(alldata, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + NoLegend()<br>
FeatureScatter(alldata, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)<br>
selected_c <- WhichCells(alldata, expression = nFeature_RNA > 200)<br>
selected_f <- rownames(alldata)[Matrix::rowSums(alldata) > 3]<br>
data.filt <- subset(alldata, features = selected_f, cells = selected_c)<br>
dim(data.filt)<br>
selected_mito <- WhichCells(data.filt, expression = percent_mito < 0.2)<br>
data.filt <- subset(data.filt, cells = selected_mito)<br>
<br>
dim(data.filt)<br>
<br>
table(data.filt$orig.ident)<br>
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")<br>
<br>
VlnPlot(data.filt, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
  NoLegend()<br>
<br>
data.filt = NormalizeData(data.filt)<br>
data.filt <- CellCycleScoring(object = data.filt, g2m.features = cc.genes$g2m.genes,
                              s.features = cc.genes$s.genes)<br>
VlnPlot(data.filt, features = c("S.Score", "G2M.Score"), group.by = "orig.ident",
        ncol = 4, pt.size = 0.1)<br>
suppressMessages(require(DoubletFinder))<br>
<br>
data.filt = FindVariableFeatures(data.filt, verbose = F)<br>
data.filt = ScaleData(data.filt, vars.to.regress = c("nFeature_RNA", "percent_mito"),
                      verbose = F)<br>
data.filt = RunPCA(data.filt, verbose = F, npcs = 20)<br>
data.filt = RunUMAP(data.filt, dims = 1:10, verbose = F)<br>
nExp <- round(ncol(data.filt) * 0.04)  # expect 4% doublets<br>
data.filt <- doubletFinder_v3(data.filt, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)<br>
<br>
DF.name = colnames(data.filt@meta.data)[grepl("DF.classification", colnames(data.filt@meta.data))]<br>
cowplot::plot_grid(ncol = 2, DimPlot(data.filt, group.by = "orig.ident") + NoAxes(),
                   DimPlot(data.filt, group.by = DF.name) + NoAxes())<br>
VlnPlot(data.filt, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)<br>
data.filt = data.filt[, data.filt@meta.data[, DF.name] == "Singlet"]<br>
dim(data.filt)<br>
<br>
umap = cbind("Barcode" = rownames(Embeddings(object = data.filt, reduction = "umap")), Embeddings(object = data.filt, reduction = "umap"))<br>
write.table(umap, file="/Users/Li/Desktop/R_scrna_seq/results/seurat_BlimpVSctl.csv", sep = ",", quote = F, row.names = F, col.names = T)<br>
save(data.filt, umap, file = "seurat_BlimpVSctl.Rdata")<br>
### to make Figure Figure S3a-c<br>
load(file = "seurat_BlimpVSctl.Rdata")<br>
batch_ids <- data.frame(barcode = rownames(data.filt@meta.data), 
                        batch_id = sample(0:2, NROW(data.filt@meta.data), replace = TRUE),
                        stringsAsFactors = FALSE)<br>
row.names(batch_ids) <- row.names(data.filt@meta.data)<br>
data.filt <- AddMetaData(object = data.filt, metadata = batch_ids, col.name = NULL)<br>
data.filt <- ScaleData(object = data.filt, vars.to.regress = 'batch_id')<br>
data.filt <- RunPCA(object = data.filt, 
                    # features = seurat@assays$RNA@var.features, 
                    ndims.print = 1:12, 
                    nfeatures.print = 12)<br>
<br>
set.seed(2020)<br>
data.filt <- FindNeighbors(object = data.filt, dims = 1:12)<br>
<br>
data.filt <- FindClusters(object = data.filt, 
                          reduction = "pca", 
                          dims = 1:12, 
                          resolution = 0.5,
                          random.seed = 2020)<br>
DimPlot(data.filt, reduction = 'umap', 
        label = TRUE, pt.size = 0.5) + NoLegend()<br>
DimPlot(data.filt, reduction = 'umap', pt.size = 0.5)<br>
DimPlot(data.filt, group.by = "orig.ident", reduction = 'umap', 
        label = TRUE, pt.size = 0.5)<br>
<br>
genes_to_check_Top = c('Cd79a','Ighd','Ms4a1','Mzb1','Jchain','Tnfrsf17','Fcgr1','Cd14','Treml4','Ace',
                       'Cd8a','Cd3d','Siglech','Ccr9','Cd200r3','Fcer1a',
                       'Tcrg-V6','Cd163l1','Tnfrsf9','Nkg7','Ccl5','Styk1','Kif13b','Mki67','Top2a')<br>
DotPlot(data.filt, group.by = 'seurat_clusters',split.by = "orig.ident",cols = c("red","blue"),
        features = unique(genes_to_check_Top)) + RotatedAxis()<br>
### to make Figure 3a-f<br>
#Subclustering of B cells <br>
B_sub_m = data.filt[,data.filt@meta.data$seurat_clusters %in% c(0,9)]<br>
sce.m = B_sub_m<br>
sce.m <- NormalizeData(sce.m, normalization.method = "LogNormalize", scale.factor = 1e4) <br>
sce.m <- FindVariableFeatures(sce.m, selection.method = 'vst', nfeatures = 2000)<br>
sce.m <- RunPCA(sce.m, features = VariableFeatures(object = sce.m)) <br>
<br>
sce.m <- FindNeighbors(sce.m, dims = 1:6)<br>
sce.m <- FindClusters(sce.m, resolution = 1 )<br>
<br>
head(Idents(sce.m), 2)<br>
table(sce.m$seurat_clusters) <br>
sce.m <- RunUMAP(sce.m, dims = 1:6)<br>
DimPlot(sce.m, reduction = 'umap') #Figure 3a<br>
DimPlot(sce. m, group.by = "orig.ident") + NoAxes() #Figure 3b<br>
> table(sce.m$orig.ident)<br>
table(Idents(sce.m))<br>
> table(Idents(sce.m),sce.m$orig.ident) #Figure 3c<br>
genes_to_check_fomzb = c('Ighd', 'Ccr7', 'Ighm', 'Mzb1')<br>
DotPlot(sce.m, group.by = 'seurat_clusters',
        features = unique(genes_to_check_fomzb)) + RotatedAxis()<br>
FeaturePlot(sce.m, features = c('Ighd', 'Ccr7', 'Ighm', 'Mzb1')) # Figure 3d<br>
<br>
df.sub.m.markers <- FindAllMarkers(object = sce.m, min.pct = 0.25, thresh.use = 0.25)<br>
write.table(df.sub.m.markers, file="/Users/Li/Desktop/R_scrna_seq/results/seurat_df_sub_m_markers_BlimpVSctl.csv", sep = ",", quote = F, row.names = F, col.names = T)<br>
df <- read.table("/Users/Li/Desktop/R_scrna_seq/results/sub.m/seurat_df_sub_m_markers_BlimpVSctl.csv", sep = ",", fill = TRUE, col.names=c("p_val","avg_log2FC","pct.1","pct.2","p_val_adj","cluster","gene"), header=T)<br>
head(df)<br>
df$label <- ifelse(df$p_val_adj<0.01,"adjust P-val<0.01","adjust P-val>=0.01")<br>
head(df)<br>
top10sig0 <- filter(df,cluster=="0") %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))<br>
head(top10sig0)<br>
top10sig1 <- filter(df,cluster=="1") %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))<br>
head(top10sig1)<br>
top10sig2 <- filter(df,cluster=="2") %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))<br>
head(top10sig2)<br>
top10sig3 <- filter(df,cluster=="3") %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))<br>
head(top10sig3)<br>
top10sig4 <- filter(df,cluster=="4") %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))<br>
head(top10sig4)<br>
top10sig5 <- filter(df,cluster=="5") %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))<br>
head(top10sig5)<br>
top10sig <- rbind(top10sig0,top10sig1,top10sig2,top10sig3,top10sig4,top10sig5)<br>
df$size <- case_when(!(df$gene %in% top10sig$gene)~ 1,
                     df$gene %in% top10sig$gene ~ 2)<br>
dt <- filter(df,size==1)<br>
head(dt)<br>
p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = cluster, y = avg_log2FC, color = label))<br>
p<br>
p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top10sig,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)<br>
p<br>
dfbar<-data.frame(x=c(0,1,2,3,4,5),
                  y=c(1.0,0.5,1.2,2.0,8.0,2.0))<br>
dfbar1<-data.frame(x=c(0,1,2,3,4,5),
                   y=c(-2.0,-1.5,-1.5,-1.8,-5.0,-1.5))<br>
p1 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)<br>
p1<br>
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
              width =0.4)<br>
p2<br>
dfcol<-data.frame(x=c(0:5),
                  y=0,
                  label=c(0:5))<br>
mycol <- c("#B2182B", "#D6604D", "#F4A582", "#D1E5F0", "#92C5DE","#4393C3")<br>
p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=0.4,
                     color = "black",
                     fill = mycol,
                     alpha = 0.6,
                     show.legend = F)<br>
p3<br>
p4 <- p3+
  geom_text_repel(
    data=top10sig,
    aes(x=cluster,y=avg_log2FC,label=gene),
    force = 1.2,
    size = 2,
    arrow = arrow(length = unit(0.008, "npc"),
                  type = "open", ends = "last")
  )<br>
p4<br>
p5 <- p4 +
  scale_color_manual(name=NULL,
                     values = c("red","black"))<br>
p5<br>
p6 <- p5+
  labs(x="Cluster",y="average logFC")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =3,
            color ="white")<br>
p6<br>
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
  )<br>
p7 #Figure 3e<br>
genes_to_check_Top = c('Fcer2a','Vpreb3','H3f3a','Tmem108','Ighd','Il4i1','Junb','Rpl35a','Rps29','Limd2',
                       'Rps15a', 'Rps16','Rps7','Rpl9','Ebf1','Rpl37a','Rps27','Fth1','Dusp2',
                       'S100a6','Apoe','Psap','Napsa','Nfatc1','Ptprj','Ms4a1','Cd9','Ighm','Cyp4f18',
                       'Ighg2c','Cxcr3','Ighg2b','Tbx21','Ass1','Fah','Cd52','Anxa2','Bhlhe41','Scimp','Cd86',
                       'Tnfrsf17','Eaf2','Jchain','Derl3','Cd28','Plpp5','Sdc1','Tent5c','Ccr10','Tigit',
                       'Ifi214','Ifit3','Ifit2','Usp18','Irf7','Ifi208','Oasl2','Irgm1','Oasl1','Ifi47','Tor3a')<br>
DotPlot(sce.m, group.by = 'seurat_clusters',
        features = unique(genes_to_check_Top)) + RotatedAxis() + theme(axis.text.x = setting) + scale_color_gradient2(low = "blue", mid = "white", high = "red") # Figure 3f<br>

### VDJ combined with scRNAseq to make Figure 4a,b & Figure S4a-d<br>
CTL <- read.csv("/Users/Li/Desktop/VDJ/CTRL_03idCTRLBCELLtyVDJdt111522/CTL_filtered_contig_annotations.csv")<br>
BLIMP1 <- read.csv("/Users/Li/Desktop/VDJ/BLIMP_01idBLIMPBCELLtyVDJdt111122/BLIMP_filtered_contig_annotations.csv")<br>
<br>
contig_list_2 <- list(CTL, BLIMP1)<br>
data("contig_list_2") #the data built into scRepertoire<br>
head(contig_list_2[[1]])<br>
combined.2 <- combineBCR(contig_list_2, samples = c( "CTL","BLIMP1"))<br>
head(combined.2$CTL)<br>
head(combined.2$BLIMP1)<br>
<br>
subset.2 <- subsetContig(combined.2, name = "sample", 
                       variables = c("CTL","BLIMP1"))<br>
quantContig(combined.2, cloneCall="strict", scale = T)<br>
abundanceContig(combined.2, cloneCall = "gene", scale = F)<br>
lengthContig(combined.2, cloneCall="aa", chain = "both") <br>
lengthContig(combined.2, cloneCall="aa", chain = "IGH") <br>
lengthContig(combined.2, cloneCall="aa", chain = "IGL")<br>
compareClonotypes(combined.2, 
                  numbers = 10, 
                  samples = c("CTL","BLIMP1"), 
                  cloneCall="aa", 
                  graph = "alluvial")<br>
scatterClonotype(combined.2, 
                 cloneCall ="aa", 
                 x.axis = "CTL", 
                 y.axis = "BLIMP1",
                 chain = "both",
                 dot.size = "total",
                 graph = "proportion",
                 exportTable = T)<br>
options(max.print=999999)<br>

scatterClonotype(combined.2, 
                 cloneCall ="aa", 
                 x.axis = "CTL", 
                 y.axis = "BLIMP1",
                 chain = "both",
                 dot.size = "total",
                 graph = "proportion")<br>

sce.m <- get(load(file = 'sce.m.B.subset.Rdata'))<br>
DimPlot(sce.m, label = T) + NoLegend()<br>
<br>
sce.m <- combineExpression(combined.2, sce.m, 
                            cloneCall="gene", 
                            group.by = "sample", 
                            proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))<br>
head(sce.m@meta.data)<br>
colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                                       "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                                       "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))<br>
DimPlot(sce.m, group.by = "cloneType") +
  scale_color_manual(values = colorblind_vector(4), na.value="grey")<br>
<br>
slot(sce.m, "meta.data")$cloneType <- factor(slot(sce.m, "meta.data")$cloneType, 
                                              levels = c("Hyperexpanded (100 < X <= 500)", 
                                                         "Large (20 < X <= 100)", 
                                                         "Medium (5 < X <= 20)", 
                                                         "Small (1 < X <= 5)", 
                                                         "Single (0 < X <= 1)", NA))<br>
DimPlot(sce.m, group.by = "cloneType", split.by = "orig.ident") +
  scale_color_manual(values = colorblind_vector(3), na.value="grey") + 
  theme(plot.title = element_blank())<br>
###to make Figure 4a,b<br>
###
data <- GetAssayData(sce.m, assay = 'RNA', layer ='counts')<br>
cell_metadata <- sce.m@meta.data<br>
gene_annotation <- data.frame(gene_short_name = rownames(data))<br>
rownames(gene_annotation) <- rownames(data)<br>
cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)<br>
<br>
cds<- preprocess_cds(cds, method="PCA", num_dim = 20)<br>
plot_pc_variance_explained(cds)<br>
<br>
cds <- reduce_dimension(cds,preprocess_method = "PCA", reduction_method ='UMAP')<br>
<br>
cds.embed<-cds@int_colData$reducedDims$UMAP<br>
int.embed <-Embeddings(sce.m,reduction="umap")<br>
int.embed <- int.embed[rownames(cds.embed),]<br>
cds@int_colData$reducedDims$UMAP <- int.embed<br>
cds<-cluster_cells(cds,reduction_method = "UMAP",k=8)<br>
<br>
cds <- learn_graph(cds, learn_graph_control=list(ncenter=197), use_partition = FALSE)<br>
cds<- order_cells(cds)<br>
plot_cells(cds)<br>
<br>
plot_cells(cds,color_cells_by = "seurat_clusters")<br>
#to make Figure S4a<br>
###<br>
df <- read.table("/Users/Li/Desktop/R_scrna_seq/results/sub.m/df_sub.m_all.csv", sep = ",", fill = TRUE, col.names=c("p_val","avg_log2FC","pct.1","pct.2","p_val_adj","cluster","gene"), header=T)<br>
dfsample <- split(df$gene,df$cluster)<br>
length(dfsample)<br>
dfsample$`0` = bitr(dfsample$`0`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")<br>
dfsample$`1` = bitr(dfsample$`1`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")<br>
dfsample$`2` = bitr(dfsample$`2`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")<br>
dfsample$`3` = bitr(dfsample$`3`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")<br>
dfsample$`4` = bitr(dfsample$`4`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")<br>
dfsample$`5` = bitr(dfsample$`5`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")<br>
genelist <- list("0" = dfsample$`0`$ENTREZID, 
                 "1" = dfsample$`1`$ENTREZID,
                 "2" = dfsample$`2`$ENTREZID,
                 "3" = dfsample$`3`$ENTREZID,
                 "4" = dfsample$`4`$ENTREZID,
                 "5" = dfsample$`5`$ENTREZID)<br>
GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO",ont = "ALL", qvalueCutoff = 0.05, pvalueCutoff = 0.05, OrgDb = "org.Mm.eg.db")<br>
dotplot(GOclusterplot)<br>
KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG", organism="mmu")<br>
dotplot(KEGGclusterplot)<br>
###to make Figure S4b,d<br>
<br>
###<br>
genes_to_check_GC = c("Irf4",  "Myc","Dkc1", "Cd52","Mif", "Hspd1","Atp5b","Ran", "Cd83","Cd72",  
                      "Bcl6","Son","Tcf3","Pax5","Bach2","Icosl", "Gadd45a",    "Fam193a" )<br>
DotPlot(sce.m, group.by = 'seurat_clusters',
        features = unique(genes_to_check_GC)) + RotatedAxis()+ scale_color_gradient2(low = "blue", mid = "white",high = "red") <br>
#to make Figure S4c<br>


