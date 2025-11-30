
library(tidyverse)
library(Seurat)
library(hdf5r)

data<- h5file("C:\\Users\\ASUS\\Downloads\\20k_NSCLC_DTC_3p_nextgem_donor_1_count_sample_feature_bc_matrix.h5",mode="r")
data$ls()
#hieu them ve file h5
nsclc.sparse.m <- Read10X_h5(filename = "C:\\Users\\ASUS\\Downloads\\20k_NSCLC_DTC_3p_nextgem_donor_1_count_sample_feature_bc_matrix.h5")
str(nsclc.sparse.m)
#Gene Exxpresion,Antibody Capture,Multiplexing Capture:
GE<-nsclc.sparse.m$`Gene Expression`
GE
#dữ liệu gene × cell 
nsclc.seurat.obj<- CreateSeuratObject(counts = GE,min.cells = 3, min.features = 200,project = 'NSCLC')
nsclc.seurat.obj
saveRDS(nsclc.seurat.obj,file='nsclc.rds')



# 1. QC -------
View(nsclc.seurat.obj@meta.data)
nsclc.seurat.obj@assays$RNA$counts
grep("^MT-", row.names(nsclc.seurat.obj), value=TRUE)
nsclc.seurat.obj[["perMT"]]<- PercentageFeatureSet(nsclc.seurat.obj,pattern = "^MT-")#gen MT
nsclc.seurat.obj[["perribo"]]<-PercentageFeatureSet(nsclc.seurat.obj,pattern="^RB[SL]")
view(nsclc.seurat.obj)
VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "perMT","perribo"), ncol = 3)
#nFeature_RNA tập hợp nhiều ở dười khoảng 2500
boxplot(nsclc.seurat.obj@meta.data[c("nFeature_RNA", "nCount_RNA", "perMT")])
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

q<- quantile(nsclc.seurat.obj$nCount_RNA,probs = c(0.25,0.98))
q[2]
IQR<-q[2]-q[1]
upper<- q[2]+1.5*IQR
upper
#Filtering
nsclc.seurat.obj<- subset(nsclc.seurat.obj,subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                            perMT< 5& nCount_RNA< q[2])
#Normalization - Sequencing depth: sl reads trung binh tren moi  cell.
#Library size: sum read hoặc counts cua cell

nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
str(nsclc.seurat.obj)

#Find Highly variable features
nsclc.seurat.obj<- FindVariableFeatures(nsclc.seurat.obj, selection.method = 'vst',nfeatures = 2000)

#Variance stabilizing transformation để t độ biến thiên
top10<-head(VariableFeatures(nsclc.seurat.obj),10)
LabelPoints(plot= VariableFeaturePlot(nsclc.seurat.obj),points=top10,repel = TRUE)
#Nếu một cell có 50,000 reads và cell khác có 10,000 reads, NormalizeData sẽ đưa chúng về cùng thang để so sánh.
# Scaling dua cac gen ve cung thang do
genes<-row.names(nsclc.seurat.obj)
genes
nsclc.seurat.obj<-ScaleData(nsclc.seurat.obj, features = genes)
#Gene A có biểu hiện rất cao, gene B thấp. Sau ScaleData, cả hai đều được đưa về thang chuẩn để đóng góp công bằng vào phân tích.
#PCA reduce dimension
nsclc.seurat.obj<- RunPCA(nsclc.seurat.obj,features = VariableFeatures(object = nsclc.seurat.obj))
print(nsclc.seurat.obj[['pca']],dims=1:5,nfeatures=5)#chon 5 thanh phan dau tien top 5 gene dong gop nhieu nhat
#Assay:du lieu tho RNA,.., slot: counts,data,scale.data, meta.data du lieu qc cho tung cell
#- counts → du lieu tho (raw UMI counts).
#data → du lieu da normalize (log-normalized).
#scale.data → du lieu da scale (mean = 0, variance = 1 cho moi gene).
nsclc.seurat.obj[["RNA"]]@layers$counts
nsclc.seurat.obj# thong tin assay,layer ...
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)#heatmap PC1

#tim elbow plot
ElbowPlot(nsclc.seurat.obj)# 15
#Find neighbors
nsclc.seurat.obj<- FindNeighbors(nsclc.seurat.obj, dims = 1:15)

# Clustering
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(nsclc.seurat.obj@meta.data)
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.1", label = TRUE,reduction = 'pca')
#gắn nhã cho từng cell
Idents(nsclc.seurat.obj)
Idents(nsclc.seurat.obj)<- "RNA_snn_res.0.1"
Idents(nsclc.seurat.obj)
#B-cell
FeaturePlot(nsclc.seurat.obj, features = c("CD19", "CD79B", "CD79A"))
#T-cell
FeaturePlot(nsclc.seurat.obj, features = c("CD3D", "CD3E", "CD3G"))#cụm 0
#NK-cell
FeaturePlot(nsclc.seurat.obj, features = c("NCAM1", "NKG7", "KLRD1"))
#Môncytes
FeaturePlot(nsclc.seurat.obj, features = c("CD14", "LYZ", "FCGR3A"))
#Enpithelial cell
FeaturePlot(nsclc.seurat.obj, features = c("EPCAM", "KRT18", "KRT19", "MUC1")) # cụm 2
#Enpedotheal cell
FeaturePlot(nsclc.seurat.obj, features = c("PECAM1", "VWF", "KDR"))
#Fibolast
FeaturePlot(nsclc.seurat.obj, features = c("COL1A1", "COL1A2", "DCN", "PDGFRA"))#cụm 1

nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)
DimPlot(nsclc.seurat.obj, reduction = "umap")
#DEA 
cluster1.marker = FindMarkers(nsclc.seurat.obj, ident.1 = 1, min.pct=0.25)
print(x = head(x = cluster1.marker, n = 5))
# find markers for every cluster compared to all remaining cells, report
# only the positive ones
nsclc.markers <- FindAllMarkers(nsclc.seurat.obj, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.2)
FeaturePlot(nsclc.seurat.obj, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), cols = c("grey", "blue"), reduction = "umap")
top10g <- nsclc.markers %>% group_by(cluster)
top10g <- top10g %>% top_n(10,avg_log2FC )
top10g
#FC = Fold Change: mức thay đổi biểu hiện gene giữa hai nhóm (ví dụ cluster A vs cluster B).

DoHeatmap(nsclc.seurat.obj, features = top10g$gene, label = TRUE)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
current.cluster.ids <- c(0, 1, 2, 3, 4)
new.cluster.ids <- c("NK cells","Enpithelial_cells","T cells","Monocyte","B cells")
names(new.cluster.ids) <- levels(nsclc.seurat.obj)
nsclc.seurat.obj <- RenameIdents(nsclc.seurat.obj, new.cluster.ids)
DimPlot(nsclc.seurat.obj,reduction = 'umap')


library(SingleR)
library(celldex)
ref <- HumanPrimaryCellAtlasData()
pred<- SingleR(test=GetAssayData(nsclc.seurat.obj,slot="data"),ref=ref, labels = ref$label.main)
nsclc.seurat.obj$SingleR.labels <- pred$labels
DimPlot(nsclc.seurat.obj, group.by = "SingleR.labels", label = TRUE)

