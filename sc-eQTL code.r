library(Seurat)
library(tidyverse)
library(dplyr)

#Data reading
data <- readRDS("data")
##QC
data <- subset(data, subset = nCount_RNA > qc_count & 
                     nFeature_RNA < qc_Feature_h & 
                     percent.mt < qc_mt & 
                     nFeature_RNA > qc_Feature_l)
##Scale data
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(data, features = rownames(data))
##Find cluster
data <- RunPCA(data, npcs = 50, verbose = FALSE)
data <- FindNeighbors(data, dims = 1:30)
data <- FindClusters(data, resolution = 0.5)

##findmarker
markers <- FindAllMarkers(data, logfc.threshold = 0.25, min.pct = 0.1, 
                          only.pos = TRUE, test.use = "wilcox")
##kegg
kegg = enrichKEGG(gene = marker, use_internal_data = T,
                         organism = "hsa",
                         pvalueCutoff = 0.9,
                         qvalueCutoff = 0.9)
##GO
GO<-enrichGO(gene = markers,
             OrgDb = "org.Hs.eg.db",
             keyType = "SYMBOL",
             ont = "ALL",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.01)    
##gsva
gsva <- gsva(expr=expression_matrix, 
            gset.idx.list=geneset, 
            kcdf="Gaussian" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
            verbose=T, 
            parallel.sz = parallel::detectCores())#调用所有核
##peer
model = PEER()
PEER_setPhenoMean(model,data2)
dim(PEER_getPhenoMean(model))
PEER_setNk(model,10)
PEER_getNk(model)
PEER_update(model)
factors = PEER_getX(model)
##fastqtl
for j in $(seq 1 200); do
    fastQTL --vcf gene_CRC.vcf.gz --bed SC.bed.gz --window 1e6 --cov cov.txt.gz --threshold 0.05 --chunk $j 200 --out resul_sc_${j}.txt.gz
  done
  zcat result_sc_*.txt.gz | gzip -c > qtl_sc_${i}.txt.gz
done
##coloc
    res <- coloc.abf(dataset1 = list(snp = data$SNP, pvalues = data$PVALUE, beta = data$BETA.x, varbeta = data$varbeta.x, type = "cc", N = 891168, s = 0.005, MAF = data$maf),
                     dataset2 = list(snp = data$SNP, pvalues = data$P, MAF = data$maf, type = "quant", N = 100, varbeta = data$varbeta.y), p1 = p1, p2 = p2, p12 = p12)
##twas
for j in $(seq 1 22);do
     Rscript FUSION.assoc_test.R \
     --sumstats ～/result.sumstats \
     --weights ./WEIGHTS/scTWAS.pos \
     --weights_dir ./WEIGHTS/ \
     --ref_ld_chr ./LDREF/1000G.EUR. \
     --chr $j \
     --out /output.txt
done

##ST
library(CellChat)
# Load data & basic processing
CRC <- Load10X_Spatial("filtered_feature_bc_matrix.h5") |> 
  AddMetaData(read.csv("clusters.csv", row.names=1))
# Cluster annotation
CRC$celltype <- ifelse(CRC$Cluster==1, "MyoFib/T", 
                      ifelse(CRC$Cluster==3, "MatrixFB", "Other"))
# CAF signature analysis
CRC <- AddModuleScore(CRC, list(
  iCAF = c("CXCL5","CXCL8","IL11"), 
  myCAF = c("SLC1A5","ACTG2")
))

##CellChat analysis
createCellChat(CRC) |> 
  computeCommunProb() |> 
  aggregateNet()
cellchat <- createCellChat(ninecell_T, DB=CellChatDB.human) |> 
  computeCommunProb() |> 
  computeCommunProbPathway() |> 
  aggregateNet()
library(Signac)
##Creat object 
atac_obj <- lapply(samples, \(s){
  CreateChromatinAssay(
    Read10X_h5(paste0(s,"/matrix.h5"))$Peaks,
    fragments = paste0(s,"/fragments.tsv.gz"),
    annotation = annotations)
}) |> merge() |> AddMetaData(cell_metadata)
##snATAC analysis
atac_core <- atac_obj |> 
  RunTFIDF() |> 
  RunSVD() |> 
  IntegrateLayers(method=HarmonyIntegration) |> 
  FindMultiModalNeighbors() |> 
  RunUMAP(reduction='wnn') |> 
  RunChromVAR(genome=BSgenome.Hsapiens.UCSC.hg38)