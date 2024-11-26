library("Seurat", lib="/home/yildizkoca/RLibrary")
library("ggplot2",lib="/home/yildizkoca/RLibrary")
library("SoupX",lib="/home/yildizkoca/RLibrary")
library("dplyr",lib="/home/yildizkoca/RLibrary")
library("fields", lib="/home/yildizkoca/RLibrary" )
library("ROCR", lib="/home/yildizkoca/RLibrary")
library("KernSmooth", lib="/home/yildizkoca/RLibrary")
library("Matrix", lib="/home/yildizkoca/RLibrary")
library("DoubletFinder",lib="/home/yildizkoca/RLibrary")


source("/home/yildizkoca/scripts/scRNA-seq_functions.R")
source("/home/yildizkoca/scripts/scRNA-seq_integration_functions.R")

setwd("/project/spott/yildizkoca/AF_TBX5KO/objects")

allouts<-readRDS("allouts_SoupX.RDS")

setwd("/project/spott/yildizkoca/AF_TBX5KO/AF_TBX5KO_SoupX/w_doublets")
ObjList<-sc_QC(allouts, species="mouse",plot=TRUE, name=NULL,version="SoupX")
rm(allouts)
ObjList_filtered<-sc_filter(ObjList, mt.perc=2.5, rb.perc=100, nfeatures = 200, plot=TRUE, name="AF_TBX5KO_SoupX")
ObjList_filtered<-sc_PCA(ObjList_filtered, nfeatures=2000, plotVF=TRUE, plotPCA=TRUE, plotJS=FALSE)
ObjList_filtered<-sc_cluster(ObjList_filtered, dims=1:20, resolution=0.5, plot=TRUE)
rm(ObjList)

#identify doublets by DoubletFinder

bcmvn<-list()
pK<-list()
for (i in 1:length(ObjList_filtered)) {
sweep.res.list<-paramSweep_v3(ObjList_filtered[[i]], PCs = 1:20, sct = FALSE)
sweep.stats<-summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn[[i]]<-find.pK(sweep.stats)
pK[[i]]<-bcmvn[[i]][which(bcmvn[[i]]["BCmetric"]==max(bcmvn[[i]]["BCmetric"])),"pK"]
ObjList_filtered[[i]]<-doubletFinder_v3(ObjList_filtered[[i]], PCs = 1:20, pN = 0.25, pK = as.numeric(as.vector(pK[[i]])), nExp =round(0.075*nrow(ObjList_filtered[[i]]@meta.data)))
}

for (i in 1:length(ObjList_filtered)) {
  ObjList_filtered[[i]][["isDoublet"]]<-ObjList_filtered[[i]]@meta.data[length(colnames(ObjList_filtered[[i]]@meta.data))]
  ObjList_filtered[[i]][["pANN"]]<-ObjList_filtered[[i]]@meta.data[length(colnames(ObjList_filtered[[i]]@meta.data))-2]
}

setwd("/project/spott/yildizkoca/AF_TBX5KO/objects")
saveRDS(ObjList_filtered[1:8], "d3_ObjList_filtered_SoupX_DF.RDS")
saveRDS(ObjList_filtered[9:16], "d7wt_ObjList_filtered_SoupX_DF.RDS")
saveRDS(ObjList_filtered[17:24], "d7ko_ObjList_filtered_SoupX_DF.RDS")

for (i in 1:length(ObjList_filtered)) {
  ObjList_filtered[[i]]<-subset(ObjList_filtered[[i]], subset= isDoublet=="Singlet")
}

setwd("/project/spott/yildizkoca/AF_TBX5KO/objects")
saveRDS(ObjList_filtered[1:8], "d3_ObjList_filtered_SoupX_db_rm.RDS")
saveRDS(ObjList_filtered[9:16], "d7wt_ObjList_filtered_SoupX_db_rm.RDS")
saveRDS(ObjList_filtered[17:24], "d7ko_ObjList_filtered_SoupX_db_rm.RDS")



