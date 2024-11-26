library("Seurat", lib="/home/yildizkoca/RLibrary")
library("ggplot2",lib="/home/yildizkoca/RLibrary")
library("SoupX",lib="/home/yildizkoca/RLibrary")
library("dplyr",lib="/home/yildizkoca/RLibrary")


#day3-wt
#d3wt1
d3wt2<-"SP3-2-WT-A" #low fraction reads in cells
d3wt3<-"SP5-4-WT-A"
d3wt4<-"IM-SP-LdV-12S-LD-03"
#day3-ko
d3ko1<-"SP4-3-KO-A" #low fraction reads in cells
d3ko2<-"SP6-5-KO-A"
d3ko3<-"SP7-1-KO-A"
d3ko4<-"SP8-6-KO-A"
d3ko5<-"IM-SP-LdV-12S-LD-04"
#day7-wt
d7wt1<-"SP2-12-WT-A" #low median gene number #REMOVE
d7wt2<-"IM-SP-LdV-12S-LD-05" #low fraction reads in cells #REMOVE
d7wt3<-"IM-SP-LdV-12S-LD-07" #low cell number
d7wt4<-"IM-SP-LdV-12S-LD-09"
d7wt5<-"IM-SP-LdV-12S-LD-11"
d7wt6<-"SP-10X-18S-LdV-20221004-MHA-WT" #low fraction reads in cells, high ambRNA in SoupX
d7wt7<-"SP-10X-18S-LdV-20221004-MHC-WT"
d7wt8<-"SP-10X-18S-LdV-20221004-MHE-WT"
d7wt9<-"SP-10X-18S-LdV-20221004-MHG-WT"
d7wt10<-"SP-10X-18S-LdV-20221004-MHI-WT" #low fraction reads in cells
#day7-ko
d7ko1<-"SP1-11-KO-A" #low median gene number #REMOVE
d7ko2<-"IM-SP-LdV-12S-LD-06" #low fraction reads in cells, high antisense map #REMOVE
d7ko3<-"IM-SP-LdV-12S-LD-08" #low cell number
d7ko4<-"IM-SP-LdV-12S-LD-10"
d7ko5<-"IM-SP-LdV-12S-LD-12"
d7ko6<-"SP-10X-18S-LdV-20221004-MHB-KO" #low fraction reads in cells
d7ko7<-"SP-10X-18S-LdV-20221004-MHD-KO"
d7ko8<-"SP-10X-18S-LdV-20221004-MHF-KO"
d7ko9<-"SP-10X-18S-LdV-20221004-MHH-KO"
d7ko10<-"SP-10X-18S-LdV-20221004-MHJ-KO" #low fraction reads in cells



setwd("/project/spott/yildizkoca/AF_TBX5KO/objects")
d3wt_counts_filtered<-readRDS("AF_TBX5KO_d3wt_counts_filtered.RDS")
d3ko_counts_filtered<-readRDS("AF_TBX5KO_d3ko_counts_filtered.RDS")
d7wt_counts_filtered<-readRDS("AF_TBX5KO_d7wt_counts_filtered.RDS")
d7ko_counts_filtered<-readRDS("AF_TBX5KO_d7ko_counts_filtered.RDS")


source("/home/yildizkoca/scripts/scRNA-seq_functions.R")
source("/home/yildizkoca/scripts/scRNA-seq_integration_functions.R")


#d3wt1 mistake, #REMOVE
#d7wt1<-"SP2-12-WT-A" #low median gene number #REMOVE
#d7wt2<-"IM-SP-LdV-12S-LD-05" #low fraction reads in cells #REMOVE
#d7ko1<-"SP1-11-KO-A" #low median gene number #REMOVE
#d7ko2<-"IM-SP-LdV-12S-LD-06" #low fraction reads in cells, high antisense map #REMOVE

filtered_list<-c(d3wt_counts_filtered[2:4],d3ko_counts_filtered,d7wt_counts_filtered[3:10],d7ko_counts_filtered[3:10])
#Create Seurat objects
ObjList_filtered<-sc_QC(filtered_list, species="mouse",plot=FALSE, name=NULL,version=NULL)
#Normalize, scale, PCA
ObjList_filtered<-sc_PCA(ObjList_filtered, nfeatures=2000, plotVF=FALSE, plotJS=FALSE, plotPCA=FALSE)
#cluster, UMAP
ObjList_filtered<-sc_cluster(ObjList_filtered, dims=1:20, resolution=0.5, plot=FALSE)

d3wt_list<-list(d3wt2,d3wt3,d3wt4)
d3ko_list<-list(d3ko1,d3ko2,d3ko3,d3ko4,d3ko5)
d7wt_list<-list(d7wt3,d7wt4,d7wt5,d7wt6,d7wt7,d7wt8,d7wt9,d7wt10)
d7ko_list<-list(d7ko3,d7ko4,d7ko5,d7ko6,d7ko7,d7ko8,d7ko9,d7ko10)

list<-c(d3wt_list,d3ko_list,d7wt_list,d7ko_list)
d=as.character("/project/spott/spott/mouse_AF")
dir_filtered<-list()
dir_raw<-list()
for (i in 1:length(list)) {
  dir_filtered[[i]]<-paste(d,list[[i]],"outs","filtered_feature_bc_matrix",sep="/")
  dir_raw[[i]]<-paste(d,list[[i]],"outs","raw_feature_bc_matrix",sep="/")
}

setwd("/project/spott/yildizkoca/AF_TBX5KO/objects")
#remove ambient RNA contamination
sc<-list()
dd<-list()
out<-list()
for (i in 1:length(ObjList_filtered)) {
  toc<-Read10X(file.path(dir_filtered[[i]]))
  tod<-Read10X(file.path(dir_raw[[i]]))
  sc[[i]]<-SoupChannel(tod,toc)
  dd[[i]]<-ObjList_filtered[[i]]@meta.data[colnames(sc[[i]]$toc),] #retrieve metadata
  dd[[i]][colnames(sc[[i]]$toc),c("UMAP_1","UMAP_2")]<-Embeddings(ObjList_filtered[[i]]$umap)[colnames(sc[[i]]$toc),] #retrieve umap coordinates
  sc[[i]]<-setClusters(sc[[i]], setNames(dd[[i]]$seurat_clusters, rownames(dd[[i]]))) #transfer metadata
  sc[[i]]<-setDR(sc[[i]],dd[[i]][colnames(sc[[i]]$toc),c("UMAP_1","UMAP_2")]) #transfer umap coordinates
  sc[[i]]<-autoEstCont(sc[[i]]) #estimate contamination fraction
  out[[i]]<-adjustCounts(sc[[i]]) #get adjusted counts
}
names(sc)<-names(ObjList_filtered)
names(dd)<-names(ObjList_filtered)
names(out)<-names(ObjList_filtered)


setwd("/project/spott/yildizkoca/AF_TBX5KO/objects")
saveRDS(sc,"allscs_SoupX.RDS")
saveRDS(dd,"alldds_SoupX.RDS")
saveRDS(out,"allouts_SoupX.RDS")
#sc<-readRDS("d3wt_sc_SoupX.RDS")
#out<-readRDS("d3wt_out_SoupX.RDS")
rm(sc)
rm(dd)
rm(ObjList_filtered)


