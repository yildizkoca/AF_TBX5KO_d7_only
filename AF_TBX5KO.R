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


#read raw matrices
d3wt_list<-list(d3wt2,d3wt2,d3wt3,d3wt4)
d=as.character("/project/spott/spott/mouse_AF")
dir_raw<-list()
d3wt_counts_raw<-list()
for (i in 1:length(d3wt_list)) {
  dir_raw[[i]]<-paste(d,d3wt_list[[i]],"outs","raw_feature_bc_matrix",sep="/")
  d3wt_counts_raw[[i]]<-as(object=Read10X(data.dir=dir_raw[[i]]), Class="CsparseMatrix")
}
names1<-c(rep("d3wt",4))
names2<-as.character(c(seq_len(4)))
names<-paste0(names1,names2)
names(d3wt_counts_raw)<-names
setwd("/project/spott/yildizkoca/AF_TBX5KO/objects")
saveRDS(d3wt_counts_raw,"AF_TBX5KO_d3wt_counts_raw.RDS")
d3wt_counts_raw<-readRDS("AF_TBX5KO_d3wt_counts_raw.RDS")
rm(d3wt_counts_raw)


d3ko_list<-list(d3ko1,d3ko2,d3ko3,d3ko4,d3ko5)
d=as.character("/project/spott/spott/mouse_AF")
dir_raw<-list()
d3ko_counts_raw<-list()
for (i in 1:length(d3ko_list)) {
  dir_raw[[i]]<-paste(d,d3ko_list[[i]],"outs","raw_feature_bc_matrix",sep="/")
  d3ko_counts_raw[[i]]<-as(object=Read10X(data.dir=dir_raw[[i]]), Class="CsparseMatrix")
}
names1<-c(rep("d3ko",5))
names2<-as.character(c(seq_len(5)))
names<-paste0(names1,names2)
names(d3ko_counts_raw)<-names
setwd("/project/spott/yildizkoca/AF_TBX5KO/objects")
saveRDS(d3ko_counts_raw,"AF_TBX5KO_d3ko_counts_raw.RDS")
#d3ko_counts_raw<-readRDS("AF_TBX5KO_d3ko_counts_raw.RDS")
rm(d3ko_counts_raw)

d7wt_list<-list(d7wt1,d7wt2,d7wt3,d7wt4,d7wt5,d7wt6,d7wt7,d7wt8,d7wt9,d7wt10)
d=as.character("/project/spott/spott/mouse_AF")
dir_raw<-list()
d7wt_counts_raw<-list()
for (i in 1:length(d7wt_list)) {
  dir_raw[[i]]<-paste(d,d7wt_list[[i]],"outs","raw_feature_bc_matrix",sep="/")
  d7wt_counts_raw[[i]]<-as(object=Read10X(data.dir=dir_raw[[i]]), Class="CsparseMatrix")
}
names1<-c(rep("d7wt",10))
names2<-as.character(c(seq_len(10)))
names<-paste0(names1,names2)
names(d7wt_counts_raw)<-names
setwd("/project/spott/yildizkoca/AF_TBX5KO/objects")
saveRDS(d7wt_counts_raw,"AF_TBX5KO_d7wt_counts_raw.RDS")
#d7wt_counts_raw<-readRDS("AF_TBX5KO_d7wt_counts_raw.RDS")
rm(d7wt_counts_raw)

d7ko_list<-list(d7ko1,d7ko2,d7ko3,d7ko4,d7ko5,d7ko6,d7ko7,d7ko8,d7ko9,d7ko10)
d=as.character("/project/spott/spott/mouse_AF")
dir_raw<-list()
d7ko_counts_raw<-list()
for (i in 1:length(d7ko_list)) {
  dir_raw[[i]]<-paste(d,d7ko_list[[i]],"outs","raw_feature_bc_matrix",sep="/")
  d7ko_counts_raw[[i]]<-as(object=Read10X(data.dir=dir_raw[[i]]), Class="CsparseMatrix")
}
names1<-c(rep("d7ko",10))
names2<-as.character(c(seq_len(10)))
names<-paste0(names1,names2)
names(d7ko_counts_raw)<-names
setwd("/project/spott/yildizkoca/AF_TBX5KO/objects")
saveRDS(d7ko_counts_raw,"AF_TBX5KO_d7ko_counts_raw.RDS")
#d7ko_counts_raw<-readRDS("AF_TBX5KO_d7ko_counts_raw.RDS")
rm(d7ko_counts_raw)

#read filtered matrices
d3wt_list<-list(d3wt2,d3wt2,d3wt3,d3wt4)
d=as.character("/project/spott/spott/mouse_AF")
dir_filtered<-list()
d3wt_counts_filtered<-list()
for (i in 1:length(d3wt_list)) {
  dir_filtered[[i]]<-paste(d,d3wt_list[[i]],"outs","filtered_feature_bc_matrix",sep="/")
  d3wt_counts_filtered[[i]]<-as(object=Read10X(data.dir=dir_filtered[[i]]), Class="CsparseMatrix")
}
names1<-c(rep("d3wt",4))
names2<-as.character(c(seq_len(4)))
names<-paste0(names1,names2)
names(d3wt_counts_filtered)<-names
setwd("/project/spott/yildizkoca/AF_TBX5KO/objects")
saveRDS(d3wt_counts_filtered,"AF_TBX5KO_d3wt_counts_filtered.RDS")
#d3wt_counts_filtered<-readRDS("AF_TBX5KO_d3wt_counts_filtered.RDS")
rm(d3wt_counts_filtered)

d3ko_list<-list(d3ko1,d3ko2,d3ko3,d3ko4,d3ko5)
d=as.character("/project/spott/spott/mouse_AF")
dir_filtered<-list()
d3ko_counts_filtered<-list()
for (i in 1:length(d3ko_list)) {
  dir_filtered[[i]]<-paste(d,d3ko_list[[i]],"outs","filtered_feature_bc_matrix",sep="/")
  d3ko_counts_filtered[[i]]<-as(object=Read10X(data.dir=dir_filtered[[i]]), Class="CsparseMatrix")
}
names1<-c(rep("d3ko",5))
names2<-as.character(c(seq_len(5)))
names<-paste0(names1,names2)
names(d3ko_counts_filtered)<-names
setwd("/project/spott/yildizkoca/AF_TBX5KO/objects")
saveRDS(d3ko_counts_filtered,"AF_TBX5KO_d3ko_counts_filtered.RDS")
#d3ko_counts_filtered<-readRDS("AF_TBX5KO_d3ko_counts_filtered.RDS")
rm(d3ko_counts_filtered)

d7wt_list<-list(d7wt1,d7wt2,d7wt3,d7wt4,d7wt5,d7wt6,d7wt7,d7wt8,d7wt9,d7wt10)
d=as.character("/project/spott/spott/mouse_AF")
dir_filtered<-list()
d7wt_counts_filtered<-list()
for (i in 1:length(d7wt_list)) {
  dir_filtered[[i]]<-paste(d,d7wt_list[[i]],"outs","filtered_feature_bc_matrix",sep="/")
  d7wt_counts_filtered[[i]]<-as(object=Read10X(data.dir=dir_filtered[[i]]), Class="CsparseMatrix")
}
names1<-c(rep("d7wt",10))
names2<-as.character(c(seq_len(10)))
names<-paste0(names1,names2)
names(d7wt_counts_filtered)<-names
setwd("/project/spott/yildizkoca/AF_TBX5KO/objects")
saveRDS(d7wt_counts_filtered,"AF_TBX5KO_d7wt_counts_filtered.RDS")
#d7wt_counts_filtered<-readRDS("AF_TBX5KO_d7wt_counts_filtered.RDS")
rm(d7wt_counts_filtered)

d7ko_list<-list(d7ko1,d7ko2,d7ko3,d7ko4,d7ko5,d7ko6,d7ko7,d7ko8,d7ko9,d7ko10)
d=as.character("/project/spott/spott/mouse_AF")
dir_filtered<-list()
d7ko_counts_filtered<-list()
for (i in 1:length(d7ko_list)) {
  dir_filtered[[i]]<-paste(d,d7ko_list[[i]],"outs","filtered_feature_bc_matrix",sep="/")
  d7ko_counts_filtered[[i]]<-as(object=Read10X(data.dir=dir_filtered[[i]]), Class="CsparseMatrix")
}
names1<-c(rep("d7ko",10))
names2<-as.character(c(seq_len(10)))
names<-paste0(names1,names2)
names(d7ko_counts_filtered)<-names
setwd("/project/spott/yildizkoca/AF_TBX5KO/objects")
saveRDS(d7ko_counts_filtered,"AF_TBX5KO_d7ko_counts_filtered.RDS")
#d7ko_counts_filtered<-readRDS("AF_TBX5KO_d7ko_counts_filtered.RDS")
rm(d7ko_counts_filtered)

setwd("/project/spott/yildizkoca/AF_TBX5KO/objects")
d3wt_counts_filtered<-readRDS("AF_TBX5KO_d3wt_counts_filtered.RDS")
d3ko_counts_filtered<-readRDS("AF_TBX5KO_d3ko_counts_filtered.RDS")
d7wt_counts_filtered<-readRDS("AF_TBX5KO_d7wt_counts_filtered.RDS")
d7ko_counts_filtered<-readRDS("AF_TBX5KO_d7ko_counts_filtered.RDS")

filtered_list<-list(d3wt_counts_filtered[2:4],d3ko_counts_filtered,d7wt_counts_filtered,d7ko_counts_filtered)
names(filtered_list)<-c("d3wt_filtered","d3ko_filtered","d7wt_filtered","d7ko_filtered")

source("/home/yildizkoca/scripts/scRNA-seq_functions.R")
source("/home/yildizkoca/scripts/scRNA-seq_integration_functions.R")


setwd("/project/spott/yildizkoca/AF_TBX5KO/AF_TBX5KO_standard")
for (i in 1:length(filtered_list)) {
  #apply standard clustering workflow
  #QC
  ObjList_filtered<-sc_QC(filtered_list[[i]], species="mouse",plot=TRUE, name=names(filtered_list)[i],version=NULL)
  #write filter stats
  ObjList_filtered_filtered<-sc_filter(ObjList_filtered, mt.perc=0.5, rb.perc=2.5, nfeatures = 200, plot=TRUE,  name=names(filtered_list)[i])
  #Normalize, scale, PCA
  ObjList_filtered<-sc_PCA(ObjList_filtered, nfeatures=2000, plotPCA=TRUE, plotJS=FALSE)
  #cluster, UMAP
  ObjList_filtered<-sc_cluster(ObjList_filtered, dims=1:20, resolution=0.5, plot=TRUE)
}


