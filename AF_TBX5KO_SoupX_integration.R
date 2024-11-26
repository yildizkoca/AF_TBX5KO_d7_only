library("Seurat", lib="/home/yildizkoca/RLibrary")
library("ggplot2",lib="/home/yildizkoca/RLibrary")
library("SoupX",lib="/home/yildizkoca/RLibrary")
library("dplyr",lib="/home/yildizkoca/RLibrary")
library("fields", lib="/home/yildizkoca/RLibrary" )
library("ROCR", lib="/home/yildizkoca/RLibrary")
library("KernSmooth", lib="/home/yildizkoca/RLibrary")
library("Matrix", lib="/home/yildizkoca/RLibrary")
library("DoubletFinder",lib="/home/yildizkoca/RLibrary")


setwd("/project/spott/yildizkoca/AF_TBX5KO/objects")
ObjList_d3<-readRDS("d3_ObjList_filtered_SoupX_db_rm.RDS")
ObjList_d7wt<-readRDS("d7wt_ObjList_filtered_SoupX_db_rm.RDS")
ObjList_d7ko<-readRDS("d7ko_ObjList_filtered_SoupX_db_rm.RDS")

source("/home/yildizkoca/scripts/scRNA-seq_functions.R")
source("/home/yildizkoca/scripts/scRNA-seq_integration_functions.R")

int_d3wt<-sc_integrate(ObjList_d3[c("d3wt2_SoupX","d3wt3_SoupX","d3wt4_SoupX")], reduction="cca")
int_d3ko<-sc_integrate(ObjList_d3[c("d3ko1_SoupX","d3ko2_SoupX","d3ko3_SoupX","d3ko4_SoupX","d3ko5_SoupX")], reduction="cca")
int_d7wt<-sc_integrate(ObjList_d7wt[c("d7wt4_SoupX","d7wt5_SoupX","d7wt6_SoupX","d7wt7_SoupX","d7wt8_SoupX","d7wt9_SoupX","d7wt10_SoupX")], reduction="cca")
int_d7ko<-sc_integrate(ObjList_d7ko[c("d7ko4_SoupX","d7ko6_SoupX","d7ko7_SoupX","d7ko8_SoupX","d7ko9_SoupX","d7ko10_SoupX")], reduction="cca")

list<-list(int_d3wt,int_d3ko,int_d7wt,int_d7ko)
names(list)<-c("int_d3wt","int_d3ko","int_d7wt","int_d7ko")

saveRDS(list, "d3_d7_wt_ko_integrated_separately.RDS")
list<-readRDS("d3_d7_wt_ko_integrated_separately.RDS")

int_wt<-sc_integrate(list[c(1,3)], reduction="cca")
saveRDS(int_wt, "wt_integrated.RDS")
int_ko<-sc_integrate(list[c(2,4)], reduction="cca")
saveRDS(int_ko, "ko_integrated.RDS")
int_all_refilt<-sc_integrate(list(int_wt,int_ko),reduction="cca")

saveRDS(int_all_refilt, "all_refilt_wt_ko_ordered_integrated.RDS")


int_d3<-sc_integrate(list(int_d3wt,int_d3ko), reduction="cca")

setwd("/project/spott/yildizkoca/AF_TBX5KO/objects")
saveRDS(int_d3, "d3_integrated_SoupX_db_rm.RDS")


int_d7<-sc_integrate(list(int_d7wt,int_d7ko), reduction="cca" )

saveRDS(int_d7, "d7_integrated_SoupX_db_rm.RDS")

list<-list(int_d3,int_d7,int_all_refilt)
names(list)<-c("integrated_d3","integrated_d7","integrated_all_refilt")

setwd("/project/spott/yildizkoca/AF_TBX5KO/AF_TBX5KO_SoupX/integrated")
for ( i in 1:length(list)) {
  list[[i]]<-sc_integrate_PCA(list[[i]], plot=TRUE, name=names(list)[i])
  list[[i]]<-sc_integrate_cluster(list[[i]], dims=1:20, idents= "orig.ident", resolution=0.5, plot=TRUE, name=names(list)[i])
  sc_integrate_cluster_markers(list[[i]], heatmap=FALSE, dotPlot=TRUE, top_n=5, user_features=NULL, name=names(list)[i])
}

setwd("/project/spott/yildizkoca/AF_TBX5KO/objects")
saveRDS(list[[1]], "d3_integrated_SoupX_db_rm_processed.RDS")
saveRDS(list[[2]], "d7_integrated_SoupX_db_rm_processed.RDS")
saveRDS(list[[3]], "all_refilt_wt_ko_ordered_integrated_SoupX_db_rm_processed.RDS")
