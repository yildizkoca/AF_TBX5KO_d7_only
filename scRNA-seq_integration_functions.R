
getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))

sc_integrate<-function(ObjList, reduction="cca", nfeatures=2000, dims=1:30) {
  if(DefaultAssay(ObjList[[1]])=="RNA") {
    ObjList<-lapply(ObjList, NormalizeData)
    ObjList<-lapply(ObjList, FindVariableFeatures, selection.method = "vst", nfeatures = nfeatures)
  }
  features<-SelectIntegrationFeatures(object.list = ObjList)
  if (reduction=="cca") {
    immune.anchors<-FindIntegrationAnchors(object.list = ObjList, anchor.features = features)
    immune.combined<-IntegrateData(anchorset = immune.anchors)
  }
  if (reduction=="rpca") {
    ObjList<-lapply(ObjList, ScaleData, features=features)
    ObjList<-lapply(ObjList, RunPCA, features=features)
    immune.anchors<-FindIntegrationAnchors(object.list = ObjList, reduction="rpca", dims=dims)
    immune.combined<-IntegrateData(anchorset = immune.anchors, dims=dims)
  }
  DefaultAssay(immune.combined)<-"integrated"
  return(immune.combined)
}


sc_integrate_PCA<-function(immune.combined, plot=FALSE, name) {
  immune.combined<-ScaleData(immune.combined)
  immune.combined<-RunPCA(immune.combined, npcs = 30)
  
  if (plot == TRUE) {
    pca_plot1<-VizDimLoadings(immune.combined, dims = 1:2, reduction = "pca") + ggtitle(name) + theme_classic()
    filename1<-paste(name,"DimLoadings.pdf",sep="_")
    ggsave(filename1, plot=pca_plot1, width=8, height = 8)
    
    pca_plot2<-DimPlot(immune.combined, reduction = "pca", group.by = "orig.ident") + ggtitle(name) + theme_classic()
    filename2<-paste(name,"DimPlotPCA.pdf",sep="_")
    ggsave(filename2, plot=pca_plot2, width=8, height = 8)
    
    pca_plot4<-ElbowPlot(immune.combined) + ggtitle(name) + theme_classic()
    filename4<-paste(name,"ElbowPlot.pdf",sep="_")
    ggsave(filename4, plot=pca_plot4, width=8, height = 8)
  }
  return(immune.combined)
}

sc_integrate_cluster<-function(immune.combined, dims, idents= "orig.ident", resolution=0.5, plot=FALSE, name) {
  immune.combined<-FindNeighbors(immune.combined, reduction = "pca", dims = dims)
  immune.combined<-FindClusters(immune.combined, resolution = resolution)
  immune.combined<-RunUMAP(immune.combined, reduction = "pca", dims = dims)
  
  if(plot == TRUE) {
    p1<-DimPlot(immune.combined, reduction = "umap", pt.size = 0.01, cols = getPalette(length(unique(immune.combined$seurat_clusters)))) + ggtitle(name)+ theme_classic()
    ggsave(paste0(name,"_DimPlotUMAPbyClusters.pdf"), plot=p1, width=8, height=8)
    for (ident in idents) {
      p2<-DimPlot(immune.combined, reduction = "umap", group.by = ident,  pt.size = 0.01, cols = getPalette(length(unique(immune.combined@meta.data[,ident])))) +  ggtitle(name) + theme_classic()
      ggsave(paste0(name,"_DimPlotUMAPby",ident,".pdf"), plot=p2, width=8, height=8)
    }
  }
  return(immune.combined)
}


sc_integrate_cluster_markers<-function(immune.combined, heatmap=FALSE, dotPlot=FALSE, top_n=5, user_features=NULL, name) {
  
  all.cluster.markers<-FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(all.cluster.markers, file= paste0(name,".all.cluster.markers.csv"))
  top<-all.cluster.markers %>% group_by(cluster) %>% top_n(n = top_n, wt = avg_log2FC) 
  
  if (heatmap==TRUE) {
    #Heatmap cluster markers for all clusters
    all.cluster.markers %>%
      group_by(cluster) %>%
      slice_max(n = 2, order_by = avg_log2FC)
    all.cluster.markers %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_log2FC) -> top10
    heatmap<-DoHeatmap(immune.combined, features = top10$gene, size=0.5, combine=TRUE) +  theme(text = element_text(size=5)) +  ggtitle(name)
    ggsave(paste0(name,"_cluster_markers_heatmap.pdf"), plot=heatmap)
  }
  
  if (dotPlot==TRUE) {
    features<-unique(top$gene)
    dp<-DotPlot(object = immune.combined, features = features, assay="RNA") + theme_classic() + theme(axis.text.x =element_text(size=16,angle=90),axis.text.y =element_text(size=16) ) + ggtitle(name) 
    filename<-paste0(name,"_top_dotPlot.pdf")
    ggsave(filename, plot=dp, width=round((length(features)+8)/4), height=round((length(unique(immune.combined@meta.data$seurat_clusters))+8)/3), limitsize = F)
    if (!is.null(user_features)) {
      dp<-DotPlot(object = immune.combined, features = user_features, assay="RNA") + theme_classic() + theme(axis.text.x =element_text(size=16, angle=90), axis.text.y =element_text(size=16)) + ggtitle(name)
      filename<-paste0(name,"_user_features_dotPlot.pdf")
      ggsave(filename, plot=dp, width=round((length(features)+8)/4), height=round((length(unique(immune.combined[[i]]@meta.data$seurat_clusters))+8)/3), limitsize = F)    
    }
  }
  return(all.cluster.markers)
}


sc_integrate_decom_norm<-function(immune.combined, name, idents, write=F) {
  dir<-getwd()
  if (write) {
    dir.create(paste("sc_integrate_decom_norm",name,sep="_"))
    setwd(paste("sc_integrate_decom_norm",name,sep="_"))
  }
  idents1<-idents
  idents2<-idents
  i_norm_list<-list()
  for (ident1 in idents1) {
    for (ident2 in idents2) {
      if (ident1 != ident2) {
        i_norm<-immune.combined@meta.data %>% group_by(.data[[ident1]], .data[[ident2]]) %>% summarise(n=n())
        i<-immune.combined@meta.data %>% group_by(.data[[ident1]]) %>% summarise(n=n())
        i_norm$prop=-1
        ids<-unique(immune.combined@meta.data[[ident1]])
        for (id in ids) {
          i_norm[i_norm[,ident1]==id,"prop"]<-i_norm[i_norm[,ident1]==id,"n"]/as.numeric(i[which(i[,ident1]==id),"n"]) 
        }
        i_norm_list[[paste0(ident1,"_into_",ident2)]]<-as.data.frame(i_norm)
        if (write) {
          write.table(i_norm, paste(name,ident1,"into",ident2,"prop.txt", sep="_"))
          ncolors<-length(unique(immune.combined@meta.data[,ident2]))
          Idents(immune.combined)<-ident1
          g<-ggplot(i_norm, aes(x=factor(.data[[ident1]], level = levels(immune.combined)), y=prop, fill=.data[[ident2]])) + geom_bar(stat="identity", position="stack") + theme_classic() + scale_fill_manual(values=getPalette(ncolors)) + ggtitle(name)
          ggsave(paste(name,ident1,"to",ident2,"distribution.pdf",sep="_"), width=12, height=8) 
        }
      }
    }
  }
  setwd(dir)
  return(i_norm_list)
}


sc_proportion_test<-function(df, ident1="orig.ident", ident2="seurat_clusters", cond="genotype", groups=c("wt","ko"), name) {
  currdir<-getwd()
  if(!dir.exists(paste("prop_test",name,sep="_"))) {
    dir.create(paste("prop_test",name,sep="_"))
  }
  setwd(paste("prop_test",name,sep="_"))
  wilcox_test_results<-list()
  perm_test_results<-list()
  wilcox_test_results_pval<-list()
  perm_test_results_pval<-list()
  df0<-list()
  idents<-as.vector(unique(df[,ident2]))
  for (id in idents) {
    df0[[id]]<-df[which(df[,ident2]== id),]
    wilcox_test_results[[id]]<-wilcox.test(as.numeric(df0[[id]][which(df0[[id]][,cond]==groups[1]),"prop"]),as.numeric(df0[[id]][which(df0[[id]][,cond]==groups[2]),"prop"]), alternative = "two.sided")
    wilcox_test_results_pval[[id]]<- wilcox_test_results[[id]]$p.value
    perm_test_results[[id]]<-perm.t.test(as.numeric(df0[[id]][which(df0[[id]][,cond]==groups[1]),"prop"]),as.numeric(df0[[id]][which(df0[[id]][,cond]==groups[2]),"prop"]), alternative = "two.sided")
    perm_test_results_pval[[id]]<-perm_test_results[[id]]$perm.p.value
  }
  sink(paste(name, ident1, "into", ident2, "prop_wilcox_test.csv", sep="_"))
  print(wilcox_test_results)
  print(wilcox_test_results_pval)
  sink()
  sink(paste(name, ident1, "into", ident2, "prop_perm.t_test.csv", sep="_"))
  print(perm_test_results)
  print(perm_test_results_pval)
  sink()
  setwd(currdir)
  return(list("wilcox"=wilcox_test_results_pval,"permutation"=perm_test_results_pval))
}






sc_dotPlot<-function(ObjList, features, font.size=12, name=NULL, assay="RNA", split.by=NULL, group.by="seurat_clusters") {
  if ( class(ObjList) != "list") {
    ObjList<-list(ObjList)
    names(ObjList)<-name
  }
  assays<-list()
  if (is.null(assay)) {
    for ( i in 1:length(ObjList)) {
      assays[[i]]<-DefaultAssay(ObjList[[i]])
    }
  }
  else {
    for ( i in 1:length(ObjList)) {
      assays[[i]]<-assay
    }
  }
  cols<-getPalette(1)
  if (!is.null(split.by)) {
    for ( i in 1:length(ObjList)) {
      cols=getPalette(length(unique(ObjList[[i]]@meta.data[,split.by])))
      dp<-DotPlot(object = ObjList[[i]], features = features, assay=assays[[i]], group.by = group.by, split.by=split.by, cols=cols) + theme_classic() + theme(axis.text.x =element_text(angle=90, font.size=font.size)) + ggtitle(name)
      filename<-paste(name,split.by,"user_features",assay,"dotPlot.pdf", sep="_")
      ggsave(filename, plot=dp, width=(length(features)/8)+4, height=round((length(unique(ObjList[[i]]@meta.data$seurat_clusters))+8)/2))  
    }
  }
  else {
    for ( i in 1:length(ObjList)) {
      dp<-DotPlot(object = ObjList[[i]], features = features, assay=assays[[i]]) + theme_classic() + theme(axis.text.x =element_text(angle=90, font.size=font.size)) + ggtitle(name)
      filename<-paste(name,split.by,"user_features",assay, "dotPlot.pdf", sep="_")
      ggsave(filename, plot=dp, width=(length(features)/8)+4, height=round((length(unique(ObjList[[i]]@meta.data$seurat_clusters))+8)/3))  
    }
  }
}

sc_DimPlot_groups<-function(ObjList, reduction="umap", group.by="seurat_clusters", split.by=NULL, name=NULL, label=F, ncol=1, order, width=12, height=8) {
  if ( class(ObjList) != "list") {
    ObjList<-list(ObjList)
    names(ObjList)<-name
  }
  for (i in 1:length(ObjList)) {
    for (ident in group.by) {
      ncolors<-length(unique(ObjList[[i]]@meta.data[,ident]))
      colors<-getPalette(ncolors)
      g<-DimPlot(ObjList[[i]], reduction=reduction, group.by=ident, split.by = split.by, pt.size = 0.01, ncol=ncol, cols=colors, label = label, raster=F, label.size = 6) + ggtitle(names(ObjList)[i]) + theme_classic() + theme(legend.text = element_text(size=20), axis.title = element_text(size=14)) 
      ggsave(paste0(names(ObjList)[i], "_DimPlot", toupper(reduction), "by_", ident,"_", split.by, ".pdf"), g, width=width , height=height)
    }
  }
}



cell_number_per_label<-function(ObjList, idents, name, position=position_dodge(), order) {
  dir<-getwd()
  if (!dir.exists(paste("cell_number_per_label",name,sep="_"))) {
    dir.create(paste("cell_number_per_label",name,sep="_"))
  }
  setwd(paste("cell_number_per_label",name,sep="_"))
  df<-list()
  for (id in idents) {
    df[[id]]<-ObjList@meta.data %>% group_by(.data[[id]]) %>% summarise(n()) %>% as.data.frame()
    colnames(df[[id]])<-c(id, "cell_no")
    for ( i in 1:nrow(df[[id]])) {
      df[[id]][i,"prop"]<- df[[id]][i,"cell_no"]/sum(df[[id]][,"cell_no"])
    }
    write.table(df[id], paste0(paste(name, "cell_number_and_proportion_per",id, sep="_"),".txt"), quote=F, row.names = F, sep="\t")
    ncolors<-length(unique(ObjList@meta.data[,id]))
    g<-ggplot(df[[id]], aes(x=factor(.data[[id]], levels=order), y=cell_no, fill=.data[[id]])) + geom_bar(stat="identity", position=position) + theme_classic() + theme(axis.text.x=element_blank(),axis.title.x = element_blank(), axis.ticks.x = element_blank()) + scale_fill_manual(values=getPalette(ncolors), breaks=order) + ggtitle(name)
    ggsave(paste(name,id,"cell_no_distribution.pdf",sep="_"), width=8, height=8)
    g<-ggplot(df[[id]], aes(x=factor(.data[[id]], levels=order), y=prop, fill=.data[[id]])) + geom_bar(stat="identity", position=position) + theme(axis.text.x=element_blank()) + theme_classic() + theme(axis.text.x=element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank()) + scale_fill_manual(values=getPalette(ncolors), breaks=order) + ggtitle(name)
    ggsave(paste(name,id,"proportion_distribution.pdf",sep="_"), width=8, height=8)
  }
  setwd(dir)
  return(df)
}













