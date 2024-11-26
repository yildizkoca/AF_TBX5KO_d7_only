
#FUNCTION
#list: a list of count matrices, type:list
#plot: whether to plot QC results, type:TRUE,FALSE
#version: version of the count matrices, type:character

getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))

sc_QC<-function(list, species, plot=F, name, version) {
  if ( class(list) != "list") {
    list<-list(list)
    names(list)="unnamed"
  }
  ObjList<-list()
  for (i in 1:length(list)) {
    o<-CreateSeuratObject(counts = list[[i]], project = names(list)[i], min.cells=3)
    ObjList[[i]]<-o
    if (is.null(version)) {
      names(ObjList)[i]<-names(list)[i]
    }
    else {
      names(ObjList)[i]<-paste(names(list)[i],version,sep="_")
    }
    
    if (species == "mouse") {
      ObjList[[i]][["percent.mt"]] <- PercentageFeatureSet(ObjList[[i]], pattern = "^mt-")
      ObjList[[i]][["percent.rb"]] <- PercentageFeatureSet(ObjList[[i]], pattern = "^Rp[sl]")
    }
    
    if (species == "human") {
      ObjList[[i]][["percent.mt"]] <- PercentageFeatureSet(ObjList[[i]], pattern = "^MT-")
      ObjList[[i]][["percent.rb"]] <- PercentageFeatureSet(ObjList[[i]], pattern = "^Rp[sl]")
    }
    
    
    if (plot == TRUE) {
      filename<-paste(names(ObjList)[i],"QC.pdf",sep="_")
      QCplot1<-VlnPlot(ObjList[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4) + theme_classic()
      ggsave(filename, plot=QCplot1, width=16, height = 8)
      
      QCplot2<-FeatureScatter(ObjList[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle(names(ObjList)[i]) + theme_classic()
      print(QCplot2)
      QCplot3<-FeatureScatter(ObjList[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle(names(ObjList)[i]) + theme_classic()
      print(QCplot3)
      QCplot4<-FeatureScatter(ObjList[[i]], feature1 = "nCount_RNA", feature2 = "percent.rb") + ggtitle(names(ObjList)[i]) + theme_classic()
      print(QCplot4)
      
      filename<-paste(names(ObjList)[i],"nCountRNA_vs_percentMT.pdf",sep="_")
      ggsave(filename, plot=QCplot2, width=8, height = 8)
      filename<-paste(names(ObjList)[i],"nCountRNA_vs_nFeatureRNA.pdf",sep="_")
      ggsave(filename, plot=QCplot3, width=8, height = 8)
      filename<-paste(names(ObjList)[i],"nCountRNA_vs_percentRb.pdf",sep="_")
      ggsave(filename, plot=QCplot4, width=8, height = 8)
      
      sink(file=paste0(name,"_allstats.csv"), append = FALSE)
      for (i in 1:length(ObjList)) {
        print(names(ObjList[[i]]))
        print(summary(ObjList[[i]]@meta.data))
      }
      sink()
    }
  }
  return(ObjList)
}


#FUNCTION
#list: a list of Seurat objects, type:list


sc_filter<-function(ObjList, nfeatures=200, mt.perc, rb.perc=100, doublets=vector(mode = "list", length = length(ObjList)), plot=FALSE, name) {
  if ( class(ObjList) != "list") {
    ObjList<-list(ObjList)
    names(ObjList)<-"unnamed"
  }
  #filter and visualize filtered cells
  ObjList_filtered<-list()
  for (i in 1:length(ObjList)) {
    ObjList[[i]][["mt.filter"]]<-ObjList[[i]]@meta.data[,"percent.mt"] < mt.perc
    ObjList[[i]][["rb.filter"]]<-ObjList[[i]]@meta.data[,"percent.rb"] < rb.perc
    ObjList[[i]][["nFeature.filter"]]<-ObjList[[i]]@meta.data[,"nFeature_RNA"] > nfeatures
    ObjList_filtered[[i]] <- subset(ObjList[[i]], subset = nFeature_RNA > nfeatures & percent.mt < mt.perc & percent.rb < rb.perc)
    if (plot==TRUE) {
      filename1<-paste(names(ObjList)[i],"_nFeature.filter_nCountRNA_vs_nFeatureRNA.pdf",sep="")
      g1<-ggplot(data=as.data.frame(ObjList[[i]]@meta.data), aes(x=nCount_RNA, y=nFeature_RNA, col=nFeature.filter)) + geom_point(aes(colour=nFeature.filter), size=0.01) + theme_classic()
      ggsave(filename=filename1, plot=g1, width = 8, height = 8)
      filename2<-paste(names(ObjList)[i],"_mt.filter_nCountRNA_vs_nFeatureRNA.pdf",sep="")
      g2<-ggplot(data=as.data.frame(ObjList[[i]]@meta.data), aes(x=nCount_RNA, y=nFeature_RNA, col=mt.filter)) + geom_point(aes(colour=mt.filter), size=0.01) + theme_classic()
      ggsave(filename=filename2, plot=g2, width = 8, height = 8)
      filename3<-paste(names(ObjList)[i],"_rb.filter_nCountRNA_vs_nFeatureRNA.pdf",sep="")
      g3<-ggplot(data=as.data.frame(ObjList[[i]]@meta.data), aes(x=nCount_RNA, y=nFeature_RNA, col=rb.filter)) + geom_point(aes(colour=rb.filter), size=0.01) + theme_classic()
      ggsave(filename=filename3, plot=g3, width = 8, height = 8)
    }
  }
  stats<-data.frame()
  for (i in 1:length(ObjList)) {
    stats[i,"sample"]<-names(ObjList)[i]
    stats[i,"features"]<-dim(ObjList[[i]])[1]
    stats[i,"cells"]<-dim(ObjList[[i]])[2]
    stats[i,"doublets"]<-length(doublets[[i]])
    stats[i,"doublets_perc"]<-length(doublets[[i]])/dim(ObjList[[i]])[2]
    stats[i,paste0("mt>",mt.perc)]<-length(which(ObjList[[i]][["mt.filter"]]==FALSE))
    stats[i,paste0("mt>",mt.perc,"_perc")]<-length(which(ObjList[[i]][["mt.filter"]]==FALSE))/dim(ObjList[[i]])[2]
    stats[i,paste0("rb>",rb.perc)]<-length(which( ObjList[[i]][["rb.filter"]]==FALSE))
    stats[i,paste0("rb>",rb.perc,"_perc")]<-length(which( ObjList[[i]][["rb.filter"]]==FALSE))/dim(ObjList[[i]])[2]
    stats[i,"features<200"]<-length(which(ObjList[[i]][["nFeature.filter"]]==FALSE))
  }
  names(ObjList_filtered)<-names(ObjList)
  write.csv(stats,paste0(name,"_stats.csv"))
  if( length(ObjList) == 1) {
    return(ObjList_filtered[[1]])
  }
  else {
    return(ObjList_filtered)
  }
}


#FUNCTION
#list: a list of Seurat objects, type:list
#plot: whether to plot PCA results, type:TRUE,FALSE
#plotJS: whether to run and plot JackStraw, type:TRUE,FALSE

sc_PCA<-function(ObjList, nfeatures=2000, plotVF=TRUE, plotPCA=TRUE, plotJS=FALSE) {
  if ( class(ObjList) != "list") {
    ObjList<-list(ObjList)
    names(ObjList)<-"unnamed"
    name="unnamed"
  }
  ObjList<-lapply(ObjList, NormalizeData)
  ObjList<-lapply(ObjList, FindVariableFeatures, selection.method = "vst", nfeatures = nfeatures)
  
  if (plotVF == TRUE) {
    for (i in 1:length(ObjList)) {
      vp<-VariableFeaturePlot(ObjList[[i]], pt.size=0.01) + theme_classic() + ggtitle(names(ObjList)[i])
      filename<-paste(names(ObjList)[i],"VarFeatures.pdf",sep="_")
      ggsave(filename, plot=vp, width=8, height = 8)
    }
  }
  
  ObjList<-mapply(ScaleData, ObjList, lapply(ObjList, rownames))
  ObjList<-mapply(RunPCA, ObjList, features=lapply(ObjList, VariableFeatures))
  
  if (plotPCA == TRUE) {
    for (i in 1:length(ObjList)) {
      pca_plot1<-VizDimLoadings(ObjList[[i]], dims = 1:2, reduction = "pca") + ggtitle(names(ObjList)[i]) + theme_classic()
      print(pca_plot1)
      filename1<-paste(names(ObjList)[i],"DimLoadings.pdf",sep="_")
      ggsave(filename1, plot=pca_plot1, width=8, height = 8)
      
      pca_plot2<-DimPlot(ObjList[[i]], reduction = "pca") + ggtitle(names(ObjList)[i]) + theme_classic()
      print(pca_plot2)
      filename2<-paste(names(ObjList)[i],"DimPlotPCA.pdf",sep="_")
      ggsave(filename2, plot=pca_plot2, width=8, height = 8)
      
      pca_plot4<-ElbowPlot(ObjList[[i]]) + ggtitle(names(ObjList)[i]) + theme_classic()
      print(pca_plot4)
      filename4<-paste(names(ObjList)[i],"ElbowPlot.pdf",sep="_")
      ggsave(filename4, plot=pca_plot4, width=8, height = 8)
    }
  }
  
  if(plotJS == TRUE) {
    for (i in 1:length(ObjList)) {
      ObjList[[i]]<-JackStraw(ObjList[[i]], num.replicate = 100)
      ObjList[[i]]<-ScoreJackStraw(ObjList[[i]], dims=1:20)
      pca_plot3<-JackStrawPlot(ObjList[[i]], dims=1:15) + theme_classic()
      print(pca_plot3)
      filename3<-paste(names(ObjList)[i],"JackStrawPlot.pdf",sep="_")
      ggsave(filename3, plot=pca_plot3, width=8, height = 8)
    }
  }
  if( length(ObjList) == 1) {
    return(ObjList[[1]])
  }
  else {
    return(ObjList)
  }
}


#FUNCTION
#list: a list of Seurat objects, type:list
#resolution
#plot: whether to plot clustering results, type:TRUE,FALSE
#heatmap: 

sc_cluster<-function(ObjList, dims, resolution=0.5, plot=FALSE) {
  if ( class(ObjList) != "list") {
    ObjList<-list(ObjList)
    names(ObjList)<-"unnamed"
    name="unnamed"
  }
  ObjList<-lapply(ObjList, FindNeighbors, dims = dims)
  ObjList<-lapply(ObjList, FindClusters, resolution = resolution)
  ObjList<-lapply(ObjList, RunUMAP, dims = dims)
  
  if(plot == TRUE) {
    for (i in 1:length(ObjList)) {
      umap<-DimPlot(ObjList[[i]], reduction = "umap", pt.size = 0.01) + ggtitle(names(ObjList)[i]) + theme_classic()
      filename<-paste(names(ObjList)[i],"DimPlotUMAP.pdf",sep="_")
      ggsave(filename, plot=umap, width=8, height = 8)
    }
  }
  if( length(ObjList) == 1) {
    return(ObjList[[1]])
  }
  else {
    return(ObjList)
  }
}


#FUNCTION
#list: a list of Seurat objects, type:list
#resolution
#plot: whether to plot clustering results, type:TRUE,FALSE
#heatmap: 

sc_cluster_markers<-function(ObjList, heatmap=FALSE, dotPlot=FALSE, top_n=5, user_features=NULL) {
  if ( class(ObjList) != "list") {
    ObjList<-list(ObjList)
    names(ObjList)<-"unnamed"
    name="unnamed"
  }
  top<-list()
  all.cluster.markers<-list()
  for (i in 1:length(ObjList)) {
    all.cluster.markers[[i]]<-FindAllMarkers(ObjList[[i]], only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    write.csv(all.cluster.markers[[i]], file= paste(names(ObjList)[i],"all.cluster.markers.csv", sep="_"))
    top[[i]]<-all.cluster.markers[[i]] %>% group_by(cluster) %>% top_n(n = top_n, wt = avg_log2FC) 
  }
  
  if (heatmap==TRUE) {
    for (i in 1:length(ObjList)) {
      all.cluster.markers[[i]] %>%
        group_by(cluster) %>%
        slice_max(n = 2, order_by = avg_log2FC)
      all.cluster.markers[[i]] %>%
        group_by(cluster) %>%
        top_n(n = 10, wt = avg_log2FC) -> top10
      heatmap<-DoHeatmap(ObjList[[i]], features = top10$gene, size=0.5, combine=TRUE) +  theme(text = element_text(size=5))
      print(heatmap)
      filename<-paste0(names(ObjList)[i],"_cluster_markers_heatmap.pdf")
      ggsave(filename, plot=heatmap, width=20, height = 15)
    }
  }
  if (dotPlot==TRUE) {
    features=list()
    for (i in 1:length(ObjList)) {
      features[[i]]<-unique(top[[i]]$gene)
      dp<-DotPlot(object = ObjList[[i]], features = features[[i]]) + theme(axis.text.x =element_text(angle=90)) + theme_classic()
      filename<-paste0(names(ObjList)[i],"_top_dotPlot.pdf")
      ggsave(filename, plot=dp, width=24, height=8)
      if (!is.null(user_features)) {
        dp<-DotPlot(object = ObjList[[i]], features = user_features) + theme(axis.text.x =element_text(angle=90)) + theme_classic()
        filename<-paste0(names(ObjList)[i],"_user_features_dotPlot.pdf")
        ggsave(filename, plot=dp, width=24, height=8)    
      }
    }
  }
  return(all.cluster.markers)
}

sc_mark_doublets<-function(ObjList, doubletList) {
  if ( class(ObjList) != "list") {
    ObjList<-list(ObjList)
    names(ObjList)<-"unnamed"
    name="unnamed"
  }
  for (i in 1:length(ObjList)) {
    db<-as.data.frame(rownames(ObjList[[i]]@meta.data))
    db[,"doublet"]<-"non-doublet"
    rownames(db)<-db[,1]
    db[intersect(doublets[[i]],colnames(ObjList[[i]])),"doublet"]<-"doublet"
    ObjList[[i]]$doublet<-db[,"doublet"]
  }
  for (i in 1:length(ObjList)) {
    umap1<-DimPlot(ObjList[[i]], reduction = "umap", pt.size = 0.01) + ggtitle(names(ObjList)[i]) + theme_classic()
    umap2<-DimPlot(ObjList[[i]], reduction = "umap", group.by = "doublet", pt.size = 0.01) + ggtitle(names(ObjList)[i]) + theme_classic()
    umap<-ggarrange(umap1,umap2)
    print(umap)
    filename<-paste(names(ObjList)[i],"_DimPlotUMAPbyDoublets.pdf",sep="")
    ggsave(filename, plot=umap, width=16, height=8)
  }
}

#identify doublets by DoubletFinder
sc_DF<-function(ObjList, dims, rate=0.075) {
  if ( class(ObjList) != "list") {
    ObjList<-list(ObjList)
    names(ObjList)<-"unnamed"
    name="unnamed"
  }
  bcmvn<-list()
  pK<-list()
  for (i in 1:length(ObjList)) {
    sweep.res.list<-paramSweep_v3(ObjList[[i]], PCs = dims, sct = FALSE)
    sweep.stats<-summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn[[i]]<-find.pK(sweep.stats)
    pK[[i]]<-bcmvn[[i]][which(bcmvn[[i]]["BCmetric"]==max(bcmvn[[i]]["BCmetric"])),"pK"]
    ObjList[[i]]<-doubletFinder_v3(ObjList[[i]], PCs = 1:12, pN = 0.25, pK = as.numeric(as.vector(pK[[i]])), nExp =round(rate*nrow(ObjList[[i]]@meta.data)))
    ObjList_filtered[[i]][["isDoublet"]]<-ObjList_filtered[[i]]@meta.data[length(colnames(ObjList_filtered[[i]]@meta.data))]
    ObjList_filtered[[i]][["pANN"]]<-ObjList_filtered[[i]]@meta.data[length(colnames(ObjList_filtered[[i]]@meta.data))-2]
  }
  return(ObjList)
}

#remove ambient RNA contamination by SoupX
sc_SoupX<-function(ObjList, dir_filtered, dir_raw) {
  if ( class(ObjList) != "list") {
    ObjList<-list(ObjList)
    names(ObjList)<-"unnamed"
    name="unnamed"
  }
  toc<-list()
  tod<-list()
  sc<-list()
  dd<-list()
  out<-list()
  for (i in 1:length(ObjList)) {
    toc[[i]]<-Read10X(file.path(dir_filtered[[i]]))
    tod[[i]]<-Read10X(file.path(dir_raw[[i]]))
    sc[[i]]<-SoupChannel(tod[[i]],toc[[i]])
    dd[[i]]<-ObjList[[i]]@meta.data[colnames(sc[[i]]$toc),] #retrieve metadata
    dd[[i]][colnames(sc[[i]]$toc),c("UMAP_1","UMAP_2")]<-Embeddings(ObjList[[i]]$umap)[colnames(sc[[i]]$toc),] #retrieve umap coordinates
    sc[[i]]<-setClusters(sc[[i]], setNames(dd[[i]]$seurat_clusters, rownames(dd[[i]]))) #transfer metadata
    sc[[i]]<-setDR(sc[[i]],dd[[i]][colnames(sc[[i]]$toc),c("UMAP_1","UMAP_2")]) #transfer umap coordinates
    sc[[i]]<-autoEstCont(sc[[i]]) #estimate contamination fraction
    out[[i]]<-adjustCounts(sc[[i]]) #get adjusted counts
  }
  names(sc)<-names(ObjList)
  names(dd)<-names(ObjList)
  names(out)<-names(ObjList)
  if( length(out) == 1) {
    return(out[[1]])
  }
  else {
    return(out)
  }
}



sc_DimPlot_metrics<-function(ObjList, metrics, name=NULL) {
  if ( class(ObjList) != "list") {
    ObjList<-list(ObjList)
    names(ObjList)<-"unnamed"
  }
  for (i in 1:length(ObjList)) {
    for (metric in metrics) {
      Idents(ObjList[[i]])<-metric
      g<-DimPlot(ObjList[[i]], reduction="umap", raster=F)
      df=as.data.frame(apply(g$data,2,as.numeric))
      g<-ggplot(df, aes(x=UMAP_1, y=UMAP_2, col=ident)) + geom_point(size=0.01) + scale_color_gradient(low="blue", high="red") + theme_classic() + ggtitle(paste(name,metric,sep="_")) 
      ggsave(paste0(name,"_DimPlotUMAPby_",metric,".pdf"), g, width=12 , height=8)
    }
  }
}

sc_jitterbox_metrics<-function(ObjList, metrics, split.by="seurat_clusters", name=NULL) {
  if ( class(ObjList) != "list") {
    ObjList<-list(ObjList)
    names(ObjList)<-name
  }
  for ( i in 1:length(ObjList)) {
    for (metric in metrics) {
      g<-ggplot(ObjList[[i]]@meta.data,aes(x=.data[[split.by]], y=.data[[metric]])) + geom_jitter(shape=16, position=position_jitter(0.2), cex=0.1) + geom_boxplot() + theme_classic() + ggtitle(paste(names(ObjList)[i],metric,sep="_"))
      ggsave(paste(names(ObjList)[i],"JitterPlot",metric,split.by,".pdf", sep="_"), g, width=12 , height=8)
    }
  }
}

sc_merge<-function(ObjList) {
  merged<-merge(ObjList[[1]], y=c(ObjList[2:length(ObjList)]), add.cell.ids = as.vector(names(ObjList)))
  return(merged)
}


#generate pseudobulk counts
generate.pseudobulk<-function(object, labels, assay="RNA", slot="counts") {
  flist<-list()
  for (i in labels) { 
    flist[[i]] <- unique(object@meta.data[,i]) 
  }
  meta<-expand.grid(flist, stringsAsFactors = FALSE)
  colnames(meta)<-labels
  rownames(meta) <- apply(meta, 1, function(x) paste0(x, collapse = '.'))
  n<-nrow(meta)
  
  #create empty matrix
  out<-matrix(nrow=dim(object[["RNA"]])[1], ncol=n, data=0)
  rownames(out)<-rownames(object[[assay]])
  colnames(out)<-rownames(meta)
  
  ncells <- c()
  ncounts <- c()
  total.cells <- dim(object[["RNA"]])[2]
  for (i in 1:n){
    cells <- 1:total.cells
    for (j in names(meta)) {
      keep  <- which(object@meta.data[[j]] == meta[i,j])
      cells <- cells[cells %in% keep]
    }
    ncells[i]<-length(cells)
    ncounts[i]<-sum(slot(object[[assay]], slot)[,cells])
    #some other thing to measure
    if (length(cells)==1) {
      out[,i] <- slot(object[[assay]], slot)[,cells]
    } 
    else {
      out[,i] <- Matrix::rowSums(slot(object[[assay]], slot)[,cells])
    }
  }
  meta$ncells<-ncells
  meta$ncounts<-ncounts/max(ncounts)
  return(list(counts=out, meta=meta))
}

#keep only pseudobulk samples with more than 'threshold' cells
filter.pseudobulk <- function(pseudobulk, threshold = 0) {
  w <- which(pseudobulk$meta$ncells > threshold)
  pseudobulk$counts <- pseudobulk$counts[,w]
  pseudobulk$meta <- pseudobulk$meta[w,]
  pseudobulk
}

# The required inputs for this function are:  seurat_input is a Seurat object;
# classification the metadata name of the cell type grouping variable in your
# Seurat object; model_formula is a formula, like above, maineffect is the
# variable name that groups the key fixed effect you want to test (used as the
# grouping for filtering genes by expression); and pseudo_factors is a vector
# with all the factors that will be used for pseudobulk grouping (ie, the terms
# that appear in the model).
de_genes <- function(seurat_input, classification, model_formula, maineffect, pseudo_factors){
  #how many levels are there of the chosen annotation?
  q <- nrow(unique(seurat_input[[classification]]))
  pseudo <- vector(mode = "list", length = q)
  names(pseudo) <- unique(seurat_input[[classification]]) %>% unlist() %>% unname()
  #generate pseudobulk of q+1 clusters
  for (k in names(pseudo)){
    cells.use <- colnames(seurat_input)[which(seurat_input[[classification]]==k)]
    pseudo[[eval(k)]] <- filter.pseudobulk(generate.pseudobulk(subset(seurat_input, cells = cells.use), labels = pseudo_factors))
  }
  print("Generated pseudobulk data. Initializing results table.")
  
  #set up output data structure
  results <- vector(mode = "list", length = q)
  names(results) <- names(pseudo)
  
  print("Results list initialized. Beginning DE testing by cluster.")
  print(results)
  #set up DE testing
  for (j in names(pseudo)){
    print(j)
    if (ncol(pseudo[[eval(j)]]$counts) < 3){
      next
    }
    # maybe not strictly necessary, but you can run into problems if there are
    # too few samples to actually work with.
    d <- DGEList(pseudo[[eval(j)]]$counts)
    d$samples <- cbind(d$samples[,c("lib.size","norm.factors")], pseudo[[eval(j)]]$meta)
    keepgenes <- filterByExpr(d$counts, group = d$samples[[maineffect]])
    d <- d[keepgenes,]
    d <- calcNormFactors(d, method = "TMM")
    v <- voomWithDreamWeights(d, model_formula, d$samples, plot=FALSE)
    modelfit <- dream(exprObj = v, formula = model_formula, data = d$samples, quiet = TRUE, suppressWarnings = TRUE)
    modelfit<-eBayes(modelfit)
    
    print(paste0("Tested cluster ", j, " for DE genes"))
    results[[eval(j)]] <- modelfit
    
    rm(d, keepgenes, v, modelfit)
    print(paste0("Finished DE testing of cluster ", j))
  }
  return(results)
}

# The resulting object is a list with a modelfit for each cell type.  You can
# get information about the test from this (e.g., the design matrix) as well as
# query it for DE results.




