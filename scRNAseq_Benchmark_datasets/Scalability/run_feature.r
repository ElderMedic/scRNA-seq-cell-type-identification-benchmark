run_CaSTLe<-function(DataPath,LabelsPath,CV_RDataPath, OutputDir, GeneOrderPath = NULL, NumGenes = NULL){
  "
  run CaSTLe
  Wrapper script to run CaSTLe on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  Data <- read.csv(DataPath,row.names = 1)
  Labels <- as.matrix(read.csv(LabelsPath))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                                CaSTLe                                     #
  #############################################################################
  library(igraph)
  library(xgboost)
  True_Labels_Castle <- list()
  Pred_Labels_Castle <- list()
  Training_Time_Castle <- list()
  Testing_Time_Castle <- list()
  
  BREAKS=c(-1, 0, 1, 6, Inf)
  nFeatures = 100
  
  for(i in c(1:n_folds)){
    # 1. Load datasets
    if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
      ds1 = Data[Train_Idx[[i]],as.vector(GenesOrder[c(1:NumGenes),i])+1]
      ds2 = Data[Test_Idx[[i]],as.vector(GenesOrder[c(1:NumGenes),i])+1]
    }
    else{
      ds1 = Data[Train_Idx[[i]],]
      ds2 = Data[Test_Idx[[i]],]
    }
    
    sourceCellTypes = as.factor(Labels[Train_Idx[[i]]])
    targetCellTypes = as.factor(Labels[Test_Idx[[i]]])
    
    start_time <- Sys.time()
    # 2. Unify sets, excluding low expressed genes
    source_n_cells_counts = apply(ds1, 2, function(x) { sum(x > 0) } )
    target_n_cells_counts = apply(ds2, 2, function(x) { sum(x > 0) } )
    common_genes = intersect( colnames(ds1)[source_n_cells_counts>10], 
                              colnames(ds2)[target_n_cells_counts>10])
    remove(source_n_cells_counts, target_n_cells_counts)
    ds1 = ds1[, colnames(ds1) %in% common_genes]
    ds2 = ds2[, colnames(ds2) %in% common_genes]
    ds = rbind(ds1[,common_genes], ds2[,common_genes])
    isSource = c(rep(TRUE,nrow(ds1)), rep(FALSE,nrow(ds2)))
    remove(ds1, ds2)
    
    # 3. Highest mean in both source and target
    topFeaturesAvg = colnames(ds)[order(apply(ds, 2, mean), decreasing = T)]
    end_time <- Sys.time()
    Training_Time_Castle[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    start_time <- Sys.time()
    # for each cell - what is the most probable classification?
    L = length(levels(sourceCellTypes))
    targetClassification = as.data.frame(matrix(rep(0,L*sum(!isSource)), nrow=L), row.names = levels(sourceCellTypes))
    
    for (cellType in levels(sourceCellTypes)) {
      
      inSourceCellType = as.factor(ifelse(sourceCellTypes == cellType, cellType, paste0("NOT",cellType)))
      
      # 4. Highest mutual information in source
      topFeaturesMi = names(sort(apply(ds[isSource,],2,function(x) { compare(cut(x,breaks=BREAKS),inSourceCellType,method = "nmi") }), decreasing = T))
      
      # 5. Top n genes that appear in both mi and avg
      selectedFeatures = union(head(topFeaturesAvg, nFeatures) , head(topFeaturesMi, nFeatures) )
      
      # 6. remove correlated features
      tmp = cor(ds[,selectedFeatures], method = "pearson")
      tmp[!lower.tri(tmp)] = 0
      selectedFeatures = selectedFeatures[apply(tmp,2,function(x) any(x < 0.9))]
      remove(tmp)
      
      # 7,8. Convert data from continous to binned dummy vars
      # break datasets to bins
      dsBins = apply(ds[, selectedFeatures], 2, cut, breaks= BREAKS)
      # use only bins with more than one value
      nUniq = apply(dsBins, 2, function(x) { length(unique(x)) })
      # convert to dummy vars
      ds0 = model.matrix(~ . , as.data.frame(dsBins[,nUniq>1]))
      remove(dsBins, nUniq)
      
      cat(paste0("<h2>Classifier for ",cellType,"</h2>"))
      
      inTypeSource = sourceCellTypes == cellType
      # 9. Classify
      xg=xgboost(data=ds0[isSource,] , 
                 label=inTypeSource,
                 objective="binary:logistic", 
                 eta=0.7 , nthread=1, nround=20, verbose=0,
                 gamma=0.001, max_depth=5, min_child_weight=10)
      
      # 10. Predict
      inTypeProb = predict(xg, ds0[!isSource, ])
      
      targetClassification[cellType,] = inTypeProb
    }
    end_time <- Sys.time()
    Testing_Time_Castle[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_Castle[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_Castle[i] <- list(rownames(targetClassification)[apply(targetClassification,2,which.max)])
  }
  True_Labels_Castle <- as.vector(unlist(True_Labels_Castle))
  Pred_Labels_Castle <- as.vector(unlist(Pred_Labels_Castle))
  Training_Time_Castle <- as.vector(unlist(Training_Time_Castle))
  Testing_Time_Castle <- as.vector(unlist(Testing_Time_Castle))
  
  dir.create(OutputDir,showWarnings = FALSE)

  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    write.csv(True_Labels_Castle,paste0(OutputDir,'/True_Labels_CaSTLe_',NumGenes,'.csv'),row.names = FALSE)
    write.csv(Pred_Labels_Castle,paste0(OutputDir,'/Pred_Labels_CaSTLe',NumGenes,'.csv'),row.names = FALSE)
    write.csv(Training_Time_Castle,paste0(OutputDir,'/Training_Time_CaSTLe',NumGenes,'.csv'),row.names = FALSE)
    write.csv(Testing_Time_Castle,paste0(OutputDir,'/Testing_Time_CaSTLe',NumGenes,'.csv'),row.names = FALSE)
  }
  else{
    write.csv(True_Labels_Castle,paste0(OutputDir,'/True_Labels_CaSTLe.csv'),row.names = FALSE)
    write.csv(Pred_Labels_Castle,paste0(OutputDir,'/Pred_Labels_CaSTLe.csv'),row.names = FALSE)
    write.csv(Training_Time_Castle,paste0(OutputDir,'/Training_Time_CaSTLe.csv'),row.names = FALSE)
    write.csv(Testing_Time_Castle,paste0(OutputDir,'/Testing_Time_CaSTLe.csv'),row.names = FALSE)
  }
  
}

run_singleCellNet<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,GeneOrderPath = NULL,NumGenes = NULL){
  "
  run singleCellNet
  Wrapper script to run singleCellNet on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  Data <- read.csv(DataPath,row.names = 1)
  colnames(Data) <- gsub('_','.',colnames(Data), fixed = TRUE)
  Labels <- as.matrix(read.csv(LabelsPath))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                              singleCellNet                                #
  #############################################################################
  library(singleCellNet)
  library(dplyr)
  True_Labels_singleCellNet <- list()
  Pred_Labels_singleCellNet <- list()
  Training_Time_singleCellNet <- list()
  Testing_Time_singleCellNet <- list()
  Data = t(as.matrix(Data))              # deals also with sparse matrix
  
  for(i in c(1:n_folds)){
    if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
      DataTrain <- Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]]
      DataTest <- Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]]
    }
    else{
      DataTrain <- Data[,Train_Idx[[i]]]
      DataTest <- Data[,Test_Idx[[i]]]
    }
    
    start_time <- Sys.time()
    cgenes2<-findClassyGenes(DataTrain, data.frame(Annotation = Labels[Train_Idx[[i]]]), "Annotation")
    cgenesA<-cgenes2[['cgenes']]
    grps<-cgenes2[['grps']]
    DataTrain<-as.matrix(DataTrain[cgenesA,])
    xpairs<-as.character(ptGetTop(DataTrain, grps,cgenes2[["cgenes_list"]]))
    pdTrain<-query_transform(DataTrain[cgenesA,], xpairs)
    rf<-sc_makeClassifier(pdTrain[xpairs,], genes=xpairs, groups=grps)
    end_time <- Sys.time()
    Training_Time_singleCellNet[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    start_time <- Sys.time()
    DataTest<-query_transform(DataTest[cgenesA,], xpairs)
    classRes <-rf_classPredict(rf, DataTest)
    end_time <- Sys.time()
    Testing_Time_singleCellNet[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_singleCellNet[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_singleCellNet[i] <- list((rownames(classRes)[apply(classRes,2,which.max)])[1:length(Test_Idx[[i]])])
  }
  True_Labels_singleCellNet <- as.vector(unlist(True_Labels_singleCellNet))
  Pred_Labels_singleCellNet <- as.vector(unlist(Pred_Labels_singleCellNet))
  Training_Time_singleCellNet <- as.vector(unlist(Training_Time_singleCellNet))
  Testing_Time_singleCellNet <- as.vector(unlist(Testing_Time_singleCellNet))
  
  
  dir.create(OutputDir,showWarnings = FALSE)
  setwd(OutputDir)
  
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    write.csv(True_Labels_singleCellNet,paste('singleCellNet_',NumGenes,'_True_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_singleCellNet,paste('singleCellNet_',NumGenes,'_Pred_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Training_Time_singleCellNet,paste('singleCellNet_',NumGenes,'_Training_Time.csv', sep = ''),row.names = FALSE)
    write.csv(Testing_Time_singleCellNet,paste('singleCellNet_',NumGenes,'_Testing_Time.csv', sep = ''),row.names = FALSE)
  }
  else{
    write.csv(True_Labels_singleCellNet,paste0(OutputDir,'/singleCellNet_True_Labels.csv'),row.names = FALSE)
    write.csv(Pred_Labels_singleCellNet,paste0(OutputDir,'/singleCellNet_Pred_Labels.csv'),row.names = FALSE)
    write.csv(Training_Time_singleCellNet,paste0(OutputDir,'/singleCellNet_Training_Time.csv'),row.names = FALSE)
    write.csv(Testing_Time_singleCellNet,paste0(OutputDir,'/singleCellNet_Testing_Time.csv'),row.names = FALSE)
  }
}

run_CHETAH<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,GeneOrderPath = NULL,NumGenes = NULL){
  "
  run CHETAH
  Wrapper script to run CHETAH on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  Data <- read.csv(DataPath,row.names = 1)
  Labels <- as.matrix(read.csv(LabelsPath))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                                CHETAH                                     #
  #############################################################################
  library(CHETAH)
  library(SingleCellExperiment)
  True_Labels_CHETAH <- list()
  Pred_Labels_CHETAH <- list()
  Total_Time_CHETAH <- list()
  Data = t(as.matrix(Data))
  
  for (i in c(1:n_folds)){
    if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
      sce <- SingleCellExperiment(assays = list(counts = Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]]), 
                                  colData = data.frame(celltypes = Labels[Train_Idx[[i]]]))
      
      sce_test <- SingleCellExperiment(assays = list(counts = Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]]), 
                                       colData = data.frame(celltypes = Labels[Test_Idx[[i]]]))
      start_time <- Sys.time()
      sce_test <- CHETAHclassifier(input = sce_test, ref_cells = sce, n_genes = NumGenes)
      end_time <- Sys.time()
    }
    else{
      sce <- SingleCellExperiment(assays = list(counts = Data[,Train_Idx[[i]]]), 
                                  colData = data.frame(celltypes = Labels[Train_Idx[[i]]]))
      
      sce_test <- SingleCellExperiment(assays = list(counts = Data[,Test_Idx[[i]]]), 
                                       colData = data.frame(celltypes = Labels[Test_Idx[[i]]]))
      start_time <- Sys.time()
      sce_test <- CHETAHclassifier(input = sce_test, ref_cells = sce)
      end_time <- Sys.time()
    }
    
    Total_Time_CHETAH[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_CHETAH[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_CHETAH[i] <- list(sce_test$celltype_CHETAH)
  }
  True_Labels_CHETAH <- as.vector(unlist(True_Labels_CHETAH))
  Pred_Labels_CHETAH <- as.vector(unlist(Pred_Labels_CHETAH))
  Total_Time_CHETAH <- as.vector(unlist(Total_Time_CHETAH))
  
  dir.create(OutputDir,showWarnings = FALSE)
  
  write.csv(True_Labels_CHETAH,paste0(OutputDir,'/CHETAH_true.csv'),row.names = FALSE)
  write.csv(Pred_Labels_CHETAH,paste0(OutputDir,'/CHETAH_pred.csv'),row.names = FALSE)
  write.csv(Total_Time_CHETAH,paste0(OutputDir,'/CHETAH_total_time.csv'),row.names = FALSE)
}



run_SingleR<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,GeneOrderPath = NULL,NumGenes = NULL){
  "
  run SingleR
  Wrapper script to run SingleR on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.

  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection,
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "

  Data <- read.csv(DataPath,row.names = 1)
  Labels <- as.matrix(read.csv(LabelsPath))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    GenesOrder = read.csv(GeneOrderPath)
  }

  #############################################################################
  #                               SingleR                                     #
  #############################################################################
  library(SingleR)
  library(Seurat)
  True_Labels_SingleR <- list()
  Pred_Labels_SingleR <- list()
  Total_Time_SingleR <- list()
  Data = t(as.matrix(Data))

  for (i in c(1:n_folds)){
    if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
      start_time <- Sys.time()
      singler = SingleR(method = "single", Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]],
                        Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]],
                        Labels[Train_Idx[[i]]])
      end_time <- Sys.time()
    }
    else{
      start_time <- Sys.time()
      singler = SingleR(method = "single", Data[,Test_Idx[[i]]], Data[,Train_Idx[[i]]], Labels[Train_Idx[[i]]])
      end_time <- Sys.time()
    }
    Total_Time_SingleR[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))

    True_Labels_SingleR[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_SingleR[i] <- list(as.vector(singler$labels))
  }
  True_Labels_SingleR <- as.vector(unlist(True_Labels_SingleR))
  Pred_Labels_SingleR <- as.vector(unlist(Pred_Labels_SingleR))
  Total_Time_SingleR <- as.vector(unlist(Total_Time_SingleR))

  dir.create(OutputDir,showWarnings = FALSE)
  
  write.csv(True_Labels_SingleR,paste0(OutputDir,'/SingleR_true.csv'),row.names = FALSE)
  write.csv(Pred_Labels_SingleR,paste0(OutputDir,'/SingleR_pred.csv'),row.names = FALSE)
  write.csv(Total_Time_SingleR,paste0(OutputDir,'/SingleR_total_time.csv'),row.names = FALSE)
}

CV_RDATA = c("D:\\scRNAseq_Benchmark-master\\scRNAseq_Benchmark_datasets\\Scalability\\Scalability_100.RData")
output = c("D:\\scRNAseq_Benchmark-master\\scRNAseq_Benchmark_datasets\\Scalability\\feature_selection\\results_25","D:\\scRNAseq_Benchmark-master\\scRNAseq_Benchmark_datasets\\Scalability\\feature_selection\\results_100","D:\\scRNAseq_Benchmark-master\\scRNAseq_Benchmark_datasets\\Scalability\\feature_selection\\results_500","D:\\scRNAseq_Benchmark-master\\scRNAseq_Benchmark_datasets\\Scalability\\feature_selection\\results_2000")
#output = c("D:\\scRNAseq_Benchmark-master\\scRNAseq_Benchmark_datasets\\Scalability\\feature_selection\\results_2000")

count <- 1 #0 for py 1 for R
features <- c(25,100,500,2000)
for (i in features){
    #run_CHETAH("D:\\scRNAseq_Benchmark-master\\scRNAseq_Benchmark_datasets\\Scalability\\Filtered_TM_data.csv", "D:\\scRNAseq_Benchmark-master\\scRNAseq_Benchmark_datasets\\Scalability\\Labels.csv", CV_RDATA[1], paste0(output[count],"CHETAH"), GeneOrderPath = "D:\\scRNAseq_Benchmark-master\\scRNAseq_Benchmark_datasets\\Scalability\\rank_genes_dropouts.csv", NumGenes = i)
    #run_CaSTLe("D:\\scRNAseq_Benchmark-master\\scRNAseq_Benchmark_datasets\\Scalability\\Filtered_TM_data.csv", "D:\\scRNAseq_Benchmark-master\\scRNAseq_Benchmark_datasets\\Scalability\\Labels.csv", CV_RDATA[1], paste0(output[count],"CaSTLe"), GeneOrderPath = "D:\\scRNAseq_Benchmark-master\\scRNAseq_Benchmark_datasets\\Scalability\\rank_genes_dropouts.csv", NumGenes = i)
    run_SingleR("D:\\scRNAseq_Benchmark-master\\scRNAseq_Benchmark_datasets\\Scalability\\Filtered_TM_data.csv", "D:\\scRNAseq_Benchmark-master\\scRNAseq_Benchmark_datasets\\Scalability\\Labels.csv", CV_RDATA[1], paste0(output[count],"SingleR"), GeneOrderPath = "D:\\scRNAseq_Benchmark-master\\scRNAseq_Benchmark_datasets\\Scalability\\rank_genes_dropouts.csv", NumGenes = i)
    #run_singleCellNet("D:\\scRNAseq_Benchmark-master\\scRNAseq_Benchmark_datasets\\Scalability\\Filtered_TM_data.csv", "D:\\scRNAseq_Benchmark-master\\scRNAseq_Benchmark_datasets\\Scalability\\Labels.csv", CV_RDATA[1], paste0(output[count],"singleCellNet"), GeneOrderPath = "D:\\scRNAseq_Benchmark-master\\scRNAseq_Benchmark_datasets\\Scalability\\rank_genes_dropouts.csv", NumGenes = i)
    print(paste0("success run",count))
    count <- count+1
}

