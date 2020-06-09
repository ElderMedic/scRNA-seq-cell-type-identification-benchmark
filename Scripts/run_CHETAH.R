run_CHETAH<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,GeneOrderPath = NULL,NumGenes = NULL){


  #Data<-DataPath
  Data <- read.csv(DataPath,row.names = 1)
  Labels <- as.matrix(read.csv(LabelsPath))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    GenesOrder = read.csv(GeneOrderPath)
  }
  
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
  
  setwd(OutputDir)
  
  if (!is.null(GeneOrderPath) & !is.null (NumGenes)){
    write.csv(True_Labels_CHETAH,paste('CHETAH_',NumGenes,'_True_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_CHETAH,paste('CHETAH_',NumGenes,'_Pred_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Total_Time_CHETAH,paste('CHETAH_',NumGenes,'_Total_Time.csv', sep = ''),row.names = FALSE)
  }
  else{
    write.csv(True_Labels_CHETAH,'CHETAH_True_Labels.csv',row.names = FALSE)
    write.csv(Pred_Labels_CHETAH,'CHETAH_Pred_Labels.csv',row.names = FALSE)
    write.csv(Total_Time_CHETAH,'CHETAH_Total_Time.csv',row.names = FALSE)
  }
}
