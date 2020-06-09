Cross_Validation <- function(LabelsPath, col_Index = 1,OutputDir){
 
  Labels <- as.matrix(read.csv(LabelsPath))
  Labels <- as.vector(Labels[,col_Index])
  
  Removed_classes <- !(table(Labels) > 10)
  Cells_to_Keep <- !(is.element(Labels,names(Removed_classes)[Removed_classes]))
  Labels <- Labels[Cells_to_Keep]
  
  random_seed<-c(1234,2345,3456,4567,5678,6789,7890,8912,9123,4321)
  times=10
  setwd(OutputDir)
  # Getting training and testing Folds
  library(rBayesianOptimization)
  n_folds = 10
  for (j in 1:times){
	Folds <- KFold(Labels,nfolds = n_folds, stratified = TRUE,seed=random_seed[j])
	Test_Folds <- c(n_folds:1)
	Train_Idx <- list()
	Test_Idx <- list()
	for (i in c(1:length(Folds))){
		Temp_Folds <- Folds
		Temp_Folds[Test_Folds[i]] <- NULL
		Train_Idx[i] <- list(unlist(Temp_Folds))
		Test_Idx[i] <- Folds[Test_Folds[i]]
		}
	remove(Temp_Folds,i,Folds)
	save(n_folds,Train_Idx,Test_Idx,col_Index,Cells_to_Keep,random_seed,times,file = paste("CV_10times10folds_",j,".RData",sep = ''))
	}
}