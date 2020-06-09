evaluate <- function(TrueLabelsPath, PredLabelsPath, time1, time2, Indices = NULL){
  
  library(aricode)
  library(sabre)
  
  true_lab <- unlist(read.csv(TrueLabelsPath))
  pred_lab <- unlist(read.csv(PredLabelsPath))
  
  Time1<- read.csv(time1)
  Time2<- read.csv(time2)
  
  if(time1==time2){
	total_time <- sum(Time1) 
  }
  else{
	total_time <- sum(Time1,Time2)
  }
  
  pred_lab_1<-pred_lab
  true_lab_1<-true_lab

  for (i in which(t(pred_lab_1) %in% c('unassigned','Unassigned','Unknown','rand','Node','ambiguous','unknown'))) {
	pred_lab_1[i]<-NA
	true_lab_1[i]<-NA
  }
  pred_lab_1<-na.omit(pred_lab_1)
  true_lab_1<-na.omit(true_lab_1)
  
  ari<-ARI(true_lab_1,pred_lab_1)
  ami<-AMI(true_lab_1,pred_lab_1)
  nvi<-NVI(true_lab_1,pred_lab_1)
  
  pearson_cor <- try(cor(true_lab_1,pred_lab_1,method = "pearson"),silent=FALSE)
  
  Vmeasure<-vmeasure(as.numeric(unlist(true_lab_1)),as.numeric(unlist(pred_lab_1)))
  
  if (! is.null(Indices)){
    true_lab <- true_lab[Indices]
    pred_lab <- pred_lab[Indices]
  }
  
  unique_true <- unlist(unique(true_lab))
  unique_pred <- unlist(unique(pred_lab))
  
  unique_all <- unique(c(unique_true,unique_pred))
  conf <- table(true_lab,pred_lab)
  pop_size <- rowSums(conf)
  
  pred_lab = gsub('Node..','Node',pred_lab)
  
  conf_F1 <- table(true_lab,pred_lab,exclude = c('unassigned','Unassigned','Unknown','rand','Node','ambiguous','unknown'))

  F1 <- vector()
  sum_acc <- 0
  
  for (i in c(1:length(unique_true))){
    findLabel = colnames(conf_F1) == row.names(conf_F1)[i]
    if(sum(findLabel)){
      prec <- conf_F1[i,findLabel] / colSums(conf_F1)[findLabel]
      rec <- conf_F1[i,findLabel] / rowSums(conf_F1)[i]
      if (prec == 0 || rec == 0){
        F1[i] = 0
      } else{
        F1[i] <- (2*prec*rec) / (prec + rec)
      }
      sum_acc <- sum_acc + conf_F1[i,findLabel]
    } else {
      F1[i] = 0
    }
  }
  
  pop_size <- pop_size[pop_size > 0]
  
  names(F1) <- names(pop_size)
  
  med_F1 <- median(F1)
  
  total <- length(pred_lab)
  num_unlab <- sum(pred_lab == 'unassigned') + sum(pred_lab == 'Unassigned') + sum(pred_lab == 'rand') + sum(pred_lab == 'Unknown') + sum(pred_lab == 'unknown') + sum(pred_lab == 'Node') + sum(pred_lab == 'ambiguous')
  per_unlab <- num_unlab / total
  
  acc <- sum_acc/sum(conf_F1)
  
  result <- list(Conf = conf, MedF1 = med_F1, F1 = F1, Acc = acc, PercUnl = per_unlab, PopSize = pop_size,ARI=ari,NVI=nvi,AMI=ami,Total_Time = total_time,pearson_corelation = pearson_cor,V_measure = Vmeasure[1])
  
  return(result)
}
