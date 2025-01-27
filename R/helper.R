
SCA_evaluate<-function(True_class, SCA_res, evaluating_class,fdr_cutoff=0.05,eval_mode='prec_recall'){
  #SCA_specific evaluation by adding fdr cutoff to set pred_clss
  #library(caret)
  true_cls<-True_class
  #pred_cls<-Pred_class
  pred_cls<-ifelse((SCA_res$qval<=fdr_cutoff & SCA_res$PCC>0),evaluating_class,'other')

  if(!is.element(evaluating_class,true_cls)){
    print('assign right class to evaluate as positive class')
  }else{
    true_cls.n<-ifelse(true_cls==evaluating_class,evaluating_class,'other')%>% factor(levels = c(evaluating_class,'other'))
    pred_cls.n<-ifelse(pred_cls==evaluating_class,evaluating_class,'other')%>% factor(levels = c(evaluating_class,'other'))
  }

  conf.m<-caret::confusionMatrix(pred_cls.n,true_cls.n, mode=eval_mode,positive=evaluating_class)
  perfm<-conf.m$byClass[c('Balanced Accuracy','Sensitivity','Specificity', 'Precision','Recall','F1')]
  #conf.m<-table(true_cls.n,pred_cls.n)
  #add auroc and auprc

  res<-list('confusion_matrix'=conf.m$table, 'performance'=perfm,'all'=conf.m)
}
get_SCA_perfm_cutoff<-function(True_class, SCA_res, evaluating_class,eval_mode='prec_recall'){

  cutoff.list<-c(0.05,0.01,0.001,0.0001,0.00001)
  eval.df<-data.frame(matrix(0,length(cutoff.list),6))
  colnames(eval.df)<-c('Balanced Accuracy','Sensitivity','Specificity','Precision','Recall','F1')
  rownames(eval.df)<-paste0('fdr',cutoff.list)
  confm<-list()
  all<-list()
  for(i in 1:length(cutoff.list)){
    res<-SCA_evaluate(True_class=True_class, SCA_res=SCA_res, evaluating_class=evaluating_class,fdr_cutoff=cutoff.list[i],eval_mode=eval_mode)
    eval.df[i,]<-res$performance
    confm[[paste0('fdr',cutoff.list[i])]]<-res$confusion_matrix
    all[[paste0('fdr',cutoff.list[i])]]<-res$all
  }

  SCA_q<-SCA_res$qval
  #SCA_cor<-SCA_res$PCC

  label_in<-ifelse(True_class==evaluating_class,1,0)
  response_in<-1-SCA_q
  if(length(unique(True_class))>1){
    pr_roc_auc<-AUPRC_AURPC(response_in,label_in)
  }else{
    pr_roc_auc<-NULL
  }


  out<-list('eval'=eval.df,'confusion_matrix'=confm,'pr_roc_auc'=pr_roc_auc)
  return(out)
}


AUPRC_AURPC<-function(response, label){
  #library(pROC)
  #library(PRROC)
  #library(caret)
  pROCroc<-pROC::roc(label,response)$auc %>% as.numeric()

  sc1=response[label==1]
  sc0=response[label==0]
  roc<-PRROC::roc.curve(scores.class0 = sc1, scores.class1 = sc0,curve=F )
  pr<-PRROC::pr.curve(scores.class0 = sc1, scores.class1 = sc0,curve=F )





  # roc$auc
  # pr$auc.integral
  #pROCroc is for sanity check
  out.df<-data.frame('pROC.AUROC'=pROCroc,'AUROC'=roc$auc,'AUPRC'=pr$auc.integral)
  return(out.df)
}


#summary function:
SCA_res_summary<-function(SCA_res){
  num_res<-length(SCA_res)
  refp_names<-names(SCA_res)
  SCA_res_summary<-list()
  for (i in 1:num_res){
    SCA_res.df<-SCA_res[[i]]$SCA_result
    lower_bound<-min(SCA_res.df$qval)
    threshold<-c(0.05)
    current<-0.01
    while(current>lower_bound){
      threshold<-c(threshold, current)
      current<-current/10

    }
    threshold<-c(threshold,lower_bound)
    thresholds <- sort(unique(threshold), decreasing = TRUE)
    summary.df<-data.frame('FDR_threshold'=thresholds,'SCA_cells'=rep(NA,length(thresholds)))
    for(j in 1:length(thresholds)){
      cells_sig<-sum(SCA_res.df$qval<=thresholds[j] & SCA_res.df$PCC>0)
      summary.df[j,2]<-cells_sig
    }
    SCA_res_summary[[refp_names[i]]]<-summary.df
  }
  return(SCA_res_summary)
}


SCA_res_assign<-function(SCA_res){
  #function to assign cell annotation as dataframe
  num_res<-length(SCA_res)
  refp_names<-names(SCA_res)
  SCA_res_assign<-list()
  for (i in 1:num_res){
    SCA_res.df<-SCA_res[[i]]$SCA_result
    lower_bound<-min(SCA_res.df$qval)
    threshold<-c(0.05)
    current<-0.01
    while(current>lower_bound){
      threshold<-c(threshold, current)
      current<-current/10

    }
    threshold<-c(threshold,lower_bound)
    thresholds <- sort(unique(threshold), decreasing = TRUE)

    assign.df<-data.frame('FDR_threshold'=thresholds,'SCA_cells'=rep(NA,length(thresholds)))
    for(j in 1:length(thresholds)){
      SCA_cells<-ifelse((SCA_res.df$qval<=thresholds[j] & SCA_res.df$PCC>0), 'SCA_cells','other')
      #summary.df[j,2]<-cells_sig
      SCA_res.df[[paste0('FDR_',thresholds[j])]]<-SCA_cells
    }
    SCA_res_assign[[refp_names[i]]]<-SCA_res.df
  }
  return(SCA_res_assign)


}
