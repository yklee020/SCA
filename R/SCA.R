

#' extract_reference_profile
#'
#' @param dat.sce.filt
#' @param label
#' @param num_signatures
#' @param up_only
#' @param rank_slot
#' @param lfc_slot
#'
#' @return
#' @export
#'
#' @examples
extract_reference_profile<-function(dat.sce.filt, label, num_signatures=500, up_only=FALSE,rank_slot='rank.AUC',lfc_slot='median.logFC.cohen'){
  #rank_slot can be c(min.AUC,mean.AUC,median.AUC)
  #I chose to use ScoreMarker function and their AUC based ranking system - take median cohen's D lfc as realative expression values
  #put class into SCA_cls
  #label<-'LT-HSC(CD38nCD45RAnCD90pCD49fp)'

  scoreMarker.refp<-scran::scoreMarkers(dat.sce.filt,colData(dat.sce.filt)$SCA_cls, assay.type='logcounts')
  refprofile<-as.data.frame(scoreMarker.refp[[label]]) %>% rownames_to_column(var = 'Gene')

  if(rank_slot!='rank.AUC'){
    print(paste0('rank by ',rank_slot))
    refprofile<-refprofile[order(refprofile[,rank_slot],decreasing = T), ]

  }else{
    print('rank by rank.AUC')
    refprofile<-refprofile[order(refprofile[,rank_slot],decreasing = F), ]
  }


  if(isTRUE(up_only)){
    Top.Bottom.n.SCA_cls<-head(refprofile,num_signatures)
  }else{
    Top.Bottom.n.SCA_cls<-rbind(head(refprofile,floor(num_signatures/2)),tail(refprofile,floor(num_signatures/2)))
  }


  #refprofile.out<-data.frame('Gene'=Top.Bottom.n.SCA_cls$Gene,'lfc'=Top.Bottom.n.SCA_cls$median.logFC.cohen)
  #'min.AUC' and 'min.logFC.cohen' can be used for cell of interest vs. all the other cell types
  refprofile.out<-data.frame('Gene'=Top.Bottom.n.SCA_cls$Gene,'lfc'=Top.Bottom.n.SCA_cls[,lfc_slot])
  return(refprofile.out)
}


prep_query_dat<-function(dat.sce,genelist.df,gene_to_use='hGene'){
  # select genes, reorder and subset query data with reference profile genes
  query_dat_DE_index<-c()
  for(i in genelist.df[[gene_to_use]]){
    ind<-which((rownames(dat.sce)) %in% i)
    query_dat_DE_index<-c(query_dat_DE_index,ind)
  }
  print(paste0('ref_profile length is:',dim(genelist.df)[1],
               ' and overlap with reference profile genes:',length(query_dat_DE_index)))

  dat.sce.srt<-dat.sce[query_dat_DE_index,]
  return(dat.sce.srt)
}

prep_refprofile<-function(dat.sce.srt,genelist.df,gene_to_use='hGene'){
  #prepare reference profile vector with the genes that only overlap with query data
  #and rearrange gene order to the query data so I can perform correlation analysis
  #query data needs to be subsetted and sorted using prep_query_dat function
  all_OL_True<-rownames(dat.sce.srt) %in% genelist.df[[gene_to_use]] %>% all()

  if(!all_OL_True){
    print('genes in the expression data are not all subset of reference profile. Try prep_query_dat function first.')
    stop()
  }else{
    #
    rownames(genelist.df)<-genelist.df[[gene_to_use]]
    ref_profile<-genelist.df[rownames(dat.sce.srt),]
    ref_profile<-ref_profile[,'lfc',drop=FALSE]

    #double check if referenec profile genes and query data genes were matched and ordered the same
    gene_ordered_True<-all.equal(rownames(ref_profile),rownames(dat.sce.srt))
    if (gene_ordered_True){
      print('reference profile genes and input expression data genes are matched and ordered')
    }else{
      print('warning!: reference profile genes and query data genes are NOT matched and ordered')
    }
  }
  return(ref_profile)
}


# cosine <- function() {
# }

SCA_obs_dist<-function(dat.sce.norm.srt,model_profile,
                       method=c('pearson','cosine','spearman'),dat.use=c('fccounts','logcounts')){
  #calculate observed correlation coefficient: reference profile vs. query data
  dat.counts<-assay(dat.sce.norm.srt,dat.use)

  #calculate correlation coefficients
  if(method =='pearson'){
    dat.CC<-apply(dat.counts, 2, function(x){cor(model_profile$lfc,x,method='pearson')})
  }else if(method =='cosine'){
    dat.CC<-apply(dat.counts, 2, function(x){cosine(model_profile$lfc,x)})
  }else if(method =='spearman'){
    dat.CC<-apply(dat.counts, 2, function(x){cor(model_profile$lfc,x,method='spearman')})
  }else{
    print('Invalid method')
    stop()
  }
  return(dat.CC)
}

SCA_random_dist<-function(dat.sce.norm,dat.sce.norm.srt,model_profile,nperm=1000,
                          method=c('pearson','cosine','spearman'),dat.use=c('fccounts','logcounts'),features.use=c('all','refprofile','across.sample.refp.gene.permutation')){
  #null dist of version01
  print(paste('number of permutations:',nperm,sep=''))

  #nsamp=1000
  ncell=dim(dat.sce.norm.srt)[2]
  N=dim(dat.sce.norm)[1]
  m=dim(dat.sce.norm.srt)[1]# number of genes in reference profile

  #make progress bar:
  pb <- txtProgressBar(min = 0, max = nperm, initial = 0)

  #PCC_random<-c()
  #place to save randomized correlation coefficient from all cells (column) for each iteration (row)
  PCC_random<-matrix(NA,nperm,ncell)
  for (itr in 1:nperm){

    #use matrix indexing to reduce the time
    if(features.use=='across.sample.refp.gene.permutation'){
      gmat<-t(replicate(m,sample.int(ncell,ncell)))
      row_indices <- matrix(rep(1:m, each = ncell), ncol = ncell, byrow = TRUE)

      # Use matrix indexing to select cells for each gene
      vec_ind<-row_indices + (gmat - 1) * m
      gmat_Gdat<- matrix(assay(dat.sce.norm.srt, dat.use)[vec_ind],nrow=m,byrow=F)
    }else{
      print('specify features methods for permutation process')
      stop()
    }


    if(method =='pearson'){
      datPCC_random<-apply(gmat_Gdat, 2, function(x){cor(model_profile$lfc,x,method='pearson')})
    }else if(method =='cosine'){
      datPCC_random<-apply(gmat_Gdat, 2, function(x){cosine(model_profile$lfc,x)})
    }else if(method =='spearman'){
      datPCC_random<-apply(gmat_Gdat, 2, function(x){cor(model_profile$lfc,x,method='spearman')})
    }else{
      print('Invalid method')
      stop()
    }

    PCC_random[itr,]<-datPCC_random

    #progress bar
    Sys.sleep(0.1)
    setTxtProgressBar(pb,itr)

  }
  merge.all<-as.vector(PCC_random)
  close(pb)
  return(merge.all)
}

calc_obs_null_dist<-function(obs_corr_coeff,null_corr_coeff){
  CC_obs_null<-list('stat' = obs_corr_coeff,'stat0'=null_corr_coeff)
  return(CC_obs_null)
}
calc_SCA_FDR<-function(PCC_obs_null){
  library(qvalue)
  #calculate positive PCC FDR
  pvalues_pos<- empPvals(stat = PCC_obs_null$stat, stat0 = PCC_obs_null$stat0)
  FDR_pos<-tryCatch({
    qvalue(p = pvalues_pos)
  },error=function(err){
    cat("An error occurred:", conditionMessage(err), "\n")
    print('trying pi0=1: BH method')
    qvalue(p = pvalues_pos,pi0 = 1)

  })

  FDR_pos.df<-data.frame('PCC'=PCC_obs_null$stat,'pval'=FDR_pos$pvalues,'qval'=FDR_pos$qvalues)
  SCA_FDR_out<-FDR_pos.df
  return(SCA_FDR_out)
}

#' SCA
#'
#' @param dat.sce.filt filtered SingleCellExperiment object
#' @param cor.method using spearman's correlation coefficient as default
#' @param ref_profile takes reference profile that either user provided or by extract_reference_profile function
#' @param num_perm
#' @param data.use
#' @param rand.features.use
#' @param ref_data_null_dist
#'
#' @return list containing 1)a dataframe with correlation coefficients (CC), p-values and FDR-q values, 2) list of a) observed CC distribution and b)null CC distribution
#' @export
#'
#' @examples
SCA<-function(dat.sce.filt,cor.method='spearman', ref_profile=NULL, num_perm=500,data.use=c('fccounts','logcounts'),
              rand.features.use='across.sample.refp.gene.permutation', ref_data_null_dist=NULL){
  #add variable for what genes to use in ref_profile
  #DEG extraction if ref_profile is NULL
  #add ref data and query data - make separate version?
  ##refdata,predefine class in colData, which gene slots to use,num genes to use (DEG) or
  #generate ref profile outside of the SCA function



  #match query data and background data
  if(!is.null(ref_data_null_dist)){
    print('find shared genes between query data and background data')
    shared_genes<-intersect(rownames(dat.sce.filt), rownames(ref_data_null_dist))
    dat.sce.filt<-dat.sce.filt[shared_genes,]
    ref_data_null_dist<-ref_data_null_dist[shared_genes,]
  }else{
    print('no external background data provided. Using query data for null distribution')
  }

  #prepare query data and reference profile
  dat.sce.filt.sort<-prep_query_dat(dat.sce=dat.sce.filt,genelist.df=ref_profile ,gene_to_use = 'Gene')
  ref_profile.n<-prep_refprofile(dat.sce.srt=dat.sce.filt.sort,genelist.df=ref_profile,gene_to_use='Gene')

  #generate observed correlation coefficient
  obs.cor<-SCA_obs_dist(dat.sce.norm.srt = dat.sce.filt.sort,model_profile = ref_profile.n,
                        method=cor.method,dat.use=data.use)
  if(!is.null(ref_data_null_dist)){
    #if reference data for null distribution is provided separately, prep ref_data_null_dist for null dist
    #ref_data_null_dist also should be sce object
    cat('\n')
    print('preprocess background data for calculating null distribution')
    ref_data_null_dist.filt<-prep_query_dat(dat.sce=ref_data_null_dist,genelist.df=ref_profile ,gene_to_use = 'Gene')
    ref_profile.ref_data<-prep_refprofile(dat.sce.srt=ref_data_null_dist.filt,genelist.df=ref_profile,gene_to_use='Gene')

    #generate null dist based on the reference profile
    null.cor<-SCA_random_dist(dat.sce.norm=ref_data_null_dist,dat.sce.norm.srt=ref_data_null_dist.filt,model_profile=ref_profile.ref_data,
                              nperm=num_perm,method=cor.method,dat.use=data.use, features.use = rand.features.use)
  }else{
    #if no reference data for null distribution is provided, generate null distribution from the query data

    null.cor<-SCA_random_dist(dat.sce.norm=dat.sce.filt,dat.sce.norm.srt=dat.sce.filt.sort,model_profile=ref_profile.n,
                              nperm=num_perm,method=cor.method,dat.use=data.use, features.use = rand.features.use)

  }

  #calculate FDR


  obs_null.cor<-calc_obs_null_dist(obs.cor,null.cor)


  SCA.FDRq<-calc_SCA_FDR(obs_null.cor)
  out.list<-list('SCA_result'=SCA.FDRq,'obs_null_SP'=obs_null.cor)

  return(out.list)
}


extract_refp_mult<-function(reference_data, cell_type_slot='SCA_cls',num_signatures=2000){
  #check if the reference data is in sce obj
  if(!class(ref_data_use)=='SingleCellExperiment'){
    print('need SCE object as input')
    stop()
  }
  cls_names<-unique(reference_data[[cell_type_slot]])

  mult.refp<-list()
  for(i in 1:length(cls_names)){
    print(paste0('get ',cls_names[i],' reference profile'))
    SCA.refP<-extract_reference_profile(dat.sce.filt=reference_data, label = cls_names[i],
                                        num_signatures =num_signatures)
    mult.refp[[cls_names[i]]]<-SCA.refP
  }

  return(mult.refp)
}
#need function that takes dataframes of genes and lfcs to make a list of reference profile
#make_refp_mult<-function(){}


SCA_multiple<-function(query_data,cor.method='spearman', ref_profiles,
                       num_perm=1000,data.use='logcounts',
                       rand.features.use='across.sample.refp.gene.permutation', ref_data_null_dist=NULL){
  if(class(ref_profiles)!='list'){
    print('need reference profiles in a list format')
    stop()
  }

  refp.names<-names(ref_profiles)
  SCA_res_mult<-list()

  for(i in 1:length(ref_profiles)){
    if (is.null(refp.names[i]) || refp.names[i] == "") {
      refp_name[i] <- paste0("reference_profile_name_", i)
    }
    print(paste0('run SCA for ',refp.names[i]))
    SCA_res<-SCA(dat.sce.filt=query_data,cor.method=cor.method, ref_profile=ref_profiles[[i]],
                 num_perm=num_perm,data.use=data.use,
                 rand.features.use=rand.features.use, ref_data_null_dist = ref_data_null_dist)
    SCA_res_mult[[refp.names[i]]]<-SCA_res

  }
  return(SCA_res_mult)

}
