#' Data for SCA
#'
#' @docType data
#' @name SCA
NULL

#' Velten et al (2017) normal HSPC data: index-sorted scRNAseq data
#'
#' @format SCE objects of normal HSPC datasets
#' \describe{
#'        \item{assays}{counts: read counts single cell gene expression matrix & logcounts: log2 counts scGEP matrix}
#'        \item{colData}{DataFrame of single cell meta data. "SCA_cls" slot contains cleaned up version of cell types annotated by original authors}
#'        }
#'        @source{add source later}
#'        @examples
#'        # example code
#'        data(Velten_HSPC)
"Velten_HSPC"

#' van Galen et al (2019) normal HSPC data (BM5 CD34p): cell types annotated by clustering and cell type specific marker gene expressions
#'
#' @format SCE objects of normal HSPC datasets
#' \describe{
#'        \item{assays}{counts: UMI counts single cell gene expression matrix & logcounts: log2 counts scGEP matrix}
#'        \item{colData}{DataFrame of single cell meta data. "SCA_cls" slot contains cell types annotated by original authors}
#'        }
#'        @source{add source later}
#'        @examples
#'        # example code
#'        data(vanGalen_HSPC)
"vanGalen_HSPC"

