#how example HSPC data was prepared

packages<-c('dplyr','SingleCellExperiment',
            'scater','scran')
for (p in packages) {
  library(p, character.only = TRUE)
}

Ind1.sce<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/SCA 2023/Data/Velten Nature Cell Biology 2017/Data used 2024/Ind1_sce_QC_20240308.RDS')
colData(Ind1.sce)$SCA_cls<-colData(Ind1.sce)$cls

column_data<-colData(Ind1.sce) %>% as.data.frame()
column_data <- column_data %>%
  mutate(SCA_cls = recode(SCA_cls, "CMP(CD38pCD10nCD45RAnCD135p)"="CMP", "GMP(CD38pCD10nCD45RApCD135p)"="GMP", "MEP(CD38pCD10nCD45RAnCD135n)"="MEP",
                          "B-NKP(CD38pCD10p)"="B-NKP","other HSC(CD38nCD45RAnCD90pCD49fn)"="other HSC","LT-HSC(CD38nCD45RAnCD90pCD49fp)"="LT-HSC","CD135pMPP"="MPP",
                          "CD135nMPP"="MPP","CD10pMLP(CD38nCD45RApCD10p)"="MLP",     "CD10nMLP(CD38nCD45RApCD10n)"="MLP","CD38nCD135n"="other", "MPP(CD38nCD45RAnCD90n)"="MPP",
                          "CD38nCD45RAp"="other", "CD38nCD45RAnCD90dim"="other", "CD38n"="other","all"="other",  "CD38pCD10n"="other"))

colData(Ind1.sce)<-DataFrame(column_data)
rownames(Ind1.sce)<-make.unique(rowData(Ind1.sce)$Final.Symbol)

BM5.34p.sce.n<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/SCA 2023/Human data/Galen_BM/Galen_BM5_34p_sce_QC_20240430.rds')
rownames(BM5.34p.sce.n)<-rowData(BM5.34p.sce.n)$HGNC.final.use
colData(BM5.34p.sce.n)$SCA_cls<-colData(BM5.34p.sce.n)$cluster
G.BM5.34p.sce01<-BM5.34p.sce.n

Velten_HSPC<-Ind1.sce
vanGalen_HSPC<-G.BM5.34p.sce01
library(usethis)
use_data(Velten_HSPC)
use_data(vanGalen_HSPC)


