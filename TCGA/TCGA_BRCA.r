

# BiocManager::install("TCGAbiolinks")
library(SummarizedExperiment)
library(TCGAbiolinks)

query.exp <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Solid Tissue Normal")
)
GDCdownload(
  query = query.exp,
  files.per.chunk = 100
)

brca.exp <- GDCprepare(
  query = query.exp, 
  save = TRUE, 
  save.filename = "brcaExp.rda"
)

# get subtype information
infomation.subtype <- TCGAquery_subtype(tumor = "BRCA")

# get clinical data
information.clinical <- GDCquery_clinic(project = "TCGA-BRCA",type = "clinical") 

# Which samples are Primary Tumor
samples.primary.tumour <- brca.exp$barcode[brca.exp$shortLetterCode == "TP"]

# which samples are solid tissue normal
samples.solid.tissue.normal <- brca.exp$barcode[brca.exp$shortLetterCode == "NT"]

dataPrep <- TCGAanalyze_Preprocessing(
  object = brca.exp, 
  cor.cut = 0.6
)                      

dataNorm <- TCGAanalyze_Normalization(
  tabDF = dataPrep,
  geneInfo = geneInfoHT,
  method = "gcContent"
)                

dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm,
  method = "quantile", 
  qnt.cut =  0.25
)   

dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[,samples.solid.tissue.normal],
  mat2 = dataFilt[,samples.primary.tumour],
  Cond1type = "Normal",
  Cond2type = "Tumor",
  fdr.cut = 0.01 ,
  logFC.cut = 2,
  method = "glmLRT",
  pipeline = "edgeR"
)  


ansEA <- TCGAanalyze_EAcomplete(
  TFname = "DEA genes Normal Vs Tumor",
  RegulonList = dataDEGs$gene_name
)  

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  GOBPTab = ansEA$ResBP,
  GOCCTab = ansEA$ResCC,
  GOMFTab = ansEA$ResMF,
  PathTab = ansEA$ResPat,
  nRGTab = dataDEGs$gene_name,
  nBar = 10
)


group1 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("NT"))
group2 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("TP"))

dataSurv <- TCGAanalyze_SurvivalKM(
  clinical_patient = dataClin,
  dataGE = dataFilt,
  Genelist = rownames(dataDEGs),
  Survresult = FALSE,
  ThreshTop = 0.67,
  ThreshDown = 0.33,
  p.cut = 0.05, 
  group1 = group1, 
  group2 = group2
)



# library(TCGAretriever)
# library(reshape2)
# library(ggplot2)
# 
# install.packages("TCGAretriever")
# # Obtain a list of cancer studies from cBio
# all_studies <- get_cancer_studies()
# 
# # Find published TCGA datasets
# keep <- grepl("tcga_pub$", all_studies[,1])
# tcga_studies <- all_studies[keep, ]
# 
# # Show results
# show_head(tcga_studies, 6, 2)
