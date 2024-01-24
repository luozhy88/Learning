#https://portal.gdc.cancer.gov/exploration?facetTab=cases&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22content%22%3A%7B%22field%22%3A%22cases.demographic.gender%22%2C%22value%22%3A%5B%22male%22%5D%7D%2C%22op%22%3A%22in%22%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.primary_site%22%2C%22value%22%3A%5B%22thyroid%20gland%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.program.name%22%2C%22value%22%3A%5B%22TCGA%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TCGA-THCA%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.samples.sample_type%22%2C%22value%22%3A%5B%22primary%20tumor%22%2C%22solid%20tissue%20normal%22%5D%7D%7D%2C%7B%22content%22%3A%7B%22field%22%3A%22genes.is_cancer_gene_census%22%2C%22value%22%3A%5B%22true%22%5D%7D%2C%22op%22%3A%22in%22%7D%5D%7D&searchTableTab=cases  样本信息从网页标题栏Exploration中查看，主要信息可以获得sample.type的信息

# BiocManager::install("TCGAbiolinks")
library(SummarizedExperiment)
library(TCGAbiolinks)

query.exp <- GDCquery(
  project = "TCGA-THCA", 
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
infomation.subtype <- TCGAquery_subtype(tumor = "THCA")

# get clinical data
information.clinical <- GDCquery_clinic(project = "TCGA-THCA",type = "clinical") 

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

# heatmap
# TCGAvisualize_Heatmap(
#   data = t(dataFilt) )
# 
# 
# TCGAvisualize_Heatmap(
#   data = t(datFilt),
#   col.metadata =  colData(brca.exp)[,
#                                    c("barcode",
#                                      "sample_type")
#   ],
#   # col.colors =  list(
#   #   groupsHC = c(
#   #     "EC1"="black",
#   #     "EC2"="red",
#   #     "EC3"="blue",
#   #     "EC4"="green3")
#   # ),
#   sortCol = "sample_type",
#   type = "expression", # sets default color
#   scale = "row", # use z-scores for better visualization. Center gene expression level around 0.
#   title = "Heatmap from concensus cluster", 
#   filename = "case2_Heatmap.png",
#   extremes = seq(-2,2,1),
#   color.levels = colorRampPalette(c("green", "black", "red"))(n = 5),
#   cluster_rows = TRUE,
#   cluster_columns = FALSE,
#   width = 1000,
#   height = 500
# )






group1 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("NT"))
group2 <- TCGAquery_SampleTypes(colnames(dataFilt

# Volcano plot
TCGAVisualize_volcano(
  x = dataDEGs$logFC,
  y = dataDEGs$FDR,
  x.cut = 1.5,
  y.cut = 0.0001,
  title = "Title example",
  xlab = expression(paste(Log[2], "FoldChange"))
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

# PCA
pca <- TCGAvisualize_PCA(dataFilt,dataDEGs, ntopgenes = 200, group1, group2)



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
