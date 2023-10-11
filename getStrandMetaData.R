#Script for getting gene meta data
library(writexl)
library(readxl)
library(writexl)
allgenes_gtf=rtracklayer::import("../Nanopore/gencode.v44.annotation.gtf.gz")

my_data <- read_excel("../Atlassian/Gene Lists/ALL/ALL.genes.xlsx")
my_data <- read_excel("../Atlassian/Gene Lists/AML/AML.genes.1-3.xlsx")

genes=my_data$`Gene name`

df = data.frame()
gene=genes[25]
for (gene in genes){
  if (gene %in% allgenes_gtf@elementMetadata@listData$gene_name){
    tmp = allgenes_gtf[allgenes_gtf@elementMetadata@listData$gene_name==gene,][1]
    tmp = as.data.frame(tmp)
    df_tmp = data.frame(Strand = tmp$strand,GeneName=tmp$gene_name)
    df = rbind(df,df_tmp)
  } else {
    message(gene," not found")
    df_tmp = data.frame(Strand = "NA",GeneName=gene)
    df = rbind(df,df_tmp)
  }
}

write_xlsx(df,"../Atlassian/Gene Lists/AML.genes.with.Strand.xlsx")
 