#Script for getting gene meta data
library(writexl)
library(readxl)
library(tidyverse)


my_data <- read_excel("Gene Lists/AML.genes.xlsx")
my_data$Start = prettyNum(my_data$Start, big.mark = ",", scientific = FALSE)
my_data$End = prettyNum(my_data$End, big.mark = ",", scientific = FALSE)

my_data_genes = my_data

my_data_svs = read_excel("Gene Lists/AML.SVs.xlsx")
my_data_svs$Start = prettyNum(my_data_svs$Start, big.mark = ",", scientific = FALSE)
my_data_svs$End = prettyNum(my_data_svs$End, big.mark = ",", scientific = FALSE)


write.table(as.data.frame(split(my_data_genes, factor(sort(rank(row.names(my_data_genes))%%2)))[[1]]),
            "Gene Lists/Text/AML.genes.1.txt",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(as.data.frame(split(my_data_genes, factor(sort(rank(row.names(my_data_genes))%%2)))[[2]]),
            "Gene Lists/Text/AML.genes.2.txt",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(my_data_svs,
            "Gene Lists/Text/AML.SVs.txt",quote = F,row.names = F,col.names = F,sep = "\t")




my_data <- read_excel("Gene Lists/ALL.genes.xlsx")
my_data$Start = prettyNum(my_data$Start, big.mark = ",", scientific = FALSE)
my_data$End = prettyNum(my_data$End, big.mark = ",", scientific = FALSE)
my_data_genes=my_data
write.table(as.data.frame(split(my_data_genes, factor(sort(rank(row.names(my_data_genes))%%2)))[[1]]),
            "Gene Lists/Text/ALL.genes.1.txt",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(as.data.frame(split(my_data_genes, factor(sort(rank(row.names(my_data_genes))%%2)))[[2]]),
            "Gene Lists/Text/ALL.genes.2.txt",quote = F,row.names = F,col.names = F,sep = "\t")




allgenes_gtf=rtracklayer::import("Gencode Release 44/gencode.v44.annotation.gtf.gz")

my_data <- read_excel("Gene Lists/AML.ALL.genes.xlsx")
genes=my_data$`Gene name`

for (gene in genes){
  if (gene %in% allgenes_gtf@elementMetadata@listData$gene_name){
    tmp = allgenes_gtf[allgenes_gtf@elementMetadata@listData$gene_name==gene,][1]
    tmp = as.data.frame(tmp)
    my_data[my_data$`Gene name`==gene,5][[1]] = tmp$strand
  } else {
    message(gene," not found")
    my_data[my_data$`Gene name`==gene,5][[1]] =  "NA"
  }
}


my_data$Start = prettyNum(my_data$Start, big.mark = ",", scientific = FALSE)
my_data$End = prettyNum(my_data$End, big.mark = ",", scientific = FALSE)
my_data_genes=my_data
write.table(as.data.frame(my_data),
            "Gene Lists/Text/AML.ALL.genes.txt",quote = F,row.names = F,col.names = F,sep = "\t")



