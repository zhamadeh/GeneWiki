#Script for getting gene meta data
library(writexl)
library(readxl)
library(tidyverse)


my_data <- read_excel("Gene Lists/AML.genes.xlsx")
my_data$Start = prettyNum(my_data$Start, big.mark = ",", scientific = FALSE)
my_data$End = prettyNum(my_data$End, big.mark = ",", scientific = FALSE)

my_data_svs = filter(my_data,Description=="Structural aberration")
my_data_genes = filter(my_data,Description!="Structural aberration")

write.table(as.data.frame(split(my_data_genes, factor(sort(rank(row.names(my_data_genes))%%2)))[[1]]),
            "Gene Lists/Text/AML.genes.1.txt",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(as.data.frame(split(my_data_genes, factor(sort(rank(row.names(my_data_genes))%%2)))[[2]]),
            "Gene Lists/Text/AML.genes.2.txt",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(my_data_svs,
            "Gene Lists/Text/AML.SVs.txt",quote = F,row.names = F,col.names = F,sep = "\t")




my_data <- read_excel("Gene Lists/ALL.genes.xlsx")
my_data$Start = prettyNum(my_data$Start, big.mark = ",", scientific = FALSE)
my_data$End = prettyNum(my_data$End, big.mark = ",", scientific = FALSE)

write.table(as.data.frame(split(my_data_genes, factor(sort(rank(row.names(my_data_genes))%%2)))[[1]]),
            "Gene Lists/Text/ALL.genes.1.txt",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(as.data.frame(split(my_data_genes, factor(sort(rank(row.names(my_data_genes))%%2)))[[2]]),
            "Gene Lists/Text/ALL.genes.2.txt",quote = F,row.names = F,col.names = F,sep = "\t")




my_data <- read_excel("Gene Lists/AML.ALL.genes.xlsx")
my_data$Start = prettyNum(my_data$Start, big.mark = ",", scientific = FALSE)
my_data$End = prettyNum(my_data$End, big.mark = ",", scientific = FALSE)

write.table(as.data.frame(split(my_data_genes, factor(sort(rank(row.names(my_data_genes))%%2)))[[1]]),
            "Gene Lists/Text/AML.ALL.genes.txt",quote = F,row.names = F,col.names = F,sep = "\t")



