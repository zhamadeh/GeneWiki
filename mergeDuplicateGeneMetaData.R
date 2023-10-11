#Script for getting gene meta data
library(writexl)
library(readxl)
library(tidyverse)
my_data <- read_excel("../Book1.xlsx")

# Assuming you have the dplyr package installed
# If not, install it using install.packages("dplyr")
library(dplyr)

# Assuming your dataframe is named my_data
# Assuming there is a column named 'Gene' that identifies each gene


collapsed_data <- my_data %>%
  group_by(`Gene name`) %>%
  summarise(
    Chrom = first(Chrom),
    Start = first(Start) ,
    End = first(End),
    Strand = first(Strand),
    `WHO classification` = paste(first(`WHO classification`),last(`WHO classification`),collapse= " ,"),
    `ICC classification` = paste(first(`ICC classification`),last(`ICC classification`),collapse= " ,"), 
    `Most common fusion` = paste(first(`Most common fusion`),last(`Most common fusion`),collapse= " ,"),  
    `Partner genes/ variants` = paste( first(`Partner genes/ variants`), last(`Partner genes/ variants`),collapse= " ,"),
    `GSC Myeloid Panel` = paste( first(`GSC Myeloid Panel`), last(`GSC Myeloid Panel`),collapse= " ,"),
   `UHN Myeloid Panel` = paste( first(`UHN Myeloid Panel`), last(`UHN Myeloid Panel`),collapse= " ,"), 
    `CGC AML` = paste( first(`CGC AML`),last(`CGC AML`),collapse= " ,"),
   `CGC B-ALL` = paste( first(`CGC B-ALL`), last(`CGC B-ALL`),collapse= " ,"), 
   `CGC T-ALL` = paste( first(`CGC T-ALL`), last(`CGC T-ALL`),collapse= " ,"), 
   `Variant Type_1` = first(`Variant Type`),
   `Variant Type_2` = last(`Variant Type`),
   Description_1 = first(Description),
   Description_2 = last(Description),
   Disease_1= first(Disease),
   Disease_2= last(Disease),
   Tier_1 = first(Tier),
   Tier_2 = last(Tier),
   `Targeted therapy_1` = first(`Targeted therapy`),
   `Targeted therapy_2` = last(`Targeted therapy`),
    Text_1 = first(Text),
    Text_2 = last(Text),
    Prognosis_1 = first(Prognosis),
    Prognosis_2 = last(Prognosis)
# Take the last 'Start' value in each group
    # Merge other columns using paste, adjust the column names accordingly
  ) %>%
  ungroup()
collapsed_data=as.data.frame(collapsed_data)

gene=collapsed_data$`Gene name`[1]
for (gene in collapsed_data$`Gene name`){
  print(gene)
  if (gene %in% allgenes_gtf@elementMetadata@listData$gene_name){
    tmp = allgenes_gtf[allgenes_gtf@elementMetadata@listData$gene_name==gene,][1]
    tmp = as.data.frame(tmp)
    if (tmp$strand == "+"){
      collapsed_data[collapsed_data$`Gene name`==gene,5][1] = "+"
    } else{
      collapsed_data[collapsed_data$`Gene name`==gene,5][1] = "-"
    }
    
    message("Found ",gene, " on strand: ",tmp$strand)
  } else {
    message(gene," not found")

  }
}


write_xlsx(collapsed_data,"../book2.xlsx")

