#Script for getting gene meta data
library(writexl)
allgenes_gtf=rtracklayer::import("../Nanopore/gencode.v44.annotation.gtf.gz")

genes=c('TAL1', 'TAL2', 'TLX1', 'TLX3', 'LMO1', 'LMO2', 'NKX2', 'SPI1', 'BHLH')

df = data.frame()
gene=genes[1]
for (gene in genes){
  tmp = allgenes_gtf[allgenes_gtf@elementMetadata@listData$gene_name==gene,][1]
  tmp = as.data.frame(tmp)
  df_tmp = data.frame(Chrom = tmp$seqnames,Start=tmp$start,End = tmp$end,Strand = tmp$strand,Tier = NA,CGL="",GSC="",ICC="",geneID=tmp$transcript_id,GeneName=tmp$gene_name)
  df = rbind(df,df_tmp)
}

write_xlsx(df,"../genes_tmp.xlsx")
