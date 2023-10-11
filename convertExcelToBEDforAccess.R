####################################################################################

# Packages
install.packages("BiocManager")
BiocManager::install("GenomicRanges")
BiocManager::install("rtracklayer")
library(GenomicRanges)
library(BiocManager)
library(rtracklayer)
library(readxl)
library(tidyverse)


####################################################################################

# All genes GTF from Bionano
#allgenes_gtf=rtracklayer::import("Gencode Release 44/gencode.v44.annotation.gtf.gz")
allgenes_gtf=rtracklayer::import("Gencode Release 44/hg38-primary_transcripts from Bionano.gtf")

####################################################################################

tier1_abnormalities=as.data.frame(read_xlsx("Gene Lists/AML.SVs.xlsx"))

tier1_abnormalities = select(tier1_abnormalities,c(Chrom,Start,End,`Gene name`))
tier1_abnormalities = tier1_abnormalities[-1,]
tier1_abnormalities$V1 = "."
tier1_abnormalities$V2 = "."
#AML - Tier 1 Structural Variants
write.table(tier1_abnormalities,"Access Files/AML.tier1.chr.abnormalities.bed",quote = F,row.names = F,col.names = F,sep = "\t")
export(tier1_abnormalities , "Access Files/AML.tier1.chr.abnormalities.bed", "bed")

####################################################################################

tier1_gene=as.data.frame(read_xlsx("Gene Lists/AML.genes.1-3.xlsx"))

tier1_gene = select(tier1_gene,c("Chrom","Start","End","Strand","Gene name"   ))
colnames(tier1_gene)[1:5]=c("chr","start","end","strand","Gene name")

tier1_gene$strand=gsub(x = tier1_gene$strand,pattern = "plus",replacement = "+")
tier1_gene$strand=gsub(x = tier1_gene$strand,pattern = "minus",replacement = "-")

####################################################################################

# 251 genes in tier 1 file
length(levels(droplevels(as.factor(tier1_gene$`Gene name`)))) 


# 61228 unique genes in All genes
length(levels(droplevels(as.factor(allgenes_gtf$gene_name))))

# 240 intersecting genes
length(intersect(levels(droplevels(as.factor(tier1_gene$`Gene name`))), 
                 levels(droplevels(as.factor(allgenes_gtf$gene_name)))))

# 3 missing genes: "CRLF2_chrX" "CRLF2_chrY" "FGFR1OP"    "GUCY1B3"    "LINC00982"  
#                 "MKL1"     "MLLT4"     "P2RY8_chrX" "P2RY8_chrY" "STAT5a" "STAT5b"   
length(setdiff(levels(droplevels(as.factor(tier1_gene$`Gene name`))), 
               levels(droplevels(as.factor(allgenes_gtf$gene_name)))))

setdiff(levels(droplevels(as.factor(tier1_gene$`Gene name`))), 
        levels(droplevels(as.factor(allgenes_gtf$gene_name))))

# 240 genes make up 89583 elements in GTF format 
tier1_gene_gtf = allgenes_gtf[allgenes_gtf$gene_name %in% intersect(levels(droplevels(as.factor(tier1_gene$`Gene name`))), 
                                                                    levels(droplevels(as.factor(allgenes_gtf$gene_name)))),]

IGH = data.frame(seqnames="chr14", start="105586437",end="106879844",width= "1293407",strand="-",
                 source="refGene",type="transcript",score="NA",phase="NA",gene_id = "IGH", transcript_id="NA",gene_name="IGH",exon_number="NA",exon_id="NA")

#IGH = data.frame(seqnames="chr14", start="105586437",end="106879844",width= "1293407",strand="-",
#                 source="refGene",type="transcript",score="NA",phase="NA",gene_id = "IGH", gene_type="NA",gene_name="IGH",level="IGH", tag ="IGH", transcript_id="IGH",transcript_type="IGH",
#                 transcript_name="IGH",transcript_support_level="NA", havana_transcript="NA", exon_number="NA",exon_id="NA",
#                 hgnc_id = "NA", havana_gene="NA", ont="NA", protein_id="NA", ccdsid="NA", artif_dupl="NA" )
#colnames(as.data.frame(tier1_gene_gtf))
#colnames(as.data.frame(IGH))

tier1_gene_gtf=GRanges(rbind(as.data.frame(tier1_gene_gtf),IGH))

#tier1_gene_gtf = GRanges(select(as.data.frame(tier1_gene_gtf),c("seqnames","start" ,"end" ,"width","strand","source","type","score" ,"phase","gene_id","transcript_id","gene_name","exon_number","exon_id")))

#tier1_gene_gtf = as.data.frame(tier1_gene_gtf)
#tier1_gene_gtf$gene_id = tier1_gene_gtf$gene_name
#tier1_gene_gtf$source = "refGene"
#tier1_gene_gtf = GRanges(tier1_gene_gtf)

         
# Ensure header is removed from this file before uploading
# "AML - Tier 1A, 1B, 2 Gene Variants"
export(tier1_gene_gtf , "Access Files/AML.tier1.gene.variants.gtf", "gtf")



####################################################################################

tier1_gene=as.data.frame(read_xlsx("Gene Lists/ALL/ALL.genes.xlsx"))

tier1_gene = select(tier1_gene,c("Chrom","Start","End","Strand","Gene"   ))
colnames(tier1_gene)[1:5]=c("chr","start","end","strand","Gene name")

tier1_gene$strand=gsub(x = tier1_gene$strand,pattern = "plus",replacement = "+")
tier1_gene$strand=gsub(x = tier1_gene$strand,pattern = "minus",replacement = "-")

####################################################################################

# 125 genes in tier 1 file
length(levels(droplevels(as.factor(tier1_gene$`Gene name`)))) 


# 61228 unique genes in All genes
length(levels(droplevels(as.factor(allgenes_gtf$gene_name))))

# 110 intersecting genes
length(intersect(levels(droplevels(as.factor(tier1_gene$`Gene name`))), 
                 levels(droplevels(as.factor(allgenes_gtf$gene_name)))))

# 15 missing genes:  [1] "11q"           "HIST1H2AD"     "HIST1H2BE"     "HIST1H2BF"    
# "HIST1H3D"      "HIST1H4D"      "IGH"           "OSTL"         
# "PAR1_deletion" "TNFRSF6"       "TRA"           "TRA-D"        
# "TRB"           "TRD"           "TRG"   
length(setdiff(levels(droplevels(as.factor(tier1_gene$`Gene name`))), 
        levels(droplevels(as.factor(allgenes_gtf$gene_name)))))

setdiff(levels(droplevels(as.factor(tier1_gene$`Gene name`))), 
        levels(droplevels(as.factor(allgenes_gtf$gene_name))))

# 240 genes make up 89583 elements in GTF format 
tier1_gene_gtf = allgenes_gtf[allgenes_gtf$gene_name %in% intersect(levels(droplevels(as.factor(tier1_gene$`Gene name`))), 
                                                                    levels(droplevels(as.factor(allgenes_gtf$gene_name)))),]

IGH = data.frame(seqnames="chr14", start="105586437",end="106879844",width= "1293407",strand="-",
                 source="refGene",type="transcript",score="NA",phase="NA",gene_id = "IGH", transcript_id="NA",gene_name="IGH",exon_number="NA",exon_id="NA")

#IGH = data.frame(seqnames="chr14", start="105586437",end="106879844",width= "1293407",strand="-",
#                 source="refGene",type="transcript",score="NA",phase="NA",gene_id = "IGH", gene_type="NA",gene_name="IGH",level="IGH", tag ="IGH", transcript_id="IGH",transcript_type="IGH",
#                 transcript_name="IGH",transcript_support_level="NA", havana_transcript="NA", exon_number="NA",exon_id="NA",
#                 hgnc_id = "NA", havana_gene="NA", ont="NA", protein_id="NA", ccdsid="NA", artif_dupl="NA" )
#colnames(as.data.frame(tier1_gene_gtf))
#colnames(as.data.frame(IGH))

tier1_gene_gtf=GRanges(rbind(as.data.frame(tier1_gene_gtf),IGH))

#tier1_gene_gtf = GRanges(select(as.data.frame(tier1_gene_gtf),c("seqnames","start" ,"end" ,"width","strand","source","type","score" ,"phase","gene_id","transcript_id","gene_name","exon_number","exon_id")))

#tier1_gene_gtf = as.data.frame(tier1_gene_gtf)
#tier1_gene_gtf$gene_id = tier1_gene_gtf$gene_name
#tier1_gene_gtf$source = "refGene"
#tier1_gene_gtf = GRanges(tier1_gene_gtf)


# Ensure header is removed from this file before uploading
# "AML - Tier 1A, 1B, 2 Gene Variants"
export(tier1_gene_gtf , "Access Files/ALL.tier1.gene.variants.gtf", "gtf")





####################################################################################

tier1_gene=as.data.frame(read_xlsx("Gene Lists/AML-ALL/AML.ALL.xlsx"))

tier1_gene = select(tier1_gene,c("Chrom","Start","End","Strand","Gene name"   ))
colnames(tier1_gene)[1:5]=c("chr","start","end","strand","Gene name")

tier1_gene$strand=gsub(x = tier1_gene$strand,pattern = "plus",replacement = "+")
tier1_gene$strand=gsub(x = tier1_gene$strand,pattern = "minus",replacement = "-")

####################################################################################

# 125 genes in tier 1 file
length(levels(droplevels(as.factor(tier1_gene$`Gene name`)))) 


# 61228 unique genes in All genes
length(levels(droplevels(as.factor(allgenes_gtf$gene_name))))

# 110 intersecting genes
length(intersect(levels(droplevels(as.factor(tier1_gene$`Gene name`))), 
                 levels(droplevels(as.factor(allgenes_gtf$gene_name)))))

# 15 missing genes: 
length(setdiff(levels(droplevels(as.factor(tier1_gene$`Gene name`))), 
               levels(droplevels(as.factor(allgenes_gtf$gene_name)))))

setdiff(levels(droplevels(as.factor(tier1_gene$`Gene name`))), 
        levels(droplevels(as.factor(allgenes_gtf$gene_name))))

# 240 genes make up 89583 elements in GTF format 
tier1_gene_gtf = allgenes_gtf[allgenes_gtf$gene_name %in% intersect(levels(droplevels(as.factor(tier1_gene$`Gene name`))), 
                                                                    levels(droplevels(as.factor(allgenes_gtf$gene_name)))),]

# Ensure header is removed from this file before uploading
# "AML - Tier 1A, 1B, 2 Gene Variants"
export(tier1_gene_gtf , "Access Files/AML.ALL.tier1.gene.variants.gtf", "gtf")
