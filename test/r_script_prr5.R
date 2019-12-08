## Script para determinar los genes dianas de un factor de transcripción
## a partir del fichero narrowPeak generado por MaCS2.

## Authors: Martin Moreno-Perez and Jose Pablo Rivas-Fernandez
## E-mail: marmorper20@alum.us.es - josrivfer1@alum.us.es
## Date: November 2019

# Loading packages
library(ChIPseeker)

library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28

library("clusterProfiler")

library("org.At.tair.db")

library("topGO")

library("pathview")


args <- commandArgs(trailingOnly = T)
picos <- args[[1]]


# Reading peaks file
prr5.peaks <- readPeakFile(peakfile = picos, header=FALSE)
prr5.peaks

# Defining the promoter region around TSS
promoter <- getPromoters(TxDb=txdb, 
                         upstream=1000, 
                         downstream=1000)
promoter

# Annotation of peaks

prr5.peakAnno <- annotatePeak(peak = prr5.peaks, 
                             tssRegion=c(-1000, 1000),
                             TxDb=txdb)
prr5.peakAnno

plotAnnoPie(prr5.peakAnno)
plotDistToTSS(prr5.peakAnno,
              title="Distribution of genomic loci relative to TSS",
              ylab = "Genomic Loci (%) (5' -> 3')")



# jpeg(filename = "plot_e_bp.jpeg")
# plot.e.bp<- plotGOgraph(x = e.bp,
#                         firstSigNodes = 10,
#                         useInfo = "all",
#                         sigForAll = TRUE,
#                         useFullNames = TRUE)
# dev.off()

# Converting annotation to data frame

prr5.annotation <- as.data.frame(prr5.peakAnno)
head(prr5.annotation)

target.genes.prr5 <- prr5.annotation$geneId[prr5.annotation$annotation == "Promoter"]
length(target.genes.prr5)

write(x = target.genes.prr5,file = "prr5_target_genes.txt")


# Reading genes

gene.set.prr5<- read.table(file = "prr5_target_genes.txt", header = F, as.is = T)[[1]]
gene.set.prr5
length(gene.set.prr5)


# Building universe
ath.genes<- as.data.frame(genes(txdb))  # probar solo genes(txdb) si no funciona
#ath.genes$seqnames
#ath.genes$seqnames == 1
ath.genes.names<- ath.genes$gene_id
#genes.chr1<- (ath.genes$gene_id)[ath.genes$seqnames == 1]
write(x = ath.genes.names,file = "ath_genes_names.txt")


# GO ENRICHMENT

idType(OrgDb = org.At.tair.db)

# CONTROL

enrich.go.prr5<- enrichGO(gene = gene.set.prr5,
                         OrgDb = org.At.tair.db,
                         keyType = "TAIR",
                         ont = "ALL",
                         pvalueCutoff = 0.05,
                         pAdjustMethod = 0.05,
                         universe = ath.genes.names)

# Error in match.arg(method) : 'arg' must be NULL or a character vector

summary(enrich.go.prr5)

enrich.go.prr5@geneSets

# ENRICHMENTS:
#       BP (biological process),
#       MF (molecular function).
#       CC (celular compartiment).

### BP
enrich.go.prr5.bp<- enrichGO(gene = gene.set.prr5,
                            OrgDb = org.At.tair.db,
                            keyType = "TAIR",
                            ont = "BP",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05,
                            universe = ath.genes.names)

head(enrich.go.prr5.bp)

enrich.go.prr5.bp.table <- as.data.frame(enrich.go.prr5.bp)

head(enrich.go.prr5.bp.table)

barplot(enrich.go.prr5.bp, showCategory=20)


### MF
enrich.go.prr5.mf<- enrichGO(gene = gene.set.prr5,
                            OrgDb = org.At.tair.db,
                            keyType = "TAIR",
                            ont = "MF",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05,
                            universe = ath.genes.names)

head(enrich.go.prr5.mf)

enrich.go.prr5.mf.table <- as.data.frame(enrich.go.prr5.mf)

head(enrich.go.prr5.mf.table)

barplot(enrich.go.prr5.mf, showCategory=20)

### CC
enrich.go.prr5.cc<- enrichGO(gene = gene.set.prr5,
                            OrgDb = org.At.tair.db,
                            keyType = "TAIR",
                            ont = "CC",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05,
                            universe = ath.genes.names)

head(enrich.go.prr5.cc)

enrich.go.prr5.cc.table <- as.data.frame(enrich.go.prr5.cc)

head(enrich.go.prr5.cc.table)

barplot(enrich.go.prr5.cc, showCategory=20)


# jpeg(filename = "plot_e_bp.jpeg")
# plot.e.bp<- plotGOgraph(x = e.bp,
#                         firstSigNodes = 10,
#                         useInfo = "all",
#                         sigForAll = TRUE,
#                         useFullNames = TRUE)
# dev.off()
# 
# 
# 
# 
# 
# e.mf<- enrichGO(gene = target.genes,
#                 OrgDb = org.At.tair.db,
#                 keyType = "TAIR",
#                 ont = "MF")
# 
# jpeg(filename = "plot_e_mf.jpeg")
# plot.e.mf<- plotGOgraph(x = e.mf,
#                         firstSigNodes = 10,
#                         useInfo = "all",
#                         sigForAll = TRUE,
#                         useFullNames = TRUE)
# dev.off()
# 
# 
# 
#   
# e.cc<- enrichGO(gene = target.genes,
#                 OrgDb = org.At.tair.db,
#                 keyType = "TAIR",
#                 ont = "CC")
# 
# jpeg(filename = "plot_e_cc.jpeg")
# plot.e.cc<- plotGOgraph(x = e.cc,
#                         firstSigNodes = 10,
#                         useInfo = "all",
#                         sigForAll = TRUE,
#                         useFullNames = TRUE)
# dev.off()


# KEGG ENRICHMENT


enrich.kegg.prr5<- enrichKEGG(gene = gene.set.prr5,
                             organism = "ath",
                             universe = ath.genes.names)

prr5.kegg.table <- as.data.frame(enrich.kegg.prr5)

head(prr5.kegg.table)





#########################################

# Preparación de los genes

my.universe.prr5<- rep(0, length(ath.genes.names))
names(my.universe.prr5) <- ath.genes.names
my.universe.prr5
my.universe.prr5[gene.set.prr5] <- 1
my.universe.prr5

# Path.id de cada función (en la función anterior)

kegg.pathview.prr5<- pathview(gene.data = my.universe, 
                             pathway.id = "ath00195",
                             species = "ath",
                             gene.idtype = "TAIR")

