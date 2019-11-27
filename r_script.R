## Script para determinar los genes dianas de un factor de transcripción
## a partir del fichero narrowPeak generado por MaCS2.

## Authors: Martin Moreno-Perez and Jose Pablo Rivas Fernandez
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
peaks.control <- args[[1]]
peaks.treatment <- args[[2]]


# Reading peaks file
col.peaks <- readPeakFile(peakfile = peaks.control, header=FALSE)
col.peaks

trat.peaks <- readPeakFile(peakfile = peaks.treatment, header=FALSE)
trat.peaks

# Defining the promoter region around TSS
promoter <- getPromoters(TxDb=txdb, 
                         upstream=1000, 
                         downstream=1000)
promoter

# Annotation of peaks
# CONTROL
col.peakAnno <- annotatePeak(peak = col.peaks, 
                             tssRegion=c(-1000, 1000),
                             TxDb=txdb)
col.peakAnno

plotAnnoPie(col.peakAnno)
plotDistToTSS(col.peakAnno,
              title="Distribution of genomic loci relative to TSS",
              ylab = "Genomic Loci (%) (5' -> 3')")

# TREATMENT
trat.peakAnno <- annotatePeak(peak = trat.peaks, 
                             tssRegion=c(-1000, 1000),
                             TxDb=txdb)
trat.peakAnno

# jpeg(filename = "plot_e_bp.jpeg")
# plot.e.bp<- plotGOgraph(x = e.bp,
#                         firstSigNodes = 10,
#                         useInfo = "all",
#                         sigForAll = TRUE,
#                         useFullNames = TRUE)
# dev.off()

plotAnnoPie(trat.peakAnno)
plotDistToTSS(trat.peakAnno,
              title="Distribution of genomic loci relative to TSS",
              ylab = "Genomic Loci (%) (5' -> 3')")

# Converting annotation to data frame
# CONTROL
col.annotation <- as.data.frame(col.peakAnno)
head(col.annotation)

target.genes.col <- col.annotation$geneId[col.annotation$annotation == "Promoter"]
length(target.genes.col)

write(x = target.genes.col,file = "col_target_genes.txt")


# TREATMENT
trat.annotation <- as.data.frame(trat.peakAnno)
head(trat.annotation)

target.genes.trat <- trat.annotation$geneId[trat.annotation$annotation == "Promoter"]
length(target.genes.trat)

write(x = target.genes.trat,file = "trat_target_genes.txt")




# Reading genes
# CONTROL
gene.set.col<- read.table(file = "col_target_genes.txt", header = F, as.is = T)[[1]]
gene.set.col
length(gene.set.col)

# TREATMENT
gene.set.trat<- read.table(file = "trat_target_genes.txt", header = F, as.is = T)[[1]]
gene.set.trat
length(gene.set.trat)


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

enrich.go.col<- enrichGO(gene = gene.set.col,
                         OrgDb = org.At.tair.db,
                         keyType = "TAIR",
                         ont = "ALL",
                         pvalueCutoff = 0.05,
                         pAdjustMethod = 0.05,
                         universe = ath.genes.names)

# Error in match.arg(method) : 'arg' must be NULL or a character vector

summary(enrich.go.col)

enrich.go.col@geneSets

# ENRICHMENTS:
#       BP (biological process),
#       MF (molecular function).
#       CC (celular compartiment).

### BP
enrich.go.col.bp<- enrichGO(gene = gene.set.col,
                            OrgDb = org.At.tair.db,
                            keyType = "TAIR",
                            ont = "BP",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05,
                            universe = ath.genes.names)

head(enrich.go.col.bp)

enrich.go.col.bp.table <- as.data.frame(enrich.go.col.bp)

head(enrich.go.col.bp.table)

barplot(enrich.go.col.bp, showCategory=20)


### MF
enrich.go.col.mf<- enrichGO(gene = gene.set.col,
                            OrgDb = org.At.tair.db,
                            keyType = "TAIR",
                            ont = "MF",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05,
                            universe = ath.genes.names)

head(enrich.go.col.mf)

enrich.go.col.mf.table <- as.data.frame(enrich.go.col.mf)

head(enrich.go.col.mf.table)

barplot(enrich.go.col.mf, showCategory=20)

### CC
enrich.go.col.cc<- enrichGO(gene = gene.set.col,
                            OrgDb = org.At.tair.db,
                            keyType = "TAIR",
                            ont = "CC",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05,
                            universe = ath.genes.names)

head(enrich.go.col.cc)

enrich.go.col.cc.table <- as.data.frame(enrich.go.col.cc)

head(enrich.go.col.cc.table)

barplot(enrich.go.col.cc, showCategory=20)


# TREATMENT

enrich.go.trat<- enrichGO(gene = gene.set.trat,
                          OrgDb = org.At.tair.db,
                          keyType = "TAIR",
                          ont = "ALL",
                          pvalueCutoff = 0.1,
                          pAdjustMethod = 0.1,
                          universe = ath.genes.names)

# Error in match.arg(method) : 'arg' must be NULL or a character vector

summary(enrich.go.trat)

enrich.go.trat@geneSets


# Para realizar los enriquecimientos por separado:
#       BP (biological process),
#       MF (molecular function).
#       CC (celular compartiment).

### BP
enrich.go.trat.bp<- enrichGO(gene = gene.set.trat,
                             OrgDb = org.At.tair.db,
                             keyType = "TAIR",
                             ont = "BP",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05,
                             universe = ath.genes.names)

head(enrich.go.trat.bp)

enrich.go.trat.bp.table <- as.data.frame(enrich.go.trat.bp)

head(enrich.go.trat.bp.table)

barplot(enrich.go.trat.bp, showCategory=20)

### MF
enrich.go.trat.mf<- enrichGO(gene = gene.set.trat,
                            OrgDb = org.At.tair.db,
                            keyType = "TAIR",
                            ont = "MF",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05,
                            universe = ath.genes.names)

head(enrich.go.trat.mf)

enrich.go.trat.mf.table <- as.data.frame(enrich.go.trat.mf)

head(enrich.go.trat.mf.table)

barplot(enrich.go.trat.mf, showCategory=20)

### CC
enrich.go.trat.cc<- enrichGO(gene = gene.set.trat,
                            OrgDb = org.At.tair.db,
                            keyType = "TAIR",
                            ont = "CC",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05,
                            universe = ath.genes.names)

head(enrich.go.trat.cc)

enrich.go.trat.cc.table <- as.data.frame(enrich.go.trat.cc)

head(enrich.go.trat.cc.table)

barplot(enrich.go.trat.cc, showCategory=20)

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

# CONTROL
enrich.kegg.col<- enrichKEGG(gene = gene.set.col,
                             organism = "ath",
                             universe = ath.genes.names)

col.kegg.table <- as.data.frame(enrich.kegg.col)

head(col.kegg.table)

# TREATMENT
enrich.kegg.trat<- enrichKEGG(gene = gene.set.trat,
                             organism = "ath",
                             universe = ath.genes.names)

trat.kegg.table <- as.data.frame(enrich.kegg.trat)

head(trat.kegg.table)



#########################################

# Preparación de los genes

## CONTROL
my.universe.col<- rep(0, length(ath.genes.names))
names(my.universe.col) <- ath.genes.names
my.universe.col
my.universe.col[gene.set.col] <- 1
my.universe.col

# Path.id de cada función (en la función anterior)

kegg.pathview.col<- pathview(gene.data = my.universe, 
                         pathway.id = "ath00195",
                         species = "ath",
                         gene.idtype = "TAIR")


## TRATAMIENTO
my.universe.trat<- rep(0, length(ath.genes.names))
names(my.universe.trat) <- ath.genes.names
my.universe.trat
my.universe.trat[gene.set.trat] <- 1
my.universe.trat

# Path.id de cada función (en la función anterior)

kegg.pathview.col<- pathview(gene.data = my.universe, 
                             pathway.id = "ath00195",
                             species = "ath",
                             gene.idtype = "TAIR")

