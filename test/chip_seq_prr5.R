## Script para determinar los genes dianas de un factor de transcripción
## a partir del fichero narrowPeak generado por MaCS2.

## Autor: Martín Moreno-Pérez y José Pablo Rivas-Fernández - marmorper20@alum.us.es y josrivfer1@alum.us.es
## Fecha: Noviembre 2019

## Instalar chipseeker y paquete de anotación de Arabidopsis thaliana

library(ChIPseeker)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28



## Leer fichero de picos

args <- commandArgs(trailingOnly = T)
peaks <- args[[1]]
directory <- args[[2]]
length_promotor <- as.numeric(args[[3]])

setwd(directory)

prr5.peaks <- readPeakFile(peakfile = peaks, header=FALSE)
prr5.peaks

## Definir la región que se considera promotor entorno al TSS
promoter <- getPromoters(TxDb=txdb, 
                         upstream=length_promotor, 
                         downstream=length_promotor)
promoter

## Anotación de los picos
prr5.peakAnno <- annotatePeak(peak = prr5.peaks, 
                             tssRegion=c(-1000, 1000),
                             TxDb=txdb)
prr5.peakAnno

plotAnnoPie(prr5.peakAnno)
plotDistToTSS(prr5.peakAnno,
              title="Distribution of genomic loci relative to TSS",
              ylab = "Genomic Loci (%) (5' -> 3')")

## Convertir la anotación a data frame
prr5.annotation <- as.data.frame(prr5.peakAnno)
head(prr5.annotation)

target.genes <- prr5.annotation$geneId[prr5.annotation$annotation == "Promoter"]
length(target.genes)
target.genes <- unique(target.genes)
length(target.genes)

write(x = target.genes,file = "prr5_target_genes.txt")

# Podemos comprobar algunos de estos genes con IGV para comprobar que realmente tienen un 
# pico asociado.

# Para realizar enriquecimientos funcionales de rutas kegg usamos clusterProfiler.

# Con HOMER encontraremos los Transcription Factors Binding Sites.

library("clusterProfiler")

library("org.At.tair.db")



# Leemos los genes desde su fichero.
gene.set<- read.table(file = "prr5_target_genes.txt", header = F, as.is = T)[[1]]
# gene.set

# Necesitaremos el universo del cual realizaremos los enriquecimientos, en nuestro caso
# serán los genes del cromosoma 1 de A. thaliana.

ath.genes<- as.data.frame(genes(txdb))
#ath.genes$seqnames
#ath.genes$seqnames == 1
#ath.genes$gene_id
genes.chr1<- (ath.genes$gene_id)[ath.genes$seqnames == 1]

library("topGO")


# ENRIQUECIMIENTO GO

idType(OrgDb = org.At.tair.db)

# e<- enrichGO(gene = gene.set,
#              OrgDb = org.At.tair.db,
#              keyType = "TAIR",
#              ont = "ALL",
#              pvalueCutoff = 0.05,
#              pAdjustMethod = 0.05,
#              universe = genes.chr1)

# summary(e)

# e@geneSets

# Para realizar los enriquecimientos por separado:
#       BP (biological process),
#       MF (molecular function).
#       CC (celular compartiment).

e.bp<- enrichGO(gene = gene.set,
             OrgDb = org.At.tair.db,
             keyType = "TAIR",
             ont = "BP",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.05,
             universe = genes.chr1)

head(e.bp)

e.bp.table <- as.data.frame(e.bp)

head(e.bp.table)

barplot(e.bp, showCategory=20)



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


# ENRIQUECIMIENTO KEGG

e.kegg<- enrichKEGG(gene = gene.set, 
                    organism = "ath", 
                    universe = genes.chr1)

e.kegg.table <- as.data.frame(e.kegg)

head(e.kegg.table)


library("pathview")

# Preparación de los genes

my.universe<- rep(0, length(genes.chr1))
names(my.universe) <- genes.chr1
#my.universe
my.universe[gene.set] <- 1
#my.universe

# Path.id de cada función (en la función anterior)

kegg.ID <- e.kegg.table$ID
kegg.ID

for (i in 1:length(kegg.ID)){
  kegg.pathview<- pathview(gene.data = my.universe, 
                           pathway.id = kegg.ID[i],
                           species = "ath",
                           gene.idtype = "TAIR")
}

