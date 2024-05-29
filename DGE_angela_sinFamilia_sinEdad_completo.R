#DGE
#Run a Differentially Expressed Genes Analysis using as input the count files 
#produced in previous step

#Load libraries
suppressPackageStartupMessages({
  library(gplots, quietly = T)
  library(ggplot2, quietly = TRUE)
  library(pheatmap, quietly = TRUE)
  library(DESeq2, quietly = TRUE)
  library('org.Rn.eg.db', quietly = TRUE, character.only = TRUE)
  library(EnhancedVolcano, quietly = TRUE)
})
require("ggrepel", quietly = TRUE)

#Paths
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

#Input
#In order to select the right folder we need to specify which group we are going to work with
divide_sex <- F
males <- F

if (divide_sex) {
  if (males) {
    suffix <- '_males'
    
  } else {
    suffix <- '_females'
  }
  
} else {
  suffix <- ''
}

configFile <- paste0('Archivo_configuracion_mPFC', suffix,'.txt')

#Outputs
resD <- paste0('DEG_results', suffix, '_sinFamilia_prueba2/')

dir.create(resD)

rawCountsF <- paste0(resD,"counts_raw.tsv")
normCountsF <- paste0(resD,"counts_normalized.tsv")
PCAF <- paste0(resD,"PCA.jpeg")
distancesF <- paste0(resD,"distances.jpeg")
dispersionF <- paste0(resD,"dispersion.tiff")
MAplotF <- paste0(resD,"maplot.jpeg")
genesTSV <- paste0(resD,"all_genes.csv")
sigTSV <- paste0(resD,"sig_pval.csv")
sigPCAF <- paste0(resD,"PCA_sig.jpeg")
volcanoF <- paste0(resD,"volcanoPlot.jpeg")
heatmapF <- paste0(resD,"heatmap.jpeg")
DESEqResultsF <- paste0(resD, 'deseq_objects.RData')

#Parameters
cutoff <- 0.05 #significancy p value adjusted
FCthres <- 4 #fold change threshold to consider in graphs

#Functions: italics for genes in heatmap
make_italics <- function(x){
  as.expression(lapply(x, function(y) bquote(italic(.(y)))))
}

#We create a function so that when the name of the gene is NA it gives back de ensembl ID
rellenar_nombres <- function(data_frame) {
  # Obtener los rownames del data frame
  row_names <- rownames(data_frame)
  
  # Obtener la columna de nombres de genes
  genes_col <- data_frame$symbol
  
  # Reemplazar los valores NA con los rownames correspondientes
  genes_col[is.na(genes_col)] <- row_names[is.na(genes_col)]
  
  # Actualizar la columna en el data frame
  data_frame$symbol <- genes_col
  
  return(data_frame)
}

#Load variables of each sample
sampleTable <- read.table(configFile, header=TRUE
                          ,colClasses= c('factor','character','factor',
                                         'factor','factor')
)

#Convert the counts into a DeSeq DataSet object
if (divide_sex) {
  data <- DESeqDataSetFromHTSeqCount(sampleTable, directory=".", 
                                     design = ~ Grupo)
} else {
  data <- DESeqDataSetFromHTSeqCount(sampleTable, directory=".", design = ~ Sexo + Grupo)
}

## Analysis
#Pre-filtering: clean some of the noise in the counts
keep <- rowSums(counts(data)) >= 10
data <- data[keep,]
#With this filter, the object goes from 30562 elements to 19538 elements

# DESeq: original DESEQ() function doesnt allow to adjust number of iterations
dds <- estimateSizeFactors(data)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, maxit = 10000)

#Save raw counts of all the samples in a single data frame to save it, as well as normalized counts
dds_raw <- counts(dds, normalized=FALSE)
dds_normalized <- counts(dds, normalized=TRUE)

#save file with counts and normalized counts
write.table(dds_raw, file=rawCountsF, quote=FALSE, 
            sep = "\t", col.names=NA)
write.table(dds_normalized, file=normCountsF, quote=FALSE, 
            sep = "\t", col.names=NA)
#Save for the GSEA
save(dds, file = DESEqResultsF)

#PCA: blind must be FALSE to take into account batch effect
#Normalization
vst <- varianceStabilizingTransformation(dds, blind = FALSE)
mat <- assay(vst)
if (!divide_sex) {
  mat <- limma::removeBatchEffect(mat, batch=vst$Sexo, 
                                  group=vst$Grupo) }
#Si no es dividido por sexo aparece la covariable del sexo

assay(vst) <- mat

jpeg(filename = PCAF, width=900, height=900, quality=300)
pca <- plotPCA(vst, intgroup = c("Grupo", "Sexo"), returnData = TRUE)
title <- "Principal Components Plot"
ggplot(pca, aes(x=PC1, y=PC2, color = Grupo, shape= Sexo)) +
  ggtitle(title) + 
  geom_point(size = 6, aes(fill=Grupo)) + scale_shape_manual(values=c(21,22)) +
  theme(plot.title = element_text(size=40, hjust = 0.5, face = "bold"), axis.title=element_text(size=20),
        legend.text=element_text(size=15),legend.title=element_text(size=15)) +
  geom_text_repel(aes(label=colnames(vst)), size=5, point.padding = 0.6)
dev.off()

#Plot distances between samples

distRL <- dist(t(mat))
distMat <- as.matrix(distRL)
hc <- hclust(distRL)
hmcol <- colorRampPalette(c("white", "blue"))(299)
title <- "Distances matrix"
jpeg(filename = distancesF, width=900, height=900, quality=300)
heatmap.2(distMat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col=rev(hmcol),
          margin=c(10, 6), main=title, key.title=NA)
invisible(dev.off())

#Plot gene dispersion
tiff(filename = dispersionF, units="in", width=5, height=5, res=300)
title <- "Per-gene dispersion estimates"
plotDispEsts(dds, main=title)
invisible(dev.off())


#Get factor levels
levels <- unique(sampleTable$Grupo)
l1 <- toString(levels[2]) #reference level has to be control
l2 <- toString(levels[1])
suffix <- paste(l1, l2, sep="_vs_")

#Get results
res <- results(dds, contrast=c("Grupo", l2, l1))
res$FoldChange <- 2^res$log2FoldChange  #have actual fold change
res <- res[colnames(res)[c(1,7,2:6)]] # order columns

# MAplot
jpeg(filename = MAplotF, width=900, height=900, quality=300)
title <- paste("MA-plot", suffix, sep=" ")
plotMA(res, alpha=cutoff, main=title, colSig = 'red', cex = 1.2)
invisible(dev.off())

#Add annotation to results
symbol <- mapIds(get('org.Rn.eg.db'), keys=row.names(res), column="SYMBOL", 
                 keytype="ENSEMBL", multiVals="first") #to obtain gene symbols

description <- mapIds(get('org.Rn.eg.db'), keys=row.names(res),
                      column="GENENAME", keytype="ENSEMBL", 
                      multiVals="first") #to obtain description

#Bind dataframe symbol and res and create a new one for description
res <- cbind(symbol, res)
res$description <- description

write.table(res, file=genesTSV, quote=FALSE, sep=",", col.names=NA)
sig_pval <- subset(res, res$pvalue < cutoff) #Select only the roles with a p-value smaller than cutoff
write.table (sig_pval, file=sigTSV, quote=FALSE, sep=",", col.names=NA)

#Volcano plot
discardNA <- !is.na(res$padj) #Remove padj which have NA due to the independent filtering
res2 <- res[discardNA,]
remove_outliers<- (res2$log2FoldChange >= -FCthres) & (res2$log2FoldChange <= FCthres) 
res2 <- res2[remove_outliers,]

# Aquí añades la función rellenar_nombres()
rellenar_nombres <- function(res) {
  for (i in seq_along(res$symbol)) {
    if (grepl("^ENS", res$symbol[i])) {
      res$symbol[i] <- "" 
    } else {
      res$symbol[i] <- res$symbol[i]
    }
  }
  return(res)
}

# Aplicas la función al dataframe res2
res2 <- rellenar_nombres(res2)

jpeg(filename = volcanoF, units="in", width=8, height=10, res=300)
EnhancedVolcano(res2, lab = res2$symbol, x = 'log2FoldChange', y = 'pvalue',
                pCutoff = 0.000976, FCcutoff= 0.3, 
                ylim = c(0, 11), 
                labSize = 5,
                legendLabSize = 9, legendIconSize = 5, drawConnectors = TRUE,
                widthConnectors = 0.5, max.overlaps = 50, title = '', arrowheads = FALSE,
                subtitle= '', gridlines.major = FALSE, gridlines.minor = FALSE)
invisible(dev.off())

# Filtrar los valores en res2 que cumplen con las condiciones dadas
filtros <- res2$log2FoldChange > 0.3 | res2$log2FoldChange < -0.3
filtros <- filtros & res2$pvalue < 0.000976

# Subconjunto de datos que cumple con los filtros
res_filtrado <- res2[filtros,]

# Número de valores que cumplen con los filtros
num_valores_filtrados <- nrow(res_filtrado)

# Imprimir el resultado
print(num_valores_filtrados)


#Get most significant genes according to cut off
significant <- subset(res, res$pvalue < cutoff)
significant <- significant[order(significant$pvalue),]
#Discard those genes with unbelievable Fold Change (outliers)
#significant <- significant[(significant$log2FoldChange >= -FCthres) & 
#(significant$log2FoldChange <= FCthres),]

#PCA only significant genes
interest_genes <- rownames(significant)
dds_sig <- dds[interest_genes,]
vst_sig <- varianceStabilizingTransformation(dds_sig, blind = FALSE)
mat_sig <- assay(vst_sig)
if (!divide_sex) {
  mat_sig <- limma::removeBatchEffect(mat_sig, batch=vst$Sexo, group=vst$Grupo)
  
}
assay(vst_sig) <- mat_sig

jpeg(filename = sigPCAF, width=900, height=900, quality=300)
pca_sig <- plotPCA(vst_sig, intgroup='Grupo')
title <- "PCA - only significant genes"
pca_sig + ggtitle(title) + 
  geom_point(size = 6) +
  theme(plot.title = element_text(size=40, hjust = 0.5, face = "bold"), 
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),legend.title=element_text(size=15)) +
  geom_text_repel(aes(label=colnames(vst_sig)), size=5, point.padding = 0.6)
invisible(dev.off())

#Heatmap
# Ordenar los genes significativos por sus valores de padj de menor a mayor
significant <- significant[order(significant$padj), ]

# Seleccionar los top 10 genes significativos
top_10_genes <- head(significant, 10)

# Subset de la matriz original utilizando solo los genes top 10
subcounts <- subset(mat, rownames(mat) %in% rownames(top_10_genes))

# Transformar la matriz de expresión a escala logarítmica
lsubcounts <- log2(subcounts + 1)

# Obtener los nombres de los genes seleccionados
top_10_genes <- rellenar_nombres(top_10_genes)
sig_symbol <- as.character(top_10_genes$symbol)

# Definir las condiciones de agrupación (l1 y l2)
conditions <- c(l1, l2)
conds <- subset(sampleTable, sampleTable$Grupo %in% conditions)
samples <- conds$Rata

# Crear un dataframe con la información de las condiciones
df <- data.frame(condition = conds$Grupo)
rownames(df) <- samples

# Definir colores para las condiciones
my_colour <- list(df = c(l1 = "skyblue", l2 = "orange"))

# Título del heatmap
title <- "Heatmap of top 10 significant genes"

# Generar el heatmap
jpeg(filename = heatmapF, units = "in", width = 8, height = 5, res = 300)
pheatmap(lsubcounts, scale = "row", cluster_rows = TRUE,
         cluster_cols = TRUE, legend = TRUE, drop_levels = TRUE, 
         labels_row = make_italics(sig_symbol), 
         main = title,
         annotation_col = df, annotation_colors = my_colour,
         treeheight_row = 30, treeheight_col = 20)
invisible(dev.off())
