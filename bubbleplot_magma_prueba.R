library(DESeq2)
library(GOplot)
library(stringr)
library(ggrepel)

input <- 'Resultados_DEG/Sin_Familia/Juntos/DEG_results_sinFamilia_sinEdad/all_genes.txt'
res <- read.csv(input, header = TRUE, sep = "\t")

colores <- c('red', 'green', 'blue', 'orange', 'yellow')

allgenes <- as.data.frame(res$log2FoldChange)
rownames(allgenes) <- res$X 
colnames(allgenes) <- 'logFC'
allgenes$gene <- res$X

#Categories validated from GSEA
circ0 <- read.csv('tabla_juntos_bubble_magma.txt', header = T, sep="\t")
#Aqui falta calcular el z score
#DAta frame intermedio 'cats'
cats <- as.data.frame(circ0$term)
cats$gene <- circ0$gene
cats$set_size <- circ0$count
colnames(cats) <- c('term', 'gene', 'set_size')
#empty vector to store z scores
z_scores <- numeric(nrow(cats))

#Hasta la línea 58, me funciona con tus datos; tuve que cambiar el símbolo por el que
#separa los genes
# Iterate through each row of 'cats' and calculate the Z-score
for (i in 1:nrow(cats)) {
  # Extract the list of genes in the current term and split it into a vector
  genes_in_term <- unlist(strsplit(cats$gene[i], "/"))
  
  # Filter 'allgenes' to get log2 fold changes for genes in the current term
  log2foldchanges_in_term <- allgenes$logFC[allgenes$gene %in% genes_in_term]
  
  # Calculate the Z-score for the term
  upregulated_count <- sum(log2foldchanges_in_term > 0)
  downregulated_count <- sum(log2foldchanges_in_term < 0)
  set_size <- length(log2foldchanges_in_term)
  
  z_score <- (upregulated_count - downregulated_count) / sqrt(set_size)
  
  # Store the calculated Z-score in the vector
  z_scores[i] <- z_score
}

# Add the Z-scores to the 'cats' data frame
cats$ZScore <- z_scores

#Relate to the circ data frame
circ0$zscore <- cats$ZScore[cats$term == circ0$term]
circ <- circ0

#Todos tus datos están en el objeto circ
#hago objeto circ2 intermedio para que no se sobreescriban los datos buenos.
circ2 <- circ
circ2$adj_pval <- circ2$magma_score

trace(GOBubble, edit = T)

jpeg(file = 'bubble_prueba1J.jpeg', units = 'in', width = 15, height = 20, res = 300)
par(mar = c(2, 2, 2, 5)) 
#Cambia 'labels' para que se vean más o menos categorías
plot <- GOBubble(circ2, labels = 2, colour = colores, ID =T, table.col=F, table.legend = F)

# Agregar la función geom_text_repel para evitar superposiciones de nombres
plot + geom_text_repel(data = circ2, aes(label = gene), size = 3, nudge_x = 0.765)

invisible(dev.off())
