## Load recount3 R package
library("recount3")

## Revisemos todos los proyectos con datos de humano en recount3
mouse_projects <- available_projects(organism = "mouse")

## Encuentra tu proyecto de interés. Aquí usaremos
## SRP151148 de ejemplo
proj_info <- subset(
  mouse_projects,
  project == "SRP151148" & project_type == "data_sources"
)
## Crea un objeto de tipo RangedSummarizedExperiment (RSE)
## con la información a nivel de genes
rse_gene_SRP151148 <- create_rse(proj_info)

## Convirtamos las cuentas por nucleotido a cuentas por lectura
## usando compute_read_counts().
## Para otras transformaciones como RPKM y TPM, revisa transform_counts().
assay(rse_gene_SRP151148, "counts") <- compute_read_counts(rse_gene_SRP151148)

## Para este estudio en específico, hagamos más fácil de usar la
## información del experimento
rse_gene_SRP151148 <- expand_sra_attributes(rse_gene_SRP151148)

colData(rse_gene_SRP151148)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP151148)))
]

#Visualisación, posiblemente se tenga que eliminar 
#iSEE::iSEE(rse_gene_SRP151148)

#Poniendo en el formato correcto las varibles de interés
#Este podría no interesar, ya que la mayoría son none
rse_gene_SRP151148$sra_attribute.in_vitro_treatment <- factor(tolower(rse_gene_SRP151148$sra_attribute.in_vitro_treatment))
rse_gene_SRP151148$sra_attribute.cell_type_short <- factor(tolower(rse_gene_SRP151148$sra_attribute.cell_type_short))
rse_gene_SRP151148$sra_attribute.tissue <- factor(tolower(rse_gene_SRP151148$sra_attribute.tissue))

summary(as.data.frame(colData(rse_gene_SRP151148)[
  ,
  grepl("^sra_attribute.[in_vitro_treatment|cell_type_short|tissue]", colnames(colData(rse_gene_SRP151148)))
]))

rse_gene_SRP151148$assigned_gene_prop <- rse_gene_SRP151148$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP151148$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP151148$assigned_gene_prop)

rse_gene_SRP151148$gene_mean

## Guardemos nuestro objeto entero por si luego cambiamos de opinión
rse_gene_SRP151148_unfiltered <- rse_gene_SRP151148

## Eliminemos a muestras malas
#hist(rse_gene_SRP151148$assigned_gene_prop)

rse_gene_SRP151148 <- rse_gene_SRP151148[, rse_gene_SRP151148$assigned_gene_prop > 0.5]

gene_means <- rowMeans(assay(rse_gene_SRP151148, "counts"))
summary(gene_means)

rse_gene_SRP151148 <- rse_gene_SRP151148[gene_means > 0.1, ]


#Porcertaje de genes que conservamos

round(nrow(rse_gene_SRP151148) / nrow(rse_gene_SRP151148_unfiltered) * 100, 2)
#49.55

## ----normalize------------------------------------------------
library("edgeR") # BiocManager::install("edgeR", update = FALSE)
dge <- DGEList(
  counts = assay(rse_gene_SRP151148, "counts"),
  genes = rowData(rse_gene_SRP151148)
)
dge <- calcNormFactors(dge)

## ----explore_gene_prop_by_age---------------------------------
library("ggplot2")
ggplot(as.data.frame(colData(rse_gene_SRP151148)), aes(y = assigned_gene_prop, x = sra_attribute.cell_type_short)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Cell_type_short")


## ----statiscal_model------------------------------------------
mod <- model.matrix(~ sra_attribute.cell_type_short + sra_attribute.in_vitro_treatment + 
                      rse_gene_SRP151148$sra_attribute.tissue + assigned_gene_prop,
                    data = colData(rse_gene_SRP151148)
)
colnames(mod)


## ----run_limma------------------------------------------------
library("limma")
vGene <- voom(dge, mod, plot = TRUE)

eb_results <- eBayes(lmFit(vGene))

de_results <- topTable(
  eb_results,
  coef = 2,
  number = nrow(rse_gene_SRP151148),
  sort.by = "none"
)
dim(de_results)
head(de_results)

## Genes diferencialmente expresados entre cell_type_short con FDR < 5%
table(de_results$adj.P.Val < 0.05)

## Visualicemos los resultados estadísticos
plotMA(eb_results, coef = 2)

volcanoplot(eb_results, coef = 2, highlight = 3, names = de_results$gene_name)
de_results[de_results$gene_name %in% c("ZSCAN2", "VASH2", "KIAA0922"), ]


## ----pheatmap-------------------------------------------------
## Extraer valores de los genes de interés
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 50, ]

## Creemos una tabla con información de las muestras
## y con nombres de columnas más amigables
df <- as.data.frame(colData(rse_gene_SRP151148)[, c("sra_attribute.cell_type_short", "sra_attribute.in_vitro_treatment",
                                                    "sra_attribute.tissue")])
colnames(df) <- c("cell_type_short", "in_vitro_treatment","tissue")

## Hagamos un heatmap
library("pheatmap")
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = df
)


## ----plot_mds-------------------------------------------------
## Para colores
library("RColorBrewer")

## Conviertiendo los grupos de tissue a colores
col.tissue <- df$tissue
levels(col.tissue) <- brewer.pal(nlevels(col.tissue), "Paired")
col.tissue <- as.character(col.tissue)

## MDS por grupos tissue
plotMDS(vGene$E, labels = df$tissue, col = col.tissue)

## Conviertiendo los valores de sex a colores
col.in_vitro_treatment <- df$in_vitro_treatment
levels(col.in_vitro_treatment) <- brewer.pal(nlevels(col.in_vitro_treatment), "Dark2")
col.in_vitro_treatment <- as.character(col.in_vitro_treatment)

## MDS por in vitro treatment
plotMDS(vGene$E, labels = df$in_vitro_treatment, col = col.in_vitro_treatment)


## ----respuesta, out.height="1100px"---------------------------
## Tenemos que usar gene_id y gene_name
rowRanges(rse_gene_SRP151148)

## Alternativamente, podriamos haber usado de_results
head(de_results, n = 3)

## Es la misma información
identical(rowRanges(rse_gene_SRP151148)$gene_id, de_results$gene_id)
identical(rowRanges(rse_gene_SRP151148)$gene_name, de_results$gene_name)

## Guardemos los IDs de nuestros genes
nombres_originales <- rownames(exprs_heatmap)

## Con match() podemos encontrar cual es cual
rownames(exprs_heatmap) <- rowRanges(rse_gene_SRP151148)$gene_name[
  match(rownames(exprs_heatmap), rowRanges(rse_gene_SRP151148)$gene_id)
]

## Vean que tambien podriamos haber usado rank()
identical(
  which(rank(de_results$adj.P.Val) <= 50),
  match(nombres_originales, rowRanges(rse_gene_SRP151148)$gene_id)
)

## Esta es otra solución
identical(
  de_results$gene_name[rank(de_results$adj.P.Val) <= 50],
  rownames(exprs_heatmap)
)

## Por último podemos cambiar el valor de show_rownames de FALSE a TRUE
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = df
)

## Guardar la imagen en un PDF largo para poder ver los nombres de los genes
pdf("pheatmap_con_nombres.pdf", height = 14, useDingbats = FALSE)
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = df
)
dev.off()


## ----"centered_and_scaled", out.height="1100px"---------------
## Versión con centering y scaling en los renglones (los genes)
pheatmap::pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = df,
  scale = "row"
)


## ----"complexheatmap", out.height="1100px"--------------------
## Misma versión pero ahora con ComplexHeatmap en vez del paquete pheatmap
ComplexHeatmap::pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = df,
  scale = "row"
)

