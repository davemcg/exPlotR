# make tiny_rse.Rdata

library(recount3)
library(geyser)
library(dplyr)
human_projects <- available_projects()
proj_info <- subset(
  human_projects,
  project == "SRP107937" & project_type == "data_sources"
)
tiny_rse <- create_rse(proj_info)
assay(tiny_rse, "counts") <- transform_counts(tiny_rse)
# first tweak that glues the gene name onto the gene id in the row names
rownames(tiny_rse) <- paste0(rowData(tiny_rse)$gene_name, ' (', row.names(tiny_rse), ')')
# creates two new metadataa fields 
colData(tiny_rse)$tissue <- colData(tiny_rse)$sra.sample_title %>% stringr::str_extract(.,'PRC|PR')
colData(tiny_rse)$disease <- colData(tiny_rse)$sra.sample_title %>% stringr::str_extract(.,'AMD|Normal')
# cut down to a small subset of genes to make the object a reasonable size for an R package
genes <- scan('inst/extdata/tiny_rse_genes.tsv.gz', what = 'character', sep = '\n')
tiny_rse <- tiny_rse[genes,]

save(tiny_rse, file = 'inst/extdata/tiny_rse.Rdata')
