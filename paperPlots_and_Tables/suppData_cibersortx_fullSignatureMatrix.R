# Load required libraries
library(dplyr)
library(tibble)
library(pheatmap)

# Set project directory and output directory
project_dir <- "D:/VM_Projects/preeclampsia_sex_chromosome_informed"
outdir <- file.path(project_dir, "cibersortx_de_output")

# Load the CIBERSORTx signature matrix
cibersortx_sigMatrix_file <- file.path(project_dir, "cibersortFullSigMatrix/CIBERSORTx_sigmatrix_Adjusted.txt")
cibersortx_sigMatrix <- read.table(file = cibersortx_sigMatrix_file, header = TRUE, sep = "\t")

# Remove rows with all zero counts
cibersortx_sigMatrix <- cibersortx_sigMatrix %>%
  filter(rowSums(across(-1)) > 0)

# Generate heatmap
pdf(
  file = file.path(outdir, "heatmap_cibersortx_signature_matrix.pdf"),
  width = 12,
  height = 8
)

cibersortx_sigMatrix %>%
  `rownames<-`(NULL) %>%
  as.data.frame() %>%
  column_to_rownames("GeneSymbol") %>%
  pheatmap(
    scale = "row",
    treeheight_row = 0,
    show_rownames = FALSE,
    main = "Cibersortx Adjusted Signature Matrix (scaled by row)"
  )

dev.off()
