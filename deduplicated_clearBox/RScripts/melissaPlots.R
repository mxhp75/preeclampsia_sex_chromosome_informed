# Plots for Melissa

# load required libraries
library(ggplot2) # for making plots
library(ggrepel) # for adding the gene annotations
library(patchwork) # for combined plots
library(dplyr) # for general data manipulation
library(tidyr)  # For pivot_longer()
library(readr) # for reading and writing data

# set project directory

projectDir <- "/home/smit1924/preeclampsia_sex_chromosome_informed"
# projectDir <- "D:/VM_Projects/preeclampsia_sex_chromosome_informed"

# set the output directory
outdir <- file.path(projectDir, "deduplicated_clearBox/output/58_sample_noCovariate")

# input files
metadata_file <- file.path(projectDir, "clearBoxRawData/cleanSampleMetadata.csv")
allTable_male_PE_file <- file.path(outdir, "allTable_male.csv")
allTable_female_PE_file <- file.path(outdir, "allTable_female.csv")


# Create list of genes to label
# labels for volcano plots
genes_to_label <- c("CP", "CDH13", "IGFBP7", "COL6A3", "VCAN", 
                    "CSHL1", "ALPP", "NR3C2", "ELOVL6", "ITGA1", "ROBO1", "HSD17B12", "ZFY")


upReg_genes <- c("CP", "GDF7",  "CDH13",  "GPC6",  "PTPN13",  "IGFBP7", "OAF", "PCOLCE",  "COL6A3",  "ANK2",  "DNAJB5",  "TRIL",
                          "SGK2",  "ASRGL1",  "ADAMTS12", "VCAN",  "SLC26A6", "TMEM150C", "NID2", "NRXN3", "OLFML1", "UNC5C", "ELOVL6", "PXDN", 
                          "CHSY3", "KCNK6", "DBN1", "PLCL1", "KREMEN1", "UBXN10", "ITGA1", "ROBO1")

downReg_genes <- c("NR3C2", "CSHL1", "ALPP")

## Metadata
metadata <- read_csv(file = metadata_file)

## scatter plot of the male (x-axis) and female (y-axis) log2 fold change ##

# Import male DE table
allTable_male_PE <- read_csv(allTable_male_PE_file)

# Import female DE table
allTable_female_PE <- read_csv(allTable_female_PE_file)

df <- left_join(allTable_female_PE %>% dplyr::select(-t, -P.Value, -B),
                allTable_male_PE %>% dplyr::select(-t, -P.Value, -B),
                by=c('hgnc_symbol', 'ensembl'),
                suffix=c("_F", "_M"))


# Calculate the maximum absolute value across both logFC columns for symmetric axis scaling
max_abs_fc <- max(abs(c(df$logFC_M, df$logFC_F)), na.rm = TRUE)
max_female_fc <- max(abs(df$logFC_F), na.rm = TRUE)

# Create the scatter plot
p <- ggplot(df, aes(x = logFC_M, y = logFC_F)) +
  # Add all points in light grey
  geom_point(color = "lightgrey", size = 2, alpha = 0.7) +
  # Add significant points colored by direction of logFC_M
  geom_point(data = subset(df, adj.P.Val_M < 0.05), 
             aes(color = ifelse(logFC_M > 0, "red3", "blue")),
             size = 2.5) +
  # Add text labels ONLY for significant points that are in genes_to_label vector
  geom_text_repel(
    data = subset(df, adj.P.Val_M < 0.05 & hgnc_symbol %in% genes_to_label),
    aes(label = hgnc_symbol),
    max.overlaps = Inf,
    size = 3,
    box.padding = 0.5,
    point.padding = 0.2,
    segment.color = "black"
  ) +
  # Customize colors for significant points
  scale_color_identity() +
  # Set identical limits for both x and y axes and ensure ticks at every 0.5
  scale_x_continuous(limits = c(-max_abs_fc, max_abs_fc), 
                     breaks = seq(floor(-max_abs_fc), ceiling(max_abs_fc), 0.5),
                     minor_breaks = seq(floor(-max_abs_fc), ceiling(max_abs_fc), 0.1)) +
  scale_y_continuous(limits = c(-max_female_fc, max_female_fc), 
                     breaks = seq(floor(-max_female_fc), ceiling(max_female_fc), 0.5),
                     minor_breaks = seq(floor(-max_female_fc), ceiling(max_female_fc), 0.1)) +
  # Add reference lines at x=0 and y=0
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  # Add labels and title
  labs(
    title = "Differential Gene Expression: Male vs Female",
    x = expression(log[2]~"Fold Change (Male)"),
    y = expression(log[2]~"Fold Change (Female)")
  ) +
  # Use a clean theme
  theme_bw() +
  # Use a fixed aspect ratio to ensure visual equality
  coord_fixed(ratio = 1)

# Display the plot
print(p)

ggsave(filename = file.path(outdir, "melissaPlots/logFC_M_v_logFC_M.pdf"),
       create.dir = TRUE,
       plot = p,
       units = "in",
       width = 10,
       height = 10,
       dpi = 300)

## violin plots of selected DE genes ##

# make violin plots for the genes Melissa wants to see
# groups (x-axis) male_control, male_PE, female_control, female_PE

# Step 1: Reshape log2CPM from wide to long format

# import the log2CPM file
log2CPM <- read_csv(file.path(projectDir, "/deduplicated_clearBox/output/58_sample_noCovariate/log2CPM.csv"))
# convert to long form
log2CPM_long <- log2CPM %>%
  as.data.frame() %>%
  tidyr::pivot_longer(cols = -ensembl, names_to = "samplename", values_to = "Expression")


# Step 2: Merge with metadata to get 'sex_outcome' group
plot_data <- log2CPM_long %>%
  dplyr::left_join(metadata[, c("samplename", "sex_outcome")], by = "samplename") %>%  # Merge using sample names
  dplyr::left_join(., allTable_male_PE[, c("ensembl", "hgnc_symbol")], by = "ensembl") # Merge to include the gene symbol

my_comparisons <- list(c("F_Control", "F_PE"), c("M_Control", "M_PE"))
# Step 3: Filter for the selected gene
selected_gene <- "CP"  # Change this for other genes
gene_data <- plot_data %>%
  filter(hgnc_symbol == selected_gene)

# Step 4: Filter for samples in "M_PE" group for labeling
label_data <- gene_data %>%
  filter(sex_outcome == "M_PE")

# Step 5: Create Violin Plot for the selected gene
p <- ggplot(gene_data, aes(x = sex_outcome,
                           y = Expression,
                           fill = sex_outcome)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test") + # added this at Seema's suggestion but I don't think it's relevant here
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.5) +  # Add jittered points
  geom_text_repel(data = label_data, aes(label = samplename), 
                  size = 3, nudge_x = 0.2, box.padding = 0.5) +  # Labels for M_PE
  theme_bw() +  # Clean theme
  scale_fill_manual(values = c("blue", "red", "blue", "red")) +
  labs(title = paste("Expression of", selected_gene, "by Sex & Outcome"),
       x = "Sex & Outcome",
       y = "Log2 CPM Expression") +
  theme(legend.position = "none")  # Hide legend if not needed

print(p)

# Step 5: Save the Violin Plot for the selected gene 
# ggsave(filename = file.path(outdir, paste0("melissaPlots/", selected_gene, "_violin_withSamplename.jpeg")),
#        create.dir = TRUE,
#        plot = p,
#        units = "in",
#        width = 10,
#        height = 10,
#        dpi = 150)

## Make the volcano plots for the publication

# Ensure consistent axis limits
max_fc <- max(abs(c(allTable_male_PE$logFC, allTable_female_PE$logFC)))
max_y <- max(-log10(c(allTable_male_PE$P.Value, allTable_female_PE$P.Value)))

# Function to create the "Differential Expression" factor
make_de_factor <- function(df) {
  factor(ifelse(df$adj.P.Val < 0.05,
                ifelse(df$logFC > 0, "Up", "Down"),
                "NS"),
         levels = c("Up", "Down", "NS"))
}

# Define a consistent scale for color
color_scale <- scale_color_manual(
  values = c("Up" = "red3", "Down" = "blue", "NS" = "grey"),
  labels = c("Up" = "Upregulated", "Down" = "Downregulated", "NS" = "Not Significant")
)

# Create the male plot
p1 <- ggplot(allTable_male_PE, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = make_de_factor(allTable_male_PE))) +
  # geom_hline(yintercept = 3.312, linetype = "dashed", color = "grey") +
  geom_text_repel(
    data = subset(allTable_male_PE, hgnc_symbol %in% genes_to_label),
    aes(
      x = logFC,
      y = -log10(P.Value),
      label = hgnc_symbol
      ),
    nudge_x = ifelse(
      subset(allTable_male_PE, hgnc_symbol %in% genes_to_label)$logFC < 0,-1, 1
      ),
    max.overlaps = Inf,
    box.padding = 1.5,
    point.padding = 0.5,
    segment.size = 0.5,
    min.segment.length = 0,
    force = 10,
    force_pull = 0.1,
    seed = 42
    ) +
  color_scale +
  scale_x_continuous(limits = c(-max_fc, max_fc),
                     breaks = seq(-max_fc, max_fc, length.out = 9),
                     # 2 decimal places
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(0, max_y)) +
  labs(x = expression(log[2]~"Fold Change"),
       y = expression(-log[10]~"P-value"),
       title = "Male",
       color = "Differential Expression") +
  theme_bw()

# Create the female plot
p2 <- ggplot(allTable_female_PE, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = make_de_factor(allTable_female_PE))) +
  # geom_hline(yintercept = 3.312, linetype = "dashed", color = "grey") +
  color_scale +
  scale_x_continuous(limits = c(-max_fc, max_fc),
                     breaks = seq(-max_fc, max_fc, length.out = 9),
                     # 2 decimal places
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(0, max_y)) +
  labs(x = expression(log[2]~"Fold Change"),
       y = expression(-log[10]~"P-value"),
       title = "Female") +  # No color legend
  guides(color = "none") + # Suppress the legend
  theme_bw()

# Combine the plots and collect legends
combinedVolcano_plot <- (p1 + p2) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Display the combined plot
print(combinedVolcano_plot)


ggsave(filename = file.path(outdir, paste0("finalManuscriptPlots/", "combinedVolcano.jpeg")),
       create.dir = TRUE,
       plot = combinedVolcano_plot,
       units = "in",
       width = 10,
       height = 10,
       dpi = 300)

# create the figure for the combined volcano plots and the log2FC scatter plot
# Define the layout using area()
design <- c(
  area(1, 1, 1, 1), area(1, 2, 1, 2),  # First row: p1 | p2
  area(2, 1, 3, 1), area(2, 2, 2, 2)   # Second row: p (wider) | empty
)

# Apply layout to plots
combinedVolcanoScatter_plot <- p1 + p2 + p + plot_spacer() +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")

print(combinedVolcanoScatter_plot)

# ggsave(filename = file.path(outdir, paste0("finalManuscriptPlots/", "combinedVolcanoScatter_plot.jpeg")),
#        create.dir = TRUE,
#        plot = combinedVolcanoScatter_plot,
#        units = "in",
#        width = 15,
#        height = 15,
#        dpi = 300)

# combine all of the violin plots into facets
## Up-Regulated

p3 <- ggplot(subset(plot_data, hgnc_symbol %in% upReg_genes), # subset the long gene expression to only include DE GOI
       aes(x = sex_outcome,
           y = Expression,
           fill = sex_outcome)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.5) +  # Add jittered points
  theme_bw() +  # Clean theme
  scale_fill_manual(values = c("blue", "red", "blue", "red")) +
  labs(title = paste("Expression of up-regulated genes by Sex & Outcome"),
       x = "Sex & Outcome",
       y = expression(log[2]~"CPM Expression")) +
  facet_wrap(~hgnc_symbol, scales = "free") +
  theme(legend.position = "none")

## Down Regulated
p4 <- ggplot(subset(plot_data, hgnc_symbol %in% downReg_genes), # subset the long gene expression to only include DE GOI
       aes(x = sex_outcome,
           y = Expression,
           fill = sex_outcome)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.5) +  # Add jittered points
  theme_bw() +  # Clean theme
  scale_fill_manual(values = c("blue", "red", "blue", "red")) +
  labs(title = paste("Expression of down-regulated genes by Sex & Outcome"),
       x = "Sex & Outcome",
       y = expression(log[2]~"CPM Expression")) +
  facet_wrap(~hgnc_symbol, scales = "free") +
  theme(legend.position = "none")

## All DE genes
# set the gene order as a factor to stop facet_wrap() plotting everything alphabetically
plot_data$hgnc_symbol <- factor(plot_data$hgnc_symbol, levels = c(upReg_genes, downReg_genes))

p5 <- ggplot(subset(plot_data, hgnc_symbol %in% c(upReg_genes, downReg_genes)), # subset the long gene expression to only include DE GOI
             aes(x = sex_outcome,
                 y = Expression,
                 fill = sex_outcome)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.5) +  # Add jittered points
  theme_bw() +  # Clean theme
  scale_fill_manual(values = c("blue", "red", "blue", "red")) +
  labs(title = paste("Expression of DE genes by Sex & Outcome"),
       x = "Sex & Outcome",
       y = expression(log[2]~"CPM Expression")) +
  facet_wrap(~hgnc_symbol, scales = "free") +
  theme(legend.position = "none")

# ggsave(filename = file.path(outdir, paste0("finalManuscriptPlots/", "all_DE_violins.jpeg")),
#        create.dir = TRUE,
#        plot = p5,
#        units = "in",
#        width = 30,
#        height = 15,
#        dpi = 300)

# co-ord flip the violins and reduce the x-axis names
## mutate the group names (ie from M_Control to M_Unc and from F_Control to F_Unc)
plot_data <- plot_data %>%
  dplyr::mutate(new_sex_outcome = dplyr::case_when(
    sex_outcome == "M_PE" ~ "M_PE",
    sex_outcome == "F_PE" ~ "F_PE",
    sex_outcome == "M_Control" ~ "M_Un",
    sex_outcome == "F_Control" ~ "F_Un",
    TRUE ~ NA_character_  # This handles any unexpected values
  )) %>%
  dplyr::mutate(., new_sex_outcome = factor(new_sex_outcome, levels = c("F_Un", "F_PE", "M_Un", "M_PE")))

# create the new plot
p6 <- ggplot(subset(plot_data, hgnc_symbol %in% c(upReg_genes, downReg_genes)), # subset the long gene expression to only include DE GOI
             aes(x = new_sex_outcome,
                 y = Expression,
                 fill = new_sex_outcome)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.5) +  # Add jittered points
  # coord_flip() +
  theme_bw() +  # Clean theme
  scale_fill_manual(values = c("palegreen3", "#004000",  "plum3", "#2A0052")) +
  labs(title = paste("Expression of DE genes by Sex & Outcome"),
       x = "Sex & Outcome",
       y = expression(log[2]~"CPM Expression")) +
  facet_wrap(~hgnc_symbol, ncol = 7, scales = "free") +
  theme(legend.position = "none")

ggsave(filename = file.path(outdir, paste0("finalManuscriptPlots/", "all_DE_violins.jpeg")),
       create.dir = TRUE,
       plot = p6,
       units = "in",
       width = 30,
       height = 15,
       dpi = 300)
