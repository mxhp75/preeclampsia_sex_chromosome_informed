# Plots for Melissa

# load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)  # For pivot_longer()

# set project directory
projectDir <- "/home/smit1924/preeclampsia_sex_chromosome_informed"
metadata_file <- file.path(projectDir, "clearBoxRawData/cleanSampleMetadata.csv")
# set the output directory
outdir <- file.path(projectDir, "deduplicated_clearBox/output/58_sample_noCovariate")
# Create list of genes to label
genes_to_label <- c("CP", "CDH13", "IGFBP7", "COL6A3", "VCAN", 
                    "CSHL1", "ALPP", "NR3C2", "ELOVL6", "ITGA1", "ROBO1")

## Metadata
metadata <- read_csv(file = metadata_file)

## scatter plot of the male (x-axis) and female (y-axis) log2 fold change ##

# Import male DE table
allTable_male_PE <- read_csv(file.path(outdir, "allTable_male.csv"))

# Import female DE table
allTable_female_PE <- read_csv(file.path(outdir, "allTable_female.csv"))

df <- left_join(allTable_female_PE %>% dplyr::select(-t, -P.Value, -B),
                allTable_male_PE %>% dplyr::select(-t, -P.Value, -B),
                by=c('hgnc_symbol', 'ensembl'),
                suffix=c("_F", "_M"))


# Calculate the maximum absolute value across both logFC columns for symmetric axis scaling
max_abs_fc <- max(abs(c(df$logFC_M, df$logFC_F)), na.rm = TRUE)

# Create the scatter plot
p <- ggplot(df, aes(x = logFC_M, y = logFC_F)) +
  # Add all points in light grey
  geom_point(color = "lightgrey", size = 2, alpha = 0.7) +
  # Add significant points colored by direction of logFC_M
  geom_point(data = subset(df, adj.P.Val_M < 0.05), 
             aes(color = ifelse(logFC_M > 0, "green", "red")),
             size = 2.5) +
  # Add text labels ONLY for significant points that are in genes_to_label vector
  geom_text_repel(
    data = subset(df, adj.P.Val_M < 0.05 & hgnc_symbol %in% genes_to_label),
    aes(label = hgnc_symbol),
    max.overlaps = Inf,
    size = 3,
    box.padding = 0.5,
    point.padding = 0.2,
    segment.color = "grey50"
  ) +
  # Customize colors for significant points
  scale_color_identity() +
  # Set identical limits for both x and y axes and ensure ticks at every 0.5
  scale_x_continuous(limits = c(-max_abs_fc, max_abs_fc), 
                     breaks = seq(floor(-max_abs_fc), ceiling(max_abs_fc), 0.5),
                     minor_breaks = seq(floor(-max_abs_fc), ceiling(max_abs_fc), 0.1)) +
  scale_y_continuous(limits = c(-max_abs_fc, max_abs_fc), 
                     breaks = seq(floor(-max_abs_fc), ceiling(max_abs_fc), 0.5),
                     minor_breaks = seq(floor(-max_abs_fc), ceiling(max_abs_fc), 0.1)) +
  # Add reference lines at x=0 and y=0
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  # Add labels and title
  labs(
    title = "Differential Gene Expression: Male vs Female",
    x = "Log2 Fold Change (Male)",
    y = "Log2 Fold Change (Female)"
  ) +
  # Use a clean theme
  theme_bw() +
  # Use a fixed aspect ratio to ensure visual equality
  coord_fixed(ratio = 1)

# Display the plot
print(p)

# ggsave(filename = file.path(outdir, "melissaPlots/logFC_M_v_logFC_M.jpeg"),
#        create.dir = TRUE,
#        plot = p,
#        units = "in",
#        width = 10,
#        height = 10,
#        dpi = 150)

## violin plots of selected DE genes ##

# make violin plots for the genes Melissa wants to see
# groups (x-axis) male_control, male_PE, female_control, female_PE

# Step 1: Reshape log2CPM from wide to long format

# import the log2CPM file
log2CPM <- read_csv("/home/smit1924/preeclampsia_sex_chromosome_informed/deduplicated_clearBox/output/58_sample_noCovariate/log2CPM.csv")
# convert to long form
log2CPM_long <- log2CPM %>%
  as.data.frame() %>%
  tidyr::pivot_longer(cols = -ensembl, names_to = "samplename", values_to = "Expression")


# Step 2: Merge with metadata to get 'sex_outcome' group
plot_data <- log2CPM_long %>%
  dplyr::left_join(metadata[, c("samplename", "sex_outcome")], by = "samplename") %>%  # Merge using sample names
  dplyr::left_join(., allTable_male_PE[, c("ensembl", "hgnc_symbol")], by = "ensembl") # Merge to include the gene symbol

# Step 3: Filter for the selected gene
selected_gene <- "VCAN"  # Change this for other genes
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

