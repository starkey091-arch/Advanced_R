
R version 4.4.2 (2024-10-31) -- "Pile of Leaves"
Copyright (C) 2024 The R Foundation for Statistical Computing
Platform: aarch64-apple-darwin20

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

[R.app GUI 1.81 (8462) aarch64-apple-darwin20]

[Workspace restored from /Users/connerstarkey/.RData]
[History restored from /Users/connerstarkey/.Rapp.history]

> ##############################################
> # Gene Expression Analysis – B473/B573
> ##############################################
> 
> library(readxl)
> library(dplyr)

Attaching package: ‘dplyr’

The following objects are masked from ‘package:stats’:

    filter, lag

The following objects are masked from ‘package:base’:

    intersect, setdiff, setequal, union

> library(tidyr)
> library(ggplot2)
> 
> # Helper for loading files safely
> load_file_or_stop <- function(path, fun, ...) {
+   if (!file.exists(path)) stop(paste("File not found:", path))
+   fun(path, ...)
+ }
> 
> expr        <- load_file_or_stop("~/Downloads/Gene_Expression_Data.xlsx", read_excel)
                                                                                                                                                                          
> gene_info   <- load_file_or_stop("~/Downloads/Gene_Information.csv", read.csv)
> sample_info <- load_file_or_stop("~/Downloads/Sample_Information.tsv", read.table,
+                                  sep = "\t", header = TRUE)
> 
> ### Rename sample columns using phenotype + patient ###
> sample_info$patient_id <- sub("patient: ", "", sample_info$patient)
> sample_info$label <- paste0(sample_info$group, "_", sample_info$patient_id)
> rename_map <- setNames(sample_info$label, rownames(sample_info))
> 
> expr_renamed <- expr
> gsm_cols <- intersect(colnames(expr_renamed), names(rename_map))
> colnames(expr_renamed)[match(gsm_cols, colnames(expr_renamed))] <- rename_map[gsm_cols]
> 
> ### Split into Tumor and Normal ###
> tumor_cols  <- sample_info$label[sample_info$group == "tumor"]
> normal_cols <- sample_info$label[sample_info$group == "normal"]
> 
> tumor_data  <- expr_renamed[, c("Probe_ID", tumor_cols)]
> normal_data <- expr_renamed[, c("Probe_ID", normal_cols)]
> 
> ### Compute means ###
> tumor_means <- tumor_data %>%
+   mutate(Tumor_mean = rowMeans(across(-Probe_ID))) %>%
+   select(Probe_ID, Tumor_mean)
> 
> normal_means <- normal_data %>%
+   mutate(Normal_mean = rowMeans(across(-Probe_ID))) %>%
+   select(Probe_ID, Normal_mean)
> 
> merged_means <- inner_join(tumor_means, normal_means, by = "Probe_ID")
> 
> ### Compute log2 fold change (assignment formula) ###
> merged_means <- merged_means %>%
+   mutate(log2FC = log2((Tumor_mean - Normal_mean) / Normal_mean)) %>%
+   filter(is.finite(log2FC))
Warning message:
There was 1 warning in `mutate()`.
ℹ In argument: `log2FC = log2((Tumor_mean - Normal_mean)/Normal_mean)`.
Caused by warning:
! NaNs produced 
> 
> ### Merge with gene info and filter DEGs ###
> full_results <- inner_join(merged_means, gene_info, by = "Probe_ID")
> DEGs <- full_results %>% filter(abs(log2FC) > 5)
> 
> ### Add direction label ###
> DEGs <- DEGs %>% mutate(Higher_in = ifelse(log2FC > 0, "Tumor", "Normal"))
> 
> ### EDA Plots ###
> ggplot(DEGs, aes(x = Chromosome)) +
+   geom_bar(fill = "steelblue") +
+   theme_minimal() +
+   labs(title = "Distribution of DEGs by Chromosome",
+        x = "Chromosome", y = "Number of DEGs")
> 
> ggplot(DEGs, aes(x = Chromosome, fill = Higher_in)) +
+   geom_bar(position = "dodge") +
+   theme_minimal() +
+   labs(title = "DEGs by Chromosome and Direction",
+        x = "Chromosome", y = "Number of DEGs")
> 
> pct_data <- DEGs %>%
+   count(Higher_in) %>%
+   mutate(Percent = n / sum(n) * 100)
> 
> ggplot(pct_data, aes(x = Higher_in, y = Percent, fill = Higher_in)) +
+   geom_col() +
+   theme_minimal() +
+   labs(title = "Percent of DEGs Upregulated in Tumor vs Normal",
+        x = "Higher Expression In", y = "Percent")
> 
> ### Heatmap + Clustermap ###
> expr_mat <- as.data.frame(expr_renamed)
> rownames(expr_mat) <- expr_mat$Probe_ID
> expr_mat$Probe_ID <- NULL
> expr_mat <- as.matrix(expr_mat)
> scaled_expr <- scale(expr_mat)
> 
> heatmap(scaled_expr, Rowv = NA, Colv = NA,
+         labRow = FALSE, main = "Gene Expression Heatmap")
> 
> heatmap(scaled_expr, labRow = FALSE,
+         main = "Gene Expression Clustermap (Clustered)")
Error: vector memory limit of 16.0 Gb reached, see mem.maxVSize()
> 
> 