> ##############################################
> # Assignment – Gene Expression Analysis in R
> ##############################################
> 
> # Load libraries -----------------------------------------------------------
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
> # 1. LOAD THE DATA --------------------------------------------------------
> 
> expr <- read_excel("~/Downloads/Gene_Expression_Data.xlsx")
                                                                                                                                                                                                                       
> gene_info <- read.csv("~/Downloads/Gene_Information.csv")
> sample_info <- read.table(
+   "~/Downloads/Sample_Information.tsv",
+   sep = "\t",
+   header = TRUE
+ )
> 
> # Quick peek
> head(expr)
# A tibble: 6 × 19
  Probe_ID     GSM820516 GSM820517 GSM820518 GSM820519 GSM820520 GSM820521 GSM820522 GSM820523 GSM820524 GSM820525 GSM820526 GSM820527 GSM820528 GSM820529 GSM820530 GSM820531 GSM820532 GSM820533
  <chr>            <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>
1 ILMN_1343291    66665.    69405.    64128.    68944.    67827.    71775.    62246.    69714.    69509.    68245.    65427.    68437.    57609.    69960.    69509.    70064.    69647.    70332.
2 ILMN_1343295    22040.    13046.    38679.    16642.    33720.    18933.    26170.     9907.    17167.    12429.    25298.    17536.    19749.    17854.    43671.    22849.    23726.    28747.
3 ILMN_1651199      227.      205.      217.      229.      226.      204.      213.      210.      230.      213.      226.      232.      209.      229.      214.      217.      196.      252.
4 ILMN_1651209      279.      254.      212.      278.      260.      265.      321.      273.      254.      310.      275.      275.      251.      256.      220.      292.      253.      238.
5 ILMN_1651210      195.      196.      175.      194.      230.      164.      245.      191.      188.      199.      221.      213.      195.      174.      185.      175.      195.      192.
6 ILMN_1651221      217.      206.      220.      205.      195.      234.      233.      222.      244.      210.      228.      256.      234.      194.      219.      271.      223.      241.
> head(gene_info)
      Probe_ID    Symbol Entrez_Gene_ID Chromosome  Cytoband
1 ILMN_1343291    EEF1A1           1915          6     6q13c
2 ILMN_1343295     GAPDH           2597         12 12p13.31d
3 ILMN_1651199 LOC643334         643334              2q37.3b
4 ILMN_1651209   SLC35E2           9906          1  1p36.33a
5 ILMN_1651210    DUSP22          56940              6p25.3b
6 ILMN_1651221 LOC642820         642820         10          
> head(sample_info)
           group    patient
GSM820516  tumor patient: 1
GSM820517 normal patient: 1
GSM820518  tumor patient: 2
GSM820519 normal patient: 2
GSM820520  tumor patient: 3
GSM820521 normal patient: 3
> # expr:   Probe_ID + GSM* columns
> # sample: rownames = GSM IDs, columns: group (tumor/normal), patient ("patient: 1")
> 
> 
> # 2. RENAME SAMPLE NAMES BASED ON PHENOTYPE -------------------------------
> 
> # Make patient IDs nice (1, 2, 3, ...)
> sample_info$patient_id <- sub("patient: ", "", sample_info$patient)
> 
> # Build labels like "tumor_1", "normal_1"
> sample_info$label <- paste0(sample_info$group, "_", sample_info$patient_id)
> 
> # Map GSM IDs (rownames) to new labels
> rename_map <- setNames(sample_info$label, rownames(sample_info))
> 
> # Apply renaming to expression matrix
> expr_renamed <- expr
> gsm_cols <- intersect(colnames(expr_renamed), names(rename_map))
> colnames(expr_renamed)[match(gsm_cols, colnames(expr_renamed))] <-
+   rename_map[gsm_cols]
> 
> # You can show this in your report as proof of renaming:
> head(expr_renamed)
# A tibble: 6 × 19
  Probe_ID     tumor_1 normal_1 tumor_2 normal_2 tumor_3 normal_3 tumor_4 normal_4 tumor_5 normal_5 tumor_6 normal_6 tumor_7 normal_7 tumor_8 normal_8 tumor_9 normal_9
  <chr>          <dbl>    <dbl>   <dbl>    <dbl>   <dbl>    <dbl>   <dbl>    <dbl>   <dbl>    <dbl>   <dbl>    <dbl>   <dbl>    <dbl>   <dbl>    <dbl>   <dbl>    <dbl>
1 ILMN_1343291  66665.   69405.  64128.   68944.  67827.   71775.  62246.   69714.  69509.   68245.  65427.   68437.  57609.   69960.  69509.   70064.  69647.   70332.
2 ILMN_1343295  22040.   13046.  38679.   16642.  33720.   18933.  26170.    9907.  17167.   12429.  25298.   17536.  19749.   17854.  43671.   22849.  23726.   28747.
3 ILMN_1651199    227.     205.    217.     229.    226.     204.    213.     210.    230.     213.    226.     232.    209.     229.    214.     217.    196.     252.
4 ILMN_1651209    279.     254.    212.     278.    260.     265.    321.     273.    254.     310.    275.     275.    251.     256.    220.     292.    253.     238.
5 ILMN_1651210    195.     196.    175.     194.    230.     164.    245.     191.    188.     199.    221.     213.    195.     174.    185.     175.    195.     192.
6 ILMN_1651221    217.     206.    220.     205.    195.     234.    233.     222.    244.     210.    228.     256.    234.     194.    219.     271.    223.     241.
> 
> 
> # 3. SPLIT INTO TUMOR AND NORMAL GROUPS -----------------------------------
> 
> # After renaming, tumor columns start with "tumor_", normal with "normal_"
> tumor_cols  <- sample_info$label[sample_info$group == "tumor"]
> normal_cols <- sample_info$label[sample_info$group == "normal"]
> 
> # Subset the renamed expression data
> tumor_data <- expr_renamed[, c("Probe_ID", tumor_cols)]
> normal_data <- expr_renamed[, c("Probe_ID", normal_cols)]
> 
> 
> # 4. COMPUTE AVERAGE EXPRESSION FOR EACH GROUP ----------------------------
> 
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
> 
> # 5. LOG2 FOLD CHANGE FOR EACH PROBE --------------------------------------
> 
> # Use standard fold change: log2(Tumor / Normal)
> merged_means <- merged_means %>%
+   mutate(log2FC = log2(Tumor_mean / Normal_mean))
> 
> # Remove non-finite values (e.g., if Normal_mean == 0)
> merged_means <- merged_means %>% filter(is.finite(log2FC))
> 
> head(merged_means)
# A tibble: 6 × 4
  Probe_ID     Tumor_mean Normal_mean  log2FC
  <chr>             <dbl>       <dbl>   <dbl>
1 ILMN_1343291     65841.      69653. -0.0812
2 ILMN_1343295     27802.      17549.  0.664 
3 ILMN_1651199       218.        221. -0.0247
4 ILMN_1651209       258.        271. -0.0716
5 ILMN_1651210       203.        189.  0.107 
6 ILMN_1651221       224.        227. -0.0200
> 
> 
> # 6. MERGE WITH GENE INFORMATION & FILTER |FC| > 5 ------------------------
> 
> full_results <- merged_means %>%
+   inner_join(gene_info, by = "Probe_ID")
> 
> DEGs <- full_results %>%
+   filter(abs(log2FC) > 5)
> 
> 
> # 7. ADD COLUMN: HIGHER IN TUMOR OR NORMAL --------------------------------
> 
> DEGs <- DEGs %>%
+   mutate(Higher_in = ifelse(log2FC > 0, "Tumor", "Normal"))
> 
> head(DEGs)
# A tibble: 6 × 9
  Probe_ID     Tumor_mean Normal_mean log2FC Symbol Entrez_Gene_ID Chromosome Cytoband Higher_in
  <chr>             <dbl>       <dbl>  <dbl> <chr>           <int> <chr>      <chr>    <chr>    
1 ILMN_1652431       545.      31267.  -5.84 CA1               759 8          8q21.2b  Normal   
2 ILMN_1681462      9065.        203.   5.48 REG1B            5968 2          2p12e    Tumor    
3 ILMN_1713462       329.      19228.  -5.87 AQP8              343 16         16p12.1b Normal   
4 ILMN_1763749       785.      25938.  -5.05 GUCA2A           2980 1          1p34.2b  Normal   
5 ILMN_1802441      9412.        213.   5.47 REG1A            5967 2          2p12e    Tumor    
6 ILMN_2192072      9960.        259.   5.27 MMP7             4316 11         11q22.2a Tumor    
> 
> 
> ###########################################################################
> # EXPLORATORY DATA ANALYSIS (EDA)
> ###########################################################################
> 
> # 8a. Histogram: number of DEGs by chromosome -----------------------------
> 
> ggplot(DEGs, aes(x = Chromosome)) +
+   geom_bar(fill = "steelblue") +
+   theme_minimal() +
+   labs(
+     title = "Distribution of DEGs by Chromosome",
+     x = "Chromosome",
+     y = "Number of DEGs"
+   )
> 
> 
> # 8b. Histogram by chromosome & direction (Tumor vs Normal) ---------------
> 
> ggplot(DEGs, aes(x = Chromosome, fill = Higher_in)) +
+   geom_bar(position = "dodge") +
+   theme_minimal() +
+   labs(
+     title = "DEGs by Chromosome and Direction of Change",
+     x = "Chromosome",
+     y = "Number of DEGs",
+     fill = "Higher Expression In"
+   )
> 
> 
> # 8c. Bar chart: % of DEGs up vs down in Tumor ----------------------------
> 
> pct_data <- DEGs %>%
+   count(Higher_in) %>%
+   mutate(Percent = n / sum(n) * 100)
> 
> ggplot(pct_data, aes(x = Higher_in, y = Percent, fill = Higher_in)) +
+   geom_col() +
+   theme_minimal() +
+   labs(
+     title = "Percent of DEGs Upregulated vs Downregulated in Tumor",
+     x = "Higher Expression In",
+     y = "Percent of DEGs"
+   )
> 
> 
> # 8d. HEATMAP OF RAW EXPRESSION (RENAMED SAMPLES) ------------------------
> 
> expr_mat <- as.data.frame(expr_renamed)
> rownames(expr_mat) <- expr_mat$Probe_ID
> expr_mat$Probe_ID <- NULL
> 
> expr_mat <- as.matrix(expr_mat)
> 
> # Scale rows so each gene shows relative high/low expression across samples
> scaled_expr <- scale(expr_mat)
> 
> # Simple heatmap (no clustering) – “gene expression by sample”
> heatmap(
+   scaled_expr,
+   Rowv = NA,
+   Colv = NA,
+   labRow = FALSE,
+   main = "Heatmap of Gene Expression (scaled by gene)"
+ )
> 
> 
> # 8e. CLUSTERMAP (clustered heatmap) --------------------------------------
> 
> heatmap(
+   scaled_expr,
+   labRow = FALSE,
+   main = "Clustermap of Gene Expression (genes & samples clustered)"
+ )
Error: vector memory limit of 16.0 Gb reached, see mem.maxVSize()
> 
> 
samples based on their expression profiles. Tumor and Normal samples tend to cluster apart,
indicating distinct global transcriptional programs between the two phenotypes.
> 

