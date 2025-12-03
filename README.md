# Advanced_R
Advanced_R_HW_Assignment 
This is the read me file for my homework documentation. 

Files Used
Gene_Expression_Data.xlsx – Probe-level raw expression matrix
Gene_Information.csv – Gene annotation including symbol and chromosome
Sample_Information.tsv – Phenotype and patient information for each sample
All files were placed in the ~/Downloads/ directory for analysis.
Methods
Loaded expression, gene annotation, and sample metadata.
Renamed sample columns using phenotype (tumor/normal) and patient ID.
Split the dataset into tumor and normal groups.
Calculated mean expression for each probe within each group.
Computed log2 fold change using the assignment formula:
log2((Tumor – Normal) / Normal).
Filtered probes with an absolute log2 fold change greater than 5.
Added a label indicating whether each DEG was higher in tumor or normal tissue.
Performed exploratory data analysis, including:
Histogram of DEGs by chromosome
Histogram of DEGs by chromosome and expression direction
Bar chart showing percentages of DEGs upregulated in tumor vs. normal
Generated:
A heatmap of scaled expression values
A clustermap with hierarchical clustering on genes and samples
Required Packages
The analysis uses the following R packages:
readxl
dplyr
tidyr
ggplot2
All packages are available through CRAN.
Script
The complete R script is included in the submission under the filename:
gene_expression_analysis.R
The script is designed to run top-to-bottom without manual intervention, assuming the data files remain in the specified directory.
Outputs
The analysis produces:
A table of differentially expressed genes (DEGs)
Visualizations for exploratory analysis
Heatmap and clustermap representing overall expression patterns
These outputs support downstream interpretation of expression differences between tumor and normal tissue.
