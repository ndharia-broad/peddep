# peddep
Code for "A First-Generation Pediatric Cancer Dependency Map" (Nature Genetics 2021)

This is the code used to run the analyses and generate the figures for Dharia et al, "A First-Generation Pediatric Cancer Dependency Map" 
published in Nature Genetics 2021. 

In order to execute the R-markdowns and other code, one must download additional data from figshare at 10.6084/m9.figshare.13516208. These 
four folders should be unzipped in the root directory of the code in order for all commands to execute properly. 

Each R-markdown performs a specific set of analyses and generates panels of figures as follows. In general, they can be run in the order
listed below to ensure all data dependencies are met prior to running a piece of code. 
1. Cell_Line_Features_Figures.Rmd: Generates heatmaps/datasets describing recurring alterations in pediatric cancer cell lines.
2. CNA_Figures.Rmd: Generates figures/datasets related to copy number alterations of cell lines.
3. Dependency_Rates_Figures.Rmd: Generates figures/datasets related to number of genetic dependencies in cell lines.
4. Dependency_Muation_Figures.Rmd: Generates figures/datasets related to estimates of mutation rates in cell lines and the relationship
to genetic dependencies as well as relationships between dependencies/mutations and confounders.
5. Dependency_Expression_Embedding_Figures.Rmd: Generates figures/datasets related to two-dimensional embedding of genetic dependencies
and gene expression in cell lines as well as comparing homogeneity of genetic dependencies and gene expression within specific tumor types.
6. Selective_Enriched_Dependency_Figures.Rmd Generates figures/datasets related to selective and/or enriched dependencies in pediatric
tumors.
7. Known_Potential_Targets_Figures.Rmd: Generates figures/datasets related to genetic dependencies on known or potential cancer therapeutic
targets in pediatric cancers.
8. Compare_MCL1_Screens_Figures.Rmd: Generates figures/datasets related to performance of MCL1 dependency in multiple CRISPR and small
molecule screens.
9. NonCERES_Selective_Enriched_Dependency_Figures.Rmd: Generates figures/datasets to identify any potential systematic bias introduced by 
CERES in the genetic dependency data.
10. Generate_Paper_Figures.Rmd: Generates draft figures from panels that the above scripts have created. 
