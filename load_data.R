### DATASETS TO USE ###
version_to_use <- '20q1'
version_url <- "https://figshare.com/articles/DepMap_20Q1_Public/11791698/2"
sample_info_url <- "https://ndownloader.figshare.com/files/21781221"
gene_dependency_url <- "https://ndownloader.figshare.com/files/21521871"
gene_effect_url <- "https://ndownloader.figshare.com/files/21521910"
gene_expression_url <- "https://ndownloader.figshare.com/files/21521940"
fusions_url <- "https://ndownloader.figshare.com/files/21521946"
mutations_url  <- "https://ndownloader.figshare.com/files/21521967"
segment_cn_url <- "https://ndownloader.figshare.com/files/21521988"
gene_cn_url <- "https://ndownloader.figshare.com/files/21521964"

`Dataset Used` <- c(version_url, sample_info_url, gene_dependency_url, gene_effect_url, 
                    gene_expression_url, fusions_url, mutations_url, segment_cn_url, 
                    gene_cn_url)
names(`Dataset Used`) <- c(version_to_use, "Sample Info", "Gene Dependency", 
                           "Gene Effect", "Gene Expression", "Fusions", "Mutations", 
                           "Segment Copy Number", "Gene Copy Number")
### END DATASETS TO USE ###

### FORMAT SAMPLE INFO ###
if(!dir.exists('data')) { dir.create('data') }

if(!file.exists("data/sample_info.csv")) # If file does not exist, download from FigShare. Otherwise load local file.
{
  sample_info <- fread(sample_info_url)
  fwrite(sample_info, "data/sample_info.csv")
} else
{
  sample_info <- fread("data/sample_info.csv")
}

# Pediatric cancers used for these analyses.
abbreviations <- c(
  "Ewing"="ES"
  ,"Synovial Sarcoma"="SS"
  ,"Osteosarcoma"="OS"
  ,"Rhabdomyosarcoma"="RMS"
  ,"Retinoblastoma"="RB"
  ,"Rhabdoid"="RT"
  ,"Medulloblastoma"="MB"
  ,"Neuroblastoma"="NB"
  ,"Hepatoblastoma"="HB"
  ,"Pediatric Glioma"="pGBM"
  ,"Pediatric CNS PNET"="pPNET"
  ,"Pediatric Germ Cell"="pGCT"
  ,"Wilms"="WT"
  ,"Pediatric Sarcoma"="pSarc"
  ,"Renal Medullary Carcinoma"="RMC"
)
pediatric_types <- names(abbreviations)

# Filter sample info for solid tumors and fix annotations.
mf <- sample_info %>%
  # Filter lines that do not have DepMap IDs
  dplyr::filter(!is.na(DepMap_ID) & DepMap_ID != '') %>% 
  # Select fields of interest
  dplyr::select(DepMap_ID, CCLE_name=CCLE_Name, Type=disease, T1=lineage, T2=lineage_subtype, T3=lineage_sub_subtype, age) %>%
  # Custom relabelling of cell lines
  # Remove hemeatology malignancies
  dplyr::filter(!grepl('HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', CCLE_name)) %>%
  dplyr::filter(!grepl('blood', T1)) %>%
  dplyr::filter(!grepl('lymphoma', T2)) %>%
  # Remove matched normal tissues (expression only, no dependency data)
  dplyr::filter(!grepl('MATCHED_NORMAL_TISSUE', CCLE_name)) %>%
  # Remove engineered cell lines that the Broad has produced for specific experiments
  dplyr::filter(!grepl('ENGINEERED', CCLE_name)) %>%
  # Remove CHLA57 as it is unclear what the origin of this cell line is; no EWS-ETS fusion and RNAseq not consistent with Ewing
  dplyr::filter(!grepl('CHLA57_BONE', CCLE_name)) %>%
  # Change type from "bone" to "Ewing" or "Osteosarcoma" etc
  dplyr::mutate(Type = case_when(
    T1 == 'bone' ~ stringr::str_to_title(T2),
    TRUE ~ Type
  )) %>% 
  # Change type from "central_nervous_system" to "Medulloblastoma" or "Glioma" etc
  dplyr::mutate(Type = case_when(
    T1 == 'central_nervous_system' ~ stringr::str_to_title(T2),
    TRUE ~ Type
  )) %>%
  # Lines that have "engineer" in their T1 classification are immortalized lines
  dplyr::mutate(Type = ifelse(
    grepl('engineer', T1),
    "Immortalized",
    Type
  )) %>%
  # Change type from "eye" to "Retinoblastoma" or "Uveal melanoma" etc
  dplyr::mutate(Type = case_when(
    T1 == 'eye' ~ stringr::str_to_title(T2),
    TRUE ~ Type
  )) %>%
  # Change type from "lung" to more specific subtypes like "NSCLC"
  dplyr::mutate(Type = case_when(
    T1 == 'lung' & T2 != "" ~ stringr::str_to_title(T2),
    TRUE ~ Type
  )) %>%
  # Change type from "soft_tissue" to more specific types like "Rhabdomyosarcoma" or "Synovial sarcoma" etc
  dplyr::mutate(Type = case_when(
    T1 == 'soft_tissue' & T2 != "" ~ stringr::str_to_title(T2),
    TRUE ~ Type
  )) %>%
  # Anything that has the word "lung" spelled out is classified as "Lung"
  dplyr::mutate(Type = ifelse(
    grepl('lung', Type),
    "Lung",
    Type
  )) %>%
  # Expand out text of "NSCLC" to full form
  dplyr::mutate(Type = ifelse(
    grepl('Nsclc', Type),
    "Non-Small Cell Lung",
    Type
  )) %>%
  # Expand out text of "SCLC" to full form
  dplyr::mutate(Type = ifelse(
    grepl('Sclc', Type),
    "Small Cell Lung",
    Type
  )) %>%
  # Classify "ATRT" as "Rhabdoid"
  dplyr::mutate(Type = ifelse(
    grepl('Atrt', Type),
    "Rhabdoid",
    Type
  )) %>%
  # Classify things that contain "rhabdoid" as "Rhabdoid" i.e. "malignant rhabdoid tumor"
  dplyr::mutate(Type = ifelse(
    grepl('[Rr]habdoid', Type),
    "Rhabdoid",
    Type
  )) %>%
  # Classify any tumor names that contain "ewing" or "Ewing" as "Ewing"
  dplyr::mutate(Type = ifelse(
    grepl("[Ee]wing",Type), 
    "Ewing", 
    Type)) %>%
  # Classify liver subgroup "hepatoblastoma" as a separate type
  dplyr::mutate(Type = ifelse(
    T2=='hepatoblastoma', 
    'Hepatoblastoma', 
    Type)) %>%
  # "PNET" tumors are only CNS tumors in our workflow so call them "CNS PNET"
  dplyr::mutate(Type = ifelse(
    T2=='PNET', 
    'CNS PNET', 
    Type)) %>%
  # Classify "glioma" as pediatric if comes from a patient <=21 yo
  dplyr::mutate(Type = ifelse(
    T2 == "glioma" & (age <= 21 & !is.na(age)), 
    'Pediatric Glioma', 
    Type)) %>%
  # Classify "CNS PNET" as pediatric if comes from a patient <=21 yo
  dplyr::mutate(Type = ifelse(
    T2 == "PNET" & (age <= 21 & !is.na(age)), 
    'Pediatric CNS PNET', 
    Type)) %>%
  # Classify "germ cell tumor" as pediatric if comes from a patient <=21 yo
  dplyr::mutate(Type = ifelse(
    (T2 == "mixed_germ_cell" | T2 == "teratoma") & (age <= 21 & !is.na(age)), 
    'Pediatric Germ Cell', 
    Type)) %>%
  # Classify any type with "wilms" in the type as "Wilms"
  dplyr::mutate(Type = ifelse(
    grepl('wilms',T2), 
    'Wilms', 
    Type)) %>%
  # Change all types to title case
  dplyr::mutate(Type = stringr::str_to_title(gsub("_"," ",Type))) %>% 
  # Classify "MPNST" as a separate type
  dplyr::mutate(Type = ifelse(
    T2 == "MPNST", 
    "MPNST", 
    Type)) %>%
  # Fix capitalization of "CNS PNET" lines
  dplyr::mutate(Type = ifelse(
    grepl("Cns Pnet",Type), 
    gsub("Cns Pnet","CNS PNET",Type), 
    Type)) %>%
  # Classify "Fibrous Histiosarcoma" as "Sarcoma"
  dplyr::mutate(Type = ifelse(
    grepl("Fibrous Histio",Type), 
    "Sarcoma", 
    Type)) %>%
  # Classify "Undifferentiated Sarcoma" as "Sarcoma"
  dplyr::mutate(Type = ifelse(
    grepl("Undifferentiated Sarcoma",Type), 
    "Sarcoma", 
    Type)) %>%
  # Classify "Pleomorphic Sarcoma" as "Sarcoma"
  dplyr::mutate(Type = ifelse(
    grepl("Pleomorphic Sarcoma",Type), 
    "Sarcoma", 
    Type)) %>%
  # Classify "Uterine Sarcoma" as "Sarcoma"
  dplyr::mutate(Type = ifelse(
    grepl("Uterine Sarcoma",Type), 
    "Sarcoma", 
    Type)) %>%
  # Classify "Thyroid Sarcoma" as "Sarcoma"
  dplyr::mutate(Type = ifelse(
    grepl("Thyroid Sarcoma",Type), 
    "Sarcoma", 
    Type)) %>%
  # Classify "Round Cell Sarcoma" as "Sarcoma"
  dplyr::mutate(Type = ifelse(
    grepl("Round Cell Sarcoma",Type), 
    "Sarcoma", 
    Type)) %>%
  # Classify "Epithelioid Sarcoma" as "Sarcoma"
  dplyr::mutate(Type = ifelse(
    grepl("Epithelioid Sarcoma",Type), 
    "Sarcoma", 
    Type)) %>%
  # Classify "Dermatofibrosarcoma Protuberans" as "Sarcoma"
  dplyr::mutate(Type = ifelse(
    grepl("Dermatofibrosarcoma Protuberans",Type), 
    "Sarcoma", 
    Type)) %>%
  # Classify "CCLFPEDS" "Sarcoma" lines as "Pediatric Sarcoma"
  dplyr::mutate(Type = ifelse(
    grepl("CCLFPEDS",CCLE_name) & grepl("Sarcoma",Type), 
    "Pediatric Sarcoma", 
    Type)) %>%
  # Classify "CCLFPEDS" "Kidney" lines as specific type
  dplyr::mutate(Type = ifelse(
    grepl("CCLFPEDS",CCLE_name) & grepl("Kidney",Type) & T3 != "", 
    stringr::str_to_title(gsub("_", " ", T3)), 
    Type)) %>%
  # Reclassify TASK1 as Ewing as it has a EWS-FLI1 fusion and RNAseq consistent with Ewing
  dplyr::mutate(Type = ifelse(CCLE_name %in% c('TASK1_CENTRAL_NERVOUS_SYSTEM'), 
                              "Ewing", 
                              Type
  )) %>%
  # Reclassify BIN67 and SCCOHT1 as SCCOHT instead of rhabdoid
  dplyr::mutate(Type = ifelse(CCLE_name %in% c('BIN67_OVARY'), 
                              "SCCOHT", 
                              Type
  )) %>%
  dplyr::mutate(Type = ifelse(CCLE_name %in% c('SCCOHT1_OVARY'), 
                              "SCCOHT", 
                              Type
  )) %>%
  # Get rid of trailing "Cancer" or "Tumor" for types
  dplyr::mutate(Type = gsub(" Cancer","",Type)) %>%
  dplyr::mutate(Type = gsub(" Tumor","",Type)) %>%
  # Get only distinct entries in sample_info
  distinct() %>%
  # Add RMS specific annotations 
  dplyr::mutate(T2=ifelse(
    CCLE_name %in% paste0(c('RC2', 'RH41', 'SJRH30', 'CW9019', 'RH28', 'RH30', 'RH4', 'RHJT'), '_SOFT_TISSUE'), 'alveolar', ifelse(
      CCLE_name %in% paste0(c('RH18DM', 'RH36', 'SCMCRM2', 'JR', 'RD', 'SMSCTR', 'TTC442'), '_SOFT_TISSUE'), 'embryonal', T2
    ))) %>%
  # Add peds vs adult annotation
  dplyr::mutate(PvA=ifelse(
    Type %in% pediatric_types, 
    'Pediatric', 
    ifelse(grepl("[fF]ibroblast", Type), 'Fibroblast', 'Adult')
  )) %>%
  # Remove lines without annotated type
  dplyr::filter(Type != "") %>%
  dplyr::filter(!grepl("Unknown", Type))
fwrite(mf, 'data/pediatric_annotations.csv')

### END FORMAT SAMPLE INFO ###

### FORMAT ADULT AND PEDIATRIC LISTS AND COLORS ###

# All adult cancer type with more than 10 lines
adult_cancers_to_use <- mf %>%
  group_by(Type) %>%
  dplyr::summarise(count=n()) %>% arrange(-count) %>%
  dplyr::filter(!(Type %in% pediatric_types)) %>%
  dplyr::filter(count > 10) %$%
  unique(Type)

# All adult cancer types
all_adult_types_remaining <- mf %>% dplyr::filter(PvA != 'Pediatric') %>% distinct(Type) %>% dplyr::filter(Type != 'Adult') %$% Type

# Subset of pediatric lines
pediatric_lines <- mf %>%
  dplyr::filter(Type %in% pediatric_types) %>%
  arrange(CCLE_name)

# Subset of adult lines
adult_lines_to_use <- mf %>% 
  dplyr::filter(Type %in% adult_cancers_to_use) %>%
  arrange(CCLE_name)

# Color palette 
color_for_subtypes_vector <- c(
  "Ewing" = "#3C5488", #darkblue 
  "Synovial Sarcoma" = "#B09C85", #lightbrown 
  "Osteosarcoma" = "#7fc97f", #lightgreen 
  "Rhabdomyosarcoma" = "#E64B35", #red 
  "Retinoblastoma" = "#7E6148", #brown 
  "Rhabdoid" = "#a6cee3", #blue 
  "Medulloblastoma" = "#F39B7F", #salmon 
  "Neuroblastoma" = "#984ea3", #purple 
  "Hepatoblastoma" = "#4DBBD5", #teal 
  "Pediatric Glioma" = "#00A087", #green
  "Pediatric CNS PNET" = "#DC0000", #darkred
  "Pediatric Germ Cell" = "orange",
  "Wilms" = "cyan",
  "Pediatric Sarcoma" = "magenta",
  "Renal Medullary Carcinoma" = "blue"
)

complete_type_colors <- c(
  color_for_subtypes_vector,
  setNames(rainbow(n=length(all_adult_types_remaining)), all_adult_types_remaining),
  c('Adult'='grey')
)

### END FORMAT ADULT AND PEDIATRIC LISTS AND COLORS ###

### LOAD DATA FUNCTIONS ###

# GENE DEPENDENCY
load_gene_dependency <- function(){
  if(!file.exists("data/gene_dependency.csv")) # If file does not exist, download from FigShare. Otherwise load local file.
  {
    gene_dependency <- fread(gene_dependency_url)
    fwrite(gene_dependency, "data/gene_dependency.csv")
  } else
  {
    gene_dependency <- fread("data/gene_dependency.csv")
  }
  # Reformat as a matrix.
  gene_dependency %<>% dplyr::filter(V1 %in% mf$DepMap_ID)
  rownames_temp <- gene_dependency$V1
  gene_dependency <- as.matrix(gene_dependency[,-1])
  rownames(gene_dependency) <- rownames_temp
  return(gene_dependency)
}

# GENE EFFECT
load_gene_effect <- function(){
  if(!file.exists("data/gene_effect.csv")) # If file does not exist, download from FigShare. Otherwise load local file.
  {
    gene_effect <- fread(gene_effect_url)
    fwrite(gene_effect, "data/gene_effect.csv")
  } else
  {
    gene_effect <- fread("data/gene_effect.csv")
  }
  # Reformat as a matrix.
  gene_effect %<>% dplyr::filter(V1 %in% mf$DepMap_ID)
  rownames_temp <- gene_effect$V1
  gene_effect <- as.matrix(gene_effect[,-1])
  rownames(gene_effect) <- rownames_temp
  return(gene_effect)
}

# GENE EXPRESSION
load_gene_expression <- function(){
  if(!file.exists("data/gene_expression.csv")) # If file does not exist, download from FigShare. Otherwise load local file.
  {
    gene_expression <- fread(gene_expression_url)
    fwrite(gene_expression, "data/gene_expression.csv")
  } else
  {
    gene_expression <- fread("data/gene_expression.csv")
  }
  # Reformat as a matrix.
  gene_expression %<>% dplyr::filter(V1 %in% mf$DepMap_ID)
  rownames_temp <- gene_expression$V1
  gene_expression <- as.matrix(gene_expression[,-1])
  rownames(gene_expression) <- rownames_temp
  return(gene_expression)
}

# FUSIONS
load_fusions <- function(){
  if(!file.exists("data/fusions.csv")) # If file does not exist, download from FigShare. Otherwise load local file.
  {
    fusions <- fread(fusions_url)
    fwrite(fusions, "data/fusions.csv")
  } else
  {
    fusions <- fread("data/fusions.csv")
  }
  fusions %<>% dplyr::filter(DepMap_ID %in% mf$DepMap_ID)
  return(fusions)
}

# MUTATIONS
load_mutations <- function(){
  if(!file.exists("data/mutations.csv")) # If file does not exist, download from FigShare. Otherwise load local file.
  {
    mutations <- fread(mutations_url)
    fwrite(mutations, "data/mutations.csv")
  } else
  {
    mutations <- fread("data/mutations.csv")
  }
  mutations %<>% dplyr::filter(DepMap_ID %in% mf$DepMap_ID)
  return(mutations)
}

# SEGMENT CN
load_segment_cn <- function(){
  if(!file.exists("data/segment_cn.csv")) # If file does not exist, download from FigShare. Otherwise load local file.
  {
    segment_cn <- fread(segment_cn_url)
    fwrite(segment_cn, "data/segment_cn.csv")
  } else
  {
    segment_cn <- fread("data/segment_cn.csv")
  }
  segment_cn %<>% dplyr::filter(DepMap_ID %in% mf$DepMap_ID)
  return(segment_cn)
}

# GENE CN
load_gene_cn <- function(){
  if(!file.exists("data/gene_cn.csv")) # If file does not exist, download from FigShare. Otherwise load local file.
  {
    gene_cn <- fread(gene_cn_url)
    fwrite(gene_cn, "data/gene_cn.csv")
  } else
  {
    gene_cn <- fread("data/gene_cn.csv")
  }
  # Reformat as a matrix.
  gene_cn %<>% dplyr::filter(V1 %in% mf$DepMap_ID)
  rownames_temp <- gene_cn$V1
  gene_cn <- as.matrix(gene_cn[,-1])
  rownames(gene_cn) <- rownames_temp
  return(gene_cn)
}

# ACHILLES COMMON ESSENTIALS
load_achilles_common_essentials <- function(){
  if(!file.exists("data/achilles_common_essentials.csv")) # If file does not exist, calculate common essentials. Otherwise load local file.
  {
    # Identify pan-dependent genes as those for whom 90% of cell lines rank the gene above a given dependency cutoff. 
    # The cutoff is determined from the central minimum in a histogram of gene ranks in their 90th percentile least dependent line.
    gene_effect <- load_gene_effect()
    gene_effect_rank <- t(colRanks(t(gene_effect)))
    rownames(gene_effect_rank) <- colnames(gene_effect)
    colnames(gene_effect_rank) <- rownames(gene_effect)
    gene_effect_rank <- t(gene_effect_rank) / colMaxs(gene_effect_rank, na.rm=T)
    achilles_common_essential_score <- apply(gene_effect_rank, 2, quantile, probs = 0.90, na.rm = T)
    d <- density(achilles_common_essential_score)
    cutoff <- optimize(approxfun(d$x,d$y),interval=c(0.1, 0.5))$minimum
    
    achilles_common_essentials <- as.data.frame(names(which(achilles_common_essential_score < cutoff)), stringsAsFactors = F) %>%
      magrittr::set_colnames("gene")
    fwrite(achilles_common_essentials, "data/achilles_common_essentials.csv", sep = "\t")
  } else
  {
    achilles_common_essentials <- fread("data/achilles_common_essentials.csv", sep = "\t")
  }
  return(achilles_common_essentials)
}

# MSIGDB GENE SETS
load_msigdb <- function(){
  if(!file.exists("data/msigdb.rds")) # If file does not exist, generate from local files. Otherwise load local rds.
  {
    msigdb_h <- read.table('input_files/msigdb/h.all.v7.1.entrez.gmt', header = F, sep = ",")
    msigdb_c2 <- fread('input_files/msigdb/c2.all.v7.1.entrez.gmt', header = F, sep = ",")
    msigdb_c5 <- fread('input_files/msigdb/c5.all.v7.1.entrez.gmt', header = F, sep = ",")
    
    msigdb <- list("H" = NULL, "C2" = NULL, "C5" = NULL)
    
    msigdb$H <- as_tibble(do.call(rbind, apply(msigdb_h, 1, function(x) 
    { 
      y <- unlist(strsplit(x, "\t")[[1]])
      return(cbind(y[1], y[3:length(y)]))
    } )) %>% set_colnames(., c("term", "gene")))
    
    msigdb$C2 <- as_tibble(do.call(rbind, apply(msigdb_c2, 1, function(x) 
    { 
      y <- unlist(strsplit(x, "\t")[[1]])
      return(cbind(y[1], y[3:length(y)]))
    } )) %>% set_colnames(., c("term", "gene")))
    
    msigdb$C5 <- as_tibble(do.call(rbind, apply(msigdb_c5, 1, function(x) 
    { 
      y <- unlist(strsplit(x, "\t")[[1]])
      return(cbind(y[1], y[3:length(y)]))
    } )) %>% set_colnames(., c("term", "gene")))
    
    saveRDS(msigdb, "data/msigdb.rds")
  } else
  {
    msigdb <- readRDS("data/msigdb.rds")
  }
  return(msigdb)
}

### END LOAD DATA FUNCTIONS ###

### CALCULATE LRT ###

load_lrt_table <- function(force_calculate = F, cl = detectCores()){
  if(!file.exists("data/lrt_table.csv") | force_calculate) # If file does not exist, calculate LRTs for the current dataset defined above.
  {
    cat("Calculating LRT values... This can take a few minutes to >1 hour depending on the number of cores for parallel processing.\n")
    
    # Load gene effect scores to calculate
    gene_effect <- load_gene_effect()
    median_gene_effect <- apply(gene_effect, 2, median, na.rm = T)
    mean_gene_effect <- apply(gene_effect, 2, mean, na.rm = T)
    
    # Run of the LRT gene calculations
    lrt_genes_data <- pbapply(gene_effect, 2, function(x) {
      source('R/normLRT_test.R')
      g = colnames(x)
      invisible(capture.output(lrt_val <- suppressMessages(suppressWarnings(normLRT_test(x)))))
      names(lrt_val) <- g
      return(lrt_val)
    }, cl = cl)

    lrt_genes_data <- data.frame(gene=names(lrt_genes_data), lrt_genes_data, stringsAsFactors = F) %>%
      dplyr::rename(lrt = lrt_genes_data) %>%
      arrange(-lrt)
    
    # As a rule, negative skewness indicates that the mean of the data values is less than the median
    lrt_genes_data %<>%
      dplyr::mutate(median=median_gene_effect[gene], mean=mean_gene_effect[gene]) %>%
      dplyr::mutate(skewed_left=mean < median)
    
    fwrite(lrt_genes_data, "data/lrt_table.csv")
  } else
  {
    cat("Loading previously calculated LRT values...\n")
    lrt_genes_data <- fread("data/lrt_table.csv")
  }
  return(lrt_genes_data)
}

### END CALCULATE LRT ###
