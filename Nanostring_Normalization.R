suppressWarnings({
  library(NACHO)
  library(tidyr)
  library(tidyverse)
  library(dplyr)
  library(stringr)
  library(psych)
})

options(scipen = 999)

##function to identify first gene column index
getStartColumn <- function(data, gene) {
  read_startCol <- unique(grep(gene, colnames(as.data.frame(data)), ignore.case = TRUE, value = FALSE))
  return(read_startCol)
}

ctrls.neg <- c("NEG_1","NEG_2","NEG_3","NEG_4","NEG_5","NEG_6","NEG_7","NEG_8")
ctrls.pos <- c("POS_1","POS_2","POS_3","POS_4","POS_5","POS_6","POS_7","POS_8")
ctrls.all <- c(ctrls.neg, ctrls.pos)

#thresholding function
threshold <- function(df, use_neg_ctrls = FALSE, cutoff = 20, kill = FALSE) {
  if (use_neg_ctrls) {
    cutoff <- max(df[, colnames(df) %in% ctrls.neg])
  }
  df[df <= cutoff | df == Inf | df == -Inf | is.na(df)] <- ifelse(kill, 0, cutoff)
  return(df)
}

build_calibrator_matrix <- function(pos_norm_calibrator_column, default_row = 'A') {
  # grab the default row that will be used for calibration
  # calibrator should be in the format with rownames(calibrator) = sample_lanes and colnames(calibrator) = gene_names
  # calibrator should not contain any control genes
  ##pattern <- paste0('\\(', default_row)
  ##basis <- pos_norm_calibrator_column[grepl(pattern, rownames(pos_norm_calibrator_column)), ]
  basis <- pos_norm_calibrator_column[rownames(pos_norm_calibrator_column) == default_row, ]
  scale <- apply(pos_norm_calibrator_column, 1, function(row) basis / row) %>% Reduce(rbind, .)
  rownames(scale) <- rownames(pos_norm_calibrator_column)
  
  return(scale)
}

technicalNormalize <- function(rawData, use_all_pos_ctrls = FALSE, calibrator.sample.names, highest_of_neg_ctrls = FALSE, read_startCol_gene = '') {
  read_startCol <- getStartColumn(rawData, gene = read_startCol_gene)
  rawData.samplenames <- rawData[, 1:(read_startCol-1)]
  rawData[, 1:(read_startCol-1)] <- NULL
  
  # positive ctrl normalization
  ctrls <- rawData[, colnames(rawData) %in% ctrls.all] 
  ctrls.mean <- ctrls %>% lapply(., mean) %>% as.data.frame(.) %>% t(.) %>% as.data.frame(.)
  colnames(ctrls.mean) <- c('reads')
  num_ctrls <- ifelse(use_all_pos_ctrls, 8, 3)
  top_ctrls <- ctrls.mean %>% arrange(desc(reads)) %>% slice(1:num_ctrls) %>% rownames(.)
  geomeans <- apply(ctrls[, colnames(ctrls) %in% top_ctrls], 1, geometric.mean)
  avg_geomean <- mean(geomeans)
  well.scale <- sapply(geomeans, function(x) avg_geomean / x)
  pos_norm_data <- rawData[, !colnames(rawData) %in% ctrls.all]
  samples <- nrow(pos_norm_data)
  genes <- ncol(pos_norm_data)
  scale <- matrix(rep(well.scale, times = genes), nrow = samples, ncol = genes)
  pos_norm_data <- scale*as.matrix(pos_norm_data) %>% as.data.frame()
  
  well.scale <- as.data.frame(well.scale)
  colnames(well.scale) <- "Positive Control Normalization Factor"
  rownames(well.scale) <- rownames(ctrls)
  
  
  # calibration matrix
  pos_norm_calibrator_column <- pos_norm_data[rownames(pos_norm_data) %in% calibrator.sample.names, ] %>% threshold(., use_neg_ctrls = highest_of_neg_ctrls)
  ##
  rownames(pos_norm_calibrator_column) <- rawData.samplenames$Row[rownames(pos_norm_data) %in% rownames(pos_norm_calibrator_column)]
  ##
  calibration_matrix <- build_calibrator_matrix(pos_norm_calibrator_column)
  ##calibratotion.rows_in_matrix <- rownames(calibration_matrix) %>% gsub("^.*Set\\s*|\\s*\\(.*$", "", .) %>% gsub(" ", "", ., fixed = TRUE)
  calibratotion.rows_in_matrix <- rownames(calibration_matrix)
  tech_norm_data <- lapply(rownames(pos_norm_data), function(sample.name) {
    ##row <- gsub("^.*Set \\s*|\\s* \\(.*$", "", sample.name)
    row <- rawData.samplenames[sample.name, 'Row']
    sample.row <- pos_norm_data[sample.name, ]
    calibration.row <- calibration_matrix[calibratotion.rows_in_matrix == row, ] %>% as.matrix()
    sample.row <- as.matrix(sample.row)
    return(as.data.frame(calibration.row * sample.row))
  }) %>% Reduce(rbind, .) %>% as.data.frame(.)
  tech_norm_data <- as.data.frame(lapply(tech_norm_data, function(x) round(x, digits = 2)))
  tech_norm_data <- threshold(tech_norm_data, kill = FALSE)
  tech_norm_data <- cbind(rawData.samplenames, tech_norm_data, ctrls, well.scale) %>% arrange(., Row)
  
  return(tech_norm_data)
}

##function for content normalization
contentNormalize <- function(TechNormData, tissue_colName, HK_genes, TOIs = 'All', read_startCol_gene = '', separate_tissues = FALSE, highest_neg_ctrl = FALSE) {
  ##subset tissue types and get tissues of interest (TOI)
  TechNormData <- as.data.frame(TechNormData)
  tissue_colName.grep <- unique(grep(paste(tissue_colName, collapse = '|'), colnames(TechNormData), ignore.case = TRUE, value = TRUE))
  if (any(TOIs == "All")) {
    TOIs <- unique(TechNormData[[tissue_colName.grep]]) %>% .[nzchar(.)] %>% .[. != 'Calibrator']
  } else {
    TOIs <- unique(grep(paste(TOIs, collapse = '|'), TechNormData[[tissue_colName.grep]], ignore.case = TRUE, value = TRUE))
  }
  TOIs_data <- TechNormData[TechNormData[[tissue_colName.grep]] %in% TOIs,]
  read_startCol <- getStartColumn(TOIs_data, gene = read_startCol_gene)
  
  ##organize data into subsets for each TOI (optional) and get HKG geomean
  if (separate_tissues) {
    calibrator <- TechNormData[TechNormData[[tissue_colName.grep]] == 'Calibrator',]
    calibrator[, "Content Normalization Factor"] <- c(rep(1, times = nrow(calibrator)))
    TOIs_data_norm <- lapply(TOIs, function(TOI) contentNormalize(TOIs_data, tissue_colName = tissue_colName.grep, 
                                                                  HK_genes = HK_genes, TOIs = TOI, 
                                                                  read_startCol_gene = read_startCol_gene, 
                                                                  separate_tissues = FALSE)) %>% 
      compact(.) %>% Reduce(rbind, .) %>% rbind(., calibrator) %>% arrange(., Row)
    return(TOIs_data_norm)
  } else {
    TOIs_ctrls <- TOIs_data[, colnames(TOIs_data) %in% ctrls.all]
    TOIs_data <- TOIs_data[, !colnames(TOIs_data) %in% ctrls.all]
    
    HKGsearch <- paste(HK_genes, collapse = '|')
    HKGs <- unique(grep(HKGsearch, colnames(TechNormData), ignore.case = TRUE, value = TRUE))
    TOIs_data_HKG <- threshold(TOIs_data, kill = TRUE)
    TOIs_data_HKG[, 1:(read_startCol-1)] <- NULL     #remove unnecessary columns
    needsNormalizing <- colnames(TOIs_data_HKG) %in% HKGs
    if (!any(needsNormalizing)) {
      out <- cbind(TOIs_data, TOIs_ctrls) %>% rbind(TechNormData[TechNormData[[tissue_colName.grep]] == 'Calibrator',])
      return(out)
    } else {
      positive_ctrl_norm <- TOIs_data[, "Positive Control Normalization Factor"] %>% as.data.frame(.)
      colnames(positive_ctrl_norm) <- "Positive Control Normalization Factor"
      rownames(positive_ctrl_norm) <- rownames(TOIs_data)
      TOIs_data[, "Positive Control Normalization Factor"] <- NULL
      
      TOIs_data_HKG <- TOIs_data_HKG[, needsNormalizing]
      TOIs_data_HKG_geomeans <- apply(TOIs_data_HKG, 1, geometric.mean)
      #TOIs_data_HKG_geomeanAVG <- mean(TOIs_data_HKG_geomeans)
      TOIs_data_HKG_geomeanAVG <- 2000 ##set at a constant value so can compare normalized counts directly across different runs and studies
      TOIs_data_HKG_scale.well <- sapply(TOIs_data_HKG_geomeans, function(x) TOIs_data_HKG_geomeanAVG / x)
      
      ##scale the data by the HKG scaling factor (content normalization)
      TOIs_data_norm <-TOIs_data
      TOIs_data_samplenames <- TOIs_data_norm[, 1:(read_startCol-1)]
      TOIs_data_norm[, 1:(read_startCol-1)] <- NULL
      TOIs_data_norm <- threshold(TOIs_data_norm, kill = TRUE)
      samples <- nrow(TOIs_data_norm)
      genes <- ncol(TOIs_data_norm)
      TOIs_data_HKG_scale <- matrix(rep(TOIs_data_HKG_scale.well, times = genes), nrow = samples, ncol = genes)
      
      TOIs_data_HKG_scale.well <- as.data.frame(TOIs_data_HKG_scale.well)
      colnames(TOIs_data_HKG_scale.well) <- "Content Normalization Factor"
      rownames(TOIs_data_HKG_scale.well) <- rownames(TOIs_data_norm)
      
      TOIs_data_norm <- TOIs_data_HKG_scale*as.matrix(TOIs_data_norm) %>% as.data.frame(.)
      TOIs_data_norm <- as.data.frame(lapply(TOIs_data_norm, function(x) round(x, digits = 2)))
      TOIs_data_norm <- threshold(TOIs_data_norm, use_neg_ctrls = highest_neg_ctrl)
      TOIs_data_norm <- cbind(TOIs_data_samplenames, TOIs_data_norm, TOIs_ctrls, positive_ctrl_norm, TOIs_data_HKG_scale.well) %>% arrange(., Row)
      
      return(TOIs_data_norm)
    }
  }
}

get_rawCounts_from_NACHO <- function(nachoCounts) {
  nacho_columns <- c('rcc.files', 'Lane_Attributes.lane_ID', 'plexset_id', 'Name', 'Count')
  counts.raw <- nachoCounts %>% as.data.frame() %>%  .[, colnames(.) %in% nacho_columns]
  colnames(counts.raw) <- c('Filename', 'Column', 'Row', 'Gene', 'Raw Counts')
  counts.raw$Row <- utf8ToInt('A') - 1 + as.integer(gsub('S', '', counts.raw$Row))
  counts.raw$Row <- lapply(counts.raw$Row, intToUtf8)
  counts.raw$Column <- str_pad(counts.raw$Column, 2, pad = "0")
  counts.raw$Readname <- paste0(rep('Set ', nrow(counts.raw)), counts.raw$Row, rep(' (', nrow(counts.raw)), 
                                counts.raw$Row, counts.raw$Column, rep(') ', nrow(counts.raw)), counts.raw$Filename)
  counts.raw <- lapply(sample.names, function (sample.name) {
    well.samples <- counts.raw[counts.raw$Readname == sample.name, colnames(counts.raw) %in% c('Gene', 'Raw Counts')] %>% as.data.frame()
    rownames(well.samples) <- well.samples[, 'Gene']
    well.samples[, 'Gene'] <- NULL
    colnames(well.samples) <- sample.name
    well.samples <- t(well.samples) %>% as.data.frame(.)
    return(well.samples)
  }) %>% Reduce(rbind, .) %>% select(sort(names(.))) %>% select(which(!colnames(.) %in% ctrls.all), ctrls.all)
  return(counts.raw)
}

## these are not case sensitive, but the spelling is important
tissue_col <- c('Organ', 'Muscle', 'Tissue')
HKGs <- c('GHITM_1', 'GTF2B_1', 'IPO13_1', 'PPIB_1', 'RER1_1', 'SDHA_1', 'SSB_1', 'TMBIM6_1', 'TXNIP_1')
TOI <- 'All'
start_gene <- 'A2M_1'
use_highest_neg_ctrl <- FALSE
separate_the_tissues <- TRUE

##The platemap should be an excel sheet with the following:
## - first cell is blank, and the first column contains the names of the samples as they would appear in nSolver
##    (i.e. "Set {Row Letter} ({Well}) {Experiment Name}_{Column Name}_{Column Number}.RCC)
## - all subsequenc columns have column names located in the same row as the first blank cell
## - The column names at least include "Row", "Column", "Group", "Muscle"/"Organ", and "Sample ID"/"Sample"
## - the calibrator is tagged as "Calibrator" in the "Group" column and "Muscle/Organ" column
rcc.dir <- file.path(
  r"{C:\Users\ayu\OneDrive - Avidity Biosciences\Documents\Data\Nanostring\CNM\2022-691-C57\2022.06.30 - 691 D21 Gas TA FOP panel\691-GasOnly-Calibrator-RCC}")
platemap.path <- file.path(
  r"{C:\Users\ayu\OneDrive - Avidity Biosciences\Documents\Data\Nanostring\CNM\2022-691-C57\2022.06.30 - 691 D21 Gas TA FOP panel\2022.06.30 - 691 FOP panel - Platemap.csv}")
expt_title <- 'Test - 691 D21 FOP Gas TA'

platemap <- read.csv(platemap.path, row.names = 1, check.names = FALSE)
sample.names <- rownames(platemap)
calibrator.sample.names <- sample.names[which(platemap$Group == 'Calibrator')]
rcc.files <- dir(rcc.dir, full.names = FALSE)
ssheet_csv <- as.data.frame(rcc.files)
rcc.nacho <- load_rcc(data_directory = rcc.dir, ssheet_csv = ssheet_csv, id_colname = 'rcc.files')

counts.raw <- get_rawCounts_from_NACHO(rcc.nacho$nacho) %>% cbind(platemap[rownames(platemap) %in% rownames(.), ], .)
counts.norm.tech <- technicalNormalize(counts.raw, calibrator.sample.names = calibrator.sample.names, 
                                       highest_of_neg_ctrls = use_highest_neg_ctrl, read_startCol_gene = start_gene)
counts.norm <- contentNormalize(counts.norm.tech, separate_tissues = separate_the_tissues, tissue_colName = tissue_col, 
                                TOIs = TOI, read_startCol_gene = start_gene, HK_genes = HKGs, highest_neg_ctrl = use_highest_neg_ctrl)

SAVE <- TRUE
output.dir <- r"{C:\Users\ayu\OneDrive - Avidity Biosciences\Documents\Data\Nanostring\CNM\2022-691-C57\2022.06.30 - 691 D21 Gas TA FOP panel}"
output.filename <- 'R_normalization_test.csv'
if (SAVE) {
  write.csv(counts.norm, file = paste(output.dir, output.filename, sep = r"{\}"))
}

##The full data file should include the following:
## - first cell is blank, and the first column contains the unique names of each sample
## - all subsequent columns have column names located in the same row as the first blank cell
## - The column names at least include "row", Group" and "Muscle"/"Organ"
## - the calibrator is tagged as "Calibrator" in the "Group" column and "Muscle/Organ" column
fulldata.path <- r"{C:\Users\ayu\Downloads\Control Charting_Raw (1).csv}"

fulldata.raw <- read.csv(fulldata.path, row.names = 1, check.names = FALSE)
fulldata.platemap <- fulldata.raw[1:(getStartColumn(fulldata.raw, start_gene) - 1)]
fulldata.sample.names <- rownames(fulldata.platemap)
fulldata.calibrator.sample.names <- fulldata.sample.names[which(fulldata.platemap$Group == 'Calibrator')]

fulldata.norm.tech <- technicalNormalize(fulldata.raw, calibrator.sample.names = fulldata.calibrator.sample.names, 
                                       highest_of_neg_ctrls = use_highest_neg_ctrl, read_startCol_gene = start_gene)
fulldata.norm <- contentNormalize(fulldata.norm.tech, separate_tissues = separate_the_tissues, tissue_colName = tissue_col, 
                                TOIs = TOI, read_startCol_gene = start_gene, HK_genes = HKGs, highest_neg_ctrl = use_highest_neg_ctrl)

SAVE <- TRUE
output.dir <- r"{C:\Users\ayu\Downloads}"
output.filename <- 'R_Control_Charting_Normalized_by_tissue.csv'
if (SAVE) {
  write.csv(fulldata.norm, file = paste(output.dir, output.filename, sep = r"{\}"))
}


