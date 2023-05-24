#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# load necessary packages
library(tidyverse)
library(readxl)
library(writexl)

# parse command line arguments
cycle_table = args[1]
sample_manifest = args[2]

# create a vector of detailed report files
results_files <- list.files(path = ".", pattern = "_Detailed_Results")

# find available cycles
available_cycles <- c()
for (i in results_files){
  cycle <- unlist(strsplit(i, split = "\\."))[1] %>%
    str_sub(-1,-1) %>%
    as.integer()
  available_cycles <- c(available_cycles, cycle)
}

# read in optimal cycle table
cycle_table <- read.csv(cycle_table)

# remove unknown samples and correct column types
# cycle_table <- cycle_table[cycle_table$cycle!="U" & cycle_table$cycle!="X" &
#                              !is.na(cycle_table$cycle),]
# cycle_table$cycle <- as.numeric(cycle_table$cycle)

# Read in the list of SNPs
SNPs <- read_excel(results_files[1], sheet = 1, skip = 15)[1] %>%
  separate(ID, into = c("sample", "SNP"), sep = "-") %>%
  subset(select = "SNP")

# make an empty data frame that will subsequently be filled for each SNP
concordance_table <- data.frame("SNP" = sort(unique(SNPs$SNP)),	
                                "Assay" = NA,
                                "Cycle" = NA, 
                                "Comments" = NA)
# determine optimal cycles
cycle_table <- cycle_table[order(cycle_table$SNP),]
concordance_table$Cycle <- as.integer(cycle_table[cycle_table$SNP==concordance_table$SNP,"cycle"])
concordance_table <- concordance_table[concordance_table$Cycle!=0,]

# Create an empty list to store the excel dataframes
df_list = list()

# Read in the Excel files and store them in the list
for (file in results_files){
  i <- as.integer(strsplit(strsplit(file, split = "\\.")[[1]][1], split = "_cycle")[[1]][2])
  df_list[[i]] <- read_excel(file, skip = 15, sheet = 1)
}

# loop through each SNP to collate calls
for (snp in concordance_table$SNP){
  
  # determine optimal cycle
  optimal_cycle <- as.integer(cycle_table[cycle_table$SNP==snp,"cycle"])
  
  if (!(optimal_cycle %in% available_cycles)){
    break
  }
  
  # read in the correct BioMark output file based on that cycle
  data <-  df_list[[optimal_cycle]]
  
  # split sample IDs and SNP IDs into their own columns
  data <- separate(data, ID, into = c("sample", "SNP"), sep = "-")
  
  snp_table <- data[data$SNP==snp,]
  new_index <- which(concordance_table$SNP==snp)
  
  concordance_table[new_index, "Assay"] <- unique(snp_table$Assay)
  concordance_table[new_index, "Cycle"] <- optimal_cycle
  
  if (snp == concordance_table$SNP[1]) {

    for (animal in sort(unique(snp_table$Name))){

      call <- snp_table[snp_table$Name==animal,"Final"][1,] %>%
        as.character()
      concordance_table <- cbind(concordance_table, as.data.frame(NA))
      concordance_table[new_index,ncol(concordance_table)] <- call
      colnames(concordance_table)[ncol(concordance_table)] <- animal

    }

  } else if (ncol(concordance_table[,5:ncol(concordance_table)]) != length(unique(data$sample))+1) {
    
    for (animal in sort(unique(snp_table$Name))){
      
      call <- snp_table[snp_table$Name==animal,"Final"][1,] %>%
        as.character()
      concordance_table[new_index,animal] <- call
      
    }
    
  } else {
    break
  }
  
}

# change assay names to SNP names
concordance_table$SNP <- str_replace(concordance_table$SNP, "A", "SNP")

# read in and format the sample manifest
samples <- read_excel(sample_manifest, skip = 1)
samples <- samples[!is.na(samples$`Animal ID`) & !is.na(samples$Relationship),]

# parse out the relationship "groups"
samples <- subset(samples, select = c("Animal ID", "Relationship"))
samples$group <- str_remove(samples$Relationship, "Dam ") %>%
  str_remove("Progeny ") %>% 
  as.numeric()

# go through each group and add them to new vectors that will eventually become
# new column headers for the concordance pivot
ordered_samples <- c()
ordered_relations <- c()
for (i in unique(samples$group)){
  
  # go group by group
  sub <- samples[samples$group==i,]
  
  # first add the dams
  ordered_samples <- c(ordered_samples, as.character(sub[grepl("Dam", sub$Relationship), "Animal ID"]))
  ordered_relations <- c(ordered_relations, as.character(sub[grepl("Dam", sub$Relationship), "Relationship"]))
  
  # then add the progeny
  ordered_samples <- c(ordered_samples, as.character(sub[grepl("Progeny", sub$Relationship), "Animal ID"]))
  ordered_relations <- c(ordered_relations, as.character(sub[grepl("Progeny", sub$Relationship), "Relationship"]))
  
}

# bind together the vectors as a data frame
# ordered_samples <- as.data.frame(ordered_samples)
# ordered_relations <- as.data.frame(ordered_relations)
# new_order <- t(cbind(ordered_samples,ordered_relations))

# define where the animals start in the concordance pivot, and sort them according
# to the new sample vector
animal_start <- which(colnames(concordance_table)=="Comments")+1
animal_cols <- concordance_table[,animal_start:ncol(concordance_table)]
concordance_table <- concordance_table[,1:(animal_start-1)]
for (i in 1:length(ordered_samples)){
  
  animal <- ordered_samples[i]
  relation <- ordered_relations[i]
  sub <- c(animal, relation, animal_cols[,animal])
  
  if (i == 1){
    
    new_cols <- data.frame(sub)
    
  } else {
    
    new_cols <- cbind(new_cols, sub)
    
  }
  
}

# prepare new rows for the concordance table
new_row1 <- c("Animal ID", NA, NA, NA)
new_row2 <- colnames(concordance_table)
colnames(concordance_table) <- NULL
concordance_table <- rbind(new_row1, new_row2, concordance_table)
concordance_table <- cbind(concordance_table, new_cols)
colnames(concordance_table) <- NULL

# write output file for viewing in excel
write_xlsx(concordance_table, "concordance-table.xlsx", col_names = F)
