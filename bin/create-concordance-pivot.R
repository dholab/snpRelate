#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# load necessary packages
library(tidyverse)
library(readxl)
library(writexl)

# parse command line arguments
cycle_table = args[1]

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

# Read in the list of SNPs
SNPs <- read_excel(results_files[1], sheet = 1, skip = 15)[1] %>%
  separate(ID, into = c("sample", "SNP"), sep = "-") %>%
  subset(select = "SNP")

# loop through each SNP to collate calls
concordance_table <- data.frame("SNP" = sort(unique(SNPs$SNP)),	
                                "Call_Map_Order" = NA,
                                "Cycle" = NA, 
                                "Comments" = NA)
for (snp in concordance_table$SNP){
  
  # determine optimal cycle
  optimal_cycle <- as.integer(cycle_table[cycle_table$SNP==snp,"cycle"])
  
  if (!(optimal_cycle %in% available_cycles)){
    optimal_cycle <- 7
  }
  
  # read in the correct BioMark output file based on that cycle
  data <- read_excel(paste("28637_BMX003_Detailed_Results_cycle", optimal_cycle, ".xlsx", sep = ""), 
                     sheet = 1, skip = 15)
  
  # split sample IDs and SNP IDs into their own columns
  data <- separate(data, ID, into = c("sample", "SNP"), sep = "-")
  
  snp_table <- data[data$SNP==snp,]
  new_index <- which(concordance_table$SNP==snp)
  
  concordance_table[new_index, "Call_Map_Order"] <- unique(snp_table$SNP)
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

# write output file for viewing in excel
write_xlsx(concordance_table, "concordance-table.xlsx")
