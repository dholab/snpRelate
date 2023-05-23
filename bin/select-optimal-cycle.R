#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Load required libraries
library(tidyverse)
library(ggplot2)
library(cluster)
library(MASS)
library(ellipse)
library(readxl)
library(writexl)

# create a vector of detailed report files
results_files <- list.files(path = ".", pattern = "_Detailed_Results")

# make a data frame that records best cycles for each assay
best_cycles <- data.frame("SNP" = as.character(NA),
                          "mean_distance" = as.numeric(0),
                          "cycle" = as.numeric(0))
best_cycles <- best_cycles[NULL,]

# create a loop that goes through each result file
for (i in results_files){
  
  # define index
  index <- which(results_files==i)
  
  # define cycle based on file name
  cycle <- unlist(strsplit(i, split = "\\."))[1] %>%
    str_sub(-1,-1) %>%
    as.integer()
  
  # Read in the data
  data <- read_excel(i, sheet = 1, skip = 15)
  
  # reduce to the columns of interest
  data <- data.frame("cycle" = cycle,
                     "ID" = data$ID,
                     "call" = data$Final,
                     "X.value" = data$`Allele X...11`,
                     "Y.value" = data$`Allele Y...12`)
  
  # split ID column
  data <- separate(data, ID, into = c("sample", "SNP"), sep = "-")
  
  # remove no-calls and NTCs
  data <- data[data$call!="No Call" &
                 data$call!="NTC",] ; rownames(data) <- NULL
  
  # loop through each SNP
  for (snp in sort(unique(data$SNP))){
    
    # add new row to best_cycles data frame
    if (!(snp %in% best_cycles$SNP)){
      best_cycles <- rbind(best_cycles, c(snp, 0, 0))
      colnames(best_cycles) <- c("SNP","mean_distance","cycle")
    }
    
    snp_table <- data[data$SNP==snp,]
    
    # Fit ellipses and find centroids
    centroids <- data.frame(Zygosity = c("XX", "XY", "YY"), X = numeric(3), Y = numeric(3))
    for (zyg in c("XX", "XY", "YY")) {
      # Subset the data by Zygosity
      subset_data <- snp_table[snp_table$call == zyg, ]
      
      if (nrow(subset_data) > 1){
        
        # Calculate the ellipse
        ell <- cov.trob(subset_data[, c("X.value", "Y.value")])
        
        # Calculate the centroid
        centroids[centroids$Zygosity == zyg, c("X", "Y")] <- ell$center
        
        if (nrow(centroids)>1){
          # Calculate pairwise distances
          distances <- dist(centroids[, c("X", "Y")])
        } else {
          distances <- c(0,0)
        }
        
      } else if (nrow(subset_data) == 0){
        
        centroids <- centroids[!centroids$Zygosity == zyg,]
        
      } else {
        
        centroids[centroids$Zygosity == zyg, c("X", "Y")] <- subset_data[, c("X.value", "Y.value")]
        
        if (nrow(centroids)>1){
          # Calculate pairwise distances
          distances <- dist(centroids[, c("X", "Y")])
        }
        
      }
      
      if (mean(centroids$X)==0){
        distances <- c(0,0)
      }
      
      # determine if this is the best-yet cycle for this SNP
      if (mean(distances) > best_cycles[best_cycles$SNP==snp, "mean_distance"] 
          | is.null(best_cycles[best_cycles$SNP==snp, "mean_distance"])){
        best_cycles[best_cycles$SNP==snp, "mean_distance"] <- mean(distances)
        best_cycles[best_cycles$SNP==snp, "cycle"] <- cycle
      }
      
    }
    
  }
  
}

# output file
write.csv(best_cycles, file = "optimal-cycles.csv", row.names = F, quote = F)
