#!/usr/bin/env python3

import sys
import os
import math
import polars as pl
import numpy
import xlsxwriter
from pysam import VariantFile

# Check if the input spreadsheet is a symlink, and if it is, follow the symlink
table_path = os.path.realpath(sys.argv[1])

# Check if the input VCF is a symlink, and if it is, follow the symlink
vcf_path = os.path.realpath(sys.argv[2])

# define a function that will produce a vector of allele frequencies based on the available data
def define_alt_freqs(table_df: pl.DataFrame, vcf_path: str):

    # parse SNP locations
    snp_names = table_df.select("Assay").to_series().to_list()
    positions = [int(s.split("_")[1]) for s in snp_names]

    # parse VCF
    vcf = VariantFile(vcf_path)

    # count samples
    sample_size = len(list((vcf.header.samples)))

    # Initialize a dictionary for SNP positions and sums of heterozygous and homozygous counts
    snp_freqs = {}

    # Loop over all variants
    for variant in vcf:

        if variant.pos not in positions:
            continue
        else:
            # Initialize counters for each genotype
            homozygous_ref_count = 0
            heterozygous_count = 0
            homozygous_alt_count = 0

            # Loop over all samples
            for sample in variant.samples:
                genotype = variant.samples[sample]['GT']
                
                # Handle missing data (denoted as -1)
                if genotype == (None, None) or genotype == (None,) or genotype == (-1, -1) or genotype == (-1,):
                    continue
                    
                # Count homozygous reference (0/0)
                if genotype == (0, 0):
                    homozygous_ref_count += 1
                # Count heterozygous (0/1 or 1/0)
                elif set(genotype) == {0, 1}:
                    heterozygous_count += 1
                # Count homozygous alternate (1/1)
                elif set(genotype) == {1} and len(genotype) > 1:
                    homozygous_alt_count += 1

            # Store the SNP position and the sum of heterozygote and homozygote alt counts in the dictionary
            snp_freqs[variant.pos] = math.sqrt(1 - ((heterozygous_count + homozygous_alt_count) / sample_size))

    # create a vector of frequencies that are in the same order as the snp names in the concordance table
    y_alt_freqs = [snp_freqs[position] if position in snp_freqs else None for position in positions]
        
    return y_alt_freqs

# define a function that will compute Hardy-Weinberg genotype probabilities for a given allele frequency
def hardy_weinberg(allele_freq: float):

    # define Hardy Weinberg equation
    parentage_freq = 1 - (allele_freq ** 2)

    return parentage_freq

# multiply Hardy-Weinberg probabilities for every SNP call to create a confidence statistic for each progeny
def composite_hw(table_df: pl.DataFrame, relations: list):
  
  composite_values = []
  freq_vector = table_df['ALT Allele Frequency'].to_list()

  for i in range(len(relations)):

    animal = relations[i]

    if "Dam" in animal:
      composite_values+=[None]
    else:
        
      hw_values = []
      calls = table_df[animal].to_list()

      for j in range(len(calls)):

        call = calls[j]

        if call == 'YY' and freq_vector[j] != None:
          new_hw = hardy_weinberg(freq_vector[j])
        elif call == 'XX' and freq_vector[j] != None:
          new_hw = hardy_weinberg(1 - freq_vector[j])
        else:
          new_hw = None
        hw_values+=[new_hw]

      composite_values+=[numpy.prod([float(item) for item in hw_values if item])]

  return composite_values

# define a function that will compute Hardy-Weinberg probabilities for each dam/progeny pair, create a new data frame, 
# and add confidence columns to the concordance table
def compute_confidences(table_path: str, vcf_path: str):

  # read in table as data frame
  table_df = pl.read_excel(table_path, read_csv_options={"skip_rows": 1})

  # define animals and their relationships
  animals = pl.read_excel(table_path, read_csv_options={"has_header": False}).row(0)[4:]
  relations = pl.read_excel(table_path, read_csv_options={"has_header": False}).row(1)[4:]

  # make Y allele frequency column
  table_df = table_df.with_columns(
      pl.Series(name="ALT Allele Frequency", values=define_alt_freqs(table_df, vcf_path))
  )

  # make X allele frequency column
  table_df = table_df.with_columns(
      (1 - pl.col("ALT Allele Frequency")).alias("REF Allele Frequency")
  )

  # construct the new row with string Polars data typing and column naming
  new_row = [None, None, None, "Confidence:"] + composite_hw(table_df, relations) + [None, None]
  new_row = [str(item) if item != None else item for item in new_row]
  row_df = pl.DataFrame({"row": new_row}).transpose(include_header=False, column_names = table_df.columns)

  for i in range(len(table_df.dtypes)):

    name = table_df.columns[i]
    row_df = row_df.with_columns(pl.col(name).cast(table_df.dtypes[i]))
  
  # stack concordance table on top of new row and write out to excel file
  table_df = table_df.vstack(row_df)
  table_df.write_excel('concordance-table-with-conf.xlsx', autofit=True, has_header=True)

  # use animal and relationship vectors to make new vectors of dam IDs and progeny IDs
  dam_IDs = []
  progeny_IDs = []
  hwe_vector = [float(i) for i in new_row if i is not None and i != "Confidence:"]
  for i in range(len(animals)):

    if "Dam" in relations[i]:
      dam_IDs+=[animals[i]]
    elif "Progeny" in relations[i]:
      progeny_IDs+=[animals[i]]
    else:
      continue

  # make a new data frame with one row for each pair
  summary_table = pl.DataFrame(
      {
          "Dam ID": dam_IDs,
          "Progeny ID": progeny_IDs,
          "Population Genotype Probability": hwe_vector,
          "Percent Confidence in Relatedness": 0.0,
          "Odds of Error": "unknown"
      }
  )

  # Loop through each row and add two more columns to conceptualize confidence differently
  summary_table = summary_table.with_columns(
      summary_table.select(
      ((1 - pl.col('Population Genotype Probability')) * 100).alias('Percent Confidence in Relatedness')
  )
  )
  odds = ["1 in " + str(1/i) for i in summary_table['Population Genotype Probability'].to_list()]
  summary_table = summary_table.with_columns(
      pl.Series(odds).alias("Odds of Error")
  )

  # return new data frame
  return summary_table

# run the function
new_table = compute_confidences(table_path, vcf_path)

# write new spreadsheet
import xlsxwriter
new_table.write_excel('final-table.xlsx', autofit=True, has_header=True)