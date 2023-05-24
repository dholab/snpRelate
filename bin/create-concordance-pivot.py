#!/usr/bin/env python3

import sys
import os
import polars as pl
import xlsxwriter

# Check if the input TSV file is a symlink, and if it is, follow the symlink
cycle_table = os.path.realpath(sys.argv[1])
sample_manifest = os.path.realpath(sys.argv[2])

# Find available cycles
results_files = [file for file in os.listdir('.') if "_Detailed_Results" in file]
available_cycles = []
for i in results_files:
    cycle = int(i.split('.')[0].split("cycle")[-1])
    available_cycles.append(cycle)

# Read in optimal cycle table and sort by SNP name
cycle_table = pl.read_csv(cycle_table).sort('SNP')

# Read in the list of SNPs
SNPs = pl.read_excel(results_files[0], read_csv_options={"skip_rows": 15})
SNPs = SNPs.select(pl.col("ID").str.split(by="-")).with_columns(
        pl.struct(
            [
                pl.col('ID').arr.get(i).alias(
                    'Sample' if i==0 else 'SNP') 
                for i in range(2)
            ]
        ).alias('id_struct')
        )
SNPs = SNPs.unnest('id_struct')
SNPs = SNPs.drop('ID')

# Create an empty data frame that will subsequently be filled for each SNP
concordance_table = pl.DataFrame({
    'SNP': SNPs['SNP'].unique().sort(),
    'Assay': "None",
    'Cycle': cycle_table['cycle'],
    'Comments': None
})

# Filter out SNPs with not optimal cycle
concordance_table = concordance_table.filter(pl.col('Cycle') != 0)

# Loop through each SNP to collate calls
for i in range(len(concordance_table)):
  
    # define SNP for this iteration
    snp = concordance_table['SNP'][i]
  
    # Determine optimal cycle
    optimal_cycle = cycle_table.filter(pl.col('SNP') == snp)['cycle'].item()
    
    if optimal_cycle not in available_cycles:
        break

    # Read in the correct BioMark output file based on that cycle
    data = pl.read_excel(f"28637_BMX003_Detailed_Results_cycle{optimal_cycle}.xlsx", read_csv_options={"skip_rows": 15})
    snp_sample = data.select(pl.col("ID").str.split(by="-")).with_columns(
        pl.struct(
            [
                pl.col('ID').arr.get(i).alias(
                    'Sample' if i==0 else 'SNP') 
                for i in range(2)
            ]
        ).alias('id_struct')
        ).unnest('id_struct').drop('ID')
    data = data.with_columns([
      snp_sample['SNP'].alias('SNP'),
      snp_sample['Sample'].alias('Sample')])

    # make table of SNPs
    snp_table = data.filter(pl.col('SNP') == snp)
    
    # Record SNP assay code in the concordance pivot
    concordance_table = concordance_table.with_columns([
      concordance_table['Assay'].set_at_idx(i,snp_table['Assay'].unique().item()).alias('Assay')
    ])

    if snp == concordance_table['SNP'][0]:
        for animal in sorted(snp_table['Name'].unique()):
            call = snp_table.filter(pl.col('Name') == animal)['Final'].unique()[0]
            concordance_table = concordance_table.with_columns([
              pl.col('Assay').alias(animal)
            ])
            concordance_table = concordance_table.with_columns([
              concordance_table[animal].set_at_idx(i,call).alias(animal)
            ])
    elif len(concordance_table[:, 5:]) != len(data['Sample'].unique()) + 1:
        for animal in sorted(snp_table['Name'].unique()):
            call = snp_table.filter(pl.col('Name') == animal)['Final'].unique()[0]
            concordance_table = concordance_table.with_columns([
              concordance_table[animal].set_at_idx(i,call).alias(animal)
            ])
    else:
        break

# Change assay names to SNP names
concordance_table = concordance_table.with_columns(concordance_table['SNP'].str.replace(r"A", "SNP").alias('SNP'))

# Read in and format the sample manifest
samples = pl.read_excel(sample_manifest, read_csv_options={"skip_rows": 1}).drop_nulls(subset=['Animal ID', 'Relationship'])

# Parse out the relationship "groups"
samples = samples.select('Animal ID', 'Relationship').with_columns(
    samples['Relationship'].str.replace('Dam ', '').str.replace('Progeny ', '').cast(pl.Int64, strict=False).alias('group')
)

# Go through each group and add them to new vectors that will eventually become 
# new column headers for the concordance pivot
ordered_samples = []
ordered_relations = []
for i in samples['group'].unique():
    sub = samples.filter(pl.col('group') == i)
    ordered_samples.append(sub.filter(
      pl.col('Relationship').str.contains('Dam')
    ).select('Animal ID').item())
    ordered_samples.append(sub.filter(
      pl.col('Relationship').str.contains('Progeny')
    ).select('Animal ID').item())
    ordered_relations.append(sub.filter(
      pl.col('Relationship').str.contains('Dam')
    ).select('Relationship').item())
    ordered_relations.append(sub.filter(
      pl.col('Relationship').str.contains('Progeny')
    ).select('Relationship').item())

# Bind together the vectors as a data frame
animal_cols = concordance_table[:,4:]
for i in range(len(ordered_samples)):
  
  animal = ordered_samples[i]
  relation = ordered_relations[i]
  sub = pl.Series([animal, relation])
  sub = sub.append(animal_cols[animal], append_chunks=False)
  
  if i == 0:
    new_cols = pl.DataFrame({"new_col_0": sub})
  else:
    col = "new_col_"
    new_cols = new_cols.with_columns([
      sub.alias("{}{}".format(col, i))
    ])

# Prepare new rows for the concordance table
new_row1 = pl.DataFrame(["Animal Name", None, None, None]).transpose(include_header=False)
new_row2 = pl.DataFrame(concordance_table[:, :4].columns).transpose(include_header=False)

# combing new rows 1 and 2 and the non-animal-column portion of the concordance
# table with a row-bind
final_table = new_row1.vstack(new_row2)
non_call_rows = concordance_table[:,:4]
for i in range(len(non_call_rows.columns)):
  old_col = non_call_rows.columns[i]
  col = "column_"
  non_call_rows = non_call_rows.rename({old_col: "{}{}".format(col, i)})
  non_call_rows = non_call_rows.with_columns([
    non_call_rows[:,i].cast(pl.Utf8).alias("{}{}".format(col, i))
    ])
final_table = final_table.vstack(non_call_rows)

# tack on the animal-column portion of the concordance table with a column-bind
final_table = final_table.hstack(new_cols)

# Write output file for viewing in Excel
final_table.write_excel('concordance-table.xlsx', autofit=True, has_header=False)
