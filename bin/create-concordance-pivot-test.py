#!/usr/bin/env python3

import sys
import os
import polars as pl

# Check if the input TSV file is a symlink, and if it is, follow the symlink
cycle_table = os.path.realpath(sys.argv[1])
sample_manifest = os.path.realpath(sys.argv[2])

# Find available cycles
results_files = os.listdir('.')
available_cycles = []
for i in results_files:
    cycle = int(i.split('.')[0][:-1])
    available_cycles.append(cycle)

# Read in optimal cycle table
cycle_table = pl.read_csv(cycle_table)

# Read in the list of SNPs
SNPs = pl.read_excel(results_files[0], sheet='Sheet1', skip=15).select('ID').str_extract(r'-(.*)').alias('SNP')

# Create an empty data frame that will subsequently be filled for each SNP
concordance_table = pl.DataFrame({
    'SNP': SNPs['SNP'].unique().sort(),
    'Assay': None,
    'Cycle': None,
    'Comments': None
})

# Determine optimal cycles
cycle_table = cycle_table.sort('SNP')
concordance_table = concordance_table.with_column(
    'Cycle', cycle_table['cycle'][cycle_table['SNP'] == concordance_table['SNP']].astype(int)
).filter(pl.col('Cycle') != 0)

# Loop through each SNP to collate calls
for snp in concordance_table['SNP']:
    # Determine optimal cycle
    optimal_cycle = int(cycle_table['cycle'][cycle_table['SNP'] == snp])
    if optimal_cycle not in available_cycles:
        break

    # Read in the correct BioMark output file based on that cycle
    data = pl.read_excel(f"28637_BMX003_Detailed_Results_cycle{optimal_cycle}.xlsx", sheet='Sheet1', skip=15)
    data = data.with_column('sample', data['ID'].str_extract(r'(.*)-.*')).with_column('SNP', data['ID'].str_extract(r'.*-(.*)'))
  
    snp_table = data[data['SNP'] == snp]
    new_index = concordance_table['SNP'] == snp
  
    concordance_table = concordance_table.with_column(
        'Assay', snp_table['Assay'].unique()
    ).with_column(
        'Cycle', optimal_cycle
    )

    if snp == concordance_table['SNP'][0]:
        for animal in sorted(snp_table['Name'].unique()):
            call = snp_table['Final'][snp_table['Name'] == animal][0]
            concordance_table = concordance_table.with_column(animal, call)
    elif len(concordance_table[:, 5:]) != len(data['sample'].unique()) + 1:
        for animal in sorted(snp_table['Name'].unique()):
            call = snp_table['Final'][snp_table['Name'] == animal][0]
            concordance_table = concordance_table.with_column(animal, call)
    else:
        break

# Change assay names to SNP names
concordance_table = concordance_table.with_column('SNP', concordance_table['SNP'].str_replace('A', 'SNP'))

# Read in and format the sample manifest
samples = pl.read_excel(sample_manifest, skip=1).drop_nulls(subset=['Animal ID', 'Relationship'])

# Parse out the relationship "groups"
samples = samples.select('Animal ID', 'Relationship').with_column(
    'group', samples['Relationship'].str_remove('Dam ').str_remove('Progeny ').astype(pl.Int)
)

# Go through each group and add them to new vectors that will eventually become new column headers for the concordance pivot
ordered_samples = []
ordered_relations = []
for i in samples['group'].unique():
    sub = samples.filter(pl.col('group') == i)
    ordered_samples.extend(sub['Animal ID'][sub['Relationship'].str_contains('Dam')])
    ordered_relations.extend(sub['Relationship'][sub['Relationship'].str_contains('Dam')])
    ordered_samples.extend(sub['Animal ID'][sub['Relationship'].str_contains('Progeny')])
    ordered_relations.extend(sub['Relationship'][sub['Relationship'].str_contains('Progeny')])

# Bind together the vectors as a data frame
animal_cols = concordance_table.columns[4:]
new_cols = pl.DataFrame({
    'Animal ID': ordered_samples,
    None: ordered_relations,
    **{animal: animal_cols[animal] for animal in animal_cols}
})

# Prepare new rows for the concordance table
new_row1 = pl.DataFrame({
    'Animal ID': ['Animal ID'],
    None: [None],
    'SNP': [None],
    'Assay': [None],
    'Cycle': [None],
    'Comments': [None]
})

new_row2 = concordance_table.columns
concordance_table = new_row1.concat(new_row2).concat(concordance_table[:, :4]).concat(new_cols)

# Write output file for viewing in Excel
concordance_table.write_excel('concordance-table.xlsx', include_index=False)
