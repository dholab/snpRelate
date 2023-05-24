#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	if ( params.cycles == "" ){

		// input channels
		ch_snp_tables = Channel
			.fromPath( "${params.input_dir}/*.xlsx" )
			.collect()

		
		// Workflow steps 
		SELECT_CYCLES (
			ch_snp_tables
		)

		PLOT_BEST_CYCLES (
			SELECT_CYCLES.out.cycles,
			ch_snp_tables
		)

		COLLATE_RESULTS (
			SELECT_CYCLES.out.cycles,
			ch_snp_tables
		)

		PREDICT_RELATEDNESS (
			COLLATE_RESULTS.out
		)

	} else {

		// input channels
		ch_snp_tables = Channel
			.fromPath( "${params.input_dir}/*.xlsx" )
			.collect()
		
		ch_cycles = Channel
			.fromPath( params.cycles )

		
		// Workflow steps
		PLOT_BEST_CYCLES (
			ch_cycles,
			ch_snp_tables
		)

		COLLATE_RESULTS (
			ch_cycles,
			ch_snp_tables
		)

		PREDICT_RELATEDNESS (
			COLLATE_RESULTS.out
		)

	}
	
}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config

// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process SELECT_CYCLES {
	
	publishDir params.results, mode: 'copy'
	
	input:
	path excel_files
	
	output:
	path "*.csv", emit: cycles

	when:
	params.cycles == ""
	
	script:
	"""
	select-optimal-cycle.R
	"""
}

process PLOT_BEST_CYCLES {

	publishDir params.results, mode: 'copy'
	
	input:
	path cycles
	path excel_files

	output:
	path "*.pdf"

	when:
	params.in_dev == false

	script:
	"""
	echo "in development"
	"""

}

process COLLATE_RESULTS {
	
	publishDir params.results, mode: 'copy'
	
	input:
	path cycles
	path excel_files
	
	output:
	path "*.xlsx"
	
	script:
	"""
	create-concordance-pivot.py ${cycles} ${params.manifest}
	"""

}

process PREDICT_RELATEDNESS {

	input:
	path concordance_table

	output:
	path "*.xlsx"

	when:
	params.in_dev == false

	script:
	"""
	echo "in development"
	"""

}

// --------------------------------------------------------------- //
