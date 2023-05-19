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

		COLLATE_RESULTS (
			SELECT_CYCLES.out.cycles,
			ch_snp_tables
		)

	} else {

		// input channels
		ch_snp_tables = Channel
			.fromPath( "${params.input_dir}/*.xlsx" )
			.collect()
		
		ch_cycles = Channel
			.fromPath( params.cycles )

		
		// Workflow steps 
		COLLATE_RESULTS (
			ch_cycles,
			ch_snp_tables
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
	
	script:
	"""
	select-optimal-cycle.R
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
	create-concordance-pivot.R ${cycles}
	"""
}

// --------------------------------------------------------------- //
