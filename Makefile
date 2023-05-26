run:
	nextflow run ./nf_workflow.nf --resume -c nextflow.config

run_transitive:
	nextflow run ./nf_workflow.nf --resume -c nextflow.config --topology=transitive