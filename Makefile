run:
	nextflow run ./nf_workflow.nf --resume -c nextflow.config

run_usi_download:
	nextflow run ./nf_workflow.nf --resume -c nextflow.config --download_usi_filename=./data/usi_files/input_files.tsv

run_transitive:
	nextflow run ./nf_workflow.nf --resume -c nextflow.config --topology=transitive