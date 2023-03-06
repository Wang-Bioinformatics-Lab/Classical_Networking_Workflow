#!/usr/bin/env nextflow

params.input_spectra = "data/spectra"

// Workflow Boiler Plate
params.OMETALINKING_YAML = "flow_filelinking.yaml"
params.OMETAPARAM_YAML = "job_parameters.yaml"

TOOL_FOLDER = "$baseDir/bin"

process filesummary {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file inputSpectra from Channel.fromPath(params.input_spectra)

    output:
    file 'result.tsv' into records_ch

    """
    python $TOOL_FOLDER/scripts/filesummary.py $inputSpectra result.tsv $TOOL_FOLDER/binaries/msaccess
    """
}