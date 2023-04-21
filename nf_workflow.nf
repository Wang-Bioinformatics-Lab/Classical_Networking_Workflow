#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_spectra = "data/spectra"

// Workflow Boiler Plate
params.OMETALINKING_YAML = "flow_filelinking.yaml"
params.OMETAPARAM_YAML = "job_parameters.yaml"

TOOL_FOLDER = "$baseDir/bin"

process filesummary {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file inputSpectra

    output:
    file 'summaryresult.tsv'

    """
    python $TOOL_FOLDER/scripts/filesummary.py $inputSpectra summaryresult.tsv $TOOL_FOLDER/binaries/msaccess
    """
}

process mscluster {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file inputSpectra

    output:
    file 'output/specs_ms.mgf'
    file 'output/clusterinfo.tsv'
    file 'output/clustersummary.tsv'

    """
    mkdir output
    python $TOOL_FOLDER/scripts/mscluster_wrapper.py $inputSpectra $TOOL_FOLDER/binaries spectra output
    """


}


workflow {
    input_spectra_ch = Channel.fromPath(params.input_spectra)
    filesummary(input_spectra_ch)
    mscluster(input_spectra_ch)
}