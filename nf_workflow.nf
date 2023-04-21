#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_spectra = "data/spectra"

params.min_cluster_size = "2"

params.pm_tolerance = "2.0"
params.fragment_tolerance = "0.5"

// Filtereing
params.min_peak_intensity = "0.0"

// Libraries
params.inputlibraries = ""

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
    python $TOOL_FOLDER/scripts/mscluster_wrapper.py \
    $inputSpectra $TOOL_FOLDER/binaries \
    spectra \
    output \
    --min_cluster_size $params.min_cluster_size \
    --pm_tolerance $params.pm_tolerance \
    --fragment_tolerance $params.fragment_tolerance \
    --min_peak_intensity $params.min_peak_intensity
    """
}

// Library Search
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
    python $TOOL_FOLDER/scripts/mscluster_wrapper.py \
    $inputSpectra $TOOL_FOLDER/binaries \
    spectra \
    output \
    --min_cluster_size $params.min_cluster_size \
    --pm_tolerance $params.pm_tolerance \
    --fragment_tolerance $params.fragment_tolerance \
    --min_peak_intensity $params.min_peak_intensity
    """
}

process librarySearchData {
    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    each file(input_library)
    each file(input_spectrum)

    output:
    file 'search_results/*' optional true

    """
    mkdir search_results
    python $TOOL_FOLDER/library_search_wrapper.py \
    $input_spectrum $input_library search_results \
    $TOOL_FOLDER/convert \
    $TOOL_FOLDER/main_execmodule.allcandidates
    """
}

process librarymergeResults {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path "results/*"

    output:
    path 'merged_results.tsv'

    """
    python $TOOL_FOLDER/tsv_merger.py \
    results merged_results.tsv
    """
}

process librarygetGNPSAnnotations {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path "merged_results.tsv"

    output:
    path 'merged_results_with_gnps.tsv'

    """
    python $TOOL_FOLDER/getGNPS_library_annotations.py \
    merged_results.tsv merged_results_with_gnps.tsv
    """
}

workflow {
    input_spectra_ch = Channel.fromPath(params.input_spectra)
    filesummary(input_spectra_ch)
    mscluster(input_spectra_ch)

    // Library Search
    libraries_ch = Channel.fromPath(params.inputlibraries + "/*.mgf" )
    search_results = librarySearchData(libraries, spectra)


}