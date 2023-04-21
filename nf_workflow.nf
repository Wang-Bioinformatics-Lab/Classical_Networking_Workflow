#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Spectra as input
params.input_spectra = "data/spectra"

// Libraries
params.inputlibraries = "data/library"

// Parameters
params.min_cluster_size = "2"

params.pm_tolerance = "2.0"
params.fragment_tolerance = "0.5"

// Filtereing
params.min_peak_intensity = "0.0"



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
process librarySearchData {
    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    each file(input_library)
    each file(input_spectrum)

    output:
    file 'search_results/*' optional true

    """
    mkdir search_results
    python $TOOL_FOLDER/scripts/library_search_wrapper.py \
    $input_spectrum $input_library search_results \
    $TOOL_FOLDER/binaries/convert \
    $TOOL_FOLDER/binaries/main_execmodule.allcandidates
    """
}

process librarymergeResults {
    publishDir "./nf_output/library_intermediate", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path "results/*"

    output:
    path 'merged_results.tsv'

    """
    python $TOOL_FOLDER/scripts/tsv_merger.py \
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
    python $TOOL_FOLDER/scripts/getGNPS_library_annotations.py \
    merged_results.tsv \
    merged_results_with_gnps.tsv
    """
}

workflow {
    input_spectra_ch = Channel.fromPath(params.input_spectra)

    filesummary(input_spectra_ch)

    (clustered_spectra_ch, clusterinfo_ch, clustersummary_ch) = mscluster(input_spectra_ch)

    // Library Search
    libraries_ch = Channel.fromPath(params.inputlibraries + "/*.mgf" )
    search_results_ch = librarySearchData(libraries_ch, clustered_spectra_ch)
    merged_results_ch = librarymergeResults(search_results_ch.collect())
    librarygetGNPSAnnotations(merged_results_ch)

}