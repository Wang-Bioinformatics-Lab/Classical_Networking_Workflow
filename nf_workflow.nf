#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Spectra as input
params.input_spectra = "data/input_spectra"

// Libraries
params.input_libraries = "data/library"

// Metadata
params.redu_metadata_integration = "No" // Yes, means we go to ReDU and grab all the relevant metadata
params.metadata_per_file_grouping = "No" // Yes means that each file can be its own group
params.metadata_filename = "data/metadata.tsv"

// Clustering Parameters
params.min_cluster_size = "2"

// Tolerance Parameters
params.pm_tolerance = "2.0"
params.fragment_tolerance = "0.5"

// Filtereing
params.min_peak_intensity = "0.0"

// Molecular Networking Options
params.similarity = "gnps"
params.topology = "classic" // or can be transitive

params.parallelism = 24
params.networking_min_matched_peaks = 6
params.networking_min_cosine = 0.7
params.networking_max_shift = 1000

// Topology Filtering
params.topology_topk = 10
params.topology_maxcomponent = 100

params.topology_cliquemincosine = 0.7

// Spectral Filtering
params.massql_filter = "None"

// Library Search Parameters
params.library_topk = 1

params.library_min_cosine = 0.7
params.library_min_matched_peaks = 6

//TODO: Implement This
params.library_filter_precursor = 1
params.library_filter_window = 1

//TODO: Implement This
params.library_analog_search = "0"
params.library_analog_max_shift = 1999

// Workflow Boiler Plate
params.OMETALINKING_YAML = "flow_filelinking.yaml"
params.OMETAPARAM_YAML = "job_parameters.yaml"

// Downloading Files
params.download_usi_filename = params.OMETAPARAM_YAML // This can be changed if you want to run locally
params.cache_directory = "data/cache"

TOOL_FOLDER = "$baseDir/bin"

process filesummary {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file inputSpectra
    val ready

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
    val ready

    output:
    file 'clustering/specs_ms.mgf'
    file 'clustering/clusterinfo.tsv'
    file 'clustering/clustersummary.tsv'

    """
    mkdir clustering
    python $TOOL_FOLDER/scripts/mscluster_wrapper.py \
    $inputSpectra $TOOL_FOLDER/binaries \
    spectra \
    clustering \
    --min_cluster_size $params.min_cluster_size \
    --pm_tolerance $params.pm_tolerance \
    --fragment_tolerance $params.fragment_tolerance \
    --min_peak_intensity $params.min_peak_intensity
    """
}

// TODO: Finish Implementing this, as this is currently an no-op
process massqlFilterSpectra {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env_massql.yml"

    input:
    file inputSpectra

    output:
    file 'specs_ms_filtered.mgf'

    """
    python $TOOL_FOLDER/scripts/massql_filter_spectra.py \
    $inputSpectra \
    specs_ms_filtered.mgf \
    --massql "$params.massql_filter"
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
    $TOOL_FOLDER/binaries/main_execmodule.allcandidates \
    --pm_tolerance $params.pm_tolerance \
    --fragment_tolerance $params.fragment_tolerance \
    --topk $params.library_topk \
    --library_min_cosine $params.library_min_cosine \
    --library_min_matched_peaks $params.library_min_matched_peaks \
    --analog_search $params.library_analog_search
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
    publishDir "./nf_output/library", mode: 'copy'

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

// Molecular Networking
process networkingGNPSPrepParams {
    publishDir "./nf_output/networking", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file spectrum_file

    output:
    file "params/*"

    """
    mkdir params
    python $TOOL_FOLDER/scripts/prep_molecular_networking_parameters.py \
        "$spectrum_file" \
        "params" \
        --parallelism "$params.parallelism" \
        --min_matched_peaks "$params.networking_min_matched_peaks" \
        --ms2_tolerance "$params.fragment_tolerance" \
        --pm_tolerance "$params.pm_tolerance" \
        --min_cosine "$params.networking_min_cosine" \
        --max_shift "$params.networking_max_shift"
    """
}

process calculatePairs {
    publishDir "./nf_output/temp_pairs", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file spectrum_file
    each file(params_file)

    output:
    file "*_aligns.tsv" optional true

    """
    $TOOL_FOLDER/binaries/main_execmodule \
        ExecMolecularParallelPairs \
        "$params_file" \
        -ccms_INPUT_SPECTRA_MS2 $spectrum_file \
        -ccms_output_aligns ${params_file}_aligns.tsv
    """
}

// Creating the metadata file
process createMetadataFile {
    publishDir "./nf_output/metadata", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file input_metadata
    file usi_information
    file input_spectra_folder

    output:
    file "merged_metadata.tsv"

    //script in case its NO_FILE
    """
    python $TOOL_FOLDER/scripts/merge_metadata.py \
    $input_metadata \
    $usi_information \
    merged_metadata.tsv \
    --include_redu $params.redu_metadata_integration \
    --per_file_grouping $params.metadata_per_file_grouping \
    --spectra_folder $input_spectra_folder
    """
}

// Calculating the groupings
process calculateGroupings {
    publishDir "./nf_output/networking", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file input_metadata
    file input_clustersummary
    file input_clusterinfo

    output:
    file "clustersummary_with_groups.tsv"

    """
    python $TOOL_FOLDER/scripts/group_abundances.py \
    $input_clusterinfo \
    $input_clustersummary \
    $input_metadata \
    clustersummary_with_groups.tsv
    """
}

// Filtering the network, this is the classic way
process filterNetwork {
    publishDir "./nf_output/networking", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file input_pairs

    output:
    file "filtered_pairs.tsv"

    """
    python $TOOL_FOLDER/scripts/filter_networking_edges.py \
    $input_pairs \
    filtered_pairs.tsv \
    filtered_pairs_old_format.tsv \
    --top_k_val $params.topology_topk \
    --max_component_size $params.topology_maxcomponent
    """
}

process filterNetworkTransitive {
    publishDir "./nf_output/networking", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    cache false

    input:
    file input_pairs
    file input_spectra

    output:
    file "filtered_pairs.tsv"

    """
    python $TOOL_FOLDER/scripts/transitive_alignment.py \
    -c $input_spectra \
    -m $input_pairs \
    -p 30 \
    -th $params.topology_cliquemincosine \
    -r filtered_pairs.tsv \
    --minimum_score $params.networking_min_cosine
    """
}

// This takes the pairs and adds the component numbers to the table
process enrichNetworkEdges {
    publishDir "./nf_output/networking", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file input_pairs
    file input_clustersummary

    output:
    file "pairs_with_components.tsv"

    """
    python $TOOL_FOLDER/scripts/enrich_network_edges.py \
    $input_clustersummary \
    $input_pairs \
    pairs_with_components.tsv
    """
}

// Enriching the Cluster Summary
process enrichClusterSummary {
    publishDir "./nf_output/networking", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file input_clustersummary
    file input_filtered_pairs
    file input_library_matches

    output:
    file "clustersummary_with_network.tsv"

    """
    python $TOOL_FOLDER/scripts/enrich_cluster_summary.py \
    $input_clustersummary \
    $input_filtered_pairs \
    $input_library_matches \
    clustersummary_with_network.tsv
    """
}

process createNetworkGraphML {
    publishDir "./nf_output/networking", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    cache false

    input:
    file input_clustersummary
    file input_filtered_pairs
    file input_library_matches

    output:
    file "network.graphml"

    """
    python $TOOL_FOLDER/scripts/create_network_graphml.py \
    $input_clustersummary \
    $input_filtered_pairs \
    $input_library_matches \
    network.graphml
    """
}

// downloading all the files
process prepInputFiles {
    publishDir "$params.input_spectra", mode: 'copyNoFollow' // Warning, this is kind of a hack, it'll copy files back to the input folder
    
    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file input_parameters
    file cache_directory

    output:
    val true
    file "*.mzML" optional true
    file "*.mzXML" optional true
    file "*.mgf" optional true 

    """
    python $TOOL_FOLDER/scripts/download_public_data_usi.py \
    $input_parameters \
    . \
    output_summary.tsv \
    --cache_directory $cache_directory
    """
}

workflow {
    // Preps input spectrum files
    input_spectra_ch = Channel.fromPath(params.input_spectra)

    // Downloads input data
    usi_download_ch = Channel.fromPath(params.download_usi_filename)
    (_download_ready, _, _, _) = prepInputFiles(usi_download_ch, Channel.fromPath(params.cache_directory))

    // File summaries
    filesummary(input_spectra_ch, _download_ready)

    // Clustering
    (clustered_spectra_intermediate_ch, clusterinfo_ch, clustersummary_ch) = mscluster(input_spectra_ch, _download_ready)

    if(params.massql_filter != "None"){
        clustered_spectra_ch = massqlFilterSpectra(clustered_spectra_intermediate_ch)
    }
    else{
        clustered_spectra_ch = clustered_spectra_intermediate_ch
    }

    // Library Search
    libraries_ch = Channel.fromPath(params.input_libraries + "/*.mgf" )
    search_results_ch = librarySearchData(libraries_ch, clustered_spectra_ch)
    
    merged_results_ch = librarymergeResults(search_results_ch.collect())
    merged_results_ch = merged_results_ch.ifEmpty(file("NO_FILE"))

    gnps_library_results_ch = librarygetGNPSAnnotations(merged_results_ch)

    // Networking
    params_ch = networkingGNPSPrepParams(clustered_spectra_ch)
    networking_results_temp_ch = calculatePairs(clustered_spectra_ch, params_ch.collect())

    merged_networking_pairs_ch = networking_results_temp_ch.collectFile(name: "merged_pairs.tsv", storeDir: "./nf_output/networking", keepHeader: true)

    // Filtering the network
    if(params.topology == "classic"){
        filtered_networking_pairs_ch = filterNetwork(merged_networking_pairs_ch)
    }
    else if (params.topology == "transitive"){
        filtered_networking_pairs_ch = filterNetworkTransitive(merged_networking_pairs_ch, clustered_spectra_ch)
    }

    filtered_networking_pairs_enriched_ch = enrichNetworkEdges(filtered_networking_pairs_ch, clustersummary_ch)

    // Handling Metadata, if we don't have one, we'll set it to be empty
    if(params.metadata_filename.length() > 0){
        if(params.metadata_filename == "NO_FILE"){
            input_metadata_ch = Channel.of(file("NO_FILE"))
        }
        else{
            input_metadata_ch = Channel.fromPath(params.metadata_filename).first()
        }
    }
    else{
        input_metadata_ch = Channel.of(file("NO_FILE"))
    }

    
    // We will also include ReDU Metadata if desired
    merged_metadata_ch = createMetadataFile(input_metadata_ch, usi_download_ch, input_spectra_ch)
    
    
    // Enriching the network with group mappings
    clustersummary_with_groups_ch = calculateGroupings(merged_metadata_ch, clustersummary_ch, clusterinfo_ch)

    // Adding component and library informaiton
    clustersummary_with_network_ch = enrichClusterSummary(clustersummary_with_groups_ch, filtered_networking_pairs_enriched_ch, gnps_library_results_ch)

    // Creating the graphml Network
    createNetworkGraphML(clustersummary_with_network_ch, filtered_networking_pairs_enriched_ch, gnps_library_results_ch)

}