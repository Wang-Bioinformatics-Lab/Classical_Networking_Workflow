import numba
import numpy as np
from numba.typed import List
import os
import concurrent.futures
import pymzml
import time
import argparse

SHIFTED_OFFSET = 6000
TOPK = 10
MINMATCHES = 6
TOPPRODUCTS = 30000
ADJACENT_BINS = np.array([-1, 0, 1], dtype=np.int64)

SpectrumTuple = numba.types.Tuple([
    numba.float32[:],  # mz
    numba.float32[:],  # intensity
    numba.float32,  # precursor_mz
    numba.int32  # precursor_charge
])

def filter_peaks_optimized(mz_array, intensity_array, precursor_mz, precursor_charge):
    """Only apply sqrt transform and L2 normalization without filtering"""
    if len(mz_array) == 0:
        return (
            np.array([], dtype=np.float32),
            np.array([], dtype=np.float32),
            precursor_mz,
            precursor_charge
        )

    # Sort by m/z first
    sorted_idx = np.argsort(mz_array)
    mz_sorted = mz_array[sorted_idx].astype(np.float32)
    int_sorted = intensity_array[sorted_idx].astype(np.float32)
    int_sorted = np.sqrt(int_sorted)

    # Apply sqrt transform and L2 normalization
    norm = np.linalg.norm(int_sorted)
    if norm > 0:
        int_sorted /= norm

    return (mz_sorted, int_sorted, np.float32(precursor_mz), np.int32(precursor_charge))

def parse_mgf_file(path):
    results=[]
    with open(path,'r') as f:
        in_ions=False
        current={"mz":[],"int":[],"pepmass":None,"charge":None}
        for line in f:
            line=line.strip()
            if not line:
                continue
            if line=="BEGIN IONS":
                in_ions=True
                current={"mz":[],"int":[],"pepmass":None,"charge":None}
            elif line=="END IONS":
                if current["pepmass"] and current["charge"]:
                    pepmass_val=float(current["pepmass"])
                    cstr = current["charge"].rstrip("+")
                    if cstr.isdigit():
                        charge_val=int(cstr)
                    else:
                        charge_val=1
                    # Always create mz and intensity arrays, even if empty
                    mz_arr = np.array(current["mz"], dtype=np.float32)
                    in_arr = np.array(current["int"], dtype=np.float32)
                    filtered = filter_peaks_optimized(mz_arr, in_arr, pepmass_val, charge_val)
                    # Append to results regardless of mz/int being empty
                    results.append(filtered)
                in_ions=False
            else:
                if in_ions:
                    if "=" in line:
                        key,val=line.split("=",1)
                        key=key.lower()
                        if key=="pepmass":
                            val_parts=val.split()
                            current["pepmass"]=val_parts[0]
                        elif key=="charge":
                            current["charge"]=val
                    else:
                        parts=line.split()
                        if len(parts)==2:
                            try:
                                mz_v=float(parts[0])
                                in_v=float(parts[1])
                                current["mz"].append(mz_v)
                                current["int"].append(in_v)
                            except:
                                pass
    return results

def parse_mzml_file(path):
    run = pymzml.run.Reader(path)
    results=[]
    for spec in run:
        if spec.ms_level==2:
            if not spec.selected_precursors:
                continue
            mz_arr=np.array(spec.mz,dtype=np.float32)
            in_arr=np.array(spec.i,dtype=np.float32)
            precursor_mz=float(spec.selected_precursors[0]['mz'])
            precursor_charge=spec.selected_precursors[0].get('charge',1)
            if not isinstance(precursor_charge,int):
                try:
                    precursor_charge=int(precursor_charge)
                except:
                    precursor_charge=1
            filtered=filter_peaks_optimized(mz_arr,in_arr,precursor_mz,precursor_charge)
            if len(filtered[0])>0:
                results.append(filtered)
    return results

def parse_one_file(path):
    ext=os.path.splitext(path)[1].lower()
    if ext==".mgf":
        return parse_mgf_file(path)
    elif ext==".mzml":
        return parse_mzml_file(path)
    else:
        # fallback
        try:
            return parse_mzml_file(path)
        except:
            return parse_mgf_file(path)

def parse_files_in_parallel(file_paths, threads=1):
    if threads<=1 or len(file_paths)<=1:
        all_spectra=[]
        for fp in file_paths:
            parsed=parse_one_file(fp)
            all_spectra.extend(parsed)
        return all_spectra
    else:
        all_spectra=[]
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as ex:
            futs=[ex.submit(parse_one_file,fp) for fp in file_paths]
            for fut in concurrent.futures.as_completed(futs):
                all_spectra.extend(fut.result())
        return all_spectra


@numba.njit
def create_index(spectra, is_shifted, tolerance, shifted_offset):
    """Create index structure for fast lookups"""
    entries = List()
    for spec_idx in range(len(spectra)):
        spec = spectra[spec_idx]
        mz_arr = spec[0]
        precursor_mz = spec[2]

        for peak_idx in range(len(mz_arr)):
            mz = mz_arr[peak_idx]
            intensity = spec[1][peak_idx]

            if is_shifted:
                bin_val = np.int64(round(
                    (precursor_mz - mz + shifted_offset) / tolerance
                ))
            else:
                bin_val = np.int64(round(mz / tolerance))

            entries.append((bin_val, spec_idx, peak_idx, mz, intensity))

    entries.sort()
    return entries


@numba.njit
def find_bin_range(entries, target_bin):
    """Binary search for bin ranges"""
    left = 0
    right = len(entries)
    while left < right:
        mid = (left + right) // 2
        if entries[mid][0] < target_bin:
            left = mid + 1
        else:
            right = mid
    start = left

    right = len(entries)
    while left < right:
        mid = (left + right) // 2
        if entries[mid][0] <= target_bin:
            left = mid + 1
        else:
            right = mid
    return start, left


@numba.njit
def compute_all_pairs(spectra, shared_entries, shifted_entries, tolerance, threshold, scoring_func):
    results = List()
    n_spectra = len(spectra)

    for query_idx in range(n_spectra):
        query_spec = spectra[query_idx]
        upper_bounds = np.zeros(n_spectra, dtype=np.float32)
        match_counts = np.zeros(n_spectra, dtype=np.int32)

        # Process both shared and shifted peaks
        for peak_idx in range(len(query_spec[0])):
            mz = query_spec[0][peak_idx]
            intensity = query_spec[1][peak_idx]
            precursor_mz = query_spec[2]

            # Shared peaks processing
            shared_bin = np.int64(round(mz / tolerance))
            shifted_bin = np.int64(round((precursor_mz - mz + SHIFTED_OFFSET) / tolerance))

            # Check both shared and shifted entries
            for entries, bin_val in [(shared_entries, shared_bin),
                                     (shifted_entries, shifted_bin)]:
                for delta in ADJACENT_BINS:
                    target_bin = bin_val + delta
                    start, end = find_bin_range(entries, target_bin)
                    pos = start
                    while pos < end and entries[pos][1] <= query_idx:
                        pos += 1
                    # Find matches in this bin
                    while pos < end and entries[pos][0] == target_bin:
                        spec_idx = entries[pos][1]
                        upper_bounds[spec_idx] += intensity * entries[pos][4]
                        match_counts[spec_idx] += 1
                        pos += 1

        # Collect candidates using threshold parameter
        candidates = List()
        for spec_idx in range(query_idx + 1, n_spectra):
            if (upper_bounds[spec_idx] >= threshold and
                    match_counts[spec_idx] >= MINMATCHES):
                candidates.append((spec_idx, upper_bounds[spec_idx]))

        # Sort candidates by upper bound score
        candidates.sort(key=lambda x: -x[1])

        # Process top candidates for exact matching
        exact_matches = List()
        for spec_idx, _ in candidates[:TOPPRODUCTS * 2]:
            target_spec = spectra[spec_idx]
            score, shared, shifted = scoring_func(spectra[query_idx], target_spec,tolerance)
            if score >= threshold:
                exact_matches.append((spec_idx, score, shared, shifted))

        # Sort and store top results
        exact_matches.sort(key=lambda x: -x[1])
        results.append((query_idx, exact_matches[:TOPPRODUCTS]))

    return results

@numba.njit(fastmath=True)
def calculate_exact_score_GNPS(query_spec, target_spec, TOLERANCE):
    """Numba-optimized cosine scoring with shift handling"""
    q_mz = query_spec[0]
    q_int = query_spec[1]
    q_pre = query_spec[2]
    q_charge = query_spec[3]

    t_mz = target_spec[0]
    t_int = target_spec[1]
    t_pre = target_spec[2]

    # Calculate precursor mass difference (assuming charge=1)
    precursor_mass_diff = (q_pre - t_pre)*q_charge
    allow_shift = True
    fragment_tol = TOLERANCE

    # Pre-allocate arrays for matches (adjust size as needed)
    max_matches = len(q_mz) * 2  # Estimate maximum possible matches
    scores_arr = np.zeros(max_matches, dtype=np.float32)
    idx_q = np.zeros(max_matches, dtype=np.int32)
    idx_t = np.zeros(max_matches, dtype=np.int32)
    match_count = 0

    # For each peak in query spectrum
    for q_idx in range(len(q_mz)):
        q_mz_val = q_mz[q_idx]
        q_int_val = q_int[q_idx]

        # For each possible shift (charge=1)
        num_shifts = 1
        if allow_shift and abs(precursor_mass_diff) >= fragment_tol:
            num_shifts += 1

        for shift_idx in range(num_shifts):
            if shift_idx == 0:
                # No shift
                adjusted_mz = q_mz_val
            else:
                # Apply shift
                adjusted_mz = q_mz_val - precursor_mass_diff

            # Find matching peaks in target using binary search
            start = np.searchsorted(t_mz, adjusted_mz - fragment_tol)
            end = np.searchsorted(t_mz, adjusted_mz + fragment_tol)

            for t_idx in range(start, end):
                if match_count >= max_matches:
                    break  # Prevent overflow
                if abs(t_mz[t_idx] - adjusted_mz)<fragment_tol:
                    # Store match information
                    scores_arr[match_count] = q_int_val * t_int[t_idx]
                    idx_q[match_count] = q_idx
                    idx_t[match_count] = t_idx
                    match_count += 1

    # Sort matches by score descending using argsort
    if match_count == 0:
        return 0.0, 0.0, 0.0

    # Get valid matches
    valid_scores = scores_arr[:match_count]
    valid_q_idx = idx_q[:match_count]
    valid_t_idx = idx_t[:match_count]

    # Argsort descending
    sort_order = np.argsort(-valid_scores)

    # Track used peaks
    q_used = np.zeros(len(q_mz), dtype=np.bool_)
    t_used = np.zeros(len(t_mz), dtype=np.bool_)
    total = 0.0
    shared = 0.0
    shifted = 0.0

    # Accumulate top matches
    for i in sort_order:
        q_idx = valid_q_idx[i]
        t_idx = valid_t_idx[i]

        if not q_used[q_idx] and not t_used[t_idx]:
            score = valid_scores[i]
            total += score

            # Determine match type
            mz_diff = abs(q_mz[q_idx] - t_mz[t_idx])
            if mz_diff <= fragment_tol:
                shared += score
            else:
                shifted += score

            q_used[q_idx] = True
            t_used[t_idx] = True

    return total, shared, shifted


@numba.njit
def calculate_exact_score_GNPS_multi_charge(query_spec, target_spec, TOLERANCE):
    """Numba-optimized cosine scoring with shift handling based on original GNPS code."""
    q_mz = query_spec[0]
    q_int = query_spec[1]
    q_pre = query_spec[2]
    q_charge = max(query_spec[3], 1)  # Ensure charge is at least 1

    t_mz = target_spec[0]
    t_int = target_spec[1]
    t_pre = target_spec[2]
    t_charge = max(target_spec[3], 1)  # Not used in calculation

    allow_shift = True
    fragment_tol = TOLERANCE

    # Compute precursor mass difference based on query's charge
    precursor_mass_diff = (q_pre - t_pre) * q_charge
    num_shifts = 1
    if allow_shift and abs(precursor_mass_diff) >= fragment_tol:
        num_shifts += q_charge

    # Pre-allocate arrays for matches
    max_matches = len(q_mz) * num_shifts * 10  # Estimate, adjust as needed
    scores_arr = np.zeros(max_matches, dtype=np.float32)
    idx_q = np.zeros(max_matches, dtype=np.int32)
    idx_t = np.zeros(max_matches, dtype=np.int32)
    is_shifted = np.zeros(max_matches, dtype=np.bool_)
    match_count = 0

    for q_idx in range(len(q_mz)):
        q_mz_val = q_mz[q_idx]
        q_int_val = q_int[q_idx]

        for cpi in range(num_shifts):
            # Calculate mass difference for this shift
            if cpi == 0:
                mass_diff_val = 0.0
            else:
                mass_diff_val = precursor_mass_diff / cpi

            adjusted_query_mz = q_mz_val - mass_diff_val
            start = np.searchsorted(t_mz, adjusted_query_mz - fragment_tol)
            end = np.searchsorted(t_mz, adjusted_query_mz + fragment_tol)

            # Collect all target peaks in this range
            for t_idx in range(start, end):
                if match_count >= max_matches:
                    break  # Prevent overflow

                scores_arr[match_count] = q_int_val * t_int[t_idx]
                idx_q[match_count] = q_idx
                idx_t[match_count] = t_idx
                is_shifted[match_count] = cpi != 0
                match_count += 1

    # Handle cases with no matches
    if match_count == 0:
        return 0.0, 0.0, 0.0

    # Extract valid entries
    valid_scores = scores_arr[:match_count]
    valid_q_idx = idx_q[:match_count]
    valid_t_idx = idx_t[:match_count]
    valid_shifted = is_shifted[:match_count]

    # Sort matches by descending score
    sort_order = np.argsort(-valid_scores)

    # Track used peaks
    q_used = np.zeros(len(q_mz), dtype=np.bool_)
    t_used = np.zeros(len(t_mz), dtype=np.bool_)
    total = 0.0
    shared = 0.0
    shifted = 0.0

    # Accumulate scores greedily
    for i in sort_order:
        q_i = valid_q_idx[i]
        t_i = valid_t_idx[i]
        if not q_used[q_i] and not t_used[t_i]:
            score = valid_scores[i]
            total += score
            if valid_shifted[i]:
                shifted += score
            else:
                shared += score
            q_used[q_i] = True
            t_used[t_i] = True

    return total, shared, shifted

@numba.njit
def calculate_exact_score(query_spec, target_spec,TOLERANCE):
    """Calculate exact cosine similarity with tolerance checks"""
    # Access tuple elements by index
    q_mz = query_spec[0]  # Already sorted by mz
    q_int = query_spec[1]
    q_pre = query_spec[2]

    t_mz = target_spec[0]  # Already sorted by mz
    t_int = target_spec[1]
    t_pre = target_spec[2]

    q_used = np.zeros(len(q_mz), dtype=np.bool_)
    t_used = np.zeros(len(t_mz), dtype=np.bool_)
    total = 0.0
    shared = 0.0
    shifted = 0.0

    # 1. Shared peak matching
    for q_idx in range(len(q_mz)):
        if q_used[q_idx]:
            continue

        mz_q = q_mz[q_idx]
        # Find matches using binary search on sorted mz
        start = np.searchsorted(t_mz, mz_q - TOLERANCE)
        end = np.searchsorted(t_mz, mz_q + TOLERANCE)

        best_score = 0.0
        best_t_idx = -1
        for t_idx in range(start, end):
            if not t_used[t_idx]:
                current_score = q_int[q_idx] * t_int[t_idx]
                if current_score > best_score:
                    best_score = current_score
                    best_t_idx = t_idx

        if best_t_idx != -1:
            total += best_score
            shared += best_score
            q_used[q_idx] = True
            t_used[best_t_idx] = True

    # 2. Shifted peak matching (critical fix)
    # Create and sort shifted mz arrays while tracking original indices
    # ---------------------------------------------------------------
    # For query
    q_shifted = q_pre - q_mz + SHIFTED_OFFSET
    q_shifted_sorted_idx = np.argsort(q_shifted)
    q_shifted_sorted = q_shifted[q_shifted_sorted_idx]

    # For target
    t_shifted = t_pre - t_mz + SHIFTED_OFFSET
    t_shifted_sorted_idx = np.argsort(t_shifted)
    t_shifted_sorted = t_shifted[t_shifted_sorted_idx]
    # ---------------------------------------------------------------

    # Match using sorted shifted mz arrays
    for q_pos in range(len(q_shifted_sorted)):
        q_orig_idx = q_shifted_sorted_idx[q_pos]
        if q_used[q_orig_idx]:
            continue

        mz_q = q_shifted_sorted[q_pos]
        start = np.searchsorted(t_shifted_sorted, mz_q - TOLERANCE)
        end = np.searchsorted(t_shifted_sorted, mz_q + TOLERANCE)

        best_score = 0.0
        best_t_pos = -1
        for t_pos in range(start, end):
            t_orig_idx = t_shifted_sorted_idx[t_pos]
            if not t_used[t_orig_idx]:
                current_score = q_int[q_orig_idx] * t_int[t_orig_idx]
                if current_score > best_score:
                    best_score = current_score
                    best_t_pos = t_pos

        if best_t_pos != -1:
            t_orig_idx = t_shifted_sorted_idx[best_t_pos]
            total += best_score
            shifted += best_score
            q_used[q_orig_idx] = True
            t_used[t_orig_idx] = True

    return total, shared, shifted


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Nextflow-compatible spectral networking"
    )
    parser.add_argument("-t", "--input_files", nargs="+", required=True,
                        help="Input MGF/mzML files")
    parser.add_argument("-o", "--output", required=True,
                        help="Output file path")
    parser.add_argument("--tolerance", type=float, default=0.01,
                        help="Fragment mass tolerance")
    parser.add_argument("--threshold", type=float, default=0.7,
                        help="Similarity score threshold")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of processing threads")
    parser.add_argument("--multi_charge", action="store_true",
                        help="Use GNPS multi-charge scoring")
    return parser.parse_args()


def main():
    args = parse_arguments()

    # Parse spectra
    spectra_list = parse_files_in_parallel(args.input_files, args.threads)
    print(f"Parsed {len(spectra_list)} spectra")

    # Convert to numba-compatible format
    numba_spectra = List()
    for spec in spectra_list:
        numba_spectra.append((
            spec[0].astype(np.float32),
            spec[1].astype(np.float32),
            np.float32(spec[2]),
            np.int32(spec[3])
        ))

    # Build indexes
    print("Building indexes...")
    shared_idx = create_index(numba_spectra, False, args.tolerance, SHIFTED_OFFSET)
    shifted_idx = create_index(numba_spectra, True, args.tolerance, SHIFTED_OFFSET)

    if args.multi_charge:
        scoring_func = calculate_exact_score_GNPS_multi_charge
    else:
        scoring_func = calculate_exact_score_GNPS
    # Compute matches
    print("Computing pairs...")
    matches = compute_all_pairs(numba_spectra, shared_idx, shifted_idx,
                                args.tolerance, args.threshold, scoring_func)

    # Write output in Nextflow-compatible format
    print("Writing results...")
    with open(args.output, 'w') as f:
        f.write("CLUSTERID1\tCLUSTERID2\tDeltaMZ\tMinRatio\tCosine\tAlignScore2\tAlignScore3\n")
        for query_idx, candidates in matches:
            q_spec = numba_spectra[query_idx]
            for cand in candidates:
                t_spec = numba_spectra[cand[0]]

                cluster1 = query_idx + 1
                cluster2 = cand[0] + 1
                delta_mz = q_spec[2] - t_spec[2]
                cosine = cand[1]

                f.write(f"{cluster1}\t{cluster2}\t"
                        f"{delta_mz:.3f}\t0.000\t{cosine:.4f}\t"
                        f"0.000\t0.000\n")

    print("Done!")


if __name__ == "__main__":
    main()