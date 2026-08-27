"""
Density estimation and Kadane maximum subarray clustering module.

Implements:
1. Overlapping interval merging
2. Local density estimation and auto-penalty coefficient calculation
3. Inter-chunk logarithmic gap penalty function
4. Stack-based peel-off Kadane maximum subarray decomposition
5. Multi-satellite parallel/vectorized clustering orchestration
"""

import math
from typing import Any, Dict, List, Optional, Set, Tuple
import numpy as np
import polars as pl

from .config import PipelineConfig


def parse_unit_size_from_name(satellite_name: str) -> int:
    """Extract monomer unit size (in bp) from satellite family name.

    E.g. 'sat-3-176' -> 176.
    If parsing fails or result is non-positive, defaults to 100 bp.
    """
    parts = satellite_name.strip().split("-")
    if parts:
        last_part = parts[-1].strip()
        if last_part.isdigit():
            val = int(last_part)
            if val > 0:
                return val
    return 100


def merge_overlapping_intervals(
    starts: np.ndarray,
    ends: np.ndarray,
    sizes: Optional[np.ndarray] = None,
    identities: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]:
    """Merge overlapping/adjacent intervals on a single chromosome.

    Assumes input intervals are already sorted by start position ascending.

    Returns:
        (m_s, m_e, m_sz, m_raw_cnt, m_ident_bp, m_raw_bp)
        where:
        - m_s: start positions of merged chunks
        - m_e: end positions of merged chunks
        - m_sz: merged chunk sizes (m_e - m_s)
        - m_raw_cnt: count of raw intervals within each merged chunk
        - m_ident_bp: sum of (raw_chunk_size * identity_score) if identities present, else None
        - m_raw_bp: sum of raw_chunk_size if identities present, else None
    """
    n = len(starts)
    if n == 0:
        empty_int = np.empty(0, dtype=np.int64)
        return empty_int, empty_int, empty_int, empty_int, None, None

    if sizes is None:
        sizes = ends - starts

    has_ident = identities is not None

    m_s = []
    m_e = []
    m_sz = []
    m_raw_cnt = []
    m_ident_bp = [] if has_ident else None
    m_raw_bp = [] if has_ident else None

    curr_s = int(starts[0])
    curr_e = int(ends[0])
    curr_raw_cnt = 1
    curr_raw_bp = int(sizes[0])
    curr_ident_bp = float(sizes[0] * identities[0]) if has_ident else 0.0

    for i in range(1, n):
        s_i = int(starts[i])
        e_i = int(ends[i])
        sz_i = int(sizes[i])
        ident_i = float(identities[i]) if has_ident else 0.0

        if s_i <= curr_e:
            if e_i > curr_e:
                curr_e = e_i
            curr_raw_cnt += 1
            if has_ident:
                curr_ident_bp += sz_i * ident_i
                curr_raw_bp += sz_i
        else:
            m_s.append(curr_s)
            m_e.append(curr_e)
            m_sz.append(curr_e - curr_s)
            m_raw_cnt.append(curr_raw_cnt)
            if has_ident:
                m_ident_bp.append(curr_ident_bp)
                m_raw_bp.append(curr_raw_bp)

            curr_s = s_i
            curr_e = e_i
            curr_raw_cnt = 1
            curr_raw_bp = sz_i
            curr_ident_bp = sz_i * ident_i if has_ident else 0.0

    m_s.append(curr_s)
    m_e.append(curr_e)
    m_sz.append(curr_e - curr_s)
    m_raw_cnt.append(curr_raw_cnt)
    if has_ident:
        m_ident_bp.append(curr_ident_bp)
        m_raw_bp.append(curr_raw_bp)

    return (
        np.array(m_s, dtype=np.int64),
        np.array(m_e, dtype=np.int64),
        np.array(m_sz, dtype=np.int64),
        np.array(m_raw_cnt, dtype=np.int64),
        np.array(m_ident_bp, dtype=np.float64) if has_ident else None,
        np.array(m_raw_bp, dtype=np.int64) if has_ident else None,
    )


def compute_auto_coefficient(
    chrom_merged_chunks: List[Tuple[np.ndarray, np.ndarray]],
    window_bp: int = 100_000,
    min_coef: float = 0.1,
    max_coef: float = 1.0,
) -> Tuple[float, int, float]:
    """Calculate the automatic penalty coefficient from median local density across the genome.

    Args:
        chrom_merged_chunks: List of (m_s, m_e) tuples per chromosome.
        window_bp: Half-width of sliding window (default: 100,000 bp).
        min_coef: Minimum penalty coefficient floor (default: 0.1).
        max_coef: Maximum penalty coefficient cap (default: 1.0).

    Returns:
        (target_coef, total_chunks, median_density)
    """
    dens_list: List[float] = []

    for m_s, m_e in chrom_merged_chunks:
        m_count = len(m_s)
        if m_count == 0:
            continue

        # Build a prefix-sum of chunk sizes once per chromosome.
        # prefix[i+1] - prefix[i]  == m_e[i] - m_s[i]  (chunk size)
        # prefix[hi] - prefix[lo]  == total unclipped bp of chunks [lo, hi)
        # This lets each window query run in O(log N) instead of O(window width).
        chunk_sizes = (m_e - m_s).astype(np.int64)
        prefix = np.empty(m_count + 1, dtype=np.int64)
        prefix[0] = 0
        np.cumsum(chunk_sizes, out=prefix[1:])

        for m in range(m_count):
            mid = (int(m_s[m]) + int(m_e[m])) // 2
            w_start = mid - window_bp
            if w_start < 0:
                w_start = 0
            w_end = mid + window_bp
            w_len = w_end - w_start

            # Binary search for the half-open slice [lo, hi) of chunks that
            # overlap [w_start, w_end):
            #   lo: first chunk whose end  > w_start  (i.e. not entirely left)
            #   hi: first chunk whose start >= w_end  (i.e. entirely right)
            lo = int(np.searchsorted(m_e, w_start, side="right"))
            hi = int(np.searchsorted(m_s, w_end, side="left"))

            if lo >= hi:
                dens_list.append(0.0)
                continue

            # Sum bp for all chunks in [lo, hi) assuming no clipping …
            cov_bp = int(prefix[hi] - prefix[lo])

            # … then correct the two boundary chunks that may stick out.
            left_clip = int(w_start - m_s[lo])
            if left_clip > 0:
                cov_bp -= left_clip

            right_clip = int(m_e[hi - 1] - w_end)
            if right_clip > 0:
                cov_bp -= right_clip

            if cov_bp < 0:
                cov_bp = 0

            dens = (float(cov_bp) / float(w_len)) if w_len > 0 else 0.0
            dens_list.append(dens)

    total_chunks = len(dens_list)
    if total_chunks == 0:
        return 0.5, 0, 0.0

    dens_arr = np.sort(np.array(dens_list, dtype=np.float64))
    if total_chunks % 2 == 1:
        median_d = float(dens_arr[(total_chunks - 1) // 2])
    else:
        median_d = float(
            (dens_arr[total_chunks // 2 - 1] + dens_arr[total_chunks // 2]) / 2.0
        )

    if median_d >= 1.0 or (1.0 - median_d) <= 1e-12:
        raw_coef = 10.0
    else:
        raw_coef = median_d / (1.0 - median_d)

    # Clamp target_coef between min_coef and max_coef
    if raw_coef > max_coef:
        target_coef = max_coef
    elif raw_coef < min_coef:
        target_coef = min_coef
    else:
        target_coef = raw_coef

    return target_coef, total_chunks, median_d


def run_kadane_chromosome(
    m_s: np.ndarray,
    m_e: np.ndarray,
    m_sz: np.ndarray,
    m_raw_cnt: np.ndarray,
    m_ident_bp: Optional[np.ndarray],
    m_raw_bp: Optional[np.ndarray],
    max_gap: float,
    min_bp: float,
    start_gap: float,
    target_gap: float,
    target_coef: float,
    has_identity: bool,
) -> List[Dict[str, Any]]:
    """Execute Kadane's maximum subarray algorithm on a single chromosome's merged chunks.

    Returns:
        List of cluster dictionaries for this chromosome.
    """
    m_count = len(m_s)
    if m_count == 0:
        return []

    a_len = 2 * m_count - 1
    a_arr = np.zeros(a_len, dtype=np.float64)
    element_type = np.zeros(a_len, dtype=np.int8)  # 0: chunk, 1: gap
    chunk_ref_start = np.zeros(a_len, dtype=np.int32)
    chunk_ref_end = np.zeros(a_len, dtype=np.int32)

    scale_log = math.log10(target_gap / start_gap) if target_gap > start_gap else 1.0
    divisor = scale_log / target_coef if target_coef > 0 else 1.0

    for m in range(m_count):
        idx = 2 * m
        a_arr[idx] = float(m_sz[m])
        element_type[idx] = 0
        chunk_ref_start[idx] = m
        chunk_ref_end[idx] = m

        if m < m_count - 1:
            gap_idx = 2 * m + 1
            g = float(m_s[m + 1] - m_e[m])
            element_type[gap_idx] = 1

            if g > max_gap:
                a_arr[gap_idx] = -1e18
            else:
                if g <= start_gap:
                    penalty = 0.0
                else:
                    coef = (math.log10(g / start_gap)) / divisor
                    if coef > target_coef:
                        coef = target_coef
                    penalty = coef * g
                a_arr[gap_idx] = -penalty

    # Peel-off recursive decomposition via stack
    stack = [(0, a_len - 1)]
    found_clusters = []

    while stack:
        curr_lo, curr_hi = stack.pop()
        if curr_lo > curr_hi:
            continue

        max_so_far = -1e18
        max_ending_here = 0.0
        best_s = curr_lo
        best_e = curr_lo
        temp_s = curr_lo

        for i in range(curr_lo, curr_hi + 1):
            max_ending_here += a_arr[i]

            if max_ending_here > max_so_far:
                max_so_far = max_ending_here
                best_s = temp_s
                best_e = i

            if max_ending_here < 0:
                max_ending_here = 0.0
                temp_s = i + 1

        if max_so_far <= 0:
            continue

        # Trim leading and trailing gap elements
        while best_s <= best_e and element_type[best_s] == 1:
            best_s += 1
        while best_e >= best_s and element_type[best_e] == 1:
            best_e -= 1

        if best_s > best_e:
            continue

        m_start_idx = int(chunk_ref_start[best_s])
        m_end_idx = int(chunk_ref_end[best_e])

        found_clusters.append((m_start_idx, m_end_idx, max_so_far, int(m_s[m_start_idx])))

        if best_s - 1 >= curr_lo:
            stack.append((curr_lo, best_s - 1))
        if best_e + 1 <= curr_hi:
            stack.append((best_e + 1, curr_hi))

    # Sort found clusters by genomic start coordinate
    found_clusters.sort(key=lambda x: x[3])

    # Filter clusters by min_bp
    valid_clusters = []
    for ms, me, score, _ in found_clusters:
        ch_bp = int(np.sum(m_sz[ms : me + 1]))
        if ch_bp >= min_bp:
            valid_clusters.append((ms, me, score, ch_bp))

    n_valid = len(valid_clusters)
    results = []

    for i in range(n_valid):
        ms, me, score, ch_bp = valid_clusters[i]
        ch_start = int(m_s[ms])
        ch_end = int(m_e[me])
        n_raw = int(np.sum(m_raw_cnt[ms : me + 1]))

        # Compute max internal gap between consecutive merged chunks
        ch_max_gap = 0
        for m in range(ms + 1, me + 1):
            g = int(m_s[m] - m_e[m - 1])
            if g < 0:
                g = 0
            if g > ch_max_gap:
                ch_max_gap = g

        span = ch_end - ch_start
        # Union coverage of merged chunks equals sum of chunk lengths since they are disjoint
        dens = (float(ch_bp) / float(span)) if span > 0 else 0.0
        if dens > 1.0:
            dens = 1.0
        sc = int(dens * 1000.0 + 0.5) if span > 0 else 0
        if sc > 1000:
            sc = 1000

        gb = str(m_s[ms] - m_e[valid_clusters[i - 1][1]]) if i > 0 else "NA"
        ga = str(m_s[valid_clusters[i + 1][0]] - m_e[me]) if i < n_valid - 1 else "NA"

        entry: Dict[str, Any] = {
            "cluster_start": ch_start,
            "cluster_end": ch_end,
            "gap_before": gb,
            "gap_after": ga,
            "n_chunks": n_raw,
            "total_chunk_bp": ch_bp,
            "max_internal_gap": ch_max_gap,
            "kadane_score": score,
            "density": dens,
            "density_score": sc,
        }

        if has_identity:
            ch_ident_bp = float(np.sum(m_ident_bp[ms : me + 1])) if m_ident_bp is not None else 0.0
            ch_raw_bp = int(np.sum(m_raw_bp[ms : me + 1])) if m_raw_bp is not None else 0
            avg_ident = int(round(ch_ident_bp / float(ch_raw_bp))) if ch_raw_bp > 0 else 0
            if avg_ident > 1000:
                avg_ident = 1000
            entry["identity_score"] = avg_ident

        results.append(entry)

    return results


def process_dataset(
    bed_df: pl.DataFrame,
    valid_sequences: Set[str],
    config: PipelineConfig,
    has_identity: bool,
) -> pl.DataFrame:
    """Execute the full Kadane density-scoring clustering pipeline across all satellite families.

    Returns:
        Polars DataFrame matching the raw clusters schema expected by classification.
    """
    if len(bed_df) == 0 or not valid_sequences:
        schema = {
            "unit_id": pl.Utf8,
            "unit_size": pl.Int64,
            "chromosome": pl.Utf8,
            "cluster_start": pl.Int64,
            "cluster_end": pl.Int64,
            "gap_before": pl.Utf8,
            "gap_after": pl.Utf8,
            "cluster_id": pl.Utf8,
            "n_chunks": pl.Int64,
            "total_chunk_bp": pl.Float64,
            "max_internal_gap": pl.Int64,
            "max_sum": pl.Float64,
            "density": pl.Float64,
            "density_score": pl.Float64,
        }
        if has_identity:
            schema["identity_score"] = pl.Float64
        return pl.DataFrame(schema=schema)

    # Filter BED to valid entries
    filtered_df = bed_df.filter(pl.col("chromosome").is_in(list(valid_sequences)))
    if len(filtered_df) == 0:
        schema = {
            "unit_id": pl.Utf8,
            "unit_size": pl.Int64,
            "chromosome": pl.Utf8,
            "cluster_start": pl.Int64,
            "cluster_end": pl.Int64,
            "gap_before": pl.Utf8,
            "gap_after": pl.Utf8,
            "cluster_id": pl.Utf8,
            "n_chunks": pl.Int64,
            "total_chunk_bp": pl.Float64,
            "max_internal_gap": pl.Int64,
            "max_sum": pl.Float64,
            "density": pl.Float64,
            "density_score": pl.Float64,
        }
        if has_identity:
            schema["identity_score"] = pl.Float64
        return pl.DataFrame(schema=schema)

    # Extract unique satellite family names and filter by min_unit_size
    sat_names = filtered_df.select("satellite_name").unique().sort("satellite_name")["satellite_name"].to_list()

    selected_units = []
    for sat in sat_names:
        usize = parse_unit_size_from_name(sat)
        if usize >= config.min_unit_size:
            selected_units.append((sat, usize))

    if not selected_units:
        schema = {
            "unit_id": pl.Utf8,
            "unit_size": pl.Int64,
            "chromosome": pl.Utf8,
            "cluster_start": pl.Int64,
            "cluster_end": pl.Int64,
            "gap_before": pl.Utf8,
            "gap_after": pl.Utf8,
            "cluster_id": pl.Utf8,
            "n_chunks": pl.Int64,
            "total_chunk_bp": pl.Float64,
            "max_internal_gap": pl.Int64,
            "max_sum": pl.Float64,
            "density": pl.Float64,
            "density_score": pl.Float64,
        }
        if has_identity:
            schema["identity_score"] = pl.Float64
        return pl.DataFrame(schema=schema)

    all_raw_cluster_rows = []
    cluster_counter = 0

    for sat_name, unit_size in selected_units:
        sat_df = filtered_df.filter(pl.col("satellite_name") == sat_name)
        if len(sat_df) == 0:
            continue

        # Sort by chromosome natural order and start position
        # Group by chromosome
        chrom_groups = sat_df.partition_by("chromosome", as_dict=True)

        # Merge chunks per chromosome
        merged_per_chrom = {}
        chunks_for_coef = []

        for chr_key, group in chrom_groups.items():
            chr_name = chr_key[0] if isinstance(chr_key, tuple) else str(chr_key)
            # Sort by start coordinate ascending
            group_sorted = group.sort("start")
            starts_np = group_sorted["start"].to_numpy()
            ends_np = group_sorted["end"].to_numpy()
            sizes_np = ends_np - starts_np
            ident_np = group_sorted["identity_score"].to_numpy() if has_identity else None

            m_s, m_e, m_sz, m_raw_cnt, m_ident_bp, m_raw_bp = merge_overlapping_intervals(
                starts_np, ends_np, sizes_np, ident_np
            )
            merged_per_chrom[chr_name] = (m_s, m_e, m_sz, m_raw_cnt, m_ident_bp, m_raw_bp)
            chunks_for_coef.append((m_s, m_e))

        # Compute auto-coefficient
        if config.target_coef_override is not None:
            target_coef = float(config.target_coef_override)
        else:
            target_coef, _, _ = compute_auto_coefficient(
                chunks_for_coef,
                window_bp=config.window_bp,
                min_coef=config.min_coef,
                max_coef=config.max_coef,
            )

        # Run Kadane on each chromosome in order
        # Preserve chromosome order as version-sorted strings
        sorted_chr_names = sorted(merged_per_chrom.keys())

        for chr_name in sorted_chr_names:
            m_s, m_e, m_sz, m_raw_cnt, m_ident_bp, m_raw_bp = merged_per_chrom[chr_name]
            chr_clusters = run_kadane_chromosome(
                m_s,
                m_e,
                m_sz,
                m_raw_cnt,
                m_ident_bp,
                m_raw_bp,
                max_gap=config.max_gap,
                min_bp=config.min_bp,
                start_gap=config.start_gap,
                target_gap=config.target_gap,
                target_coef=target_coef,
                has_identity=has_identity,
            )

            for cl in chr_clusters:
                cluster_counter += 1
                row = {
                    "unit_id": sat_name,
                    "unit_size": unit_size,
                    "chromosome": str(chr_name),
                    "cluster_start": cl["cluster_start"],
                    "cluster_end": cl["cluster_end"],
                    "gap_before": cl["gap_before"],
                    "gap_after": cl["gap_after"],
                    "cluster_id": f"cluster_{cluster_counter}",
                    "n_chunks": cl["n_chunks"],
                    "total_chunk_bp": float(cl["total_chunk_bp"]),
                    "max_internal_gap": cl["max_internal_gap"],
                    "max_sum": float(cl["kadane_score"]),
                    "density": float(cl["density"]),
                    "density_score": float(cl["density_score"]),
                }
                if has_identity:
                    row["identity_score"] = float(cl.get("identity_score", 0))
                all_raw_cluster_rows.append(row)

    if not all_raw_cluster_rows:
        schema = {
            "unit_id": pl.Utf8,
            "unit_size": pl.Int64,
            "chromosome": pl.Utf8,
            "cluster_start": pl.Int64,
            "cluster_end": pl.Int64,
            "gap_before": pl.Utf8,
            "gap_after": pl.Utf8,
            "cluster_id": pl.Utf8,
            "n_chunks": pl.Int64,
            "total_chunk_bp": pl.Float64,
            "max_internal_gap": pl.Int64,
            "max_sum": pl.Float64,
            "density": pl.Float64,
            "density_score": pl.Float64,
        }
        if has_identity:
            schema["identity_score"] = pl.Float64
        return pl.DataFrame(schema=schema)

    return pl.DataFrame(all_raw_cluster_rows)
