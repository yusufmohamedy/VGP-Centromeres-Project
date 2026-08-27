"""
Unit tests for density_clustering module.
"""

import numpy as np
import pytest

from centromere_pipeline.config import PipelineConfig
from centromere_pipeline.density_clustering import (
    compute_auto_coefficient,
    merge_overlapping_intervals,
    parse_unit_size_from_name,
    run_kadane_chromosome,
)


def test_parse_unit_size():
    assert parse_unit_size_from_name("sat-3-176") == 176
    assert parse_unit_size_from_name("sat-1-240") == 240
    assert parse_unit_size_from_name("sat-13-40") == 40
    assert parse_unit_size_from_name("invalid_name") == 100
    assert parse_unit_size_from_name("sat-0-0") == 100


def test_merge_overlapping_intervals():
    starts = np.array([100, 150, 300, 350, 400], dtype=np.int64)
    ends = np.array([180, 220, 320, 410, 450], dtype=np.int64)
    # Expected merges:
    # 1. [100, 180] + [150, 220] -> [100, 220], size 120, count 2
    # 2. [300, 320] -> [300, 320], size 20, count 1
    # 3. [350, 410] + [400, 450] -> [350, 450], size 100, count 2

    m_s, m_e, m_sz, m_raw_cnt, _, _ = merge_overlapping_intervals(starts, ends)

    np.testing.assert_array_equal(m_s, [100, 300, 350])
    np.testing.assert_array_equal(m_e, [220, 320, 450])
    np.testing.assert_array_equal(m_sz, [120, 20, 100])
    np.testing.assert_array_equal(m_raw_cnt, [2, 1, 2])


def test_merge_with_identities():
    starts = np.array([100, 150], dtype=np.int64)
    ends = np.array([160, 200], dtype=np.int64)
    sizes = ends - starts  # 60, 50
    identities = np.array([90.0, 80.0], dtype=np.float64)

    m_s, m_e, m_sz, m_raw_cnt, m_ident_bp, m_raw_bp = merge_overlapping_intervals(
        starts, ends, sizes, identities
    )

    np.testing.assert_array_equal(m_s, [100])
    np.testing.assert_array_equal(m_e, [200])
    np.testing.assert_array_equal(m_sz, [100])
    np.testing.assert_array_equal(m_raw_cnt, [2])
    # Expected ident_bp = 60 * 90.0 + 50 * 80.0 = 5400 + 4000 = 9400.0
    # raw_bp = 60 + 50 = 110
    np.testing.assert_allclose(m_ident_bp, [9400.0])
    np.testing.assert_array_equal(m_raw_bp, [110])


def test_compute_auto_coefficient_bounds():
    # Sparse intervals
    chunks = [(np.array([1000]), np.array([1100]))]
    coef, count, median_d = compute_auto_coefficient(chunks, window_bp=10000, min_coef=0.1, max_coef=1.0)
    assert count == 1
    assert 0.1 <= coef <= 1.0

    # Empty chunks
    coef, count, median_d = compute_auto_coefficient([])
    assert count == 0
    assert coef == 0.5


def test_run_kadane_single_chunk():
    m_s = np.array([1000], dtype=np.int64)
    m_e = np.array([2000], dtype=np.int64)
    m_sz = np.array([1000], dtype=np.int64)
    m_raw_cnt = np.array([1], dtype=np.int64)

    clusters = run_kadane_chromosome(
        m_s, m_e, m_sz, m_raw_cnt, None, None,
        max_gap=1e7, min_bp=0, start_gap=1e3, target_gap=1e6, target_coef=0.5, has_identity=False
    )
    assert len(clusters) == 1
    cl = clusters[0]
    assert cl["cluster_start"] == 1000
    assert cl["cluster_end"] == 2000
    assert cl["total_chunk_bp"] == 1000
    assert cl["density_score"] == 1000  # 100% density
    assert cl["kadane_score"] == 1000.0


def test_run_kadane_hard_gap_split():
    # Two chunks separated by a massive gap exceeding max_gap (10 Mb)
    m_s = np.array([1000, 20_000_000], dtype=np.int64)
    m_e = np.array([2000, 20_001_000], dtype=np.int64)
    m_sz = np.array([1000, 1000], dtype=np.int64)
    m_raw_cnt = np.array([1, 1], dtype=np.int64)

    clusters = run_kadane_chromosome(
        m_s, m_e, m_sz, m_raw_cnt, None, None,
        max_gap=1e7, min_bp=0, start_gap=1e3, target_gap=1e6, target_coef=0.5, has_identity=False
    )
    # Must form two separate clusters, not merge across the hard gap
    assert len(clusters) == 2
    assert clusters[0]["cluster_start"] == 1000
    assert clusters[1]["cluster_start"] == 20_000_000


# ---------------------------------------------------------------------------
# Reference implementation of the OLD O(N²) density loop for comparison.
# Used only by the exactness tests below – do not call from production code.
# ---------------------------------------------------------------------------

def _compute_auto_coefficient_reference(chrom_merged_chunks, window_bp=100_000,
                                        min_coef=0.1, max_coef=1.0):
    """Verbatim copy of the original nested-loop implementation."""
    import numpy as np
    dens_list = []
    for m_s, m_e in chrom_merged_chunks:
        m_count = len(m_s)
        if m_count == 0:
            continue
        k_start = 0
        for m in range(m_count):
            mid = (int(m_s[m]) + int(m_e[m])) // 2
            w_start = mid - window_bp
            if w_start < 0:
                w_start = 0
            w_end = mid + window_bp
            w_len = w_end - w_start
            while k_start < m_count and m_e[k_start] <= w_start:
                k_start += 1
            cov_bp = 0
            for k in range(k_start, m_count):
                if m_s[k] >= w_end:
                    break
                o_s = max(m_s[k], w_start)
                o_e = min(m_e[k], w_end)
                if o_e > o_s:
                    cov_bp += int(o_e - o_s)
            dens = (float(cov_bp) / float(w_len)) if w_len > 0 else 0.0
            dens_list.append(dens)
    total_chunks = len(dens_list)
    if total_chunks == 0:
        return 0.5, 0, 0.0
    dens_arr = np.sort(np.array(dens_list, dtype=np.float64))
    if total_chunks % 2 == 1:
        median_d = float(dens_arr[(total_chunks - 1) // 2])
    else:
        median_d = float((dens_arr[total_chunks // 2 - 1] + dens_arr[total_chunks // 2]) / 2.0)
    if median_d >= 1.0 or (1.0 - median_d) <= 1e-12:
        raw_coef = 10.0
    else:
        raw_coef = median_d / (1.0 - median_d)
    if raw_coef > max_coef:
        target_coef = max_coef
    elif raw_coef < min_coef:
        target_coef = min_coef
    else:
        target_coef = raw_coef
    return target_coef, total_chunks, median_d


def _make_dense_chunks(n=500, spacing=500, chunk_len=200, start=0):
    """Generate n chunks spaced `spacing` bp apart, each of length `chunk_len`."""
    starts = np.array([start + i * spacing for i in range(n)], dtype=np.int64)
    ends = starts + chunk_len
    return starts, ends


def test_compute_auto_coefficient_exact_match_dense():
    """New O(N log N) implementation must return the same result as the old O(N²) one."""
    m_s, m_e = _make_dense_chunks(n=300, spacing=800, chunk_len=400)
    chunks = [(m_s, m_e)]

    ref_coef, ref_count, ref_med = _compute_auto_coefficient_reference(chunks, window_bp=100_000)
    new_coef, new_count, new_med = compute_auto_coefficient(chunks, window_bp=100_000)

    assert new_count == ref_count
    assert new_med == ref_med, f"median mismatch: {new_med} vs {ref_med}"
    assert new_coef == ref_coef, f"coef mismatch: {new_coef} vs {ref_coef}"


def test_compute_auto_coefficient_exact_match_sparse():
    """Sparse chunks (window mostly empty) – exactness check."""
    m_s, m_e = _make_dense_chunks(n=50, spacing=500_000, chunk_len=1_000)
    chunks = [(m_s, m_e)]

    ref_coef, ref_count, ref_med = _compute_auto_coefficient_reference(chunks, window_bp=100_000)
    new_coef, new_count, new_med = compute_auto_coefficient(chunks, window_bp=100_000)

    assert new_count == ref_count
    assert new_med == ref_med
    assert new_coef == ref_coef


def test_compute_auto_coefficient_exact_match_near_zero():
    """Chunks at position 0 (w_start clamped to 0) – exactness check."""
    m_s = np.array([0, 200, 500], dtype=np.int64)
    m_e = np.array([150, 350, 700], dtype=np.int64)
    chunks = [(m_s, m_e)]

    ref_coef, ref_count, ref_med = _compute_auto_coefficient_reference(chunks, window_bp=100_000)
    new_coef, new_count, new_med = compute_auto_coefficient(chunks, window_bp=100_000)

    assert new_count == ref_count
    assert new_med == ref_med
    assert new_coef == ref_coef


def test_compute_auto_coefficient_multi_chrom_exact():
    """Multiple chromosomes – exactness check across the aggregated density list."""
    chunks = [
        _make_dense_chunks(n=100, spacing=1_000, chunk_len=500, start=0),
        _make_dense_chunks(n=80,  spacing=2_000, chunk_len=300, start=10_000_000),
        _make_dense_chunks(n=20,  spacing=50_000, chunk_len=100, start=5_000_000),
    ]
    chunks_as_pairs = [(s, e) for s, e in chunks]

    ref = _compute_auto_coefficient_reference(chunks_as_pairs, window_bp=100_000)
    new = compute_auto_coefficient(chunks_as_pairs, window_bp=100_000)

    assert new[1] == ref[1]            # total count
    assert new[2] == ref[2]            # median density
    assert new[0] == ref[0]            # coefficient

