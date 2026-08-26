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
