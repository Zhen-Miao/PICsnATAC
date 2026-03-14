"""Tests for theoretical distributions and insertion rate estimation."""

import numpy as np
import pytest

from picsnATAC.insertion_rates import (
    _get_theoretical_c1,
    _get_theoretical_c12,
    _insertion_to_c1,
)


class TestTheoreticalC1:
    def test_probabilities_sum_to_one(self):
        """Distribution should sum to approximately 1."""
        for rate in [0.1, 0.5, 1.0, 5.0, 10.0]:
            p = _get_theoretical_c1(rate, peak_length=500)
            assert abs(p.sum() - 1.0) < 1e-6, f"rate={rate}: sum={p.sum()}"

    def test_zero_rate_all_zero_counts(self):
        """With very low insertion rate, almost all probability at count=0."""
        p = _get_theoretical_c1(0.001, peak_length=500)
        assert p[0] > 0.999

    def test_output_length(self):
        p = _get_theoretical_c1(1.0, peak_length=500)
        assert len(p) == 6


class TestTheoreticalC12:
    def test_probabilities_reasonable(self):
        """Sum of probabilities for counts >= 1 should be <= 1."""
        p = _get_theoretical_c12(1.0, peak_length=500)
        assert p.sum() <= 1.0

    def test_output_length(self):
        p = _get_theoretical_c12(1.0, peak_length=500, cap_insertion=20)
        assert len(p) == 10  # cap_insertion // 2


class TestInsertionToC1:
    def test_matrix_shape(self):
        rates = np.array([0.1, 0.5, 1.0])
        lengths = np.array([200, 500, 1000])
        mat = _insertion_to_c1(rates, lengths)
        assert mat.shape == (3, 3)

    def test_monotonic_in_rate(self):
        """Expected count should increase with insertion rate."""
        rates = np.arange(0.1, 5.1, 0.1)
        lengths = np.array([500])
        mat = _insertion_to_c1(rates, lengths)
        assert np.all(np.diff(mat[:, 0]) >= 0)
