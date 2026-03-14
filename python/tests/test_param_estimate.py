"""Tests for the EM parameter estimation."""

import numpy as np
import pytest
from scipy import sparse

from picsnATAC.param_estimate import get_r_by_ct_mat_pq


class TestGetRByCTMatPQ:
    def test_basic_em(self):
        """Basic EM should converge and return reasonable values."""
        np.random.seed(42)
        # Create a simple binary matrix: 100 peaks x 20 cells
        data = np.random.binomial(1, 0.3, size=(100, 20)).astype(float)
        mat = sparse.csc_matrix(data)

        labels = np.array(["A"] * 10 + ["B"] * 10)

        result = get_r_by_ct_mat_pq(
            cell_type_set=["A", "B"],
            r_by_c=mat,
            cell_type_labels=labels,
            n_features_per_cell=100,
            verbose=False,
        )

        assert "p_by_t" in result
        assert "q_vec" in result
        assert result["p_by_t"].shape == (100, 2)
        assert result["q_vec"].shape == (20,)
        # All probabilities should be between 0 and 1
        assert np.all(result["p_by_t"] >= 0)
        assert np.all(result["p_by_t"] <= 1)
        assert np.all(result["q_vec"] >= 0)
        assert np.all(result["q_vec"] <= 1)
