"""Tests for DAR_by_LRT."""

import numpy as np
import pytest
from scipy import sparse

from picsnATAC.dar_test import DAR_by_LRT


class TestDARByLRT:
    def test_input_validation(self):
        """Should raise on mismatched inputs."""
        mat = sparse.random(10, 5, format="csc")
        rates = np.ones(4)  # wrong length
        labels = np.array(["A", "A", "B", "B", "B"])

        with pytest.raises(ValueError, match="Number of cells"):
            DAR_by_LRT(mat, rates, labels)

    def test_estimation_approach_validation(self):
        mat = sparse.random(10, 4, format="csc")
        rates = np.ones(4) * 0.5
        labels = np.array(["A", "A", "B", "B"])
        plen = np.ones(10) * 500

        with pytest.raises(ValueError, match="estimation_approach"):
            DAR_by_LRT(mat, rates, labels, plen=plen, estimation_approach="INVALID")
