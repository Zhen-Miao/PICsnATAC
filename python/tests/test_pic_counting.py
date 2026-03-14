"""Tests for PIC counting - validates against known R output."""

import numpy as np
import pytest
from scipy import sparse

from picsnATAC import PIC_counting, load_fragments


class TestLoadFragments:
    def test_basic_loading(self, fragment_file, cells):
        df = load_fragments(fragment_file, cells, verbose=False)
        assert len(df) > 0
        assert set(df.columns) == {"chrom", "start", "end", "cell_barcode"}
        assert all(df["cell_barcode"].isin(set(cells)))

    def test_no_cells_found(self, fragment_file):
        with pytest.raises(ValueError, match="Cell barcodes not found"):
            load_fragments(fragment_file, ["NONEXISTENT-1"], verbose=False)


class TestPICCounting:
    def test_pic_counting_matches_r(self, fragment_file, cells, peaks):
        """Validate against known R output.

        R test expects:
            pic_mat@i == c(8, 8, 3, 1, 4, 6, 3, 7, 8, 4)  (0-based row indices)
            pic_mat@x == c(1, 3, 1, 1, 1, 1, 1, 1, 1, 1)   (values)
            colnames(pic_mat) == cells
        """
        result = PIC_counting(
            cells=cells,
            fragment_tsv_gz_file_location=fragment_file,
            peak_sets=peaks,
            verbose=False,
        )

        assert isinstance(result.matrix, sparse.csc_matrix)
        assert result.cell_names == cells
        assert result.matrix.shape[0] == len(peaks)
        assert result.matrix.shape[1] == len(cells)

        # R's dgCMatrix stores 0-based row indices in @i and values in @x
        # These are stored column by column in CSC format
        expected_i = np.array([8, 8, 3, 1, 4, 6, 3, 7, 8, 4])
        expected_x = np.array([1, 3, 1, 1, 1, 1, 1, 1, 1, 1])

        np.testing.assert_array_equal(result.matrix.indices, expected_i)
        np.testing.assert_array_equal(result.matrix.data, expected_x)

    def test_extend_size_negative(self, fragment_file, cells, peaks):
        with pytest.raises(ValueError, match="extend_size"):
            PIC_counting(
                cells=cells,
                fragment_tsv_gz_file_location=fragment_file,
                peak_sets=peaks,
                extend_size=-1,
            )

    def test_empty_cells(self, fragment_file, peaks):
        with pytest.raises(ValueError, match="cell names are empty"):
            PIC_counting(
                cells=[],
                fragment_tsv_gz_file_location=fragment_file,
                peak_sets=peaks,
            )
