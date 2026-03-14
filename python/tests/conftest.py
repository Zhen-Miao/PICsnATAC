"""Shared test fixtures for PICsnATAC tests."""

from pathlib import Path

import pandas as pd
import pytest

TEST_DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture
def cells():
    meta = pd.read_csv(TEST_DATA_DIR / "atac_singlecell_test_sub.csv")
    filtered = meta.loc[meta["is__cell_barcode"] == 1]
    return filtered["barcode"].tolist()


@pytest.fixture
def peaks():
    df = pd.read_csv(
        TEST_DATA_DIR / "atac_pbmc_5k_nextgem_peaks.bed",
        sep="\t",
        header=None,
        names=["chrom", "start", "end"],
        comment="#",
    )
    return df


@pytest.fixture
def fragment_file():
    return str(TEST_DATA_DIR / "atac_fragment_test.tsv.gz")
