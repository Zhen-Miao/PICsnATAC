"""PICsnATAC: Paired Insertion Counting for single-nucleus ATAC-seq data.

Python port of the R package PICsnATAC (Nature Methods 21.1, 2024: 32-36).
"""

from picsnATAC.pic_counting import (
    PIC_counting,
    PICResult,
    count_peaks,
    load_fragments,
)
from picsnATAC.insertion_rates import (
    obs_to_insertion_ME,
    obs_to_insertion_MLE_lam,
    obs_to_insertion_MLE_obj,
)
from picsnATAC.dar_test import DAR_by_LRT
from picsnATAC.param_estimate import get_r_by_ct_mat_pq

__version__ = "0.3.1"

__all__ = [
    "PIC_counting",
    "PICResult",
    "count_peaks",
    "load_fragments",
    "obs_to_insertion_ME",
    "obs_to_insertion_MLE_lam",
    "obs_to_insertion_MLE_obj",
    "DAR_by_LRT",
    "get_r_by_ct_mat_pq",
]
