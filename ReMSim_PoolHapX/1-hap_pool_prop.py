import csv
from pathlib import Path
import sys
from typing import List

import numpy as np


def main(num_haps: int, num_pools: int, output_fpath: Path) -> None:
    """
    Simulate haplotype proportions for each pool.

    Parameters
    ----------
    num_haps : int
        Number of haplotypes.

    num_pools : int
        Number of pools.

    output_fpath : pathlib.Path
        Output file.
    """

    header: List[str] = ["Hap_ID"] + ["h" + str(x) for x in range(0, num_haps)]
    with output_fpath.open("w") as o_f:
        wtr = csv.writer(o_f, delimiter="\t")
        wtr.writerow(header)

        pool: int
        for pool in range(0, num_pools):
            row: list = ["p" + str(pool)] \
                + np.random.dirichlet(np.ones(num_haps), size=1)[0].round(3).tolist()

            wtr.writerow(row)


if __name__ == "__main__":
    # CLI input parameters.
    num_haps: int = int(sys.argv[1])
    num_pools: int = int(sys.argv[2])
    output_fpath: Path = Path(sys.argv[3])

    main(num_haps, num_pools, output_fpath)
