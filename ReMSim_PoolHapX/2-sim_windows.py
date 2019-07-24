import csv
from pathlib import Path
import random
import sys
import time
from typing import Dict, List


def _create_chrom_dict(chrom_len_fpath: Path) -> Dict[str, int]:
    """
    Create a dictionary of chromosome names and their lengths.

    Parameters
    ----------
    chrom_len_fpath : pathlib.Path
        The tab delimited file containing a header, chromosome names in the first column, and the
        corresponding lengths in the second column.

    Returns
    -------
    dict of {str, int}
        A dictionary of chromosome names and lengths.
    """

    chrom_dict: Dict[str, int] = {}
    with chrom_len_fpath.open() as c_f:
        next(c_f)  # skip header
        line: str
        for line in c_f:
            chrom_dict[line.split("\t")[0]] = int(line.split("\t")[1])

    return chrom_dict


def main(num_sims: int, chrom_len_fpath: Path, region_len: int, output_fpath: Path) -> None:
    """
    Simulate ReMSim simulation windows.

    Parameters
    ----------
    num_sims : int
        Number of windows to simulate.

    chrom_len_fpath : pathlib.Path
        File containing information on chromosome name and chromosome length of the source organism
        for simulation.

    region_len : int
        The size of the simulation window.

    output_fpath : pathlib.Path
        Output file.
    """

    random.seed(time.time())  # seed the random number generator with current time
    chrom_dict: Dict[str, int] = _create_chrom_dict(chrom_len_fpath)
    header: List[str] = ["#sim", "chrom", "0-start", "0-end"]
    with output_fpath.open("w") as o_f:
        wtr = csv.writer(o_f, delimiter="\t")
        wtr.writerow(header)

        sim: int
        for sim in range(0, num_sims):
            chrom: str = random.choice(list(chrom_dict.keys()))
            chrom_len: int = chrom_dict[chrom]
            start: int = random.sample(range(0, chrom_len - region_len), 1)[0]
            wtr.writerow([sim, chrom, start, start + region_len - 1])


if __name__ == "__main__":
    # CLI input parameters.
    num_sims: int = int(sys.argv[1])
    chrom_len_fpath: Path = Path(sys.argv[2])
    region_len: int = int(sys.argv[3])
    output_fpath: Path = Path(sys.argv[4])

    main(num_sims, chrom_len_fpath, region_len, output_fpath)
