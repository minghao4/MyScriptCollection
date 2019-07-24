import csv
from multiprocessing import Pool as Pool_f
from multiprocessing.pool import Pool
import os
from pathlib import Path
import sys
from typing import Dict, List, Tuple


def _create_sim_region(sim: str, regions_fpath: Path) -> Tuple[str, int]:
    """
    Create a tuple of the chromosome name and the 0-based simulation window start position.

    Parameters
    ----------
    sim : str
        The simulation number/name.

    regions_fpath : pathlib.Path
        File containing simulation windows.

    Returns
    -------
    tuple of {str, int}
        A tuple of the chromosome name and the 0-based start of the simulation window.
    """

    with regions_fpath.open() as r_f:
        rdr = csv.reader(r_f, delimiter="\t")
        next(rdr)  # skip header
        line: List[str]
        for line in rdr:
            if line[0] == sim[3:]:
                return (line[1], int(line[2]))


def _create_chrom_files_dict(chrom_fpaths_fpath: Path) -> Dict[int, Path]:
    """
    Create a dictionary of chromosome names and the location of their fasta files.

    Parameters
    ----------
    chrom_fpaths_fpath : pathlib.Path
        A tab delimited file containing a header, chromosome names in the first column, and the
        system filepaths in the second column.

    Returns
    -------
    dict of {int, pathlib.Path}
        A dictionary of the chromosome name and the fasta system filepath.
    """

    chrom_files_dict: Dict[int, Path] = {}
    with chrom_fpaths_fpath.open() as cf_f:
        next(cf_f)  # skip header
        for line in cf_f:
            chrom_files_dict[line.split("\t")[0]] = Path(line.split("\t")[1])

    return chrom_files_dict


def _remsim(
    remsim_fpath: Path,
    config_fpath: Path,
    ref_fpath: Path,
    source: str,
    read_count: int,
    output_dir: Path,
    start: int,
    prefix: str,
) -> None:
    """
    Construct the command for running ReMSim with the given parameters and execute the system call.

    Parameters
    ----------
    remsim_fpath : pathlib.Path
        ReMSim.

    config_fpath : pathlib.Path
        The ReMSim configuration JSON file.

    ref_fpath : pathlib.Path
        The chromosome fasta file for ReMSim to simulate from.

    source : str
        The source organism name.

    read_count : int
        The number of reads for ReMSim to simulate.

    output_dir : Path
        The output directory for ReMSim.

    start : int
        The 0-based start position of the ReMSim simulation window.

    prefix : str
        The file prefix of the ReMSim output files.
    """

    cmd: str = "python {} -c {} -r {} -f simulate -sr {} -rc {} -o {} -s {} -p {}".format(
        remsim_fpath, config_fpath, ref_fpath, source, read_count, output_dir, start, prefix
    )
    os.system(cmd)


def _pool_processes(
    mp_pool: Pool,
    pool_hap_props: List[str],
    total_read_count: int,
    configs_dpath: Path,
    ref_fpath: Path,
    output_dir: Path,
    start: int,
    prefix: str,
) -> None:
    """
    Given a multiprocessing pool and the haplotype proportions for a particular pool, simulate reads
    using ReMSim.

    Parameters
    ----------
    mp_pool : multiprocessing.pool.Pool
        A multiprocessing pool.

    pool_hap_props : list of str
        List of haplotype proportions for the current pool.

    total_read_count : int
        Total read count of all haplotypes combined.

    configs_dpath : pathlib.Path
        Directory holding all haplotype configuration JSONs for ReMSim.

    ref_fpath : pathlib.Path
        Chromosome fasta file used by ReMSim.

    output_dir : pathlib.Path
        ReMSim output directory.

    start : int
        0-based start position of ReMSim simulation window.

    prefix : str
        ReMSim output file prefix.
    """

    j: int
    hap_prop: str
    for j, hap_prop in enumerate(pool_hap_props):
        config_fpath: Path = list(configs_dpath.glob("{}*.json".format(j)))[0]
        source: str = config_fpath.name[:-5]
        read_count: int = int(float(hap_prop) * total_read_count)
        mp_pool.apply_async(
            _remsim, args=(config_fpath, ref_fpath, source, read_count, output_dir, start, prefix)
        )


def main(
    sim_prop_fpath: Path,
    regions_fpath: Path,
    chrom_fpaths_fpath: Path,
    total_read_count: int,
    configs_dpath: Path,
    remsim_fpath: Path,
    outputs_dpath: Path,
    NUM_CPUS: int,
) -> None:
    """
    For each pool, simulate a proportional number of reads from each haplotype adding up to a total
    read count.

    Parameters
    ----------
    sim_prop_fpath : pathlib.Path
        File containing haplotype proportions for each pool for the current simulation number.

    regions_fpath : pathlib.Path
        File containing simulation window start and end positions.

    chrom_fpaths_fpath : pathlib.Path

    total_read_count : int
        Total number of reads for ReMSim to simulate.

    configs_dpath : pathlib.Path
        Directory containing ReMSim configuration JSON files for all haplotypes.

    remsim_fpath : pathlib.Path
        System filepath of remsim.py

    outputs_dpath : pathlib.Path
        Directory containing output directories for ReMSim.

    NUM_CPUS : int
        Number of workers available for the multiprocess pool.
    """

    sim: str = sim_prop_fpath.name.split("_")[0]
    sim_chrom_start: Tuple[int, int] = _create_sim_region(sim, regions_fpath)
    chrom_files_dict: Dict[str, Path] = _create_chrom_files_dict(chrom_fpaths_fpath)

    with sim_prop_fpath.open() as sp_f:
        rdr = csv.reader(sp_f, delimiter="\t")
        next(rdr)  # skip header

        i: int
        pool_hap_props: List[str]
        for i, pool_hap_props in enumerate(rdr):
            ref_fpath: Path = chrom_files_dict[sim_chrom_start[0]]
            start: int = sim_chrom_start[1]
            prefix: str = "{}p{}".format(sim, i)
            output_dir: Path = Path(outputs_dpath, sim)
            mp_pool: Pool = Pool_f(NUM_CPUS)
            _pool_processes(
                mp_pool,
                pool_hap_props[1:],
                total_read_count,
                configs_dpath,
                ref_fpath,
                output_dir,
                start,
                prefix,
            )
            mp_pool.close()
            mp_pool.join()


if __name__ == "__main__":
    # CLI input parameters.
    sim_prop_fpath: Path = Path(sys.argv[1])
    regions_fpath: Path = Path(sys.argv[2])
    chrom_fpaths_fpath: Path = Path(sys.argv[3])
    total_read_count: int = int(sys.argv[4])
    configs_dpath: Path = Path(sys.argv[5])
    remsim_fpath: Path = Path(sys.argv[6])
    outputs_dpath: Path = Path(sys.argv[7])
    NUM_CPUS: int = int(sys.argv[8])

    main(
        sim_prop_fpath,
        regions_fpath,
        chrom_fpaths_fpath,
        total_read_count,
        configs_dpath,
        remsim_fpath,
        outputs_dpath,
        NUM_CPUS,
    )
