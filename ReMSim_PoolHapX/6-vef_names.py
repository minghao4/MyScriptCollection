from pathlib import Path
import sys


def main(vef_dpath: Path, project_name: str, output_dpath: Path) -> None:
    """
    Create the sample names TXT file for a given project; to be used in PoolHapX.

    Parameters
    ----------
    vef_dpath : pathlib.Path
        Directory holding VEF files for each pool.

    project_name : str
        Project name.

    output_dpath : pathlib.Path
        Directory to write the output file to.
    """

    output_file: Path = Path(output_dpath, "{}_sample_names.txt".format(project_name))
    with output_file.open("w") as out_f:
        vef_file: Path
        for vef_file in vef_dpath.glob("*.vef"):
            out_f.write(vef_file.name[:-4] + "\n")


if __name__ == "__main__":
    # CLI input parameters.
    vef_dpath: Path = Path(sys.argv[1])
    project_name: str = sys.argv[2]
    output_dpath: Path = Path(sys.argv[3])

    main(vef_dpath, project_name, output_dpath)
