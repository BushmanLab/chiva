import os
import sys
import argparse
import subprocess

from ruamel.yaml import YAML
from pathlib import Path

def main( argv = sys.argv ):
    """Initiate an cHIVa project run in Snakemake."""

    try:
        conda_prefix = os.environ.get("CONDA_PREFIX")
    except (KeyError, IndexError):
        raise SystemExit(
            "Could not determine Conda prefix. Activate your cHIVa "
            "environment and try this command again.")

    usage_str = "chiva %(prog)s <path/to/config.file> <options> -- <snakemake.options>"

    description_str = (
        "Initiate the processing of an cHIVa project givin a configuration "
        "file. Arguments after '--' are passed to Snakemake asis.")
    
    parser = argparse.ArgumentParser(
        prog = "run", 
        usage = usage_str,
        description = description_str
    )

    parser.add_argument(
        "config", 
        help = ("name of config file (%(default)s)"),
        metavar = "CONFIG_FILE"
    )

    parser.add_argument(
        "-i", "--chiva_dir", 
        default = os.getenv("CHIVA_DIR", os.getcwd()),
        help = "Path to cHIVa installation."
    )

    # The remaining args (after --) are passed to Snakemake
    args, remaining = parser.parse_known_args(argv)

    snakefile = Path(args.chiva_dir)/"Snakefile"
    if not snakefile.exists():
        sys.stderr.write(
            "Error: could not find a Snakefile in directory '{}'\n".format(
                args.chiva_dir))
        sys.exit(1)

    snakemake_args = ['snakemake',
                      '--snakefile', str(snakefile),
                      '--configfile', str(args.config),
                      '--dir', str(args.chiva_dir)] + remaining
    print("Running: "+" ".join(snakemake_args))

    cmd = subprocess.run(snakemake_args)
    
    sys.exit(cmd.returncode)
    
        
def check_existing(path, force = False):
    if path.is_dir():
        raise SystemExit(
            "Error: specified file '{}' exists and is a directory".format(path))
    if path.is_file() and not force:
        raise SystemExit(
            "Error: specified file '{}' exists. Use --force to "
            "overwrite.".format(path))
    return path
