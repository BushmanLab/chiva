import os
import sys
import argparse
import subprocess

from ruamel.yaml import YAML
from pathlib import Path

def main( argv = sys.argv ):
    """Generate a custom report from cHIVa output files.."""

    try:
        conda_prefix = os.environ.get("CONDA_PREFIX")
    except (KeyError, IndexError):
        raise SystemExit(
            "Could not determine Conda prefix. Activate your cHIVa "
            "environment and try this command again.")

    root_dir = os.getenv("CHIVA_DIR")
    r_script = Path(root_dir + "/tools/rscripts/generate_cHIVa_report.R")
    
    if not r_script.is_file():
        sys.stderr.write(
            "Error: Could not find a {0} in directory '{1}'\n".format(
                "generate_cHIVa_report.R", root_dir + "/tools/rscripts/"
            )
        )
        sys.exit(1)
    
    r_comps = ["Rscript", str(r_script)] + argv

    cmd = subprocess.run(r_comps)

    sys.exit(cmd.returncode)
