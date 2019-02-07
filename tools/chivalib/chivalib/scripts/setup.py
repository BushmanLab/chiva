import os
import sys
import argparse
import subprocess

from ruamel.yaml import YAML
from pathlib import Path

def main( argv = sys.argv ):
    """Create a new project directory with necessary subdirectories."""

    try:
        conda_prefix = os.environ.get("CONDA_PREFIX")
    except (KeyError, IndexError):
        raise SystemExit(
            "Could not determine Conda prefix. Activate your cHIVa "
            "environment and try this command again.")

    usage_str = "chiva %(prog)s <path/to/config.file>"
    
    description_str = (
        "Setup a new cHIVa project given a project configuration file."
    )
    
    parser = argparse.ArgumentParser(
        prog = "setup", 
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
        help = "Path to cHIVa installation")

    # Set arguments
    args, remaining = parser.parse_known_args(argv)

    bash_script = Path(args.chiva_dir) / "tools/bashscripts/make_proc_dir.sh"
    
    if not bash_script.exists():
        sys.stderr.write(
            "Error: could not find make_proc_dir.sh in directory:'{}'\n".format(
                args.chiva_dir + "/tools/bashscripts/"
            )
        )
        sys.exit(1)
        
    
    config_file = Path(args.config)
    
    if not config_file.exists():
        abs_config_file = Path(args.chiva_dir) / args.config
        if not abs_config_file.exists():
            sys.stderr.write(
                "Error: could not find config file:'{}'\n".format(
                    args.config
                )
            )
            sys.exit(1)
        else:
            config_file=abs_config_file
        
    cmd_comps = ["bash", str(bash_script), str(config_file)]
    cmd = subprocess.run(cmd_comps)
  
    sys.exit(cmd.returncode)
