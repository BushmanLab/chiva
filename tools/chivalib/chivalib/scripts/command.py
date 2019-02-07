import sys
import argparse
import subprocess

from chivalib import __version__
from chivalib.scripts.run import main as Run
from chivalib.scripts.setup import main as Setup
#from chivalib.scripts.config import main as Config
from chivalib.scripts.list_samples import main as ListSamples
from chivalib.scripts.report import main as Report
from chivalib.scripts.clean import main as Clean

def main():

    usage_str = "%(prog)s [-h/--help,-v/--version] <subcommand> <path/to/config.file> <options> -- <snakemake.options>"
    description_str = (
        "subcommands:\n"
        "  setup        \tCreate a new config file for a project using local data.\n"
        "  run          \tExecute the cHIVa pipeline.\n"
        "  report       \tGenerate a custom report from cHIVa output files.\n"
        "  list_samples \tOutput a list of samples from a project.\n"
        "  config       \t[inDev] Modify or update cHIVa config files.\n"
        "  clean        \tCleanup project directory to reduce size. Keeps terminal files."
    ).format(version=__version__)

    parser = argparse.ArgumentParser(
        prog = "chiva",
        usage = usage_str,
        description = description_str,
        epilog = "For more help, see the docs at http://chiva.readthedocs.io.",
        formatter_class = argparse.RawDescriptionHelpFormatter,
        add_help = False
    )

    parser.add_argument(
        "command", default = "None", help = argparse.SUPPRESS, nargs = "?"
    )
    
    parser.add_argument(
        "-v", "--version", action = "version",
        version = "%(prog)s {}".format(__version__)
    )

    args, remaining = parser.parse_known_args()
    
    sub_cmds = ["setup", "run", "config", "list_samples", "report", "clean"]
    
    if not args.command in sub_cmds:
        parser.print_help()
        if not args.command in ['None']:
            sys.stderr.write("  Unrecognized subcommand, '{}'.\n".format(
                args.command
            ))
        sys.exit(1)

    if args.command == "setup":
        Setup(remaining)
    elif args.command == "run":
        Run(remaining)
    elif args.command == "config":
        raise SystemExit(
          print("  'chiva config' subcommand is currently under development.\n"
                "  Checkout https://github.com/cnobles/cHIVa/ for updates   \n"
                "  and announcements. Thanks for using cHIVa!               \n"
          )
        )
        #Config(remaining)
    elif args.command == "list_samples":
        ListSamples(remaining)
    elif args.command == "report":
        Report(remaining)
    elif args.command == "clean":
        Clean(remaining)
    else:
        parser.print_help()
