from os import getenv, getcwd
from subprocess import run, PIPE
from setuptools import setup, find_packages

def get_chiva_version(with_hash = False):
    chiva_version_path = getenv("CHIVA_DIR", getcwd()) + "/.version"
    chiva_version = open(chiva_version_path, "r").readlines()[0].rstrip()
    commit_hash = run(
      ["git", "rev-parse", "--short", "HEAD"], stdout=PIPE
      )
    commit_str = commit_hash.stdout.decode('utf-8').rstrip()
    if with_hash:
        return chiva_version + "+" + commit_str
    else:
        return chiva_version

setup(
    name = "chiva",
    version = get_chiva_version(),
    packages = find_packages(),
    entry_points = { 'console_scripts': [
        'chiva = chivalib.scripts.command:main'
    ] }
)
