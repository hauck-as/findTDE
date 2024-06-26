# -*- coding: utf-8 -*-
# Copyright (c) 2024 Alexander Hauck. Distributed under the terms of the MIT License.
import sys
import argparse
import warnings
from pathlib import Path

from findtde import __version__
from findtde.cli.main_functions import *
import pydefect.cli.main as pyd_parse


description = """findTDE comprises a set of scripts to facilitate easy,
high-throughput calculations of threshold displacement energies for materials
in VASP/LAMMPS using ab initio/classical molecular dynamics."""


epilog = f'Author: Alexander Hauck Version: {__version__}'


def parse_args_main(args):
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    subparsers = parser.add_subparsers()

    return parser.parse_args(args)


def main():
    args = parse_args_main(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()