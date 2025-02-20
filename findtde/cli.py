"""findTDE package containing command line interface (CLI) functions."""

from pathlib import Path
import click

from findtde import io


@click.group("findtde", no_args_is_help=True)
def findtde():
    """findTDE CLI group."""
    pass


@findtde.command()
def run():
    pass


@findtde.command()
def multi():
    pass