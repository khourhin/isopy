# -*- coding: utf-8 -*-

"""Console script for isopy."""
import sys
import click


@click.command()
def main(args=None):
    """Console script for isopy."""
    click.echo("Replace this message by putting your code into "
               "isopy.cli.main")
    click.echo("See click documentation at http://click.pocoo.org/")
    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
