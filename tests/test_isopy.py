#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `isopy` package."""

import pytest
import pandas as pd
from click.testing import CliRunner
import shutil

from isopy.isopy import Transcripts
from isopy import cli


def test_command_line_interface():
    """Test the CLI."""
    runner = CliRunner()
    result = runner.invoke(cli.main)
    assert result.exit_code == 0
    assert "isopy.cli.main" in result.output
    help_result = runner.invoke(cli.main, ["--help"])
    assert help_result.exit_code == 0
    assert "--help  Show this message and exit." in help_result.output


def test_Transcripts(tmpdir, datadir):
    """ Test the Transcripts object
    """
    # isopy.Transcripts()
    transcripts = Transcripts(
        (datadir / "mock_exon_composition.csv"),
        (datadir / "mock_exons.fas"),
        fasta_out=(tmpdir / "mock_transcripts.fasta"),
        exon_bed_out=(tmpdir / "mock_exons.bed"),
        junction_bed_out=(tmpdir / "mock_junction.bed"),
    )

    print(transcripts)
