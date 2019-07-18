#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `isopy` package."""

from easydev import md5
from click.testing import CliRunner

from isopy.isopy import Transcripts, ExonIdentifier
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

    transcripts = Transcripts(
        (datadir / "mock_exon_composition.csv"),
        (datadir / "mock_exons.fas"),
        fasta_out=(tmpdir / "mock_transcripts.fas"),
        exon_bed_out=(tmpdir / "mock_exons.bed"),
        junction_bed_out=(tmpdir / "mock_junctions.bed"),
    )

    assert md5(transcripts.fasta) == md5(datadir / "ref_mock_transcripts.fas")
    assert md5(transcripts.exon_bed.path) == md5(datadir / "ref_mock_exons.bed")
    assert md5(transcripts.junction_bed.path) == md5(datadir / "ref_mock_junctions.bed")


# def test_ExonIdentifier(tmpdir, datadir):
#     """ Test the ExonIdentifier object
#     """
#     transcripts = Transcripts(
#         (datadir / "mock_exon_composition.csv"),
#         (datadir / "mock_exons.fas"),
#         fasta_out=(tmpdir / "mock_transcripts.fasta"),
#         exon_bed_out=(tmpdir / "mock_exons.bed"),
#         junction_bed_out=(tmpdir / "mock_junction.bed"),
#     )

#     exons = ExonIdentifier(
#         reads_fas=transcripts.fasta,
#         genome_fas=(datadir / "mock_genome.fas"),
#         out_dir=(tmpdir / "ali"),
#     )

#     assert md5(datadir / "mock_exons.fas") == md5(exons.fasta)
