#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `isopy` package."""

from easydev import md5
from click.testing import CliRunner
import pandas as pd

from isopy.isopy import Transcripts, ExonIdentifier, TranscriptCluster
from isopy import cli

THREADS = 4


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


def test_ExonIdentifier(tmpdir, datadir):
    """ Test the ExonIdentifier object
    """
    transcripts = Transcripts(
        (datadir / "mock_exon_composition.csv"),
        (datadir / "mock_exons.fas"),
        fasta_out=(tmpdir / "mock_transcripts.fasta"),
        exon_bed_out=(tmpdir / "mock_exons.bed"),
        junction_bed_out=(tmpdir / "mock_junction.bed"),
    )

    exons = ExonIdentifier(
        read_files_fas=[transcripts.fasta],
        # Has to be converted to str for bedtools to use it properly for chromSizes
        genome_fas=str(datadir / "mock_genome.fas"),
        out_dir=(tmpdir / "ali"),
    )

    assert md5(datadir / "mock_exons.fas") == md5(exons.exon_fas_fn)

    # Test the TranscriptCluster object

    clusters = TranscriptCluster(
        exons.exon_fas_fn, (tmpdir / "mock_transcripts.fasta"), exons.out_dir, THREADS
    )

    # Test if exon initial and infered exon composition are the same
    pd.testing.assert_frame_equal(
        clusters.exon_composition_df,
        pd.read_csv(datadir / "mock_exon_composition.csv", index_col=0),
    )


def test_multi_fasta(tmpdir, datadir):
    """ Test the implementation with multiple fasta files
    """

    fastas = [
        (datadir / "mock_multi_fasta_1.fas"),
        (datadir / "mock_multi_fasta_2.fas"),
    ]

    exons = ExonIdentifier(
        read_files_fas=fastas,
        genome_fas=str(datadir / "mock_genome.fas"),
        out_dir=(tmpdir / "ali_multifasta"),
    )

    clusters_1 = TranscriptCluster(exons.exon_fas_fn, fastas[0], out_dir=tmpdir)
    clusters_2 = TranscriptCluster(exons.exon_fas_fn, fastas[1], out_dir=tmpdir)

    pd.testing.assert_frame_equal(
        clusters_1.exon_composition_df,
        pd.read_csv(datadir / "ref_multifas_fasta1_exon_composition.csv", index_col=0),
    )

    pd.testing.assert_frame_equal(
        clusters_2.exon_composition_df,
        pd.read_csv(datadir / "ref_multifas_fasta2_exon_composition.csv", index_col=0),
    )
