#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `isopy` package."""

from easydev import md5
from click.testing import CliRunner
import pandas as pd

import isopy.isopy as iso
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

    exon_bed, junction_bed, fasta_out = iso.build_fasta_bed_from_exon_composition(
        (datadir / "mock_exons.fas"),
        pd.read_csv(datadir / "mock_exon_composition.csv", index_col=0),
        fasta_out=(tmpdir / "mock_transcripts.fas"),
        exon_bed_out=(tmpdir / "mock_exons.bed"),
        junction_bed_out=(tmpdir / "mock_junctions.bed"),
    )

    assert md5(fasta_out) == md5(datadir / "ref_mock_transcripts.fas")
    assert md5(exon_bed.path) == md5(datadir / "ref_mock_exons.bed")
    assert md5(junction_bed.path) == md5(datadir / "ref_mock_junctions.bed")


def test_ExonIdentifier(tmpdir, datadir):
    """ Test the ExonIdentifier object
    """
    exon_bed, junction_bed, fasta_out = iso.build_fasta_bed_from_exon_composition(
        (datadir / "mock_exons.fas"),
        pd.read_csv(datadir / "mock_exon_composition.csv", index_col=0),
        fasta_out=(tmpdir / "mock_transcripts.fasta"),
        exon_bed_out=(tmpdir / "mock_exons.bed"),
        junction_bed_out=(tmpdir / "mock_junction.bed"),
    )

    exons = iso.ExonIdentifier(
        read_files_fas=[fasta_out],
        # Has to be converted to str for bedtools to use it properly for chromSizes
        genome_fas=str(datadir / "mock_genome.fas"),
        out_dir=(tmpdir / "ali"),
    )

    assert md5(datadir / "mock_exons.fas") == md5(exons.exon_fas_fn)


def test_multi_fasta(tmpdir, datadir):
    """ Test the implementation with multiple fasta files
    """

    fastas = [
        (datadir / "mock_multi_fasta_1.fas"),
        (datadir / "mock_multi_fasta_2.fas"),
    ]

    exons = iso.ExonIdentifier(
        read_files_fas=fastas,
        genome_fas=str(datadir / "mock_genome.fas"),
        out_dir=(tmpdir / "ali_multifasta"),
    )

    clusters_1 = iso.TranscriptCluster(exons.exon_fas_fn, [fastas[0]], out_dir=tmpdir)
    clusters_2 = iso.TranscriptCluster(exons.exon_fas_fn, [fastas[1]], out_dir=tmpdir)

    # Formatting necessary to get correct index and index names
    cluster_1_ref = pd.read_csv(
        datadir / "ref_multifas_fasta1_exon_composition.csv",
        index_col=["library", "query"],
    )

    cluster_1_ref.rename_axis(["subject"], axis=1, inplace=True)

    pd.testing.assert_frame_equal(clusters_1.exon_composition_df, cluster_1_ref)

    # Formatting necessary to get correct index and index names
    cluster_2_ref = pd.read_csv(
        datadir / "ref_multifas_fasta2_exon_composition.csv",
        index_col=["library", "query"],
    )
    cluster_2_ref.rename_axis(["subject"], axis=1, inplace=True)

    pd.testing.assert_frame_equal(clusters_2.exon_composition_df, cluster_2_ref)

    cluster_all = iso.TranscriptCluster(exons.exon_fas_fn, fastas, out_dir=tmpdir)

    cluster_all_ref = pd.read_csv(
        datadir / "ref_multifas_fasta_all_exon_composition.csv",
        index_col=["library", "query"],
    )

    cluster_all_ref.rename_axis(["subject"], axis=1, inplace=True)

    pd.testing.assert_frame_equal(cluster_all.exon_composition_df, cluster_all_ref)


def test_IsoAnalysis(tmpdir, datadir):
    """ Test the whole pipeline IsoAnalysis

    TOFIX So far only check the IsoAnalysis.read_transcript_key table
    """

    fastas = [
        (datadir / "mock_multi_fasta_1.fas"),
        (datadir / "mock_multi_fasta_2.fas"),
        (datadir / "mock_multi_fasta_3.fas"),
    ]

    genome_fas = str(datadir / "mock_genome.fas")

    isoseq = iso.IsoAnalysis(fastas, genome_fas, out_dir=(tmpdir / "isopy_out"))

    pd.testing.assert_frame_equal(
        isoseq.read_transcript_key,
        pd.read_csv(
            datadir / "ref_read_transcript_key.csv", index_col=["library", "read_id"]
        ),
    )
