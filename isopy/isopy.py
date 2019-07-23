"""A library for alternative splicing analysis using long read technology """

from pynextgen.basics_fasta import Fasta
from pybedtools import BedTool
import pynextgen.bed as bed
import pandas as pd
import numpy as np
import subprocess
import os
import logging
from natsort import natsorted

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

# Number of basis in 5' and 3' of a junction to make it specific in our dataset
OVERHANG_5 = 9
OVERHANG_3 = 6
OUT_DIR = "isopy_out"


def _exec_command(cmd, silent=False):
    """ Wrapper for proper execution of subprocesses
    """

    try:
        proc = subprocess.run(
            cmd,
            shell=True,
            check=True,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            universal_newlines=True,
        )
        if not silent and proc.stdout:
            print(proc.stdout.strip())
            return proc.stdout

    except subprocess.CalledProcessError as exc:
        raise ChildProcessError(exc.stderr)


class ExonIdentifier(object):
    """ A tool to identify exons from BAM files

        :param reads_fas: Path to the fasta file with long reads.
        :param genome_fas: Path to the fasta file of the genome to map to.
        :param out_dir: Path to the output directory.
        :param min_read_count: The minimum number of read supporting an exon for the the exon to be considered.
        """

    def __init__(self, reads_fas, genome_fas, out_dir=OUT_DIR, min_read_count=0):
        self.genome = genome_fas
        self.reads = reads_fas
        self.out_dir = out_dir
        os.makedirs(self.out_dir)

        self.bam = self._map_reads_to_genome()
        self.exons_df = self._extract_raw_exons()
        self._filter_exons_by_canonical_splice_sites()
        self._filter_exons_by_read_support(min_read_count)

        self.exons_fas_fn, self.exons_bed_fn = self.export()

    def _map_reads_to_genome(self):
        """ Map long reads to the genome. Using Minimap2.
        """

        bam = os.path.basename(os.path.splitext(self.reads)[0] + ".bam")
        bam = os.path.join(self.out_dir, bam)

        minimap_cmd = f"minimap2 -ax splice -uf {self.genome} {self.reads} \
        | samtools view -bh \
        | samtools sort \
        > {bam}"

        logging.debug("Minimap2 START")
        _exec_command(minimap_cmd)
        logging.debug("Minimap2 DONE")

        return bam

    def _extract_raw_exons(self):
        """ From the bam alignment, get a dataframe with the regions of putative exons
        """

        def get_sequence_and_splice_sites(row):
            """ Function to apply to an exon dataframe to obtain exon sequence and splice sites
            """

            # Get 2 bases upstream and downstream of the exon
            sequence = BedTool.seq(
                (row["chrom"], row["start"] - 2, row["end"] + 2), self.genome
            )
            donor = sequence[0:2]
            acceptor = sequence[-2:]

            # Series are used to cast the result of apply in 2 columns
            return pd.Series([sequence[2:-2], (donor, acceptor)])

        logging.debug("Raw exon extraction START")
        # Extract mapped intervals from the bam alignement file
        exons_df = BedTool(self.bam).bam_to_bed(split=True).sort().to_dataframe()

        # Keep only interval information
        exons_df = exons_df.loc[:, ["chrom", "start", "end"]]

        # Add number of reads found to support a particular exon
        exons_df = (
            exons_df.groupby(["chrom", "start", "end"])
            .size()
            .reset_index(name="frequency")
        )

        # Add sequence and splice sites information
        exons_df[["sequence", "splice_sites"]] = exons_df.apply(
            get_sequence_and_splice_sites, axis=1
        )

        logging.debug("Raw exon extraction DONE")
        return exons_df

    def _filter_exons_by_canonical_splice_sites(self):
        """ Keep only exons with canonical splices sites
        (ie GT-AG, AT-AC, GT-AC, AT-AG)
        """

        logging.debug("Filter exons by canonical splice sites START")

        canonical_junctions = (("AG", "GT"), ("AC", "AT"), ("AC", "GT"), ("AG", "AT"))
        self.exons_df = self.exons_df.loc[
            self.exons_df.splice_sites.isin(canonical_junctions)
        ]

        logging.debug("Filter exons by canonical splice sites DONE")

    def _filter_exons_by_read_support(self, n):
        """ Filter only exons supported by at least n reads.
        """

        logging.debug("Filter exons by read support START")

        self.exons_df = self.exons_df.loc[self.exons_df.frequency >= n]

        logging.debug("Filter exons by read support DONE")

    def export(self):
        """ Convert the dataframe self.exons_df to a bed file and fasta file
        and name the exons.
        """

        exon_fas_fn = os.path.join(self.out_dir, "test.fas")
        exon_bed_fn = os.path.join(self.out_dir, "test.bed")

        export_df = self.exons_df.copy()
        export_df.insert(3, "exon_id", [f"E{n}" for n in range(1, len(export_df) + 1)])

        with open(exon_fas_fn, "w") as fas_out:
            for index, exon in export_df.iterrows():
                fas_out.write(f">{exon.exon_id}\n{exon.sequence}\n")

        export_df.drop(["sequence", "splice_sites"], axis=1, inplace=True)
        BedTool.from_dataframe(export_df).saveas(exon_bed_fn)

        return (exon_fas_fn, exon_bed_fn)


class Transcripts(object):
    """Representation of a set of transcripts. Transcripts are herein first defined by
    their exon_composition

    :param exon_composition_matrix: An dataframe in csv format with as index the transcripts IDs and as columns the exon names. IMPORTANT: The exons in columns should be ordered in order of appearance in the trancript !!
        
    :param exon_fasta: A fasta file with the sequence of each exons with IDs corresponding to the IDs in exon_composition_matrix
    """

    def __init__(
        self,
        exon_composition_csv,
        exon_fasta,
        fasta_out="artificial_transcripts.fas",
        exon_bed_out="artificial_transcripts_exons.bed",
        junction_bed_out="artificial_transcripts_junctions.bed",
    ):
        self.exon_composition_df = pd.read_csv(exon_composition_csv, index_col=0)
        self.exon_bed, self.junction_bed, self.fasta = self._build_fasta_bed_from_exon_composition(
            exon_fasta, fasta_out, exon_bed_out, junction_bed_out
        )
        self.exon_bed.to_saf()
        self.junction_bed.to_saf()

        self.junction_composition_df = self._build_junction_composition()

    def __repr__(self):

        # Since row are transcripts_id and columns exons
        n_transcripts, n_exons = self.exon_composition_df.shape
        return f"<Collection Object of {n_transcripts} transcripts with {n_exons} potential exons>"

    def _get_exons(self, row):
        """
        Extract a list of exons from a line of the exon composition matrix
        """

        exons = list(row[row == True].index)
        return exons

    def _get_junctions(self, row):
        """
        Extract a string with all junctions (e.g 'E1_E2') separated by ','.
        """
        exons = list(row[row == True].index)
        junctions = ",".join(
            ["_".join(exons[i : i + 2]) for i in range(len(exons) - 1)]
        )
        return junctions

    def _build_junction_composition(self):
        """ Returns a DataFrame with summary information on the population of transcripts.
        """

        junctions = self.exon_composition_df.apply(self._get_junctions, axis=1)
        junctions.index.name = "transcript_id"

        # This produce the junction composition table
        junctions_df = (
            junctions.str.split(",", expand=True)
            .stack()
            .reset_index()
            .pivot_table(index="transcript_id", columns=0, aggfunc="size")
        ).fillna(0)

        index_order = natsorted(list(junctions_df.index))
        column_order = natsorted(list(junctions_df.columns))

        junctions_df = junctions_df.loc[index_order, column_order]
        junctions_df = junctions_df.astype(bool)

        return junctions_df

    def _build_fasta_bed_from_exon_composition(
        self, exon_fasta, fasta_out, exon_bed_out, junction_bed_out
    ):
        """Build the sequence exon by exon depending on the exon_composition_matrix and

        :param exon_fasta: Path to the fasta file with exon nucleotide sequence.
        :param fasta_out: Path to recomposed transcripts fasta file.
        """
        transcript_exon_list = self.exon_composition_df.apply(self._get_exons, axis=1)
        exon_seqs = Fasta(exon_fasta).to_dict()

        seq_dict = {}

        exon_bed_df = pd.DataFrame()
        junction_bed_df = pd.DataFrame()

        for transcript_id, exon_compo in transcript_exon_list.iteritems():
            seq_dict[transcript_id] = "".join(
                [exon_seqs[exon_id] for exon_id in exon_compo]
            )

            # Define exon intervals on transcripts, using length of each exons summed
            exon_annot = pd.DataFrame(
                {
                    "chr": [transcript_id for _ in exon_compo],
                    "start": [0]
                    + list(
                        np.cumsum([len(exon_seqs[exon_id]) for exon_id in exon_compo])[
                            :-1
                        ]
                    ),
                    "end": np.cumsum(
                        [len(exon_seqs[exon_id]) for exon_id in exon_compo]
                    ),
                    "name": [exon_id for exon_id in exon_compo],
                    "score": "0",
                    "strand": ".",
                }
            )

            junction_annot = exon_annot.copy()
            junction_annot.end = junction_annot.start + OVERHANG_3
            junction_annot.start = junction_annot.start - OVERHANG_5
            junction_names = [
                f"junction_{exon_1}_{exon_2}"
                for exon_1, exon_2 in zip(
                    junction_annot.name[:-1], junction_annot.name[1:]
                )
            ]
            junction_annot = junction_annot[1:]
            junction_annot.name = junction_names

            # Build iteratively the bed formatted annotation table
            exon_bed_df = pd.concat([exon_bed_df, exon_annot], axis=0)
            junction_bed_df = pd.concat([junction_bed_df, junction_annot], axis=0)

        exon_bed = bed.df_to_bed(exon_bed_df, exon_bed_out)
        junction_bed = bed.df_to_bed(junction_bed_df, junction_bed_out)

        with open(fasta_out, "w") as f:
            for transcript in seq_dict:
                f.write(f">{transcript}\n{seq_dict[transcript]}\n")

        return exon_bed, junction_bed, fasta_out
