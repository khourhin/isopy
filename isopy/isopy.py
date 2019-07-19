"""A library for alternative splicing analysis using long read technology """

from pynextgen.basics_fasta import Fasta
from pybedtools import BedTool
import pynextgen.bed as bed
import pandas as pd
import numpy as np
import subprocess
import os

# Number of basis in 5' and 3' of a junction to make it specific in our dataset
OVERHANG_5 = 9
OVERHANG_3 = 6
OUT_DIR = "isopy_out"
GENOME = "mm10"


def exec_command(cmd, silent=False):
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
    def __init__(self, reads_fas, genome_fas, out_dir=OUT_DIR, min_read_count=0):
        """ A tool to identify exons from BAM files

        :param reads_fas: Path to fasta files with long reads.
        :param genome_fas: Path to fasta file of the genome to map to.
        :param out_dir: Path to the output directory.
        :param min_read_count: The minimum number of read supporting an exon for the the exon to be considered.
        """
        self.genome = genome_fas
        self.reads = reads_fas
        self.out_dir = out_dir
        os.makedirs(self.out_dir)

        self.bam = self._map_reads_to_genome()
        self.exons_df = self._extract_raw_exons()
        self._filter_exons_by_canonical_splice_sites()
        self._filter_exons_by_read_support(min_read_count)

        self.export_to_bed()

    def _map_reads_to_genome(self):
        """ Map long reads to the genome. Using Minimap2.
        """

        bam = os.path.basename(os.path.splitext(self.reads)[0] + ".bam")
        bam = os.path.join(self.out_dir, bam)

        minimap_cmd = f"minimap2 -ax splice -uf {self.genome} {self.reads} \
        | samtools view -bh \
        | samtools sort \
        > {bam}"

        exec_command(minimap_cmd)

        return bam

    def _extract_raw_exons(self):
        """ From the bam alignment, get a dataframe with the regions of putative exons
        """

        def get_splice_sites(row):
            """ Function to apply to an exon dataframe to obtain splice sites
            """

            # Get 2 bases upstream and downstream of the exon
            sequence = BedTool.seq(
                (row["chrom"], row["start"] - 2, row["end"] + 2), self.genome
            )
            donor = sequence[0:2]
            acceptor = sequence[-2:]

            return (donor, acceptor)

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

        # Add splice sites information
        exons_df.loc[:, "splice_sites"] = exons_df.apply(get_splice_sites, axis=1)

        return exons_df

    def _filter_exons_by_canonical_splice_sites(self):
        """ Keep only exons with canonical splices sites
        (ie GT-AG, AT-AC, GT-AC, AT-AG)
        """

        canonical_junctions = (("AG", "GT"), ("AC", "AT"), ("AC", "GT"), ("AG", "AT"))
        self.exons_df = self.exons_df.loc[
            self.exons_df.splice_sites.isin(canonical_junctions)
        ]

    def _filter_exons_by_read_support(self, n):
        """ Filter only exons supported by at least n reads.
        """

        self.exons_df = self.exons_df.loc[self.exons_df.frequency >= n]

    def export_to_bed(self):
        """ Convert the dataframe self.exons_df to a bed file
        """

        bed_df = self.exons_df.copy()
        bed_df.drop(["splice_sites"], axis=1, inplace=True)
        bed_df.insert(3, "name", [f"E{n}" for n in range(1, len(bed_df) + 1)])

        BedTool.from_dataframe(bed_df).saveas(os.path.join(self.out_dir, "test.bed"))


class Transcripts(object):
    """Representation of a set of transcripts. Transcripts are herein first defined by
    their exon_composition
    """

    def __init__(
        self,
        exon_composition_csv,
        exon_fasta,
        fasta_out="artificial_transcripts.fas",
        exon_bed_out="artificial_transcripts_exons.bed",
        junction_bed_out="artificial_transcripts_junctions.bed",
    ):
        """
        :param exon_composition_matrix: An dataframe in Excel format with as index
        the transcripts IDs and as columns the exon names. IMPORTANT: The exons
        in columns should be ordered in order of appearance in the trancript !!
        
        :param exon_fas_ta: A fasta file with the sequence of each exons with IDs
        corresponding to the IDs in exon_composition_matrix

        """
        self.exon_composition_df = pd.read_csv(exon_composition_csv, index_col=0)
        self.exon_bed, self.junction_bed, self.fasta = self._build_fasta_bed_from_exon_composition(
            exon_fasta, fasta_out, exon_bed_out, junction_bed_out
        )
        self.exon_bed.to_saf()
        self.junction_bed.to_saf()

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
