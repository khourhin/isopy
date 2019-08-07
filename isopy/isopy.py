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
from pathlib import Path

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

# Number of basis in 5' and 3' of a junction to make it specific in our dataset
OVERHANG_5 = 9
OVERHANG_3 = 6
OUT_DIR = Path("isopy_out")


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


def build_fasta_bed_from_exon_composition(
    exon_fas, cluster_df, fasta_out, exon_bed_out, junction_bed_out
):
    """ Build fasta and bed files with pseudo-transcripts. These pseudo-transcript are obtained by a concatenation of the exon sequences (found in exon_fas) following an exon composition matrix (cluster_df)
    """

    def _get_exons(row):
        """
        Extract a list of exons from a line of the exon composition matrix
        """

        exons = list(row[row == True].index)
        return exons

    transcript_exon_list = cluster_df.drop(["read_ids", "frequency"], axis=1).apply(
        _get_exons, axis=1
    )
    exon_seqs = Fasta(exon_fas).to_dict()

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
                    np.cumsum([len(exon_seqs[exon_id]) for exon_id in exon_compo])[:-1]
                ),
                "end": np.cumsum([len(exon_seqs[exon_id]) for exon_id in exon_compo]),
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
            for exon_1, exon_2 in zip(junction_annot.name[:-1], junction_annot.name[1:])
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


class IsoAnalysis(object):
    """ The complete isopy alternative splicing analysis
    """

    def __init__(
        self, read_files_fas, genome_fas, out_dir=OUT_DIR, min_read_count=0, threads=1
    ):

        self.out_dir = Path(out_dir)

        self.exons = ExonIdentifier(
            read_files_fas,
            genome_fas,
            out_dir=self.out_dir,
            min_read_count=min_read_count,
        )
        self.clusters = TranscriptCluster(
            self.exons.exon_fas_fn, read_files_fas, out_dir=self.out_dir
        )
        build_fasta_bed_from_exon_composition(
            self.exons.exon_fas_fn,
            self.clusters.cluster_df,
            "test.fas",
            "exons.bed",
            "junctions.bed",
        )

        self.exons_df = self.exons.exons_df
        self.clusters_df = self.clusters.cluster_df
        self.exons_composition_df = self.clusters.exon_composition_df
        self.junction_pseudotrans = self.clusters.junction_pseudtrans
        self.junction_libraries = self.clusters.junction_libraries

        self.blast_df = self.clusters.blast_df


class ExonIdentifier(object):
    """ A tool to identify exons from long reads and a reference genome

        :param read_files_fas: List of paths to the fasta files with long reads.
        :param genome_fas: Path to the fasta file of the genome to map to.
        :param out_dir: Path to the output directory.
        :param min_read_count: The minimum number of read supporting an exon for the the exon to be considered.
        """

    def __init__(self, read_files_fas, genome_fas, out_dir, min_read_count=0):
        self.genome = genome_fas
        self.read_files_fas = read_files_fas
        self.out_dir = Path(out_dir)
        os.makedirs(self.out_dir)

        self.merged_bam = self.out_dir / "merged.bam"
        self._map_reads_to_genome()
        self.exons_df = self._extract_raw_exons()
        self._filter_exons_by_canonical_splice_sites()
        self._filter_exons_by_read_support(min_read_count)

        self.exon_fas_fn, self.exon_bed_fn = self.export()

    def _map_reads_to_genome(self):
        """ Map long reads to the genome. Using Minimap2.
        """
        bams = []

        for read_file in self.read_files_fas:

            bam = os.path.basename(os.path.splitext(read_file)[0] + ".bam")
            bam = self.out_dir / bam

            minimap_cmd = f"minimap2 -ax splice -uf {self.genome} {read_file} \
            | samtools view -bh \
            | samtools sort \
            > {bam}"

            logging.debug("Minimap2 START")
            _exec_command(minimap_cmd)
            logging.debug("Minimap2 DONE")

            bams.append(str(bam))

        samtools_merge_cmd = f"samtools merge {self.merged_bam} {' '.join(bams)}"

        logging.debug("Samtools merge START")
        _exec_command(samtools_merge_cmd)
        logging.debug("Samtools merge DONE")

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
        exons_df = (
            BedTool(str(self.merged_bam)).bam_to_bed(split=True).sort().to_dataframe()
        )

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


class TranscriptCluster(object):
    """ Cluster long reads based on their exon composition. 
    
    :param exon_fas: Path to the fasta file with the exon sequences
    :param read_files_fas: List of path to the fasta files with reads
    :param out_dir: Path to the output directory
    :param threads: Number of threads to use
    """

    def __init__(
        self, exon_fas, read_files_fas, out_dir=OUT_DIR, prefix=None, threads=1
    ):

        self.exon_fas = exon_fas
        self.read_files_fas = read_files_fas
        self.out_dir = out_dir

        self.threads = threads

        self._format_blast_db()
        self.blast_df = self._run_blast_versus_exons()
        self.exon_composition_df = self._get_exon_composition_from_blast()
        self.cluster_df = self._cluster_reads()
        self.junction_pseudtrans = self._get_junction_composition(all_reads=False)
        self.junction_libraries = self._get_junction_composition(all_reads=True)
        # self.export()

    def __repr__(self):

        return "<TranscriptCluster object for Iso-seq analyis>"

    def _format_blast_db(self):
        """ Generate the blast db for exons using makeblastdb
        """

        blastdb_cmd = f"makeblastdb -dbtype nucl -in {self.exon_fas}"
        _exec_command(blastdb_cmd)

    def _run_blast_versus_exons(self):
        """ Execute blast of transcripts against exons identified by ExonIdentifier
        """

        blast_dfs = []
        blast_out_fn = "blast_out.tab"

        for fasta in self.read_files_fas:

            blast_cmd = f"blastn -num_threads {self.threads} -outfmt 6 -task 'blastn-short' \
            -evalue 0.001 -penalty -2 -query {fasta} -db {self.exon_fas} > {(self.out_dir / blast_out_fn)}"

            _exec_command(blast_cmd)

            fields = [
                "query",
                "subject",
                "identity",
                "alignment_length",
                "mismatches",
                "gap_opens",
                "query_start",
                "query_end",
                "subject_start",
                "subject_end",
                "evalue",
                "bit_score",
            ]

            blast_df = pd.read_csv(
                (self.out_dir / blast_out_fn), sep="\t", names=fields
            )
            blast_df.loc[:, "library"] = os.path.basename(fasta)

            blast_dfs.append(blast_df)

        blast_df = pd.concat(blast_dfs)
        return blast_df

    def _get_exon_composition_from_blast(self, id_threshold=90):
        """ Filter blast results by identity percent threshold
        """

        exon_composition_df = (
            pd.pivot_table(
                self.blast_df,
                values="identity",
                columns="subject",
                index=["library", "query"],
            )
            > id_threshold
        )

        column_order = natsorted(exon_composition_df.columns)
        exon_composition_df = exon_composition_df.loc[:, column_order]

        return exon_composition_df

    def _cluster_reads(self, prefix="alternative_transcript_"):
        """ Cluster reads based on their exon composition

        :param prefix: The prefix to put for transcript cluster names
        """

        # Group by exon composition and create tuple of reads id corresponding to an exon composition
        cluster_df = self.exon_composition_df.groupby(
            self.exon_composition_df.columns.to_list()
        ).apply(lambda x: tuple(x.index))
        cluster_df = cluster_df.reset_index(name="read_ids")

        cluster_df["frequency"] = cluster_df["read_ids"].apply(len)
        cluster_df.sort_values("frequency", ascending=False, inplace=True)

        cluster_df["transcript_id"] = [
            f"{prefix}{n}" for n in range(1, cluster_df.index.size + 1)
        ]
        cluster_df.set_index("transcript_id", inplace=True)

        return cluster_df

    def _get_junctions(self, row):
        """
        Extract a string with all junctions (e.g 'E1_E2') separated by ','.
        """
        exons = list(row[row == True].index)
        junctions = ",".join(
            ["_".join(exons[i : i + 2]) for i in range(len(exons) - 1)]
        )
        return junctions

    def _get_junction_composition(self, all_reads=True):
        """ Returns a DataFrame with summary information on the population of transcripts.
        """

        if all_reads:
            junctions = self.exon_composition_df.apply(self._get_junctions, axis=1)
            # This produce the junction composition table
            junctions_df = (
                junctions.str.split(",", expand=True)
                .stack()
                .reset_index()
                .pivot_table(index=["library", "query"], columns=0, aggfunc="size")
            ).fillna(0)

        else:
            junctions = self.cluster_df.drop(["read_ids", "frequency"], axis=1).apply(
                self._get_junctions, axis=1
            )
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

    def export(self):
        """ Export the results of the clustering to self.out_dir:
        - Exon composition csv
        - Reads clustered in various fasta file corresponding to their cluster name (transcript/isoform)
        """

        self.exon_composition_df.to_csv(self.out_dir / f"exon_composition.csv")

        # for k, reads in self.cluster_df["read_ids"].iteritems():
        #     with open(self.out_dir / f"{k}.fas", "w") as f:
        #         for seq in Fasta(self.read_files_fas).parse():
        #             if seq.name_simple in reads:
        #                 f.write(str(seq))
