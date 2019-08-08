import argparse
import isopy.isopy as iso


def parse_arguments():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument(
        "-r",
        "--read_files_fas",
        type=list,
        help="List of path to the fasta files with reads",
    )
    parser.add_argument(
        "-g", "--genome_fas", help="Path to the fasta genome file to map to"
    )
    parser.add_argument(
        "-o", "--out_dir", help="Path to the output directory (which will be created)"
    )
    parser.add_argument(
        "-m",
        "--min_read_count",
        help="Minimum number of long reads which supports an exon to consider this exon",
    )
    parser.add_argument("-t", "--threads", help="Number of threads to use")

    return parser.parse_args()


def main():

    args = parse_arguments()

    isoseq = iso.IsoAnalysis(
        args.read_files_fas,
        args.genome_fas,
        args.out_dir,
        args.min_read_count,
        args.threads,
    )


if __name__ == "__main__":
    main()
