from __future__ import annotations
from parsers import parse_fasta, parse_fastq
from fm import fm_index
import argparse
import sys
import pickle
import os


def main():
    argparser = argparse.ArgumentParser(
        description="Readmapper",
        usage="\n\treadmap -p genome\n\treadmap -d dist genome reads"
    )
    argparser.add_argument(
        "-p", action="store_true",
        help="preprocess the genome."
    )
    argparser.add_argument(
        "-d", type=int, metavar="integer",
        default=1, help="max edit distance."
    )
    argparser.add_argument(
        "genome",
        help="Simple-FASTA file containing the genome.",
        type=argparse.FileType('r')
    )
    argparser.add_argument(
        "reads", nargs="?",
        help="Simple-FASTQ file containing the reads.",
        type=argparse.FileType('r')
    )
    args = argparser.parse_args()

    if args.p:
        print(f"Preprocess {args.genome}")
        genome = parse_fasta(args.genome)
        processed_genome = {chr: fm_index(genome[chr]) for chr in genome}
        with open(f"{args.genome.name}.bin", "wb") as file:
            pickle.dump(processed_genome, file)

    else:
        # here we need the optional argument reads
        if args.reads is None:
            argparser.print_help()
            sys.exit(1)
        try:
            with open(f"{args.genome.name}.bin", "rb") as file:
                genome = pickle.load(file)
        except FileNotFoundError:
            print(
                f"{args.genome.name} has not been preprocessed, you should probably do that... see how here:")
            argparser.print_help()
            sys.exit(1)

        reads = parse_fastq(args.reads)

        for chr in genome:
            for read in reads:
                print(genome[chr].approx_match(reads[read], args.d))


if __name__ == '__main__':
    main()
