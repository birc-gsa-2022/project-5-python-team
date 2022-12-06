from __future__ import annotations
from parsers import parse_fasta, parse_fastq
from tree import SuffixTree
from cigar import edits_to_cigar
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
        processed_genome = {chr: SuffixTree(genome[chr]) for chr in genome}
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

        out=[]
        for i, chr in enumerate(genome):
            for j, read in enumerate(reads):
                gg=genome[chr].search_approx(reads[read],args.d)
                if gg:
                    for edit in gg:
                        for position in edit[0]:
                            out.append(
                                f'{read}\t{chr}\t{int(position)+1}\t{edits_to_cigar(edit[1])}\t{reads[read]}'
                            )
        out.sort()
        print('\n'.join(out))


if __name__ == '__main__':
    main()
