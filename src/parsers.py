'''Module for parsing simple fasta & simple fastq files'''
from __future__ import annotations
import argparse
from typing import TextIO


def parse_fasta(fasta_file: TextIO) -> dict[str, str]:
    seqs = {}
    name = ""
    tmp: list[str] = []
    for line in fasta_file:
        if line and line[0] == '>':
            if name:
                seq = ''.join(tmp)
                seqs[name] = seq
            name = line[1::].strip()
            tmp: list[str] = []
        else:
            clean = ''.join(line.split())
            tmp.append(clean.strip())
    if name:
        seq = ''.join(tmp)
        seqs[name] = seq
    return seqs


def parse_fastq(fastq_file: TextIO) -> dict[str, str]:
    out: dict[str, str] = {}
    tmp = ''.join(fastq_file.readlines())
    seqs = tmp.split('@')
    for line in seqs:
        if line:
            seq = line.split('\n')
            out[seq[0]] = seq[1]
    return out


def main():
    argparser = argparse.ArgumentParser(
        description="Extract Simple-FASTA records"
    )
    argparser.add_argument(
        "fasta",
        type=argparse.FileType('r')
    )
    args = argparser.parse_args()

    lib = parse_fastq(args.fasta)
    for key in lib:
        print(f'{key}\t{lib[key]}')

    args.fasta.close()


if __name__ == '__main__':
    main()
