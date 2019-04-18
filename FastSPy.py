#!/usr/bin/env python3
from gzip import open as gopen
from numpy import uintc,zeros
from sys import stderr,stdout
from time import time
import argparse
EMPTY_FASTA = "Empty FASTA"
BAD_FASTA = "Malformed FASTA"
DIFFERENT_LENGTH = "Sequences are different lengths"

# parse user arguments
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r', '--reference', required=True, type=str, help="Reference Alignment (FASTA)")
    parser.add_argument('-e', '--estimated', required=True, type=str, help="Estimated Alignment (FASTA)")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output File")
    parser.add_argument('-v', '--verbose', action='store_true', help="Verbose")
    args = parser.parse_args()
    global VERBOSE; VERBOSE = args.verbose
    if VERBOSE:
        print("Output: %s" % args.output, file=stderr)
        print("Reference: %s" % args.reference, file=stderr)
        print("Estimated: %s" % args.estimated, file=stderr)
        stderr.flush()
    if args.output == 'stdout':
        outfile = stdout
    else:
        outfile = open(args.output,'w')
    if args.reference.lower().endswith('.gz'):
        ref_stream = gopen(args.reference,'r')
    else:
        ref_stream = open(args.reference,'r')
    if args.estimated.lower().endswith('.gz'):
        est_stream = gopen(args.estimated,'r')
    else:
        est_stream = open(args.estimated,'r')
    return ref_stream, est_stream, outfile

# read sequences from a FASTA stream
def read_seqs_FASTA(stream):
    seqs = list()
    for line in stream:
        if isinstance(line,bytes):
            l = line.decode().strip()
        else:
            l = line.strip()
        if l[0] == '>':
            seqs.append('')
        else:
            seqs[-1] += l
    return seqs

# build S matrix, where S(i,j) = (a,b) denoting that the j-th nucleotide of the i-th sequence (no gaps) is in site a of reference and site b of estimated
def build_S(ref_aln, est_aln):
    num_rows = len(ref_aln); num_cols = max(len(seq)-seq.count('-') for seq in ref_aln)
    Sa = zeros((num_rows,num_cols), dtype=uintc); Sb = zeros((num_rows,num_cols), dtype=uintc)
    for S,aln in ((Sa,ref_aln), (Sb,est_aln)):
        for i in range(len(aln)):
            S_row = S[i]; aln_row = aln[i]; j = 0
            for a in range(len(aln_row)):
                if aln_row[a] != '-':
                    S_row[j] = a+1; j += 1 # store values using 1-based indexing so index = 0 implies sequence was shorter than that
    return Sa,Sb

# Compute the SP-FN error
def SPFN(ref_aln,est_aln):
    if VERBOSE:
        print("Computing S matrix...", file=stderr)
        stderr.flush()
    Sa,Sb = build_S(ref_aln,est_aln)
    A = [int(j*(j-1)/2) for j in range(1,len(ref_aln)+1)]
    ref_num_hom_pairs = sum(A[sum(row[j] != '-' for row in ref_aln)-1] for j in range(len(ref_aln[0])))
    List = [list() for _ in range(len(ref_aln[0]))]
    for i in range(len(ref_aln)):
        Sb_row = Sb[i]; ref_row = ref_aln[i]; j = 0
        for a in range(len(ref_row)):
            if ref_row[a] != '-':
                b = Sb_row[j]; List[a].append(b); j += 1
    for x in range(len(List)):
        m = zeros(len(est_aln[0]), dtype=uintc)
        for y in List[x]:
            m[y-1] += 1
    Nsum = 0; seen_y = set()
    for x in range(len(List)):
        for y in List[x]:
            if y not in seen_y:
                Nsum += A[m[y-1]]; seen_y.add(y)
    return ref_num_hom_pairs-Nsum

# run FastSPy
if __name__ == "__main__":
    ref_stream, est_stream, outfile = parse_args()
    if VERBOSE:
        START_TIME = time()
        print("Reading reference alignment...", file=stderr)
        stderr.flush()
    ref_aln = read_seqs_FASTA(ref_stream); ref_stream.close()
    if VERBOSE:
        print("Successfully read %d sequences of length %d from reference alignment" % (len(ref_aln), len(ref_aln[0])), file=stderr)
        print("Reading estimated alignment...", file=stderr)
        stderr.flush()
    est_aln = read_seqs_FASTA(est_stream); est_stream.close()
    if VERBOSE:
        print("Successfully read %d sequences of length %d from estimated alignment" % (len(est_aln), len(est_aln[0])), file=stderr)
        stderr.flush()
    error = SPFN(ref_aln,est_aln)
    print(error)
    if VERBOSE:
        END_TIME = time(); print("Total time: %f seconds" % (END_TIME-START_TIME), file=stderr)
