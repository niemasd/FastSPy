#!/usr/bin/env python3

# read FASTA file
def read_FASTA(stream):
    seqs = {}; name = None; seq = ''
    for line in stream:
        if isinstance(line,bytes):
            l = line.decode().strip()
        else:
            l = line.strip()
        if len(l) == 0:
            continue
        if l[0] == '>':
            if name is not None:
                if len(seq) == 0:
                    raise RuntimeError("Malformed FASTA")
                seqs[name] = seq
            name = l[1:]
            if name in seqs:
                raise RuntimeError("Duplicate sequence ID: %s" % name)
            seq = ''
        else:
            seq += l
    if name is None or len(seq) == 0:
        raise RuntimeError("Malformed FASTA")
    seqs[name] = seq
    return seqs

# run FastSPy
if __name__ == "__main__":
    import argparse; from gzip import open as gopen; from sys import stderr
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r', '--reference', required=True, type=str, help="Reference Alignment (FASTA)")
    parser.add_argument('-e', '--estimated', required=True, type=str, help="Estimated Alignment (FASTA)")
    parser.add_argument('-v', '--verbose', action='store_true', help="Verbose")
    args = parser.parse_args()
    global VERBOSE; VERBOSE = args.verbose
    if VERBOSE:
        print("Reading reference alignment from: %s" % args.reference, file=stderr); stderr.flush()
    if args.reference.lower().endswith('.gz'):
        ref_stream = gopen(args.reference)
    else:
        ref_stream = open(args.reference)
    ref = read_FASTA(ref_stream); ref_stream.close()
    if VERBOSE:
        print("Successfully read %d sequences from reference alignment" % len(ref), file=stderr)
        print("Reading estimated alignment from: %s" % args.estimated, file=stderr); stderr.flush()
    if args.estimated.lower().endswith('.gz'):
        est_stream = gopen(args.estimated)
    else:
        est_stream = open(args.estimated)
    est = read_FASTA(est_stream); est_stream.close()
    if VERBOSE:
        print("Successfully read %d sequences from estimated alignment" % len(est), file=stderr); stderr.flush()
