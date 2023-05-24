#! /usr/bin/env python3

import re


def main():
    import argparse

    ap = argparse.ArgumentParser(description="""
        Prints the frames of an xyz file as given in the index file.
        Can be used to reorder, cut down or duplicate frames.""")
    ap.add_argument('xyz_file', type=argparse.FileType('r'), help="Concatenated xyz frames")
    ap.add_argument('index_file', type=argparse.FileType('r'), help="One frame number per line (1-based)")
    args = ap.parse_args()

    lines = args.xyz_file.readlines()
    args.xyz_file.close()

    indices = [int(l) - 1 for l in args.index_file]
    args.index_file.close()

    offsets = []
    for i, line in enumerate(lines):
        m = re.match(r'^\d+$', line)
        if m:
            offsets.append(i)
    offsets.append(i+1)

    for i in indices:
        print("".join(lines[offsets[i]:offsets[i+1]]), end='')

if __name__ == '__main__':
    main()
