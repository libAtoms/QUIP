#! /usr/bin/env python3

import argparse
from pathlib import Path
from xml.etree import ElementTree


def main():
    ap = argparse.ArgumentParser(description="Convert GAP sparseX XML output to input for use with sparseX_method=FILE.")
    ap.add_argument('file', type=argparse.FileType('r'), help="GAP main XML output file")
    args = ap.parse_args()

    text = args.file.read()
    args.file.close()

    root = ElementTree.fromstring(text)
    gp_coords = root.findall('GAP_params/gpSparse/gpCoordinates')
    for gp_coord in gp_coords:
        n_dims = int(gp_coord.attrib['dimensions'])
        sparseXs = gp_coord.findall('sparseX')
        cutoffs = {entry.attrib['i']: entry.attrib['sparseCutoff'] for entry in sparseXs}
        file = gp_coord.attrib['sparseX_filename']
        in_path = Path() / Path(file)
        out_path = in_path.parent / (in_path.name + '.input')
        cutoff_gen = iter(cutoffs.values())
        with open(in_path) as source, open(out_path, 'w') as target:
            for i, line in enumerate(source):
                if (i % n_dims == 0):
                    target.write(next(cutoff_gen) + '\n')
                target.write(line)


if __name__ == '__main__':
    main()
