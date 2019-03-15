#!/usr/bin/env python

"""
Combine GAP XML

Take a number of GAP xml files and produce a simgle xml file with
all the descriptors combined in a single potential.

"""

import argparse
from os import path

from quippy.qpxml import combine_xml


def main():
    """
    Run the processing.
    """

    parser = argparse.ArgumentParser(description='Tool to combine GAP xml')
    parser.add_argument('-f', '--force', action='store_true',
                        help='Overwrite existing output file.')
    parser.add_argument('xml_filenames', nargs='*',
                        help='Input files to be combined, first one is used '
                             'as the base file.')
    parser.add_argument('output_filename',
                        help='Name of file to output combined xml')
    parser.add_argument('-k', '--keep-xyz', action='store_true',
                        help='Keep xyz data in output potentials, '
                             'only keeps data from the first potential.')
    parser.add_argument('-l', '--label',
                        help='Use the custom label for the output potential.')

    args = parser.parse_args()

    # No trashing existing files by accident
    if path.exists(args.output_filename) and not args.force:
        raise IOError('Output file exists: {}'.format(args.output_filename))

    base_filename = args.xml_filenames[0]
    extra_filenames = args.xml_filenames[1:]
    remove_xyz = not args.keep_xyz
    label = args.label

    with open(args.output_filename, 'w') as output_file:
        output_file.write(combine_xml(base_filename, extra_filenames,
                                      remove_xyz, label))


if __name__ == "__main__":
    main()
