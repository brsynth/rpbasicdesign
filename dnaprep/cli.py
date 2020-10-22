#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""CLI for dnaprep"""

__author__ = 'Thomas Duigou'
__license__ = 'MIT'
__date__ = '2020.10.21'

import sys
import logging
import argparse

from dnaprep import BASICDesigner


def __cli():
    """CLI for dnaprep."""

    help = "Convert rpSBML enzyme info in to BASIC construct. UniProt IDs corresponding enzyme variants are extracted" \
           " from rpSBMl files. Promoters and RBSs are randomly chosen from a default list (see TODO). CDSs, in other" \
           " words gene variants, of enzymes are randomly chosen from amongst the UniProt IDs extracted. Constructs" \
           " generated can be stored as (i) a CSV file ready to be used by DNA-Bot, (ii) as SBOL files. "
    parser = argparse.ArgumentParser(description=help, prog='python -m dnaprep.cli')
    parser.add_argument('--monocistronic_design', help='Build monocistronic constructs. Default to false, ie'
                                                       'polycistronic constructs will be generated.',
                        default=True, type=lambda x: (str(x).lower() == 'false'))
    # parser.add_argument('--verbose, -v', help="Let's the tool express itself to its full potential ;)",
    #                     default=False)
    parser.add_argument('--rpsbml_file', help='rpSBML file from which enzymes UniProt IDs will be collected',
                        required=True)
    parser.add_argument('--sample_size', help='Number of construct to generate.', default=6, type=int)
    parser.add_argument('--o_dnabot_file', help='Output file ready to be used by DNA-Bot. Existing file will be'
                                                'overwritten.')
    parser.add_argument('--o_sbol_dir', help='Output folder to write SBOL depictions of constructs. It will be '
                                             'created if it does not exist yet. Existing files will be overwritten.')

    # Logging
    logging.basicConfig(
            stream=sys.stderr, level=logging.INFO,
            datefmt='%d/%m/%Y %H:%M:%S',
            format='%(asctime)s -- %(levelname)s -- %(message)s'
            )

    # Compute
    args = parser.parse_args()
    o = BASICDesigner.BASICDesigner(monocistronic_design=args.monocistronic_design)
    o.enzyme_from_rpsbml(rpsbml_file=args.rpsbml_file)
    nb_constructs = o.combine(sample_size=args.sample_size)
    logging.info(f'{nb_constructs} generated.')
    if args.o_dnabot_file:
        nb_constructs = o.write_dnabot_input(out_file=args.o_dnabot_file)
        logging.info(f'{nb_constructs} constructs written in {args.o_dnabot_file} file.')
    if args.o_sbol_dir:
        o.write_sbol(out_dir=args.o_sbol_dir)
        logging.info(f'{nb_constructs} constructs written as SBOL in {args.o_dnabot_file} folder.')


if __name__ == "__main__":
    __cli()
