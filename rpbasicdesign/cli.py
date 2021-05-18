#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""CLI for dnaprep"""

__author__ = 'Thomas Duigou'
__license__ = 'MIT'
__date__ = '2020.10.21'

import sys
import logging
import argparse

from rpbasicdesign.Designer import Designer


def __cli():
    """CLI for dnaprep."""

    desc = "Convert rpSBML enzyme info in to BASIC construct. UniProt IDs corresponding enzyme variants are extracted" \
           " from rpSBMl files. Promoters and RBSs are randomly chosen from a default list. CDSs, in other" \
           " words gene variants, of enzymes are randomly chosen from amongst the UniProt IDs extracted. Constructs" \
           " generated can be stored as (i) a CSV file ready to be used by DNA-Bot, (ii) as SBOL files. "
    parser = argparse.ArgumentParser(description=desc, prog='python -m dnaprep.cli')
    parser.add_argument('--rpsbml_file',
                        help='rpSBML file from which enzymes UniProt IDs will be collected',
                        required=True)
    parser.add_argument('--linker_parts_file',
                        help='File listing available linkers for constructs.', type=str)
    parser.add_argument('--user_parts_file',
                        help='File listing user parts (eg backbone, promoters) available for constructs.', type=str)
    parser.add_argument('--monocistronic',
                        help='Build monocistronic constructs. Default to false, ie polycistronic constructs will be '
                             'generated.',
                        default=False, type=lambda x: (str(x).lower() == 'true'))
    parser.add_argument('--lms_id', help='part ID to be used as the LMS methylated linker', default='LMS')
    parser.add_argument('--lmp_id', help='part ID to be used as the LMP methylated linker', default='LMP')
    parser.add_argument('--backbone_id', help='part ID to be used as the backbone', default='BASIC_SEVA_37_CmR-p15A.1')
    parser.add_argument('--sample_size', help='Number of construct to generate.', default=3, type=int)
    parser.add_argument('--o_dnabot_dir',
                        help='Output folder to write construct and coord part files. It will be created if it does '
                             'not exist yet. Existing files will be overwritten. Default: out/dnabot_in',
                        default='out/dnabot_in')
    parser.add_argument('--o_sbol_dir',
                        help='Output folder to write SBOL depictions of constructs. It will be created if it does not '
                             'exist yet. Existing files will be overwritten. Default: out/sbol_export',
                        default='out/sbol_export')

    # Logging
    logging.basicConfig(
            stream=sys.stderr, level=logging.INFO,
            datefmt='%d/%m/%Y %H:%M:%S',
            format='%(asctime)s -- %(levelname)s -- %(message)s'
            )

    # Compute
    args = parser.parse_args()
    o = Designer(monocistronic=args.monocistronic,
                 lms_id=args.lms_id, lmp_id=args.lmp_id, backbone_id=args.backbone_id,
                 linker_parts_file=args.linker_parts_file,
                 user_parts_file=args.user_parts_file,
                 )
    o.enzyme_from_rpsbml(rpsbml_file=args.rpsbml_file)
    nb_constructs = o.combine(sample_size=args.sample_size)
    logging.info(f'{nb_constructs} generated.')
    if args.o_dnabot_dir:
        nb_constructs = o.write_dnabot_inputs(out_dir=args.o_dnabot_dir)
        logging.info(f'{nb_constructs} constructs written in {args.o_dnabot_dir} folder.')
    if args.o_sbol_dir:
        o.write_sbol(out_dir=args.o_sbol_dir)
        logging.info(f'{nb_constructs} constructs written as SBOL in {args.o_sbol_dir} folder.')


if __name__ == "__main__":
    __cli()
