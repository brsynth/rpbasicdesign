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

     desc = "Convert rpSBML enzyme info in to BASIC construct. UniProt IDs"  \
               " corresponding enzyme variants are extracted rpSBMl files."     \
               " Promoters and RBSs are randomly chosen from a default list."   \
               " CDSs, in other words gene variants, of enzymes are randomly"   \
               " chosen from amongst the UniProt IDs extracted. Constructs"     \
               " generated can be stored as (i) a CSV file ready to be used by" \
               " DNA-Bot, (ii) as SBOL files."
     parser = argparse.ArgumentParser(
          description=desc,
          prog='python -m rpbasicdesign.cli'
          )
     parser.add_argument(
          '--rpsbml_file',
               help='rpSBML file from which enzymes UniProt IDs will be collected.',
               required=True
               )
     parser.add_argument(
          '--parts_files',
               help='List of files providing available linkers and user parts '
                    '(backbone, promoters, ...) for constructs. '
                    'Default: [data/biolegio_parts.csv, user_parts.csv]',
               type=str,
               nargs='+'
               )
     parser.add_argument(
          '--lms_id',
          help='part ID to be used as the LMS methylated linker. '
               'Default: LMS',
          default='LMS'
          )
     parser.add_argument(
          '--lmp_id',
          help='part ID to be used as the LMP methylated linker. '
               'Default: LMP',
          default='LMP'
          )
     parser.add_argument(
          '--backbone_id',
          help='part ID to be used as the backbone. '
               'Default: BASIC_SEVA_37_CmR-p15A.1',
          default='BASIC_SEVA_37_CmR-p15A.1'
          )
     parser.add_argument(
          '--sample_size',
          help='Number of construct to generate.'
               'Default: 88',
          default=88,
          type=int
          )
     parser.add_argument(
          '--cds_permutation',
          help='Whether all combinations of CDS permutation should be built '
               'Default: true',
          default=True,
          type=lambda x: (not str(x).lower() == 'false')
          )
     parser.add_argument(
          '--max_enz_per_rxn',
          help='Maximum number of enyzme to consider per reaction. '
               'If more enzymes are available for a given reaction, '
               'then only the last one listed in the MIRIAM annotation '
               'section will be kept.',
          default=1,
          type=int
          )
     parser.add_argument(
          '--o_dnabot_dir',
          help='Output folder to write construct and plate files. '
               'It will be created if it does not exist yet. Existing '
               'files will be overwritten. '
               'Default: out/dnabot_in',
          default='out/dnabot_in'
          )
     parser.add_argument(
          '--o_sbol_dir',
          help='Output folder to write SBOL depictions of constructs. '
               'It will be created if it does not exist yet. Existing '
               'files will be overwritten. '
               'Default: out/sbol_export',
          default='out/sbol_export'
          )

     # Logging
     logging.basicConfig(
          stream=sys.stderr, level=logging.INFO,
          datefmt='%d/%m/%Y %H:%M:%S',
          format='%(asctime)s -- %(levelname)s -- %(message)s'
          )

     # Compute
     args = parser.parse_args()
     o = Designer(
          # polycistronic=args.polycistronic,
          lms_id=args.lms_id,
          lmp_id=args.lmp_id,
          backbone_id=args.backbone_id,
          parts_files=args.parts_files,
          max_enz_per_rxn=args.max_enz_per_rxn
          )
     o.get_selenzyme_annotation(rpsbml_path=args.rpsbml_file)
     nb_constructs = o.combine(
          sample_size=args.sample_size,
          cds_permutation=args.cds_permutation
          )
     logging.info(f'{nb_constructs} constructs generated.')
     # Write
     if args.o_dnabot_dir:
          nb_constructs = o.write_dnabot_inputs(out_dir=args.o_dnabot_dir)
          logging.info(f'{nb_constructs} constructs written as CSV in {args.o_dnabot_dir} folder.')
     if args.o_sbol_dir:
          o.write_sbol(out_dir=args.o_sbol_dir)
          logging.info(f'{nb_constructs} constructs written as SBOL in {args.o_sbol_dir} folder.')


if __name__ == "__main__":
     __cli()
