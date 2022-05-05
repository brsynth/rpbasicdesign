#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Methods for conversion"""

__author__ = 'Thomas Duigou'
__license__ = 'MIT'
__date__ = '2020.10.21'

import itertools
import logging
import os
import pprint
import random
import sys
from csv import DictReader, DictWriter
from copy import deepcopy
from math import factorial
from pathlib import Path
from typing import Dict
from operator import getitem

from rptools.rplibs import rpSBML, rpPathway
from rpbasicdesign import DNABOT_PART_HEADER
from rpbasicdesign.Construct import Construct
from rpbasicdesign.Part import Part


# CONSTANTS
MAX_ENZ_PER_RXN = 1
MAX_RXN_PER_CONSTRUCT = 3


def _gen_plate_coords(nb_row=8, nb_col=12):
    """Generator that generates the label coordinates for plate

    :param nb_row: number of rows in the plate (max: 26)
    :type nb_row: int
    :param nb_col: number of columns in the plate (max: 99)
    :type nb_col: int
    :returns: well coordinate
    :rtype: str
    """
    assert nb_row <= 26 and nb_col <= 99
    row_names = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    col_names = [str(_) for _ in range(1, 100)]
    for icol in range(nb_col):
        for jrow in range(nb_row):
            yield ''.join([row_names[jrow], col_names[icol]])


class Designer:
    """Convert rpSBML enzyme info in to BASIC construct.

    WARNING: promoters and RBSs are randomly chosen from amongst all available

    :param polycistronic: True if polycistronic design should be prepared
    :type polycistronic: bool
    :param verbose: True to increase log verbosity
    :type verbose: bool
    :param lms_id: part ID that corresponds to the LMS methylated linker
    :type lms_id: str
    :param lmp_id: part ID that corresponds to the LMP methylated linker
    :type lmp_id: str
    :param backbone_id: part ID that corresponds to the backbone
    :type backbone_id: str
    :param parts_files: files listing available linkers and user parts (backbone, promoters, ...)for constructs
    :type parts_files: list of str

    :return: Designer object
    :rtype: <Designer>
    """
    def __init__(self,
        # polycistronic=True,
        verbose=False,
        lms_id='LMS', lmp_id='LMP',
        backbone_id='BASIC_SEVA_37_CmR-p15A.1',
        parts_files=None,
        max_enz_per_rxn=MAX_ENZ_PER_RXN,
        max_rxn_per_construct=MAX_RXN_PER_CONSTRUCT
        ):
        # Default settings
        self._DATA_PATH = Path(__file__).resolve().parent / 'data'
        self._DEFAULT_DATA = [
            self._DATA_PATH / 'biolegio_parts.csv',
            self._DATA_PATH / 'user_parts.csv'
        ]

        # File to copy-paste, for the user convenience
        self._BIOLEGIO_PLATE_FILE = self._DATA_PATH / 'biolegio_plate.csv'

        self._verbose = verbose
        self._polycistronic_design = True
        self._lms_id = lms_id
        self._lmp_id = lmp_id
        self._backbone_id = backbone_id
        self._max_enz_per_rxn = max_enz_per_rxn
        self._max_rxn_per_construct = max_rxn_per_construct

        # Data files
        if parts_files is None:
            self._parts_files = self._DEFAULT_DATA
        else:
            self._parts_files = parts_files

        # Storage
        self._parts = {}
        self.constructs = []

        # Get resources
        self._get_linkers_and_parts_from_file()
        self._check_linkers_and_parts()

    def _get_backbone_part(self) -> Part:
        return self._parts[self._backbone_id]

    def _get_lms_part(self) -> Part:
        return self._parts[self._lms_id]
    
    def _get_lmp_part(self) -> Part:
        return self._parts[self._lmp_id]

    def _get_neutral_linker_parts(self) -> list:
        parts = []
        for part in self._parts.values():
            if part.linker_class == 'neutral linker':
                parts.append(part)
        return parts

    def _get_rbs_ids(self) -> list:
        """Return RBS part IDs"""
        part_ids = []
        for part in self._parts.values():
            if part.biological_role == 'rbs':
                part_ids.append(part.id)
        return sorted(part_ids)

    def _get_promoter_ids(self) -> list:
        """Return promoter parts"""
        part_ids = []
        for part in self._parts.values():
            if part.biological_role == 'promoter':
                part_ids.append(part.id)
        return sorted(part_ids)

    def _get_cds_ids(self) -> list:
        """Return CDS parts"""
        part_ids = []
        for part in self._parts.values():
            if part.biological_role == 'cds':
                part_ids.append(part.id)
        return sorted(part_ids)

    def _get_cds_ids_at_step(self, cds_step: str) -> list:
        """Return the CDS IDs for the specified step"""
        part_ids = []
        for part in self._parts.values():
            if part.biological_role == 'cds':
                if cds_step in part.cds_steps:
                    part_ids.append(part.id)
        return sorted(part_ids)

    def _get_parts_from_ids(self, part_ids: list) -> list:
        parts = []
        for part_id in part_ids:
            parts.append(self._get_part_from_id(part_id))
        return parts

    def _part_exists(self, part_id: str) -> bool:
        if part_id in self._parts:
            return True
        return False

    def _get_part_from_id(self, part_id: str) -> Part:
        return self._parts.get(part_id, None)

    def _get_rbs_ortho_id(self, part_id: str) -> str:
        assert len(part_id.split('-')) == 2
        return part_id.split('-')[0]

    def _get_rbs_seq_id(self, part_id: str) -> dict:
        assert len(part_id.split('-')) == 2
        return part_id.split('-')[1]
    
    def _gen_rbs_id(self, rbs_ortho_id: str, rbs_seq_id: str) -> str:
        return f'{rbs_ortho_id}-{rbs_seq_id}'

    def _get_linkers_and_parts_from_file(self):
        """Collect linkers and parts from the parts files."""
        for parts_file in self._parts_files:
            with open(parts_file) as ifh:
                for item in DictReader(ifh):
                    if item['id'].startswith('#'):  # Skip if commented
                        continue
                    elif item['id'] in self._parts:
                        logging.warning(f'Warning, part {item["id"]} duplicated, only the last definition kept.')
                    elif item['type'].lower() in ['neutral linker', 'methylated linker', 'peptide fusion linker']:
                        self._parts[item['id']] = Part(id=item['id'], basic_role='linker', biological_role='misc',
                                                    linker_class=item['type'].lower(), seq=item['sequence'])
                    elif item['type'].lower() == 'rbs linker':
                        self._parts[item['id']] = Part(id=item['id'], basic_role='linker', biological_role='rbs',
                                                    linker_class=item['type'].lower(), seq=item['sequence'])
                    elif item['type'].lower() == 'backbone':
                        self._parts[item['id']] = Part(id=item['id'], basic_role='backbone', biological_role='ori',
                                                    seq=item['sequence'])
                    elif item['type'].lower() == 'constitutive promoter':
                        self._parts[item['id']] = Part(id=item['id'], basic_role='part', biological_role='promoter',
                                                    seq=item['sequence'])
                    else:
                        logging.warning(f'Part "{item["id"]}" not imported because it does not fall any supported part '
                                        f'type.')

    def _check_linkers_and_parts(self):
        """Check that a minimal set of ressources have been stored."""
        _EXPECTED_ROLES = {'linker': 0, 'part': 0}
        for item in self._parts.values():
            if item.basic_role in _EXPECTED_ROLES.keys():
                _EXPECTED_ROLES[item.basic_role] += 1
        for role, count in _EXPECTED_ROLES.items():
            if count == 0:
                raise BaseException(f'No linker or part of type "{role}" provided. Does any have been provided? Exit.')

    def get_selenzyme_annotation(
        self, 
        rpsbml_path: str
        ) -> Dict:
        rpsbml = rpSBML(str(rpsbml_path))
        pathway = rpPathway.from_rpSBML(rpsbml=rpsbml)
        for idx_rxn, rxn_id in enumerate(pathway.get_reactions_ids()):
            # Stop if too many reactions
            if idx_rxn > self._max_rxn_per_construct:
                raise ValueError(
                    f'Number of reactions exceed the defined allowed number of ',
                    f'enzymes : {self._max_rxn_per_construct}. Execution cancelled.')
            # 
            rxn = pathway.get_reaction(rxn_id)
            enzymes = rxn.get_selenzy()
            # Stop if no enzyme available
            if len(enzymes) == 0:
                raise ValueError(
                    f'Missing UniProt IDs from selenzyme annotation for '
                    f'reaction {rxn_id}. Execution cancelled.')
            # Collect enzyme ordered by score, the first is the best
            for idx_enz, enz in enumerate(
                    sorted(enzymes.items(), key=lambda x: getitem(x[1], 'score'), reverse=True),
                    start=1):
                # Skip worst enzyme if too many
                if idx_enz > self._max_enz_per_rxn:
                    logging.warning(
                        f'Max number of enzyme per reaction reached ({self._max_enz_per_rxn}) '
                        f'for reaction {rxn_id}. Only the best one(s) are kept.')
                    break
                uniprot_id, _ = enz
                if uniprot_id in self._parts:
                    self._parts[uniprot_id].cds_steps.append(rxn_id)
                else:
                    self._parts[uniprot_id] = Part(id=uniprot_id, basic_role='part', 
                            biological_role='cds', cds_steps=[rxn_id],
                            seq='atgc')


    def combine(self, sample_size: int, random_seed: int =42, cds_permutation: bool =True) -> int:
        """Generate random constructs

        NOTICE:
        - special attention is made to prevent combination that would contain the same promoter, rbs and CDS.
        - special attention is made to prevent reusing the same RBS suffix in a given construct, to this end RBS
            linker IDs are expected to be in the form AAA-BBB with "AAA" being the linker suffix ID.
        - randomness is reset at the begining of each execution by `random.seed(random_seed)`

        :param sample_size: expected number of distinct constructs
        :type sample_size: int
        :param random_seed: seed to be used for random-related operations
        :type random_seed: int
        :param cds_permutation: whether all combinations of CDS permutation should be built 
        :type cds_permutation: bool
        :return: number of constructs generated
        :rtype: int
        """
        random.seed(random_seed)
        # Collect CDS ===============================================
        cds_ids = self._get_cds_ids()
        cds_steps = list()
        for cds in self._get_parts_from_ids(cds_ids):
            cds_steps += cds.cds_steps
        cds_steps = sorted(list(set(cds_steps)))
        # cds_steps = sorted(list(set([cds.cds_steps for cds in self._get_parts_from_ids(cds_ids)])))
        if len(cds_steps) == 0:
            logging.error(f'No CDS registered so far, generation of combination aborted')
            return 0
        # Collect promoters =========================================
        promoter_ids = self._get_promoter_ids()
        # Collect RBS ===============================================
        # Reminder: a RBS ID is composed of (i) the ID of the orthogonal RBS linker
        # and (ii) the ID of the RBS sequence.
        rbs_ids = self._get_rbs_ids()
        rbs_seq_ids = set()
        rbs_ortho_ids = set()
        for rbs_id in rbs_ids:
            rbs_seq_ids.add(self._get_rbs_seq_id(rbs_id))
            rbs_ortho_ids.add(self._get_rbs_ortho_id(rbs_id))
        rbs_seq_ids = sorted(list(rbs_seq_ids))
        rbs_ortho_ids = sorted(list(rbs_ortho_ids))
        # Checks ====================================================
        nb_required_promoters = 1 if self._polycistronic_design else len(cds_steps)
        if len(promoter_ids) < nb_required_promoters:
            logging.error(
                f'Not enough promoter parts provided: '
                f'{nb_required_promoters} required but '
                f'{len(promoter_ids)} provided. Exit')
            sys.exit()
        if len(rbs_ortho_ids) < len(cds_steps):
            logging.error(
                f'Not enough distinct ortholog sequences '
                'for building RBS linkers. Exit')
            sys.exit()
        # Print info ================================================
        max_combinations = len(promoter_ids) * (len(rbs_seq_ids) ** len(cds_steps))
        logging.info(f'Requested sample size: {sample_size}')
        logging.info(f'Min required promoter part: {nb_required_promoters}')
        logging.info(f'RBS seq IDs: {rbs_seq_ids}')
        logging.info(f'RBS ortho IDs: {rbs_ortho_ids}')
        logging.info(f'CDS Steps: {cds_steps}')
        logging.info(f'Perform CDS permutation? {cds_permutation}')
        logging.info(f"Max combinations -- without permutation: {max_combinations}")
        logging.info(f"Max combinations -- with    permutation: {max_combinations * factorial(len(cds_steps))}")
        # Argument checking =========================================
        if sample_size > 88:
            sample_size = 88
            logging.warning('Sample size is limited to 88. 88 constructs '
                            'will be return at the most.')
        # Build combinations ========================================
        # Promoters will be considered later
        distinct_blocks = {}
        for cds_step in cds_steps:
            if cds_step not in distinct_blocks:
                distinct_blocks[cds_step] = []
            for rbs_seq_id in rbs_seq_ids:
                for cds_id in self._get_cds_ids_at_step(cds_step):
                    block = {
                        'cds_step': cds_step,
                        'rbs_seq_id': rbs_seq_id,
                        'cds_id': cds_id
                        }
                    distinct_blocks[cds_step].append(block)
        # Build all combinations between RBS sequence and CDS
        rbs_cds_combis = list(itertools.product(*[l for l in distinct_blocks.values()]))
        # Add RBS ortho seq IDs (ie the BASIC adapters used in the RBS linkers)
        #   One construct should not use several time a given ortho seq
        tmp_combis = []
        for combi in rbs_cds_combis:
            tmp_combi = []
            for block_idx, block in enumerate(combi):
                tmp_block = deepcopy(block)
                tmp_block['rbs_ortho_id'] = rbs_ortho_ids[block_idx]
                tmp_combi.append(tmp_block)
            tmp_combis.append(tmp_combi)
        rbs_cds_combis = tmp_combis
        # TODO: remove combination that does not have corresponding part in self._part
        for combi in rbs_cds_combis:
            for block in combi:
                rbs_ortho_id = block['rbs_ortho_id']
                rbs_seq_id = block['rbs_seq_id']
                rbs_id = self._gen_rbs_id(rbs_ortho_id, rbs_seq_id)
                if not self._part_exists(rbs_id):
                    logging.error(f'RBS part {rbs_id} does not exist.')
        # Perform permutation of CDS if requested
        if cds_permutation:
            tmp = []
            for draft_construct in rbs_cds_combis:
                tmp.extend(list(itertools.permutations(draft_construct)))
            rbs_cds_combis = tmp
        else:
            # rbs_cds_combinations stay as they are
            pass  
        # Append promoter(s)
        pro_rbs_cds_combis = []
        if not self._polycistronic_design:
            # monocistronic: not yet implemented
            raise NotImplementedError(
                'Monocistronic constructs are not implemented.'
            )
        else:
            # polycistronic: associate one promoter to the first CDS only
            for promoter_id, combi in list(itertools.product(promoter_ids, rbs_cds_combis)):
                tmp = deepcopy(combi)
                tmp[0]['promoter_id'] = promoter_id
                pro_rbs_cds_combis.append(tmp)
        # Perform the sampling
        if sample_size > len(pro_rbs_cds_combis):
            sample_size = len(pro_rbs_cds_combis)
            logging.warning(
                'Requested sample size bigger than actual number of '
                'possible constructs. All constructs will be generated '
                'instead.'
                )
        combis_sample = random.sample(pro_rbs_cds_combis, k=sample_size)  # Remember the seed
        logging.info(f'Nb generated combinations: {len(combis_sample)}')
        # Instanciate the construct objects =========================
        for combi in combis_sample:
            blocks = []
            for item in combi:
                block = {}
                # Promoter if any defined
                if 'promoter_id' in item.keys():
                    promoter_id = item['promoter_id']
                    block['promoter'] = self._get_part_from_id(promoter_id)
                # RBS, combination of the orthologous BASIC seq and the RBS itself
                rbs_id = self._gen_rbs_id(item['rbs_ortho_id'], item['rbs_seq_id'])
                block['rbs'] = self._get_part_from_id(rbs_id)
                # CDS, ie the enzyme sequence
                cds_id = item['cds_id']
                block['cds'] = self._get_part_from_id(cds_id)            
                blocks.append(block)
            self.constructs.append(Construct(
                backbone=self._get_backbone_part(),
                lms=self._get_lms_part(),
                lmp=self._get_lmp_part(),
                blocks=blocks,
                nlinkers=self._get_neutral_linker_parts(),
                polycistronic_design=self._polycistronic_design
            ))

        return len(self.constructs)

    def write_dnabot_inputs(self, out_dir):
        """Write constructs in CSV format expected by DNA-Bot

        :param out_dir: folder path where construct and plate files be written
        :type out_dir: str
        :return: number of constructs written to file
        :rtype: int
        """
        __CONSTRUCT_FILE = 'constructs.csv'
        __USER_PLATE_FILE = 'user_parts_plate.csv'
        __BIOLEGIO_PLATE_FILE = 'biolegio_plate.csv'
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        # Construct file
        with open(os.path.join(out_dir, __CONSTRUCT_FILE), 'w') as ofh:
            plate_coords = _gen_plate_coords(nb_row=8, nb_col=12)
            writer = DictWriter(
                f=ofh,
                fieldnames=Construct.get_construct_file_header(),
                delimiter=',',
                restval=''
            )
            writer.writeheader()
            nb_constructs = 0
            for construct in self.constructs:
                nb_constructs += 1
                writer.writerow(construct.get_construct_file_row(coord=next(plate_coords)))
        # Custom parts (ie not biolegio / linker)
        with open(os.path.join(out_dir, __USER_PLATE_FILE), 'w') as ofh:
            plate_coords = _gen_plate_coords(nb_row=8, nb_col=12)
            writer = DictWriter(
                f=ofh,
                fieldnames=DNABOT_PART_HEADER,
                delimiter=',',
                restval=''
            )
            writer.writeheader()
            part_ids = set()
            for construct in self.constructs:  # Collect parts that are used
                part_ids |= set(construct.get_part_ids(basic_roles=['part', 'backbone']))
            for _ in sorted(part_ids):
                writer.writerow({'Part/linker': _, 'Well': next(plate_coords), 'Part concentration (ng/uL)': ''})
        # Biolegio plate
        from_file = self._BIOLEGIO_PLATE_FILE
        to_file = Path(out_dir) / __BIOLEGIO_PLATE_FILE
        import shutil
        shutil.copy(from_file, to_file)
        return nb_constructs

    def write_sbol(self, out_dir):
        """Write constructs as SBOL files

        NOTICE: out_dir will be created if not existing yet.

        :param out_dir: out dir path
        :type out_dir: str
        :return: number of SBOL files written
        :rtype: int
        """
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        plate_coords = _gen_plate_coords(nb_row=8, nb_col=12)
        nb_constructs = 0
        for construct in self.constructs:
            nb_constructs += 1
            coord = next(plate_coords)
            fname = 'BASIC_construct_{}.xml'.format(coord)
            construct_id = 'BASIC_construct_{}'.format(coord)
            with open(os.path.join(out_dir, fname), 'w') as ofh:
                ofh.write(construct.get_sbol(construct_id=construct_id, validate=False).writeString())
        return nb_constructs

    # def _iter_sample_fast_custom(self, iterable, sample_size):
    #     """Reservoir sampling.
    #
    #     NOTICE: this approach is not adapted for large combinatorics (billions of combination) as it
    #     still require to iterate other each combination, hence it takes a lot of time.
    #
    #     :param iterable: iterable variable to be sampled
    #     :param sample_size: expected size of the sample
    #     :return: list, the final sample
    #
    #     Adapted from: https://stackoverflow.com/questions/12581437/python-random-sample-with-a-generator-iterable-iterator
    #     """
    #     results = []
    #     iterator = iter(iterable)
    #     # Fill in the first sample_size elements:
    #     idx = 0
    #     while len(results) != sample_size:
    #         try:
    #             item = next(iterator)
    #         except StopIteration:
    #             raise ValueError("Sample larger than population.")
    #         if not part_duplicated(item, types=['promoter', 'rbs', 'cds']):
    #             results.append(item)
    #             idx += 1
    #     random.shuffle(results)  # Randomize their positions
    #     for item in iterator:
    #         if not part_duplicated(item, types=['promoter', 'rbs', 'cds']):
    #             r = random.randint(0, idx)  # at a decreasing rate based on idx, replace random items
    #             if r < sample_size:
    #                 results[r] = item
    #             idx += 1  # idx only incremented for "valid" (not duplicates) combination
    #
    #     if len(results) < sample_size:
    #         raise ValueError("Sample larger than population.")
    #     return results
