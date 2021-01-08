#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Methods for conversion"""

__author__ = 'Thomas Duigou'
__license__ = 'MIT'
__date__ = '2020.10.21'

import os
import logging
import pkgutil
import random
from csv import DictReader, DictWriter

from libsbml import SBMLReader
from sbol import setHomespace, Document, ComponentDefinition, Sequence
from sbol import SO_PROMOTER, SO_CDS, SO_RBS, SO_MISC, SO_PLASMID, SO_CIRCULAR

logging.basicConfig(level=logging.DEBUG)

_DNABOT_CONSTRUCT_HEADER = ['Well',
                            'Linker 1', 'Part 1', 'Linker 2', 'Part 2', 'Linker 3', 'Part 3', 'Linker 4', 'Part 4',
                            'Linker 5', 'Part 5', 'Linker 6', 'Part 6', 'Linker 7', 'Part 7', 'Linker 8', 'Part 8',
                            'Linker 9', 'Part 9', 'Linker 10', 'Part 10',
                            ]
_DNABOT_PART_HEADER = ['Part/linker', 'Well', 'Part concentration (ng/uL)']


def _gen_plate_coords(nb_row=8, nb_col=12, by_row=True):
    """Generator that generates the label coordinates for plate

    :param nb_row: int, number of rows in the plate (max: 26)
    :param nb_col: int, number of columns in the plate (max: 99)
    :param by_row, bool, True to work by row
    """
    assert nb_row <= 26 and nb_col <= 99
    row_names = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    col_names = [str(_) for _ in range(1, 100)]
    if by_row:
        i_dim = row_names
        i_nb = nb_row
        j_dim = col_names
        j_nb = nb_col
    else:
        i_dim = col_names
        i_nb = nb_col
        j_dim = row_names
        j_nb = nb_row
    for i in range(i_nb):
        for j in range(j_nb):
            yield ''.join([i_dim[i], j_dim[j]])


class Part:
    """Handle information on a single part.

    :param id: str, part ID
    :param role: str, part role, (eg promoter, rbs, cds)
    :param type: str, part type (eg linker, ..)
    :param seq: str, part DNA sequence
    :return: <Part>
    """

    def __init__(self, id, role, type=None, seq=''):
        self.id = id
        self.role = role
        self.type = type
        self.seq = seq.lower()

    def __repr__(self):
        return str({
            'id': self.id,
            'role': self.role,
            'type': self.type,
            'seq': (self.seq[:10] + '...') if len(self.seq) > 11 else self.seq
        })

    def get_sbol_id(self):
        """Get a part ID compliant with SBOL rules.

        NOTICE: Invalid characters are replaced by "_"
        """
        valid_characters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_'
        i = 0
        pid = self.id
        while i < len(pid):
            if pid[i] not in valid_characters:
                pid = pid.replace(pid[i], '_')
            i += 1
        return pid


class BASICConstruct:
    """Everything about construct.

    :param backbone: <Part>, backbone
    :param LMS: <Part>, suffix methylated linker
    :param LMP: <Part>, prefix methylated linker
    :param blocks: [{k: <Part>}, .. ], blocks to be included
    :param nlinkers: [<Part>, ..], available neutral linker
    :param monocistronic_design: bool, True if construct should be monocistronic
    :return: <BASICConstruct>
    """

    def __init__(self, backbone, LMS, LMP, blocks, nlinkers, monocistronic_design=True):
        self._nlinkers = nlinkers.copy()  # Prevent any side effect
        self._monocistronic_design = monocistronic_design

        self._parts = []
        self._parts.append(LMS)
        self._parts.append(backbone)
        self._parts.append(LMP)
        for block_idx, block in enumerate(blocks):
            if monocistronic_design or block_idx == 0:
                self._parts.append(block['promoter'])
            self._parts.append(block['rbs'])
            self._parts.append(block['cds'])
            if monocistronic_design and block_idx != len(blocks) - 1:
                self._parts.append(self._nlinkers.pop(0))

    def __repr__(self):
        return str(self.get_dnabot_row())

    def __eq__(self, other):
        if not isinstance(other, BASICConstruct):
            return NotImplemented
        if len(self._parts) != len(other._parts):
            return False
        for idx in range(len(self._parts)):
            if self._parts[idx].id != other._parts[idx].id:
                return False
        return True

    def get_dnabot_row(self, coord=None):
        """Get the row format expected by DNA-Bot.

        :param coord: str, plate coordinate (opt)
        :return: [str, ..], row
        """
        header_iter = iter(_DNABOT_CONSTRUCT_HEADER)
        row = {next(header_iter): coord}
        for part in self._parts:
            row[next(header_iter)] = part.id
        return row

    def has_duplicates(self):
        """Check if some parts are used several times.

        :return: bool, True if any parts is duplicated
        """
        if len(self._parts) != len(set([part.id for part in self._parts])):
            return True
        return False

    def get_part_ids(self, role=[]):
        """Get the list of part IDs involved in the construct

        :param role: [str, ...]: list of part role to include. If empty, parts of any role are included
        :return: [str, ...]: list of part IDs
        """
        ids = []
        for part in self._parts:
            if len(role) == 0 or part.role in role:
                ids.append(part.id)
        return ids

    def get_sbol(self, construct_id='BASIC_construct', validate=False):
        """Get the SBOL string representation of the construct.

        The object outputted is SBOL document which can be written
        to a file using the "writeString" method.

        WARNING: validation needs internet connexion.

        :param construct_id: str, BASICConstruct object ID
        :param validate: bool, perform online SBOL validation
        :return: <sbol.Document>, SBOL object
        """

        _SBOL_ROLE_ASSOC = {
            'misc': SO_MISC,
            'promoter': SO_PROMOTER,
            'rbs': SO_RBS,
            'cds': SO_CDS,
            'backbone': SO_CIRCULAR
        }

        setHomespace('https://localhost')
        doc = Document()

        components = []
        for part in self._parts:
            component = ComponentDefinition(part.get_sbol_id())
            component.roles = _SBOL_ROLE_ASSOC[part.role]
            component.sequence = Sequence(part.get_sbol_id(), part.seq)
            doc.addComponentDefinition(component)
            components.append(component)

        plasmid = ComponentDefinition(construct_id)
        doc.addComponentDefinition(plasmid)
        plasmid.assemblePrimaryStructure(components)

        if validate:
            logging.info(doc.validate())

        return doc


class BASICDesigner:
    """Convert rpSBML enzyme info in to BASIC construct.

    WARNING: promoters and RBSs are randomly chosen from amongst all available

    :return: <BASICDesigner>
    """

    def __init__(self, monocistronic_design=True, verbose=False):
        self._default_part_file = 'data/default_parts.tsv'
        self._ref_coord_part_file = 'data/default_coords.csv'
        self._verbose = verbose
        self._monocistronic_design = monocistronic_design
        self._parts = {}
        self.constructs = []

        self._MAX_ENZ = 5
        self._SEED = 42

        # Load default data
        content = pkgutil.get_data(__name__, self._default_part_file).decode().splitlines()
        for item in DictReader(content, delimiter='\t'):
            if item['id'].startswith('#'):  # Skip if commented
                continue
            if item['id'] in self._parts:
                logging.warning(f'Warning, part {id} duplicated, only the last definition kept.')
            if item['type'].lower() == 'backbone':
                self._parts[item['id']] = Part(id=item['id'], role='backbone',
                                               type=item['type'].lower(), seq=item['sequence']
                                               )
            elif item['type'].lower() == 'constitutive promoter':
                self._parts[item['id']] = Part(id=item['id'], role='promoter',
                                               type=item['type'].lower(), seq=item['sequence']
                                               )
            elif item['type'].lower() in ['neutral linker', 'methylated linker']:
                self._parts[item['id']] = Part(id=item['id'], role='misc',
                                               type=item['type'].lower(), seq=item['sequence']
                                               )
            elif item['type'].lower() == 'rbs linker':
                self._parts[item['id']] = Part(id=item['id'], role='rbs',
                                               type=item['type'].lower(), seq=item['sequence']
                                               )
            else:
                logging.warning(f'Part "{id}" not imported because it does not fall any supported part type.')

    def _read_MIRIAM_annotation(self, annot):
        """Return the MIRIAM annotations of species.

        :param annot: SBML annotation block
        :return: dict of annotations
        """
        try:
            to_keep = {}
            bag = annot.getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
            for i in range(bag.getNumChildren()):
                str_annot = bag.getChild(i).getAttrValue(0)
                if str_annot == '':
                    logging.warning('This contains no attributes: ' + str(bag.getChild(i).toXMLString()))
                    continue
                dbid = str_annot.split('/')[-2].split('.')[0]
                if len(str_annot.split('/')[-1].split(':')) == 2:
                    cid = str_annot.split('/')[-1].split(':')[1]
                else:
                    cid = str_annot.split('/')[-1]
                if dbid not in to_keep:
                    to_keep[dbid] = []
                to_keep[dbid].append(cid)
            return to_keep
        except AttributeError:
            return {}

    def enzyme_from_rpsbml(self, rpsbml_file):
        """Extract enzyme from rpSBML annotation

        WARNING: the rpSBML file is expected to follow the specific schema of annotations used for rpSBMLs

        :param rpsbml_file: str, path to rpSBML from which enzyme will be extracted.
        """
        reader = SBMLReader()
        document = reader.readSBML(rpsbml_file)
        model = document.getModel()
        for nb_rxn, reaction in enumerate(model.getListOfReactions()):
            if nb_rxn > self._MAX_ENZ:
                logging.warning(f'Number of reactions exceed the defined allowed number of enzymes {self._MAX_ENZ}')
            if not reaction.id.endswith('_sink'):
                annot = reaction.getAnnotation()
                for uni_id in self._read_MIRIAM_annotation(annot)['uniprot']:
                    if uni_id in self._parts and self._verbose:
                        logging.warning(f'Warning, part {uni_id} duplicated, only the last definition kept.')
                    self._parts[uni_id] = Part(id=uni_id, role='cds', type=f'enzyme {reaction.id}', seq='tata')

    def combine(self, sample_size, lms_id='LMS', lmp_id='LMP', backbone_id='BASIC_SEVA_37_CmR-p15A.1'):
        """Generate random constructs

        NOTICE: special attention is made to prevent combination that would contain the same promoter, rbs and CDS.

        :param sample_size: int, expected number of distinct constructs
        :param lms_id: str, part ID that corresponds to the LMS methylated linker
        :param lmp_id: str, part ID that corresponds to the LMP methylated linker
        :param backbone_id: str, part ID that corresponds to the backbone
        :return: number of constructs generated
        """
        random.seed(self._SEED)
        cds_levels = set([part.type for part in self._parts.values() if part.role == 'cds'])
        if len(cds_levels) == 0:
            logging.error(f'No CDS registered so far, method aborted')
            return []
        nb_promoter = len(cds_levels) if self._monocistronic_design else 1
        nb_dupe_parts = 0
        nb_dupe_constructs = 0
        while len(self.constructs) < sample_size:
            # Get variants
            promoters = random.choices(
                population=[part.id for part in self._parts.values() if part.role == 'promoter'],
                k=nb_promoter
            )
            rbss = random.choices(
                population=[part.id for part in self._parts.values() if part.role == 'rbs'],
                k=len(cds_levels)
            )
            cdss = []
            for cds_lvl in cds_levels:
                variants = [part.id for part in self._parts.values() if part.type == cds_lvl]
                cdss.append(random.choice(variants))
            # Build blocks
            blocks = []
            for idx in range(len(cds_levels)):
                block = {'rbs': self._parts[rbss[idx]], 'cds': self._parts[cdss[idx]]}
                if idx == 0 or self._monocistronic_design:
                    block['promoter'] = self._parts[promoters[idx]]
                blocks.append(block)
            # Instantiate
            construct = BASICConstruct(
                backbone=self._parts['BASIC_SEVA_37_CmR-p15A.1'],
                LMS=self._parts['LMS'],
                LMP=self._parts['LMP'],
                blocks=blocks,
                nlinkers=[part for part in self._parts.values() if part.type == 'neutral linker'],
                monocistronic_design=self._monocistronic_design
            )
            if construct.has_duplicates():
                nb_dupe_parts += 1
                continue
            elif construct in self.constructs:
                nb_dupe_constructs += 1
                continue
            else:
                self.constructs.append(construct)
        # Logs
        if nb_dupe_parts:
            logging.debug(f'{nb_dupe_parts} constructs skipped because of duplicated parts within a construct.')
        if nb_dupe_constructs:
            logging.debug(f'{nb_dupe_constructs} constructs skipped because the construct already exists.')
        return len(self.constructs)

    def write_dnabot_inputs(self, out_dir):
        """Write constructs in CSV format expected by DNA-Bot

        :param out_dir: str, folder path where construct and coord files be written
        :return: int, number of constructs written to file
        """
        __CONSTRUCT_FILE = 'constructs.csv'
        __COORD_REF_FILE = 'ref_parts_coord.csv'
        __COORD_CUSTOM_FILE = 'parts_coord.csv'
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        with open(os.path.join(out_dir, __CONSTRUCT_FILE), 'w') as ofh:
            plate_coords = _gen_plate_coords(nb_row=8, nb_col=12, by_row=True)
            writer = DictWriter(f=ofh, fieldnames=_DNABOT_CONSTRUCT_HEADER, delimiter=',', restval='')
            writer.writeheader()
            nb_constructs = 0
            for construct in self.constructs:
                nb_constructs += 1
                writer.writerow(construct.get_dnabot_row(coord=next(plate_coords)))
        with open(os.path.join(out_dir, __COORD_CUSTOM_FILE), 'w') as ofh:
            plate_coords = _gen_plate_coords(nb_row=8, nb_col=12, by_row=True)
            writer = DictWriter(f=ofh, fieldnames=_DNABOT_PART_HEADER, delimiter=',', restval='')
            writer.writeheader()
            part_ids = set()
            for construct in self.constructs:  # Collect parts that are used
                part_ids |= set(construct.get_part_ids(role=['cds']))
            for _ in sorted(part_ids):
                writer.writerow({'Part/linker': _, 'Well': next(plate_coords), 'Part concentration (ng/uL)': ''})
        with open(os.path.join(out_dir, __COORD_REF_FILE), 'wb') as ofh:
            ofh.write(pkgutil.get_data(__name__, self._ref_coord_part_file))
        return nb_constructs

    def write_sbol(self, out_dir):
        """Write constructs as SBOL files

        NOTICE: out_dir will be created if not existing yet.

        :param out_dir: str, out dir path
        :return: int, number of SBOL files written
        """
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        plate_coords = _gen_plate_coords(nb_row=8, nb_col=12, by_row=True)
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

    # def combine_old(self, sample_size=3, lms_id='LMS', lmp_id='LMP', backbone_id='BASIC_SEVA_37_CmR-p15A.1'):
    #     """Generate the combinatorics.
    #
    #     NOTICE: special attention is made to prevent combination that would contain the same promoter, rbs and CDS.
    #
    #     :param lms_id: str, part ID corresponding to the LMS methylated BASIC linker
    #     :param lmp_id: str, part ID corresponding to the LMP methylated BASIC linker
    #     :param backbone_id: str, part ID corresponding to the backbone
    #     :return constructs: [<BASICConstruct>, ...]
    #
    #     For curious people, here an example of how to know the total number of combination. Let's consider the
    #     following number of possible parts for each promoter, rbs, cds:
    #     - 5 promoters
    #     - 9 RBSs
    #     - 3 enzymes, 1 sequence for each
    #     ... then the number of combination that does not reuse promoter nor rbs
    #     and provide 3 monocistronic units: (5*4*3)*(9*8*7) = 30240
    #     """
    #     # Building unit blocks
    #     to_combine = [
    #         [{'promoter': part.id} for part in self._parts.values() if part.role == 'promoter'],
    #         [{'rbs': part.id} for part in self._parts.values() if part.role == 'rbs'],
    #         [{'cds': part.id, 'step': part.type} for part in self._parts.values() if part.role == 'cds']
    #     ]
    #     # Reshape
    #     logging.debug('Reshape')
    #     unit_blocks = []
    #     for item in product(*to_combine):
    #         block = {}
    #         for d in item:
    #             block.update(d)
    #         unit_blocks.append(block)
    #     # Split by CDS step
    #     logging.debug('Split by CDS step')
    #     by_step_blocks = {}
    #     for block in unit_blocks:
    #         step = block['step']
    #         if step not in by_step_blocks:
    #             by_step_blocks[step] = []
    #         by_step_blocks[step].append(block)
    #     # Sampling unit blocks combination using reservoir sampling
    #     logging.debug('Sampling unit blocks combination using reservoir sampling')
    #     by_step_blocks_combinatorics = product(*[block for block in by_step_blocks.values()])
    #     logging.debug('Sampling unit blocks combination using reservoir sampling, second step')
    #     block_combinations = self._iter_sample_fast_custom(iterable=by_step_blocks_combinatorics, sample_size=sample_size)
    #     # Switch ID with full description
    #     logging.debug('Switch ID with full description')
    #     for blocks in block_combinations:
    #         for block in blocks:
    #             for part_type, pid in block.items():
    #                 if pid in self._parts:
    #                     block[part_type] = self._parts[pid]
    #     # Build the constructs
    #     constructs = []
    #     for blocks in block_combinations:
    #         constructs.append(
    #             BASICConstruct(
    #                 backbone=self._parts[backbone_id],
    #                 LMS=self._parts[lms_id],
    #                 LMP=self._parts[lmp_id],
    #                 blocks=blocks,
    #                 nlinkers=[part for part in self._parts.values() if part.type == 'neutral linker'],
    #                 monocistronic_design=self._monocistronic_design
    #             )
    #         )
    #     return constructs
