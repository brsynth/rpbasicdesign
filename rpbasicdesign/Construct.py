#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging


from sbol2 import setHomespace, Document, ComponentDefinition, Sequence
from sbol2 import SO_PROMOTER, SO_CDS, SO_RBS, SO_MISC, SO_CIRCULAR
from rpbasicdesign import DNABOT_CONSTRUCT_HEADER


class Construct:
    """Everything about construct.

    :param backbone: backbone
    :type backbone: <Part>
    :param lms: suffix methylated linker
    :type lms: <Part>
    :param lmp: prefix methylated linker
    :type lmp: <Part>
    :param blocks: blocks embedding parts to be used
    :type blocks: list of dict of <Part>, [{k: <Part>}, .. ]
    :param nlinkers: available neutral linker (if junction necessary, used for monocistronic designs)
    :type nlinkers: list of <Part>
    :param monocistronic_design: True if construct should be monocistronic
    :type monocistronic_design: bool
    :return: a Construct object
    :rtype: <Construct>
    """
    def __init__(self, backbone, lms, lmp, blocks, nlinkers, monocistronic_design=True):
        self._nlinkers = nlinkers.copy()  # Prevent any side effect
        self._monocistronic_design = monocistronic_design

        self._parts = []
        self._parts.append(lms)
        self._parts.append(backbone)
        self._parts.append(lmp)
        for block_idx, block in enumerate(blocks):
            if monocistronic_design or block_idx == 0:
                self._parts.append(block['promoter'])
            self._parts.append(block['rbs'])
            self._parts.append(block['cds'])
            if monocistronic_design and block_idx != len(blocks) - 1:
                self._parts.append(self._nlinkers.pop(0))

    def __repr__(self):
        return str(self.get_construct_file_row())

    def __eq__(self, other):
        if not isinstance(other, Construct):
            return NotImplemented
        if len(self._parts) != len(other._parts):
            return False
        for idx in range(len(self._parts)):
            if self._parts[idx].id != other._parts[idx].id:
                return False
        return True

    @staticmethod
    def get_construct_file_header():
        """Return the header for DNABOT construct file

        :return: construct file header
        :rtype: list of str
        """
        return DNABOT_CONSTRUCT_HEADER

    def get_construct_file_row(self, coord=None):
        """Get the row format expected by DNA-Bot.

        :param coord: plate coordinate (opt)
        :type coord: str
        :return: formatted row to ready to be written
        :rtype: list of str
        """
        header_iter = iter(self.get_construct_file_header())
        row = {next(header_iter): coord}
        for part in self._parts:
            row[next(header_iter)] = part.id
        return row

    def has_duplicate_part(self):
        """Check if some parts are used several times.

        :return: True if any part is duplicated
        :rtype: bool
        """
        if len(self._parts) != len(set([part.id for part in self._parts])):
            return True
        return False

    def has_duplicate_suffix(self):
        """Check if some linkers are used several times

        :return: True if any linker is duplicated
        :rtype: bool
        """
        suffixes = [part.get_prefix_suffix()[1] for part in self._parts if part.basic_role == 'linker']
        if len(suffixes) != len(set(suffixes)):
            return True
        return False

    def get_part_ids(self, basic_roles=[]):
        """Get the list of part IDs involved in the construct

        :param basic_roles: part roles to include. If empty, parts of any basic_role are included
        :type basic_roles: list of str
        :return: part IDs
        :rtype: list of str
        """
        ids = []
        for part in self._parts:
            if len(basic_roles) == 0 or part.basic_role in basic_roles:
                ids.append(part.id)
        return ids

    def get_sbol(self, construct_id='BASIC_construct', validate=False):
        """Get the SBOL string representation of the construct.

        The object outputted is SBOL document which can be written
        to a file using the "writeString" method.

        WARNING: validation needs internet connexion.

        :param construct_id: Construct object ID
        :type construct_id: str
        :param validate: perform online SBOL validation
        :type validate: bool
        :return: SBOL object
        :rtype: <sbol.Document>
        """

        _SBOL_ROLE_ASSOC = {
            'misc': SO_MISC,
            'promoter': SO_PROMOTER,
            'rbs': SO_RBS,
            'cds': SO_CDS,
            'ori': SO_CIRCULAR
        }

        setHomespace('https://localhost')
        doc = Document()

        components = []
        for part in self._parts:
            component = ComponentDefinition(part.get_sbol_id())
            component.roles = _SBOL_ROLE_ASSOC[part.biological_role]
            component.sequence = Sequence(part.get_sbol_id(), part.seq)
            doc.addComponentDefinition(component)
            components.append(component)

        plasmid = ComponentDefinition(construct_id)
        doc.addComponentDefinition(plasmid)
        plasmid.assemblePrimaryStructure(components)

        if validate:
            logging.info(doc.validate())

        return doc