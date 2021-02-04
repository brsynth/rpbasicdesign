#!/usr/bin/env python
# -*- coding: utf-8 -*-


class Part:
    """Handle information on a single part.

    :param id: part ID
    :type id: str
    :param basic_role: structural role of the part according to BASIC method, (eg linker, part, backbone, ...)
    :type basic_role: str
    :param biological_role: biological role of the part, (eg promoter, rbs, cds)
    :type biological_role: str
    :param linker_class: linker class (only for linkers) or None
    :type linker_class: str or None
    :param cds_step: "step" of the CDS (eg first, second, third enzyme of a pathway...)
    :type cds_step: int or None
    :param seq: part DNA sequence
    :type seq: str
    :return: Part object
    :rtype: <Part>
    """

    def __init__(self, id, basic_role, biological_role,
                 linker_class=None, cds_step=None, seq=''):
        self.id = id
        self.basic_role = basic_role
        self.biological_role = biological_role
        self.linker_class = linker_class
        self.cds_step = cds_step
        self.seq = seq.lower()

    def __repr__(self):
        return str({
            'id': self.id,
            'basic_role': self.basic_role,
            'biological_role': self.biological_role,
            'linker_class': self.linker_class,
            'cds_step': self.cds_step,
            'seq': (self.seq[:10] + '...') if len(self.seq) > 11 else self.seq
        })

    def get_prefix_suffix(self):
        """Return the expected prefix suffix according to naming convention. Caution this should be used only
        for linkers.

        Notice:
        Prefix and suffix are named so-called when regarding the surrounded part, ie
        [linker A suffix]-[A first DNA part]-[linker B prefix] that would combine with
        [linker B suffix]-[Another DNA part]-[linker C prefix] which would give in the final construct
        [A first DNA part]-[linker B]-[Another DNA part]

        This means, when parsing a "full" linker:
         - the first half is the suffix
         - the second half is the prefix

         By convention:
         - prefix IDs end by "-P"
         - suffix IDs end by "-S"

        :return: prefix and suffix
        :rtype: tuple(str, str)
        """
        if self.basic_role == 'linker':
            if self.linker_class == 'rbs linker':
                suffix = f"{self.id.split('-')[0]}-S"
            else:
                suffix = f"{self.id}-S"
            prefix = f"{self.id}-P"
            return prefix, suffix
        return None

    def get_sbol_id(self):
        """Get a part ID compliant with SBOL rules.

        NOTICE: Invalid characters are replaced by "_"

        :returns: SBOL-ready part ID
        :rtype: str
        """
        valid_characters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_'
        i = 0
        pid = self.id
        while i < len(pid):
            if pid[i] not in valid_characters:
                pid = pid.replace(pid[i], '_')
            i += 1
        return pid
