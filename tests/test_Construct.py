"""
Test Construct
"""
import pytest
from hashlib import sha1
from rpbasicdesign.Part import Part
from rpbasicdesign.Construct import Construct

_LR1 = Part(id='LR1', basic_role='linker', biological_role='rbs', linker_class='rbs linker')
_LR2 = Part(id='LR2', basic_role='linker', biological_role='rbs', linker_class='rbs linker')
_CDS1 = Part(id='CDS1', basic_role='part', biological_role='cds')
_CDS2 = Part(id='CDS2', basic_role='part', biological_role='cds')
_P1 = Part(id='P1', basic_role='part', biological_role='promoter')
_P2 = Part(id='P2', basic_role='part', biological_role='promoter')
_LN1 = Part(id='LN1', basic_role='linker', biological_role='misc', linker_class='neutral linker')
_LN2 = Part(id='LN2', basic_role='linker', biological_role='misc', linker_class='neutral linker')
_BB = Part(id='BB', basic_role='backbone', biological_role='misc')
_LMS = Part(id='LMS', basic_role='methylated linker', biological_role='misc')
_LMP = Part(id='LMP', basic_role='methylated linker', biological_role='misc')

__BLOCKS = [
    {'promoter': _P1, 'rbs': _LR1, 'cds': _CDS1},
    {'promoter': _P2, 'rbs': _LR2, 'cds': _CDS2}
]


def test_init_1():
    # Empty is not OK
    with pytest.raises(TypeError):
        Construct()


def test_init_2():
    # Keywords are not controlled
    Construct(
        backbone='cubitus', lms='lms', lmp='lmp',
        blocks=[], nlinkers=[]
    )


def test_str():
    o = Construct(
        backbone=_BB, lms=_LMS, lmp=_LMP,
        blocks=__BLOCKS, nlinkers=[_LN1, _LN2]
    )
    str(o)


def test_eq():
    C1 = Construct(
        backbone=_BB, lms=_LMS, lmp=_LMP,
        blocks=__BLOCKS, nlinkers=[_LN1, _LN2],
        monocistronic_design=True
    )
    C2 = Construct(
        backbone=_BB, lms=_LMS, lmp=_LMP,
        blocks=__BLOCKS, nlinkers=[_LN1, _LN2],

        monocistronic_design=False
    )
    assert C1 == C1
    assert C1 != C2


def test_construct_file_header():
    # Expect a list, to be used as header file
    assert isinstance(Construct.get_construct_file_header(), list)


def test_has_duplicate_part():
    # CDS are considered as parts
    o = Construct(
        backbone=_BB, lms=_LMS, lmp=_LMP,
        blocks=[
            {'promoter': _P1, 'rbs': _LR1, 'cds': _CDS1},
            {'promoter': _P2, 'rbs': _LR2, 'cds': _CDS2}
            ],
        nlinkers=[_LN1, _LN2]
    )
    assert not o.has_duplicate_part()
    o = Construct(
        backbone=_BB, lms=_LMS, lmp=_LMP,
        blocks=[
            {'promoter': _P1, 'rbs': _LR1, 'cds': _CDS1},
            {'promoter': _P2, 'rbs': _LR2, 'cds': _CDS1}
            ],
        nlinkers=[_LN1, _LN2]
    )
    assert o.has_duplicate_part()
    # Linker are also considered as part here
    o = Construct(
        backbone=_BB, lms=_LMS, lmp=_LMP,
        blocks=[
            {'promoter': _P1, 'rbs': _LR1, 'cds': _CDS1},
            {'promoter': _P2, 'rbs': _LR1, 'cds': _CDS2}
            ],
        nlinkers=[_LN1, _LN2]
    )
    assert o.has_duplicate_part()


def test_has_duplicate_suffix():
    o = Construct(
        backbone=_BB, lms=_LMS, lmp=_LMP,
        blocks=[
            {'promoter': _P1, 'rbs': _LR1, 'cds': _CDS1},
            {'promoter': _P2, 'rbs': _LR2, 'cds': _CDS2}
            ],
        nlinkers=[_LN1, _LN2]
    )
    assert not o.has_duplicate_suffix()
    o = Construct(
        backbone=_BB, lms=_LMS, lmp=_LMP,
        blocks=[
            {'promoter': _P1, 'rbs': _LR1, 'cds': _CDS1},
            {'promoter': _P2, 'rbs': _LR1, 'cds': _CDS2}
            ],
        nlinkers=[_LN1, _LN2]
    )
    assert o.has_duplicate_suffix()


def test_get_part_ids():
    o = Construct(
        backbone=_BB, lms=_LMS, lmp=_LMP,
        blocks=[
            {'promoter': _P1, 'rbs': _LR1, 'cds': _CDS1},
            {'promoter': _P2, 'rbs': _LR2, 'cds': _CDS2}
            ],
        nlinkers=[_LN1, _LN2]
    )
    assert o.get_part_ids() == ['LMS', 'BB', 'LMP', 'P1', 'LR1', 'CDS1', 'LN1', 'P2', 'LR2', 'CDS2']


def test_get_sbol():
    o = Construct(
        backbone=_BB, lms=_LMS, lmp=_LMP,
        blocks=[
            {'promoter': _P1, 'rbs': _LR1, 'cds': _CDS1},
            {'promoter': _P2, 'rbs': _LR2, 'cds': _CDS2}
            ],
        nlinkers=[_LN1, _LN2]
    )
    doc = o.get_sbol(construct_id='TEST')
    # Item order are random in SBOL file
    query = '\n'.join(sorted(doc.writeString().split('\n')))
    assert sha1(query.encode()).hexdigest() == '59ea3b9f26e29546b62475799a401f09cfa0e6f8'
