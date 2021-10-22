"""
Test Designer
"""

import pytest

from pathlib import Path
from rpbasicdesign.Designer import Designer


_RPSBML_PATH = Path(__file__).resolve().parent / 'input' / 'muconate_example.xml'
_DATA_PATH = Path(__file__).resolve().parent / 'data'
_PARTS_FILES = [ _DATA_PATH / 'biolegio_parts.csv', _DATA_PATH / 'user_parts.csv' ]


def test_init_1():
    # Empty is OK
    Designer()
    # With some default values
    Designer(lms_id='LMS', lmp_id='LMP',
             backbone_id='BASIC_SEVA_37_CmR-p15A.1')

def test_init_2():
    # Not enough linkers / parts is not ok
    with pytest.raises(BaseException):
        Designer(lms_id='LMS', lmp_id='LMP',
                backbone_id='BASIC_SEVA_37_CmR-p15A.1',
                parts_files=_PARTS_FILES[0]
        )

def test_extract_from_rpsbml():
    o = Designer(lms_id='LMS', lmp_id='LMP',
                 backbone_id='BASIC_SEVA_37_CmR-p15A.1')
    o.enzyme_from_rpsbml(rpsbml_file=_RPSBML_PATH)


def test_combine():
    o = Designer(lms_id='LMS', lmp_id='LMP',
                 backbone_id='BASIC_SEVA_37_CmR-p15A.1')
    o.enzyme_from_rpsbml(rpsbml_file=_RPSBML_PATH)
    assert o.combine(sample_size=0, random_seed=42) == 0
    assert o.combine(sample_size=12, random_seed=42, cds_permutation=False) == 12
    assert isinstance(o.constructs, list)
    print(o.constructs[0].get_part_ids())
    assert o.constructs[0].get_part_ids() == [
        'LMS', 'BASIC_SEVA_37_CmR-p15A.1',
        'LMP', 'PJ23104_BASIC',
        'U1-RBS1', 'Q8XEC3',
        'U2-RBS1', 'P23262',
        'U3-RBS2', 'P95607']


def test_write_dnabot_inputs(tmp_path):
    o = Designer(lms_id='LMS', lmp_id='LMP',
                 backbone_id='BASIC_SEVA_37_CmR-p15A.1')
    o.enzyme_from_rpsbml(rpsbml_file=_RPSBML_PATH)
    o.combine(sample_size=10, random_seed=42)
    o.write_dnabot_inputs(out_dir=tmp_path)
    files = list(tmp_path.iterdir())
    assert tmp_path / 'constructs.csv' in files
    assert tmp_path / 'user_parts_plate.csv' in files


def test_write_sbol(tmp_path):
    o = Designer(lms_id='LMS', lmp_id='LMP',
                 backbone_id='BASIC_SEVA_37_CmR-p15A.1')
    o.enzyme_from_rpsbml(rpsbml_file=_RPSBML_PATH)
    o.combine(sample_size=10, random_seed=42)
    o.write_sbol(out_dir=tmp_path)
    files = list(tmp_path.iterdir())
    assert len(files) == 10
