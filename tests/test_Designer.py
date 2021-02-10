"""
Test Designer
"""

from pathlib import Path
from rpbasicdesign.Designer import Designer


_RPSBML_PATH = Path(__file__).resolve().parent / 'input' / 'rp_1_3.sbml.xml'


def test_init_1():
    # Empty is OK
    Designer()
    # With some default values
    Designer(monocistronic=True, lms_id='LMS', lmp_id='LMP',
             backbone_id='BASIC_SEVA_37_CmR-p15A.1')


def test_extract_from_rpsbml():
    o = Designer(monocistronic=True, lms_id='LMS', lmp_id='LMP',
                 backbone_id='BASIC_SEVA_37_CmR-p15A.1')
    o.enzyme_from_rpsbml(rpsbml_file=_RPSBML_PATH)


def test_combine():
    o = Designer(monocistronic=True, lms_id='LMS', lmp_id='LMP',
                 backbone_id='BASIC_SEVA_37_CmR-p15A.1')
    o.enzyme_from_rpsbml(rpsbml_file=_RPSBML_PATH)
    assert o.combine(sample_size=0, random_seed=42) == 0
    assert o.combine(sample_size=12, random_seed=42) == 12
    assert isinstance(o.constructs, list)
    assert o.constructs[0].get_part_ids() == [
        'LMS', 'BASIC_SEVA_37_CmR-p15A.1', 'LMP', 'PJ23119_BASIC',
        'U1-RBS1', 'D2WKD9', 'L1', 'PJ23104_BASIC',
        'U2-RBS3', 'O48935', 'L2', 'PJ23101_BASIC',
        'U3-RBS1', 'O66952']


def test_write_dnabot_inputs(tmp_path):
    o = Designer(monocistronic=True, lms_id='LMS', lmp_id='LMP',
                 backbone_id='BASIC_SEVA_37_CmR-p15A.1')
    o.enzyme_from_rpsbml(rpsbml_file=_RPSBML_PATH)
    o.combine(sample_size=10, random_seed=42)
    o.write_dnabot_inputs(out_dir=tmp_path)
    files = list(tmp_path.iterdir())
    assert tmp_path / 'constructs.csv' in files
    assert tmp_path / 'user_parts_coord.csv' in files


def test_write_sbol(tmp_path):
    o = Designer(monocistronic=True, lms_id='LMS', lmp_id='LMP',
                 backbone_id='BASIC_SEVA_37_CmR-p15A.1')
    o.enzyme_from_rpsbml(rpsbml_file=_RPSBML_PATH)
    o.combine(sample_size=10, random_seed=42)
    o.write_sbol(out_dir=tmp_path)
    files = list(tmp_path.iterdir())
    assert len(files) == 10
