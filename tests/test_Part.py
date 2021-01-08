"""
Test BASIC Part
"""
import pytest
from rpbasicdesign.BASICDesigner import Part


def test_init_1():
    # Empty is not OK
    with pytest.raises(TypeError):
        o = Part()


def test_init_2():
    # Args are just strings, keywords are not controlled
    o = Part(id='ID', basic_role='misc', biological_role='unknown')


def test_repr():
    o = Part(id='ID', basic_role='misc', biological_role='unknown')
    o.__repr__()


def test_get_prefix_suffix():
    # general case for linkers
    o = Part(id='IDx', basic_role='linker', biological_role='misc', linker_class='neutral linker')
    assert o.get_prefix_suffix() == ('IDx-P', 'IDx-S')
    # special case for RBS suffix
    o = Part(id='IDx-RBSy', basic_role='linker', biological_role='rbs', linker_class='rbs linker')
    assert o.get_prefix_suffix() == ('IDx-RBSy-P', 'IDx-S')
    # prefix and suffix are only of linkers
    o = Part(id='IDx', basic_role='misc', biological_role='unknown')
    assert o.get_prefix_suffix() is None


def test_get_sbol_id():
    o = Part(id='IDx', basic_role='misc', biological_role='unknown')
    assert o.get_sbol_id() == 'IDx'
    # non alphanumeric character forbidden
    o = Part(id='I-D?x!', basic_role='misc', biological_role='unknown')
    assert o.get_sbol_id() == 'I_D_x_'
