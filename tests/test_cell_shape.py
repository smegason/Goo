import pytest
import bpy
from goo import Cell, CellType
from goo.force import * 
from goo.division import BisectDivisionLogic, BooleanDivisionLogic


@pytest.fixture
def setup_blender():
    bpy.ops.wm.read_factory_settings(use_empty=True)  # Reset to empty scene
    cellsA = CellType("A")
    cellsA.homo_adhesion_strength = 1000
    oriented_cell1 = cellsA.create_cell(name="A1", 
                                        loc=(10, 0, 0), 
                                        color=(0.5, 0.5, 0), 
                                        scale=(2, 1.5, 1.5))
    oriented_cell2 = cellsA.create_cell(name="A2", 
                                        loc=(-10, 0, 0), 
                                        color=(0.5, 0.5, 0), 
                                        scale=(1.5, 2, 1.5))
    oriented_cell3 = cellsA.create_cell(name="A3", 
                                        loc=(0, 0, 0), 
                                        color=(0.5, 0.5, 0), 
                                        scale=(1.5, 1.5, 2))
    
    oriented_cell4 = cellsA.create_cell(name="A4", 
                                        loc=(0, 10, 0), 
                                        color=(0.5, 0.5, 0), 
                                        scale=(1, 1, 1), 
                                        rotation=(45, 0, 0))
    
    oriented_cell5 = cellsA.create_cell(name="A5", 
                                        loc=(0, -10, 0), 
                                        color=(0.5, 0.5, 0), 
                                        scale=(1, 1, 1), 
                                        rotation=(0, 45, 0))
    
    yield oriented_cell1, oriented_cell2, oriented_cell3, oriented_cell4, oriented_cell5


def test_cell_major_axis(setup_blender): 
    cell1, cell2, cell3 = setup_blender[0:3]

    major1 = cell1.major_axis()
    major2 = cell2.major_axis()
    major3 = cell3.major_axis()

    assert (pytest.approx(major1._axis, abs=1e-6) == (-1, 0, 0) or (1, 0, 0))
    assert (pytest.approx(major2._axis, abs=1e-6) == (0, -1, 0) or (0, 1, 0))
    assert (pytest.approx(major3._axis, abs=1e-6) == (0, 0, -1) or (0, 0, 1))


def test_cell_major_cell_rotation(setup_blender):
    cell4, cell5 = setup_blender[3:5]
    major4 = cell4.major_axis()
    major5 = cell5.major_axis()

    assert (pytest.approx(major4._axis, abs=1e-6) == 
            (0.7071067811865476, 0, 0.7071067811865476) 
            or (-0.7071067811865476, 0, -0.7071067811865476))
    assert (pytest.approx(major5._axis, abs=1e-6) == 
            (0, 0.7071067811865476, 0.7071067811865476) 
            or (0, -0.7071067811865476, -0.7071067811865476))